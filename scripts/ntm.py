#!/usr/bin/env python3
# Lanxin Zhang

import numpy as np
import pandas as pd
import sys
from .sample import Sample
from .utils import *


def solvePE_ntm(sample, t_start, t_end, m_dose, count_doses,
                    adherence, timestep, timestep_p3, recall_rate,
                    if_ll=False, t_extend=100):

    """
     compute the extinction probabilities for V, T1 and T2
     t_start, t_end: time inteval of interest,
                    can be negative if PEP is required. [hr]
     timestep: time step for p1 and p2 [min]
     timestep_p3: time step for p3 [min]
     m_dose: mass of oral taken DTG of one dose [mg]
     count_doses: count of doses in a regimen
     adherence: % of pills taken
     recall_rate: delta
     if_ll: if long-lived and latent infected cells are considered
    """

    m_dose *= 1e6     # convert the mass into ng
    P_La5 = 8e-6 if if_ll else 0    # probability of latent infected cell
    # reset the sample, in case the concentration profile already exists
    sample.reset()
    t_end += t_extend

    def computeSteps(recall_rate, timestep):
        """
        compute the look-forward horizons according to delta
         k1, k2, k3: timespan for integration [min]
        """
        a = getPropensities([1,1,1], 0)
        # compute the first term and reuse it
        k1 = math.ceil(-math.log(1 - recall_rate) / (a[1] + a[4]) / timestep*60)
        k2 = math.ceil(-math.log(1 - recall_rate) / a[2] / timestep * 60)
        k3 = math.ceil(-math.log(1 - recall_rate) / (a[6] + a[3]) / timestep*60)
        return k1, k2, k3

    k1, k2, k3 = computeSteps(recall_rate, timestep)    # look-forward horizons
    if t_start < 0:     # if PEP is required
        pre_duration= int(-t_start * 60 / timestep)
        t_start = 0
    else:
        pre_duration = 0

    # if the timespan of interest is not long enough for the given recall rate
    if k2 > (t_end - t_start) * 60 / timestep:
        t_end = int(t_start + k2 / 60 * timestep)

    # the count of cycle to update p3 in one minute
    if timestep_p3 != timestep:
        cycle_p3 = int(1 / timestep_p3)
    else:
        cycle_p3 = 1

    # an array of drug concentrations in time-order, each step represents 1 min
    sample_c = [0] * pre_duration + [sample.getConcentration(m_dose,
                                     t / 60, adherence, count_doses)
                                     for t in np.arange(t_start*60,
                                                        t_end*60+k2*timestep,
                                                        timestep)]
    # propensities
    sample_a0 = getPropensities([1, 1, 1], 0, if_ll=if_ll)
    # extinction probabilities without drug, for initialization
    PE_v_0, PE_t1_0, p30 = getExtinctionProbability([1, 0, 0], 0, if_ll=if_ll)[1]
    # initialize the output arrays
    p1 = [PE_v_0] * (len(sample_c) + 1)
    p2 = [PE_t1_0] * (len(sample_c) + 1)
    p3 = [p30] * (len(sample_c) + 1)

    # an array of a5 values, computed from sample_c
    a5 =[getPropensities([1, 1, 1], max(i, 0))[5] for i in sample_c]

    # compute an array of d1, with same length of k1:
    sum_a147 = sample_a0[1] + sample_a0[4] + sample_a0[7]
    d1 = np.array([sample_a0[4] / sum_a147 * np.exp(- sum_a147 * i/60*timestep)
                   * (1 - np.exp(- sum_a147 / 60*timestep)) for i in range(k1)])

    # compute an array of d3, with same length of k3*cycle_p3 :
    sum_a36 = sample_a0[3] + sample_a0[6]
    d3 = np.array([sample_a0[6] / sum_a36 * np.exp(- sum_a36 * i/60*timestep_p3)
                   * (1 - np.exp(- sum_a36 / 60 * timestep_p3))
                   for i in range(k3 * cycle_p3)])

    # integral of a5 with same length of sample_c
    a5_integral = []
    for i in a5:
        try:
            a5_integral.append(a5_integral[-1] + i / 60 * timestep)
        except:
            a5_integral.append(i / 60 * timestep)
    a5_integral = np.array(a5_integral) / (1-P_La5)

    # the term involves a5 in d2: a5(t + i*dt)*exp(-integral(0, t+i*dt) a5(x)dx)
    # unit of a5: [1/min]
    a5_term = np.array(a5) / 60 * (np.exp(- np.array(a5_integral)))

    for i in range(len(sample_c)-k2-1, -1, -1):
        p1[i] = sample_a0[1] / (sample_a0[1] + sample_a0[4] + sample_a0[7]) + \
                np.dot(d1, p2[i + 1: i + 1 + k1])
        d2 = a5_term[i: i + k2] * np.exp(a5_integral[i]) * timestep * \
             (np.exp(- (sample_a0[2] + sample_a0[8]) *
                     np.arange(0, k2) / 60 * timestep))
        p2[i] = (1 - sum(d2) / (1-P_La5)) + np.dot(d2, p3[i + 1: i + 1 + k2])

        # updating p3: since the value of a6 is significantly larger than
        # other propensities, the delta t for p3 should be smaller, otherwise
        # the error of p3 will be relatively large
        tmp = p3[i+1]
        for j in range(cycle_p3):
            tmp = sample_a0[3] / (sample_a0[3] + sample_a0[6]) + \
                  np.sum(d3 * p1[i + 1] * tmp)
        p3[i] = tmp
    return p1, p2, p3


def run_ntm(mdose, cdose, adh, inputfile, ts, te, timesteps, rate, if_ll):
    df_sample = pd.read_csv(inputfile)
    if df_sample.shape[0] != 1:
        sys.stderr.write('This method takes only one sample, more than one '
                         'are given \n')
        exit(1)
    sample = Sample(*df_sample.loc[0, :].values)
    p1, p2, p3 = solvePE_ntm(sample, ts, te, mdose, cdose, adh,
                             timesteps[0], timesteps[2], rate, if_ll)
    return p1, p2, p3
