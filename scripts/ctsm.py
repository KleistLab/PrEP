#!/usr/bin/env python3
# Lanxin Zhang

import numpy as np
import pandas as pd
import sys
from .sample import Sample
from .utils import *


def solvePE_ctsm(sample, t_start, t_end, m_dose, count_doses, adherence,
                 timestep, timestep_p3, if_ll=False, p10=0.9017, p20=0.5212, p30=0.01496, t_extend=100):

    """solve extinction probability using constant time step method

    return three lists containing extinction probabilities for V, T1, T2

    sample: a Sample object which contains PK parameters
    t_start, t_end: time span of interest [hr],
                can be negative if PEP is required. [hr]
    m_dose: mass of oral taken DTG of one dose [mg]
    count_doses: count of doses in a regimen
    adherence: % of pills taken
    timestep: time step for p1 and p2 [min]
    timestep_p3: time step for p3 [min]
    if_ll: if long-lived and latent infected cells are considered
    p10, p20, p30: initial values
    t_extend: duration to extend the end time [hr]
    """

    m_dose *= 1e6     # convert the mass into ng
    P_La5 = 8e-6 if if_ll else 0    # probability of latent infected cell

    # reset the sample, in case the concentration profile already exists
    sample.reset()

    a0 = getPropensities([1,1,1], 0, if_ll)        # propensities with absence of DTG

    # list of concentrations, from time point t_start to t_end with timestep
    if timestep < 1:
        cycle = int(1 / timestep)
        sample_c = list()
        for i in range(t_start*60, t_end*60):
            for j in range(cycle):
                sample_c.append(sample.getConcentration(m_dose,
                                                        (i+j/60)/60,
                                                        adherence,
                                                        count_doses))
    else:
        sample_c = [sample.getConcentration(m_dose,i/60,adherence,count_doses)\
                    for i in np.arange(t_start*60, (t_end+t_extend)*60, timestep)]
    # list of a5 values, calculated from the concentration list
    a5 =[getPropensities([1,1,1], max(i, 0), if_ll)[5] for i in sample_c]
    # extinction probabilities with absence of DTG
    # initialize three lists for PE with  P10, P20, P30
    p1 = [p10] * len(sample_c)
    p2 = [p20] * len(sample_c)
    p3 = [p30] * len(sample_c)

    # some constants used in DP
    count_t = 60 / timestep     # count of steps in an hour
    sum_a147 = a0[1] + a0[4] + a0[7]
    e147 = np.exp(-sum_a147/count_t)    # divied by count_t = multiplied by dt
    sum_a36 = a0[3] + a0[6]

    # DP backwards
    for i in range(len(p1)-2, -1, -1):
        p1[i] = a0[1] / sum_a147 * (1 - e147) + e147 * p1[i + 1] + \
                (1 - e147) * (a0[4] / sum_a147) * p2[i + 1]
        sum_a258 = a0[2] + a5[i] / (1-P_La5)
        e258 = np.exp(-sum_a258/count_t)
        p2[i] = a0[2] / sum_a258 * (1 - e258) + e258 * p2[i + 1] + \
                (1 - e258) * (a5[i] / sum_a258) * p3[i + 1]

        if timestep_p3:
        # for p3: since the value of a6 is significantly larger than
        # other propensities, time step for p3 should be smaller,
        # otherwise the error of p3 will be relatively large
            cycle_p3 = int(1/timestep_p3)
            e36 = np.exp(-sum_a36/(60*cycle_p3))
            tmp = p3[i + 1]
            for _ in range(cycle_p3):
                tmp = a0[3] / sum_a36 * (1 - e36) + e36 * tmp + \
                      (1 - e36) * (a0[6] / sum_a36) * tmp * p1[i + 1]
            p3[i] = tmp
        else:
            e36 = np.exp(-sum_a36/count_t)
            p3[i] = a0[3] / sum_a36 * (1 - e36) + e36 * p3[i + 1] + \
                    (1 - e36) * (a0[6] / sum_a36) * p3[i + 1] * p1[i + 1]

    return p1, p2, p3


def run_ctsm(mdose, cdose, adh, inputfile, ts, te, timesteps, if_ll):
    df_sample = pd.read_csv(inputfile)
    if df_sample.shape[0] != 1:
        sys.stderr.write('This method takes only one sample, more than one '
                         'are given \n')
        exit(1)
    sample = Sample(*df_sample.loc[0, :].values)
    p1, p2, p3 = solvePE_ctsm(sample, ts, te, mdose, cdose, adh, timesteps[0],
                              timesteps[2], if_ll)
    return p1, p2, p3
