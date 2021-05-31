#!/usr/bin/env python3
# Lanxin Zhang

import numpy as np
import pandas as pd
import math, sys
from functools import partial
import multiprocessing
from enum import Enum
from .sample import Sample
from .utils import *


class reaction(Enum):
    VCLEAR = 1
    T1CLEAR = 2
    T2CLEAR = 3
    VTOT1 = 4
    T1TOT2 = 5
    T2TOV = 6
    LONGLIVED = 7
    LATENT = 8


def extrande(sample, m_dose, Y, adherence, t, count_doses, threshold = 1e-4):
    """run EXTRANDE with initial state

    return True if extinction occurs, return False if infection occurs
    Args:
    sanple: a Sample object which contains PK parameters
    m_dose: mass of oral taken DTG of one dose [mg]
    Y: initial state [V, T1,T2] (np array)
    adherence: possibility to take drug
    t: time point at which the simulation begins (exposure occurs),
          relative to the time 0, on which the first dose is taken[hr]
          (can be negative, before first dose is taken)
    count_doses: count of doses in a regimen
    threshold: stopping criteria (out of simplex)
    """
    
    np.random.seed()
    sample.reset()
    m_dose *= 1e6     # convert the mass into ng
    P_Ma4 = 1.25e-4   # probability of long-lived infected cell if a4 fires
    P_La5 = 8e-6      # probability of latent infected cell if a5 fires

    # stoichiometric matrix of the viral compartment
    #     r1 r2 r3 r4 r5 r6
    S = np.array([
           [-1, 0, 0, -1, 0, 1],  # Viruses
           [0, -1, 0, 1, -1, 0],  # T1 cells
           [0, 0, -1, 0, 1, 0]    #T2 cells
            ])

    while True:
        B = sum(getPropensities(Y, 0))  # upper bound of propensities
        # time to the next reaction event is exponentially distributed with B
        tau = np.random.exponential(1/B)
        t0 = t
        t = t + tau     # update time

        # plasma concentration of DTG at time t
        if m_dose > 0:
            Dt = sample.getConcentration(m_dose, t, adherence, count_doses)
        else:
            Dt = 0

        a0 = getPropensities(Y, Dt)     # reaction propensities at time t
        r = np.random.random()      # random number from uniform distribution
        if sum(a0) >= B * r:       # one of the six reaction fires
            tmp, j = 0, 0
            for i in range(1, 7):       # choose the fired reaction
                tmp += a0[i]
                if tmp >= B * r:
                    j = i
                    break

            if j == reaction.T1TOT2.value:      # reaction 5 is chosen
                r_latent = np.random.random()
                if r_latent <= P_La5:   # latently infected cells emerged
                    return False

            elif j == reaction.VTOT1.value:      # reaction 4 is chosen
                r_M1 = np.random.random()
                if r_M1 <= P_Ma4:       # long lived infected cell emerged
                    return False

            Y0 = Y.copy()
            Y = Y + S[:, j-1]           # update the state of system

            if not np.any(Y):             # Extinction event occurs
                return True

            # compute the extinction probability of state Y and Dmax
            else:
                Dmax = sample.getMaximum(m_dose, t, adherence, count_doses)
                PE = getExtinctionProbability(Y, Dmax)[0]
                if PE < threshold:      # Infection event occurs
                    return False
                else:   # trajectory at time t within the extinction simplex
                    pass
        else:       # thinning reaction fires without changing the system
            pass


def runSingleSimulation(sample, tps, m, count_doses, adherence, vload):
    # run extrande over sample for timepoins in tps,
    # help function for multiprocessing

    try:
        v0 = int(vload)
        Y = np.array([v0, 0, 0])
    except:
        ishomo = True if vload == 'homo' else False
        Y = None

    result = {t:[0, 0] for t in tps}
    for t in tps:
        if Y is None:
            v0 = transmittedVirusCount(ishomo)
            if v0:
                Y = np.array([v0, 0, 0])
            else:
                result[t][0] += 1
                continue

        if extrande(sample, m, Y, adherence, t, count_doses):
            result[t][0] += 1
        else:
            result[t][1] += 1

    return result

def  runSimulations(cycle, samples, **kwargs):
    # run simulations with multiprocessing to speed up
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=cores)
    helpfunc = partial(runSingleSimulation, **kwargs)
    result = {t:[0,0] for t in kwargs['tps']}
    for cyc in range(cycle):
        cnt = 0
        for i in pool.imap_unordered(helpfunc, samples):
            for key in i:
                result[key][0] += i[key][0]
                result[key][1] += i[key][1]
            cnt += 1
            sys.stdout.write('done {}/{}\r'.format(cnt, len(samples)))

        print('Cycle {} done'.format(cyc+1))
    print('All simulations done')
    res = pd.DataFrame(result)
    return res

def run_extrande(mdose, cdose, adh, csimul,tps, vload, inputfile):
    param = {'tps':tps, 'm':mdose, 'count_doses':cdose,
             'adherence':adh, 'vload':vload}
    df_sample = pd.read_csv(inputfile)
    samples = [Sample(*df_sample.loc[i,:].values)
               for i in range(df_sample.shape[0])]
    c_samples = len(samples)
    if c_samples < 1000:        # build large sample set for multiprocessing
        samples += samples * (1000 // c_samples - 1)

    cycle = csimul // len(samples)
    res_df = runSimulations(cycle, samples, **param)
    return res_df
