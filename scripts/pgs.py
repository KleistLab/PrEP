#!/usr/bin/env python3
# Lanxin Zhang

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import sys
from .sample import Sample
from .utils import *


def get_a5(sample, t, m_dose, count_doses, adherence):

    """
    function to get the value of a5, for solving ODEs

    return the value of a5 at time point t
    """

    D = sample.getConcentration(m_dose, t, adherence, count_doses)
    return getPropensities([1,1,1], D)[5]


def peModel(t, y, a, sample, m_dose, count_doses, adherence, if_ll=False):

    """
    ODE of Extinction probability model, for solving ODE with scipy

    return a vector of the right hand side of the ODEs
    t: time
    y: vector [P1, P2, P3]
    a: propensities
    sample: a Sample object which contains PK parameters
    m_dose: mass of oral taken DTG of one dose [mg]
    count_doses: count of doses in a regimen
    adherence: % of pills taken
    if_ll: if long-lived and latent infected cells are considered
    """

    P_Ma4 = 1.25e-4 if if_ll else 0
    P_La5 = 8e-6 if if_ll else 0    # probability of latent infected cell
    a1 = a[1]
    a2 = a[2]
    a3 = a[3]
    a4 = a[4]
    #a5 = 0.0145
    a5 = get_a5(sample, t, m_dose, count_doses, adherence)
    a6 = a[6]
    P1, P2, P3 = y[0], y[1], y[2]
    dydt = [-a1*(1 - P1) - a4*(P2 - P1 / (1 - P_Ma4)),
                    - a2*(1 - P2) - a5*(P3 - P2 / (1 - P_La5)),
                    - a3*(1 - P3) - a6*(P3*P1 - P3)]

    return dydt


def solvePE_pgs(sample, t_start, t_end, m_dose, count_doses, adherence,
                if_ll=False, p10=0.9017, p20=0.5212, p30=0.01496, t_extend=100):

    """
    solve ODE for extinction probability based on a concentration-time curve

    return bunch object containing solutions

    sample: a Sample object which contains PK parameters
    t_start, t_end: time span of interest [hr],
            can be negative if PEP is required. [hr]
    m_dose: mass of oral taken DTG of one dose [mg]
    count_doses: count of doses in a regimen
    adherence: % of pills taken
    if_ll: if long-lived and latent infected cells are considered
    p10, p20, p30: initial condition (extinction probabilities) for solving ODE
    (here constants are used, but I just keep these as input for the further
    functional extension)
    """

    # reset the sample, in case the concentration profile already exists
    sample.reset()
    m_dose *= 1e6     # convert the mass into ng
    t_end += t_extend
    a = getPropensities([1,1,1], 0)
    t = [t_end, t_start]
    y0 = [p10, p20, p30]
    sol = solve_ivp(peModel, t, y0, method='LSODA',
                    args=[a, sample, m_dose, count_doses, adherence, if_ll],
                    dense_output=True)
    return sol


def run_pgs(mdose, cdose, adh, inputfile, ts, te, tspanres, if_ll):
    df_sample = pd.read_csv(inputfile)
    if df_sample.shape[0] != 1:
        sys.stderr.write('This method takes only one sample, more than one '
                         'are given \n')
        exit(1)
    sample = Sample(*df_sample.loc[0, :].values)
    solution = solvePE_pgs(sample, ts, te+100, mdose, cdose, adh, if_ll).sol
    p1 = [solution(i/60)[0] for i in np.arange(ts*60, te*60, int(tspanres*60))]
    p2 = [solution(i/60)[1] for i in np.arange(ts*60, te*60, int(tspanres*60))]
    p3 = [solution(i/60)[2] for i in np.arange(ts*60, te*60, int(tspanres*60))]
    return p1, p2, p3
