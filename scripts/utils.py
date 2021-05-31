#!/usr/bin/env python3
# Lanxin Zhang

import math
import numpy as np
import configparser
import os
import scipy.stats as st


def getPropensities(Y, D, if_ll=False):
    """reaction propensities, DTG specific

    return value a: dict of propensities, i.e. a_k [1/hr]

    Y: viral dynamical state, np array, [V, T1, T2]
    D: plasma concentration of DTG [nM]
    if_ll: if long-lived and latent infected cells are considered
    """

    config = configparser.ConfigParser()
    filepath = os.path.abspath(os.path.join(os.path.dirname( __file__ ),
                               os.pardir, 'config', 'config.ini'))
    config.read(filepath)



    P_Ma4 = float(config['viraldynamics']['P_Ma4'])
    P_La5 = float(config['viraldynamics']['P_La5'])
    IC50 = float(config['viraldynamics']['IC50'])
    m = float(config['viraldynamics']['m'])
    CL = float(config['viraldynamics']['CL'])
    rho = float(config['viraldynamics']['rho'])
    beta = float(config['viraldynamics']['beta'])
    lambda_T = float(config['viraldynamics']['lambda_T'])
    delta_T = float(config['viraldynamics']['delta_T'])
    delta_T1 = float(config['viraldynamics']['delta_T1'])
    delta_T2 = float(config['viraldynamics']['delta_T2'])
    delta_PIC = float(config['viraldynamics']['delta_PIC'])
    k = float(config['viraldynamics']['k'])
    N_T = float(config['viraldynamics']['N_T'])

    eta = D ** m / (IC50 ** m + D ** m)     # direct effect of DTG
    CL_T = (1 / rho - 1) * beta # clearence rate of unsuccessful infected virus
    T_u = lambda_T / delta_T    # steady state level of uninfected T-cells
    a = [0] * 9                 # list of propensities
    a[1] = (CL + CL_T * T_u) * Y[0] / 24        # V -> *
    a[2] = (delta_PIC + delta_T1) * Y[1] / 24   # T1 -> *
    a[3] = delta_T2 * Y[2] / 24                 # T2 -> *
    a[4] = beta * T_u * Y[0] / 24               # V -> T1
    a[5] = (1 - eta) * k * Y[1] / 24            # T1 -> T2
    a[6] = N_T * Y[2] / 24                      # T2 -> T2 + V
    a[7] = 0                                    # V -> long-lived infected cell
    a[8] = 0                                    # T1 -> latent infected cell

    # if long-lived and latent infected cells are considered
    if if_ll:
        a[7] = 1e-4 * Y[0] / 24
        a[8] = a[5] * P_La5
        #a[4] = a[4] * (1 - P_Ma4)
        a[5] = a[5] * (1 - P_La5)

    return a


def getExtinctionProbability(Y, D, if_ll=False):
    """Calculate the extinction probability (DTG specific).

    return value PE: extinction probability

    Y: viral dynamical state, np array, [V, T1, T2]
    D: plasma concentration of DTG [nM]
    if_ll: if long-lived and latent infected cells are considered
    """

    def solution_quadratic_equation(a, b, c):
        # solve ax^2 + bx + c = 0
        tmp = b * b - 4 * a * c
        assert tmp >= 0
        tmp_sqr = tmp ** 0.5
        return (-b - tmp_sqr) / (2 * a)

    a = getPropensities([1, 1, 1], D, if_ll)
    alpha1 = a[1] / (a[1] + a[4] + a[7])
    alpha2 = a[7] / (a[1] + a[4] + a[7])
    alpha3 = a[4] / (a[1] + a[4] + a[7])
    beta1 = a[2] / (a[2] + a[5] + a[8])
    beta2 = a[8] / (a[2] + a[5] + a[8])
    beta3 = a[5] / (a[2] + a[5] + a[8])
    gamma = a[6] / (a[3] + a[6])
    e_a = alpha3 * beta3 * gamma
    e_b = gamma - alpha2 * gamma - beta2 * gamma - 1 - e_a
    e_c = 1 - gamma
    PE_t2 = solution_quadratic_equation(e_a, e_b, e_c)
    PE_t1 = beta1 + beta3 * PE_t2
    PE_v = alpha1 + alpha3 * PE_t1

    # extinction probability
    PE = PE_v ** Y[0] * PE_t1 ** Y[1] * PE_t2 ** Y[2]

    return PE, [PE_v, PE_t1, PE_t2]

def transmittedVirusCount(is_homo=True):
    """Get the number virus entering a replication-competent compartment
    according to the transmission mode.

    return value: number of virus entering a replication-competent compartment

    is_homo: boolean, if transmission mode is homosexual
    """

    config = configparser.ConfigParser()
    filepath = os.path.abspath(os.path.join(os.path.dirname( __file__ ),
                               os.pardir, 'config', 'config.ini'))
    config.read(filepath)

    VL_mean = float(config['inoculumsize']['VL_mean'])
    VL_std = float(config['inoculumsize']['VL_std'])
    c = float(config['inoculumsize']['c'])
    VL = 10 ** (np.random.normal(VL_mean, VL_std))  # VL:log-normal distributed
    f_VL = math.ceil(VL ** c)   # n in binomial distribution

    # transmission probability according to different transmission mode
    # (p in binomial distribution)
    if is_homo:
        r0 = float(config['inoculumsize']['r_homo'])
    else:
        r0 = float(config['inoculumsize']['r_hetero'])

    n_virus = 0
    # simulation of Bernoulli Experiment
    for i in range(f_VL):
        r = np.random.random()
        if r <= r0:
            n_virus += 1
    return n_virus


def calculate_CI(n_extinction, n_infection, p_value=0.95, Y=[1,0,0]):
    """Compute the confidential interval using Greenwoods formula for Extrande

           return: confidential interval for PE and efficacy
    """
    p_value_1tile = 1 - (1 - p_value) / 2
    z_score = st.norm.ppf(p_value_1tile)
    variance = (n_extinction / (n_extinction + n_infection)) **2 * (1 / n_extinction - 1 / (n_extinction + n_infection))
    ci = z_score * variance ** 0.5
    p0 = getExtinctionProbability(Y, 0)[0]
    ci_phi = ci / (1 - p0)
    return ci, ci_phi
