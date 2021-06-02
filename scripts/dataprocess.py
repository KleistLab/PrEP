#!/usr/bin/env python3
# Lanxin Zhang

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os.path
from .utils import *


def extrande_output(df, file):
    # process the output dataframe df from run_extrande and output a csv file
    # with name file

    p0 = getExtinctionProbability([1,0,0], 0, True)[1][0]
    n_total = df.iloc[0,0] + df.iloc[1, 0]
    df.index = ['n_extinction', 'n_infection']
    df.loc['PE'] = df.loc['n_extinction', :] / (df.loc['n_extinction', :] + \
                                                df.loc['n_infection',:])
    df.loc['efficacy'] = (1 - (1 - df.loc['PE', :])/(1 - p0))*100
    df.loc['variance', :] = df.loc['PE', :]**2 * \
                            (1/df.loc['n_extinction',:] - 1/n_total)
    df.loc['95CL', :]  = 1.96 * df.loc['variance', :]**0.5 / (1 - p0)*100
    df.to_csv(file)
    return df


def extrande_plot(df, file, mdose):
    # plot the data df which is processed by output_extrande
    # and output a pdf file

    t = [float(i) for i in df.columns]
    phi = df.loc['efficacy', :]
    errorbar = df.loc['95CL', :]
    fig, ax = plt.subplots()

    ax.errorbar(t, phi, yerr=errorbar, capsize=2, lw=1, c='crimson')
    ax.set_xticks(t)
    ax.set_xlabel('Timing of viral exposure relative to first dose [hr]')
    ax.set_ylabel('Prophylactic efficacy [%]')
    ax.set_title('Prophylactic efficacy for {} mg DTG'.format(mdose))
    fig.savefig('{}'.format(file))


def deterministics_output(p1, p2, p3, file, tspanres, timesteps, ts, te):
    # process the output arrays from run_ctsm and output a csv file
    # with name file

    t = np.arange(ts, te, tspanres)
    step = int(tspanres/timesteps[0]*60)
    p0 = getExtinctionProbability([1,0,0], 0, True)[1][0]
    efficacy = [1 - (1 - p1[i])/(1-p0) for i in range(0, len(p1), step)]
    p1 = p1[::step]
    p2 = p2[::step]
    p3 = p3[::step]
    if len(efficacy) > len(t):
        efficacy = efficacy[:len(t)]
        p1 = p1[:len(t)]
        p2 = p2[:len(t)]
        p3 = p3[:len(t)]
    df = pd.DataFrame([t, p1, p2, p3, efficacy])
    df.index = ['time point', 'PE_V', 'PE_T1', 'PE_T2', 'efficacy']
    df.to_csv(file)
    return df


def deterministics_plot(df, file, mdose):
    t = df.loc['time point', :] / 24
    days = int(t[len(t)-1] - t[0]) + 1
    if days <= 1:
        ticks = 2
        t *= 24
    elif days < 10:
        ticks = 1
    else:
        ticks = days // 10
    fig, axs = plt.subplots(1, 2, figsize=[15,6])
    axs[0].plot(t, df.loc['efficacy', :]*100, lw=1, c='crimson')
    axs[0].set_ylabel('Prophylactic efficacy [%]')
    axs[0].set_title('Prophylactic efficacy for {} mg DTG'.format(mdose))
    axs[0].set_xticks(np.arange(int(t[0]), int(t[len(t)-1])+2, ticks))
    if days <= 1:
        axs[0].set_xlabel('Timing of viral exposure relative to first dose [hr]')
    else:
        axs[0].set_xlabel('Timing of viral exposure relative to first dose [day]')

    axs[1].plot(t, df.loc['PE_V', :]*100, lw=1, c='crimson')
    axs[1].plot(t, df.loc['PE_T1', :]*100, lw=1, c='forestgreen')
    axs[1].plot(t, df.loc['PE_T2', :]*100, lw=1, c='dodgerblue')
    axs[1].set_xticks(np.arange(int(t[0]), int(t[len(t)-1])+2, ticks))
    axs[1].set_ylabel('Extinction probability [%]')
    axs[1].set_title('Extinction probability for {} mg DTG'.format(mdose))
    axs[1].legend(('$P_{E\_V}$', '$P_{E\_T_1}$', '$P_{E\_T_2}$'))
    if days <= 1:
        axs[1].set_xlabel('Timing of viral exposure relative to first dose [hr]')
    else:
        axs[1].set_xlabel('Timing of viral exposure relative to first dose [day]')
    fig.savefig('{}'.format(file))
