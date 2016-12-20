#!/usr/bin/env python
#
# tools.py
#
# (C) The James Hutton Institute 2016
# Author: Leighton Pritchard

"""
tools.py

This module provides helper functions used in the supplementary information
notebooks and scripts for the Holmes et al. (2017) paper.
"""

from matplotlib import pyplot as plt

import numpy as np
import pandas as pd
import scipy
import seaborn as sns


# PRNG seed
SEED = 123456789


def corrfunc(x, y, **kws):
    """Return a matplotlib axis with text describing the Spearman
    correlation coefficient for x and y

    This function is written to support plot_correlation
    """
    coeff, _ = scipy.stats.spearmanr(x, y)
    ax = plt.gca()
    ax.annotate("r = {:.3f}".format(coeff),
                xy=(.3, .5), size=15,
                xycoords=ax.transAxes)


def plot_correlation(df, title=None):
    """Render Seaborn PairGrid of columns in df, with Pearson correlation
    coefficients in the upper triangle, and KDE plots on the diagonal.
    """
    g = sns.PairGrid(df)
    g.map_lower(plt.scatter)
    g.map_diag(sns.kdeplot, legend=False)
    g.map_upper(corrfunc)
    g.set(xticklabels=[])
    g.set(title=title or '')
    return g


def quantile_norm(df, columns=None):
    """Normalise the columns of df to each have the same distribution"""
    df_matrix = df.as_matrix(columns=columns)
    quantiles = np.mean(np.sort(df_matrix, axis=0), axis=1)
    ranks = scipy.stats.mstats.rankdata(df_matrix, axis=0).astype(int) - 1
    norm_matrix = quantiles[ranks]
    return(pd.DataFrame(data=norm_matrix, index=df.index,
                        columns=columns or df.columns))


def plot_normalised(ctl_in, ctl_out, trt_in, trt_out):
    """Return violin plots of input/output control/treatment distributions"""
    fig, axes = plt.subplots(2, 2, figsize=(12,6))
    fig.subplots_adjust(hspace=.25)
    axes = axes.ravel()
    for ttl, arr, ax in zip(("control input", "control output",
                             "treatment input", "treatment output"),
                            (ctl_in, ctl_out, trt_in, trt_out),
                            axes):
        ax.set_title(ttl)
        sns.violinplot(np.log(arr), ax=ax)


def wide_to_long_df(df, stage):
    """Convert wide dataframe to long

    This function is brittle, and only for Holmes et al SI
    """
    if not stage:
        stagestr = 'input'
    else:
        stagestr = 'output'

    df.reset_index(level=0, inplace=True)  # make probes a column
    df = pd.melt(df, id_vars=['Systematic'],
                 value_vars=['{0}.1'.format(stagestr),
                             '{0}.2'.format(stagestr),
                             '{0}.3'.format(stagestr)])
    df.columns = ['probe', 'class', stagestr]
    df.loc[:, 'replicate'] = df['class'].astype(str).str[-1].astype(np.int64)
    df = df[['probe', 'replicate', stagestr]]
    df.set_index(['probe', 'replicate'], inplace=True)
    return df
        

def wide_to_long_join(df_in, df_out, treatment):
    """Convert two wide dataframes to long and join on common index

    This function is brittle and only for Holmes et al SI
    """
    if treatment:
        treatval = 1
    else:
        treatval = 0            
    df = pd.merge(wide_to_long_df(df_in, 0), wide_to_long_df(df_out, 1),
                  left_index=True, right_index=True)
    df['treatment'] = treatval
    df.reset_index(inplace=True)
    return df


def wide_to_long(ctl_in, ctl_out, trt_in, trt_out):
    """Convert four dataframes from wide to long format

    This function returns a dataframe with columns:
    
    * probe
    * replicate
    * treatment
    * input
    * output
    * log_input
    * log_output
    """
    ctl_long = wide_to_long_join(ctl_in, ctl_out, treatment=False)
    trt_long = wide_to_long_join(trt_in, trt_out, treatment=False)
    data = ctl_long.append(trt_long, ignore_index=True)
    data['log_input'] = np.log(data['input'])
    data['log_output'] = np.log(data['output'])
    data = data[['probe', 'replicate', 'treatment',
                 'input', 'output',
                 'log_input', 'log_output']]
    return data


def plot_input_output_violin(data):
    """Plot Seaborn violin plot of log input and output data"""
    input_v_output = pd.melt(data,
                             id_vars=['probe', 'replicate', 'treatment'],
                             value_vars=['log_input', 'log_output'])
    input_v_output.columns = ['probe', 'replicate', 'treatment',
                              'stage', 'log_intensity']
    g = sns.violinplot(data=input_v_output, x="treatment", y="log_intensity",
                       hue="stage", split=True)
    g.set_xticklabels(['control', 'treatment'])
    g.set_ylabel("log(intensity)")
    g.set_xlabel("")
    g.set_title("log(intensity) distribution by treatment and input/output")
