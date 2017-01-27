#!/usr/bin/env python3.5
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
import pickle
import random
import scipy
import seaborn as sns

from collections import defaultdict

from Bio import SeqIO

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
    * repXtrt (combination of replicate and treatment)
    * input
    * output
    * log_input
    * log_output
    """
    ctl_long = wide_to_long_join(ctl_in, ctl_out, treatment=False)
    trt_long = wide_to_long_join(trt_in, trt_out, treatment=True)
    data = ctl_long.append(trt_long, ignore_index=True)
    data['log_input'] = np.log(data['input'])
    data['log_output'] = np.log(data['output'])
    data['repXtrt'] = 'rep' + data['replicate'].map(str) +\
                      'trt' + data['treatment'].map(str)
    data = data[['probe',
                 'replicate', 'treatment', 'repXtrt',
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


def unique_probe_matches(blastfiles):
    """Returns a dataframe of unique queries and their unique matches"""
    # Columns in a BLASTN+ -outfmt 6 file
    blast_columns = ['probe', 'match', 'identity', 'length', 'mismatch',
                     'gapopen', 'qstart', 'qend', 'sstart', 'send',
                     'evalue', 'bitscore']
    df = None
    for bfile in blastfiles:
        if df is None:
            df = pd.read_csv(bfile, sep="\t", names=blast_columns)
        else:
            df = df.append(pd.read_csv(bfile, sep="\t",
                                       names=blast_columns))
    df = df.drop_duplicates('probe')  # Drop rows with repeated probes
    return df


def annotate_seqdata(df, seqfiles):
    """Returns the passed dataframe, annotated with locus tags"""
    ids = []
    locus_tags = []
    for seqfile in seqfiles:
        for seq in SeqIO.parse(seqfile, 'fasta'):
            labels = seq.description.split(' ')
            for label in labels:
                if label.startswith('[locus_tag'):
                    ids.append(seq.id)
                    locus_tags.append(label.split('=')[-1][:-1])
    seqdf = pd.DataFrame({'match': ids, 'locus_tag': locus_tags})
    return pd.merge(df, seqdf, 'inner', ['match'])    


def index_column(df, colname):
    """Return the dataframe, with an index column for 'probe's"""
    col_ids = df[colname].unique()
    nvals = len(col_ids)
    col_lookup = dict(zip(col_ids, range(nvals)))
    df['{0}_index'.format(colname)] = df[colname].replace(col_lookup).values
    return df
    

def reduce_dataset(df, colname, n=2000, seed=True):
    """Returns the passed dataframe, with a reduced set of rows"""
    if seed:
        random.seed(SEED)  # for reproducibility of random choice

    col_ids = df[colname].unique()
    nvals = len(col_ids)

    indices = [random.randint(0, nvals) for i in range(n)]  
    reduced = df.loc[df['{0}_index'.format(colname)].isin(indices)]

    # create indices and values for probes
    new_ids = reduced[colname].unique()
    nvals = len(new_ids)
    new_lookup = dict(zip(new_ids, range(nvals)))

    # add data column with probe index from probe_lookup
    reduced['{0}_index'.format(colname)] =\
        reduced[colname].replace(new_lookup).values

    return reduced


def reduce_dataset_by_column_value(df, colname, values):
    """Returns the passed dataframe, with only the passed column values"""
    col_ids = df[colname].unique()
    nvals = len(col_ids)

    # Reduce dataset
    reduced = df.loc[df['locus_tag'].isin(values)]

    # create indices and values for probes
    new_ids = reduced[colname].unique()
    nvals = len(new_ids)
    new_lookup = dict(zip(new_ids, range(nvals)))

    # add data column with probe index from probe_lookup
    reduced['{0}_index'.format(colname)] =\
        reduced[colname].replace(new_lookup).values

    return reduced


def extract_fit_variable_summary(fit, varname, index=None):
    """Returns summary information for a variable in the passed Stan fit object

    Calculates mean, std, median, and 5%, 25%, 75% and 95% percentiles
    for the passed variable, returning them as a dataframe.
    """
    # Using Pandas methods
    mean = pd.Series(fit[varname].mean(0), index=index)
    se = pd.Series(fit[varname].std(0), index=index)

    # Need to use numpy functions
    median = pd.Series(np.median(fit[varname], 0), index=index)
    perc_2_5 = pd.Series(np.percentile(fit[varname], 2.5, 0), index=index)
    perc_25 = pd.Series(np.percentile(fit[varname], 25, 0), index=index)
    perc_75 = pd.Series(np.percentile(fit[varname], 75, 0), index=index)
    perc_97_5 = pd.Series(np.percentile(fit[varname], 97.5, 0), index=index)

    return pd.DataFrame({'%s_mean' % varname: mean,
                         '%s_sem' % varname: se,
                         '%s_median' % varname: median,
                         '%s_2.5pc' % varname: perc_2_5,
                         '%s_97.5pc' % varname: perc_97_5,
                         '%s_25pc' % varname: perc_25,
                         '%s_75pc' % varname: perc_75})



def extract_df_variable_summary(df, varname, index=None):
    """Returns summary information for a variable in the passed datframe object

    This function expects a dataframe of pickled fit information

    Calculates mean, std, median, and 5%, 25%, 75% and 95% percentiles
    for the passed variable, returning them as a dataframe.
    """
    # Using Pandas methods
    mean = pd.Series(df[varname][0].mean(0), index=index)
    se = pd.Series(df[varname][0].std(0), index=index)

    # Need to use numpy functions
    median = pd.Series(np.median(df[varname][0], 0), index=index)
    perc_2_5 = pd.Series(np.percentile(df[varname][0], 2.5, 0), index=index)
    perc_25 = pd.Series(np.percentile(df[varname][0], 25, 0), index=index)
    perc_75 = pd.Series(np.percentile(df[varname][0], 75, 0), index=index)
    perc_97_5 = pd.Series(np.percentile(df[varname][0], 97.5, 0), index=index)

    return pd.DataFrame({'%s_mean' % varname: mean,
                         '%s_sem' % varname: se,
                         '%s_median' % varname: median,
                         '%s_2.5pc' % varname: perc_2_5,
                         '%s_97.5pc' % varname: perc_97_5,
                         '%s_25pc' % varname: perc_25,
                         '%s_75pc' % varname: perc_75})


def extract_variable_summaries(obj, otype='fit',
                               varnames=['a', 'b', 'g', 'd'],
                               indices=None,
                               data=None):
    """Return dataframe of parameter estimate summaries

    For this modelling there is a specific issue with estimating variables on
    arrays (length 6), and estimating them on probes (length around 6000),
    and having to combine them.

    The calls to extract_*_variable_summary() return a dataframe for each
    variable. We broadcast values for a and g across the probe dataset, and
    join values for b and d directly.
    """
    # Choice of function depends on object being passed
    functions = {'fit': extract_fit_variable_summary,
                 'df': extract_df_variable_summary}

    # Get dataframes for each fitted variable summary, keyed by variable name
    dfdict = defaultdict()
    for varname, index in zip(varnames, indices):
        dfdict[varname] = functions[otype](obj, varname, index)
        dfdict[varname].reset_index(inplace=True)

    # Broadcast parameter estimates across probes
    df = pd.merge(data, dfdict['a'],
                  left_on='repXtrt', right_on='index')
    df = pd.merge(df, dfdict['b'],
                  left_on='locus_tag', right_on='index')
    df = pd.merge(df, dfdict['g'],
                  left_on='repXtrt', right_on='index')
    df = pd.merge(df, dfdict['d'],
                  left_on='locus_tag', right_on='index')
    
    # Broadcast parameter estimates across locus tags
    lt = pd.DataFrame(data['locus_tag'].unique())
    lt.columns = ['locus_tag']
    lt = pd.merge(lt, dfdict['b'],
                  left_on='locus_tag', right_on='index')
    lt = pd.merge(lt, dfdict['d'],
                  left_on='locus_tag', right_on='index')
    
    df.drop('index_x', 1, inplace=True)
    df.drop('index_y', 1, inplace=True)    
    lt.drop('index_x', 1, inplace=True)
    lt.drop('index_y', 1, inplace=True)    

    lt.sort('locus_tag', inplace=True)

    return df, lt


def boxplot_medians(estimates, varnames=['a', 'b', 'g', 'd']):
    """Plot 2x2 boxplot of parameter median estimates"""
    fig, axes = plt.subplots(int(len(varnames)/2), 2,
                             figsize=(12, 2 * len(varnames)))
    axes = axes.ravel()
    fig.subplots_adjust(hspace=0.3)

    for idx, varname in enumerate(varnames):
        sns.boxplot(estimates['{0}_median'.format(varname)],
                    ax=axes[idx])
        axes[idx].set_title("Median {0}".format(varname))


def split_estimates(df, org):
    """Split the passed dataframe into either Sakai or DH10B subsets"""
    if org == 'dh10b':
        subset = df.loc[df['locus_tag'].str.startswith('ECDH10B')]
    else:
        subset = df.loc[~df['locus_tag'].str.startswith('ECDH10B')]
    return subset


def plot_treatment_vs_control(df):
    """Plot median treatment vs control parameters"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 8))
    axes = axes.ravel()
    fig.subplots_adjust(hspace=0.3)
    
    for idx, xvar, yvar, ax in zip(range(2),
                                   ['a_median', 'a_median',
                                    'b_median', 'b_median'],
                                   ['g_median', 'd_median',
                                    'g_median', 'd_median'],
                                   axes):
        ax.scatter(df[xvar], df[yvar], alpha=0.2)
        ax.set_xlabel(xvar)
        ax.set_ylabel(yvar)


def label_positive_effects(df):
    """Label the locus tags as having positive effects on treatment, control,
    or both.
    """
    df['trt_pos'] = df['d_25pc'] > 0
    df['ctl_pos'] = df['b_25pc'] > np.percentile(df['b_median'], 97.5)
    df['combined'] = df['trt_pos'] & df['ctl_pos']
    return df


def plot_parameter(df, varname, title='', thresh=0):
    """Plot the estimated parameter median, and 50% CI, in locus tag order
    
    credibility intervals are coloured blue if they include the threshold,
    red (value below threshold) or green (value above threshold) otherwise
    """
    vals = df['{0}_median'.format(varname)]
    cilo = df['{0}_25pc'.format(varname)]
    cihi = df['{0}_75pc'.format(varname)]
    
    plt.figure(figsize=(20,8))
    #plt.scatter(range(len(test_data)), test_data['beta'], alpha=1, color='k')
    plt.scatter(range(len(df)), vals, c='k', marker='.')
    for idx, val, vlo, vhi in zip(range(len(df)),
                                  vals, cilo, cihi):
        if vlo < thresh < vhi:
            color = 'b-'
        elif val < thresh:
            color = 'r-'
        elif val > thresh:
            color = 'g-'
        else:
            color = 'k-'
        plt.plot([idx, idx], [vlo, vhi], color, alpha=0.4)
    plt.title("{0} [threshold: {1:.2f}]".format(title, thresh))
    plt.xlim(0, len(df));


def get_annotation(tag, anndict):
    try:
        return anndict[tag]
    except KeyError:
        return None


def annotate_locus_tags(df, gbfilepath):
    """Add gene product annotations from gbfiles to passed dataframe

    The annotations are added/placed in a column called "annotation", and are
    identified on the basis of the "locus_tag" column
    """
    products = dict()
    startpos = defaultdict(int)
    for record in SeqIO.parse(gbfilepath, 'genbank'):
        products.update({ft.qualifiers['locus_tag'][0]:ft.qualifiers['product'][0]
                         for ft in record.features if
                         (ft.type == 'CDS' and
                          'product' in ft.qualifiers)})
        startpos.update({ft.qualifiers['locus_tag'][0]:
                         int(ft.location.nofuzzy_start)
                         for ft in record.features if
                         ft.type == 'gene'})
    df['annotation'] = df['locus_tag'].apply(get_annotation,
                                             args=(products,))
    df['startpos'] = df['locus_tag'].apply(get_annotation,
                                           args=(startpos,))
    return df


def parse_full_fit(picklefilename, datafilename):
    """Parses the full model fit into a Pandas dataframe which is returned

    The returned dataframe has columns for mean, SEM, median, and 2.5, 25,
    75, 97.5 percentiles
    """
    # Load fit
    with open(picklefilename, 'rb') as ifh:
        fit = pickle.load(ifh)
    indata = pd.read_csv(datafilename, sep="\t")
    locus_tags = indata['locus_tag'].unique()

    # Get dataframes for each fitted variable summary, and join them
    dflist = []
    for varname in ['a', 'b', 'g', 'd']:
        dflist.append(extract_variable_summaries(fit, varname, locus_tags))

    return pd.concat(dflist, axis=1)


def plot_errors(df):
    """Plot distributions of absolute and relative error in crossvalidation"""
    fig, axes = plt.subplots(1, 2, figsize=(12,4))
    fig.subplots_adjust(hspace=.25)
    axes = axes.ravel()
    for ttl, col, ax in zip(("absolute error", "relative error"),
                            ("y_pred_abs_error", "y_pred_rel_error"),
                            axes):
        ax.set_title(ttl)
        sns.boxplot(df[col], ax=ax)    

def plot_error_vs_column(df, colname):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes = axes.ravel()
    for ttl, col, ax in zip(("absolute error", "relative error"),
                             ("y_pred_abs_error", "y_pred_rel_error"),
                             axes):
        ax.set_title("{0} v {1}".format(ttl, colname))
        ax.set_xlabel(colname)
        ax.set_ylabel(ttl)
        ax.scatter(df[colname], df[col], alpha=0.05)
