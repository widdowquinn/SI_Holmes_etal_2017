from matplotlib import pyplot as plt

import numpy as np
import pandas as pd
import pickle
import scipy
import seaborn as sns

from Bio import SeqIO


def plot_fit_params(fit, params):
    nparams = len(params)
    fig, axes = plt.subplots(nparams, 2,
                             figsize=(10, nparams * 2.5))
    axes = axes.ravel()
    fig.subplots_adjust(hspace=0.3)
    
    for idx, vbl, title in zip(range(nparams), params,
                               params):
        # density plot
        sns.kdeplot(fit[vbl], ax=axes[idx * 2])
        axes[idx * 2].set_xlim(fit[vbl].min(),
                              fit[vbl].max())
        # scatterplot
        axes[idx * 2 + 1].plot(fit[vbl], 'o', alpha=0.3)
    
        # labels
        axes[idx * 2].set_title(title)
        axes[idx * 2 + 1].set_title(title)
    axes[nparams-1 * 2].set_xlabel("value")
    axes[nparams-1 * 2 + 1].set_xlabel("iteration")


def print_intervals(fit, param, percentile):
    assert 0 <= percentile <= 100, "We need percentile in [0, 100]"
    lq, uq = 50 - percentile/2., 50 + percentile/2.
    qs = tuple(np.percentile(fit[param], [lq, uq]))
    return '{0} {1}%CI: {2:.3f}..{3:.3f}'.format(param, percentile, *qs)


def plot_threshold_errors(means, errors, thresh, upper=True, param=''):
    # apply threshold and order values
    if upper:
        means_t = means.loc[means >= thresh]
    else:
        means_t = means.loc[means <= thresh]
    errors_t = errors.loc[errors.index.isin(means_t.index)]
    order = means_t.sort_values().index
    # plot ordered means and errors
    plt.scatter(range(len(means_t)), means_t[order])
    for idx, mn, se in zip(range(len(means_t)),
                           means_t[order], errors_t[order]):
        if se > abs(mn):
            color = 'r-'
        else:
            color = 'b-'
        plt.plot([idx, idx], [mn - se, mn + se], color)
    plt.xlabel("ordered probe ID")
    plt.ylabel("mean {0}".format(param))
    plt.xlim(-1, len(means_t))
    if upper:
        modetxt = 'upper'
    else:
        modetxt = 'lower'
    plt.title("Variation in mean estimates ({0})".format(modetxt))


def plot_input_output_violin(data):
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


def corrfunc(x, y, **kws):
    r, _ = scipy.stats.pearsonr(x, y)
    ax = plt.gca()
    ax.annotate("r = {:.3f}".format(r),
                xy=(.3, .5), size=15,
                xycoords=ax.transAxes)
    

def plot_correlation(data, title=None):
    g = sns.PairGrid(data)
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


def wide_to_long(df, stage):
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
    if treatment:
        treatval = 1
    else:
        treatval = 0            
    df = pd.merge(wide_to_long(df_in, 0), wide_to_long(df_out, 1),
                  left_index=True, right_index=True)
    df['treatment'] = treatval
    df.reset_index(inplace=True)
    return df


def save_model(model, filename):
    with open(filename, 'wb') as f:
        pickle.dump(model, f)


def load_model(filename):
    return pickle.load(open(filename, 'rb'))


blast_columns = ['probe', 'match', 'identity', 'length', 'mismatch', 'gapopen',
                 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']


def subset_blast(df, blastfile, seqfile, prefix):
    hits = pd.read_csv(blastfile, sep="\t", names=blast_columns)
    df['{0}_match'.format(prefix)] = \
        df['probe'].isin(hits['probe'].values).astype(int)
    blastdf = pd.merge(df, hits[['probe', 'match']], 'inner', ['probe'])
    ids = []
    locus_tags = []
    for seq in SeqIO.parse(seqfile, 'fasta'):
        labels = seq.description.split(' ')
        for label in labels:
            if label.startswith('[locus_tag'):
                ids.append(seq.id)
                locus_tags.append(label.split('=')[-1][:-1])
    seqdf = pd.DataFrame({'match': ids, 'locus_tag': locus_tags})
    return pd.merge(blastdf, seqdf, 'inner', ['match'])

