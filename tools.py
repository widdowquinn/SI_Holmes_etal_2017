from matplotlib import pyplot as plt

import numpy as np
import seaborn as sns

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
