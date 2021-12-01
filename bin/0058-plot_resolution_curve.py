#!/usr/bin/env python

__date__ = '2020-05-01'
__version__ = '0.0.1'

import argparse
from distutils.version import LooseVersion
import os
import random
import numpy as np
import pandas as pd
import csv
import math

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec
from matplotlib import colors
from matplotlib import cm

from sklearn import metrics

plt.rcParams.update({'font.size': 18})

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)


# This function is based off of scanpy:
# https://github.com/theislab/scanpy/blob/master/scanpy/plotting/_tools/scatterplots.py
def panel_grid(hspace, wspace, ncols, num_panels):
    """Init plot."""
    n_panels_x = min(ncols, num_panels)
    n_panels_y = np.ceil(num_panels / n_panels_x).astype(int)
    if wspace is None:
        #  try to set a wspace that is not too large or too small given the
        #  current figure size
        wspace = 0.75 / rcParams['figure.figsize'][0] + 0.02
    # each panel will have the size of rcParams['figure.figsize']
    fig = plt.figure(
        figsize=(
            n_panels_x * rcParams['figure.figsize'][0] * (1 + wspace),
            n_panels_y * rcParams['figure.figsize'][1],
        )
    )
    left = 0.2 / n_panels_x
    bottom = 0.13 / n_panels_y
    gs = gridspec.GridSpec(
        nrows=n_panels_y,
        ncols=n_panels_x,
        left=left,
        right=1 - (n_panels_x - 1) * left - 0.01 / n_panels_x,
        bottom=bottom,
        top=1 - (n_panels_y - 1) * bottom - 0.1 / n_panels_y,
        hspace=hspace,
        wspace=wspace
    )
    return fig, gs


def _create_colors(classes):
    n_cts = len(classes)
    color_norm = colors.Normalize(vmin=-n_cts / 3, vmax=n_cts)
    ct_arr = np.arange(n_cts)
    ct_colors = cm.YlGnBu(color_norm(ct_arr))

    return ct_colors


def plot_roc(ax, y_prob, y_test, classes):
    """Plot ROC curve. Based off of NaiveDE library."""
    ct_colors = _create_colors(classes)

    for i, cell_type in enumerate(classes):
        fpr, tpr, _ = metrics.roc_curve(y_test == cell_type, y_prob[:, i])
        ax.plot(fpr, tpr, c=ct_colors[i], lw=2)

    ax.plot([0, 1], [0, 1], color='k', ls=':')


def plot_prec_recall(ax, y_prob, y_test, classes):
    """Plot ROC curve. Based off of NaiveDE library."""
    ct_colors = _create_colors(classes)

    for i, cell_type in enumerate(classes):
        precision, recall, _ = metrics.precision_recall_curve(
            y_test == cell_type, y_prob[:, i]
        )
        ax.plot(recall, precision, c=ct_colors[i], lw=2)

    # TODO add balanced result
    # ax.plot([0, 1], [0, 1], color='k', ls=':')


def _make_plots(
    df,
    out_file_base
):
    # Get a sorted list of our resolutions
    resolutions = np.sort(df['cluster__resolution'].unique())
    # Make sure resolution is set to a categorical variable.
    df['resolution'] = df['cluster__resolution'].astype(
        'category'
    )

    # Get a sorted list of the number of neighbors used in clustering graph
    n_neighbors = np.sort(df['neighbors__n_neighbors'].unique())
    # Make sure resolution is set to a categorical variable.
    df['n_neighbors'] = df['neighbors__n_neighbors'].astype(
        'category'
    )

    # get all pairwise combinations of resolutions and neighbors
    # iters = [(x, y) for x in resolutions for y in n_neighbors]
    iters = [(x, y) for x in n_neighbors for y in resolutions]

    hspace = 0.5
    # Make a grid of ROC plots for each resolution.
    ncols = 5
    nrows = math.ceil(len(iters) / ncols)
    fig, grid = panel_grid(
        hspace=hspace,
        wspace=None,
        ncols=ncols,
        num_panels=len(iters)
    )
    # fig, axs = plt.subplots(
    #     len(df['resolution'].unique()),
    #     sharex=True,
    #     sharey=True
    # )
    i__ax = 0
    i__row = 0
    # for i__res, i__n_neighbors in iters:
    for i__n_neighbors, i__res in iters:
        # Get the subset of data
        df_tmp = df[df['resolution'] == i__res]
        df_tmp = df_tmp[df_tmp['n_neighbors'] == i__n_neighbors]
        # Drop extra cell type columns for clusters not found at this
        # resolution.
        df_tmp = df_tmp.dropna(axis='columns')
        # Get classification / cell type columns
        df_class_cols = [i for i in df_tmp.columns if 'class__' in i]
        df_nonclass_cols = [
            i for i in df_tmp.columns if i not in df_class_cols
        ]

        # For dev if multiple files
        # df_tmp = df_tmp[df_tmp['file'] == df_tmp['file'].values[-1]]
        # Check to make sure we have the right data ... for every non-class
        # column, there should be only one value. Otherwise, we are not
        # plotting a unique run.
        for index, value in df_tmp[df_nonclass_cols].nunique().items():
            if 'cell_label' not in index and value > 1:
                print(index, value)
                raise Exception('ERROR: multiple resolution files combined.')

        # Set the facet title.
        plt_title = '{} resolution'.format(i__res)
        if len(n_neighbors) > 1:
            plt_title = '{}\n{} n_neighbors'.format(
                plt_title,
                i__n_neighbors
            )
        plt_title = '{}\n{} clusters'.format(plt_title, len(df_class_cols))
        plt_title = plt_title.rstrip()

        # Get the proper axis for this plot.
        ax = fig.add_subplot(grid[i__ax])
        # Add ROC plot
        plot_roc(
            # axs[i__ax],  # If use subplots
            ax,
            df_tmp[df_class_cols].values,
            df_tmp['cell_label_true'].values,
            df_class_cols
        )
        ax.set_title(plt_title)

        if i__ax % ncols == 0:
            i__row += 1
            plt.ylabel('TPR')
        if i__row == nrows:
            plt.xlabel('FPR')

        i__ax += 1

    fig.savefig(
        '{}-resolution_roc.png'.format(out_file_base),
        # dpi=250,
        bbox_inches='tight'
    )
    plt.close(fig)

    # Make a grid of precision-recall plots for each resolution.
    fig, grid = panel_grid(
        hspace=hspace,
        wspace=None,
        ncols=ncols,
        num_panels=len(iters)
    )
    i__ax = 0
    i__row = 0
    # for i__res, i__n_neighbors in iters:
    for i__n_neighbors, i__res in iters:
        # Get the subset of data
        df_tmp = df[df['resolution'] == i__res]
        df_tmp = df_tmp[df_tmp['n_neighbors'] == i__n_neighbors]
        # Drop extra cell type columns for clusters not found at this
        # resolution.
        df_tmp = df_tmp.dropna(axis='columns')
        # Get classification / cell type columns
        df_class_cols = [i for i in df_tmp.columns if 'class__' in i]
        df_nonclass_cols = [
            i for i in df_tmp.columns if i not in df_class_cols
        ]

        # Check to make sure we have the right data ... for every non-class
        # column, there should be only one value. Otherwise, we are not
        # plotting a unique run.
        for index, value in df_tmp[df_nonclass_cols].nunique().items():
            if 'cell_label' not in index and value > 1:
                print(index, value)
                raise Exception('ERROR: multiple resolution files combined.')

        # Set the facet title.
        plt_title = '{} resolution'.format(i__res)
        if len(n_neighbors) > 1:
            plt_title = '{}\n{} n_neighbors'.format(
                plt_title,
                i__n_neighbors
            )
        plt_title = '{}\n{} clusters'.format(plt_title, len(df_class_cols))
        plt_title = plt_title.rstrip()

        # Get the proper axis for this plot.
        ax = fig.add_subplot(grid[i__ax])
        # Add plot
        plot_prec_recall(
            ax,
            df_tmp[df_class_cols].values,
            df_tmp['cell_label_true'].values,
            df_class_cols
        )
        ax.set_title(plt_title)

        if i__ax % ncols == 0:
            i__row += 1
            plt.ylabel('Precision')
        if i__row == nrows:
            plt.xlabel('Recall')

        i__ax += 1

    fig.savefig(
        '{}-resolution_precision_recall.png'.format(out_file_base),
        # dpi=300,
        bbox_inches='tight'
    )
    plt.close(fig)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Plots AUC from model_report.tsv.gz
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '--y_prob_dfs',
        action='store',
        dest='model_reports',
        required=True,
        help='List of tab-delimited files of y_prob_dfs.\
            List should be split by "::" (e.g. file1.tsv.gz::file2.tsv.gz).'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='resolution_auc',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # Get the out file base.
    out_file_base = options.of

    # Read in a list of model_report.tsv.gz files.
    # NOTE: cluster id will not necissarily be comparable across reports.
    files = options.model_reports.split('::')

    list_modelreport = []
    for i in range(len(files)):
        _tmp = pd.read_csv(
            files[i],
            sep='\t'
        )
        _tmp['file'] = files[i]
        list_modelreport.append(_tmp)

    # Make one long dataframe.
    df = pd.concat(list_modelreport)

    # Save the results
    df.to_csv(
        '{}-merged_test_result.tsv.gz'.format(out_file_base),
        sep='\t',
        index=False,
        quoting=csv.QUOTE_NONNUMERIC,
        na_rep='',
        compression=compression_opts
    )

    # Make a combined plot for all n_neighbors values
    _make_plots(
        df,
        out_file_base
    )

    # Make separate plots for each n_neighbors values
    n_neighbors = df['neighbors__n_neighbors'].unique()
    for i__n_neighbors in n_neighbors:
        df_tmp = df.loc[df['n_neighbors'] == i__n_neighbors, :].copy()
        _make_plots(
            df_tmp,
            '{}-n_neighbors={}'.format(out_file_base, i__n_neighbors)
        )


if __name__ == '__main__':
    main()
