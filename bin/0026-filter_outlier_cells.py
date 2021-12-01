#!/usr/bin/env python

__date__ = '2021-01-20'
__version__ = '0.0.1'

import argparse
import random
import os
from distutils.version import LooseVersion
import numpy as np
import scanpy as sc
import pandas as pd
from sklearn.svm import OneClassSVM
from sklearn.covariance import EllipticEnvelope
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
import matplotlib
matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# matplotlib.style.use('ggplot')
import seaborn as sns

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Performs automatic outlier detection over an anndata dataset,
            plotting the outlier cells.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '--outliers_fraction',
        action='store',
        dest='outliers_fraction',
        default=0.0,
        type=float,
        help='Anticipated fraction of outlier cells. If 0.0, then runs sklearn \
            methods with "auto" as the anticipated number of outlier cells.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '--metadata_columns',
        action='store',
        dest='metadata_columns',
        default='pct_counts_gene_group__mito_transcript,log1p_total_counts,log1p_n_genes_by_counts',
        help='Columns to use for outliers.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '--cell_qc_column',
        action='store',
        dest='cell_qc_column',
        default='cell_passes_qc',
        help='If column exists, cells are first filtered for those that \
            evaluate to true in this column, then this column is updated \
            based on the QC filtering. If this column does not exist then \
            it is added. \
            (default: %(default)s)'
    )

    parser.add_argument(
        '--method',
        action='store',
        dest='method',
        default='IsolationForest',
        help='Method for outlier detection. \
            (default: %(default)s)'
    )

    parser.add_argument(
        '--max_samples',
        action='store',
        dest='max_samples',
        default=0.1,
        type=float,
        help='The fraction of cells to draw from X to train each estimator. \
            Only valid if method == IsolationForest. \
            (default: %(default)s)'
    )
    # TODO: add option where user can say "run per each experiment id"

    parser.add_argument(
        '--cell_filtered_per_experiment_file',
        action='store',
        dest='cell_filtered_per_experiment',
        default='None',
        help='File detailing samples filtered per experiment. \
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='outliers',
        help='Basename of output files. \
            (default: %(default)s)'
    )

    parser.add_argument(
        '--anndata_compression_opts',
        action='store',
        dest='anndata_compression_opts',
        default=4,
        type=int,
        help='Compression level in anndata. A larger value decreases disk \
            space requirements at the cost of compression time. \
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=options.h5)
    adata.obs['cell_id'] = adata.obs.index
    adata_original = adata.copy()

    # Drop out previous QCed cells
    cell_qc_column = options.cell_qc_column
    if cell_qc_column in adata.obs.columns:
        n_cells_original = adata.shape[0]
        adata = adata[adata.obs[cell_qc_column], :]
        print('Filtered out {} previously flagged cells using {}'.format(
            n_cells_original - adata.shape[0],
            cell_qc_column
        ))
    else:
        adata.obs[cell_qc_column] = True

    # Get ballpark number of outliers
    outliers_fraction = options.outliers_fraction
    n_cells = adata.shape[0]
    n_outliers = outliers_fraction * n_cells
    if n_outliers == 0:
        outliers_fraction = 'auto'

    # Get a list of the data to use to calculate outliers
    metadata_columns = options.metadata_columns.split(',')
    # metadata_columns = [
    #     'pct_counts_gene_group__mito_transcript',
    #     #'total_counts',
    #     'log1p_total_counts',
    #     #'n_genes_by_counts',
    #     'log1p_n_genes_by_counts'
    # ]

    method = options.method
    if method == 'LocalOutlierFactor':
        # fit the model for outlier detection (default)
        clf = LocalOutlierFactor(
            #n_neighbors=100,
            contamination=outliers_fraction
        )
    elif method == 'IsolationForest':
        max_samples = options.max_samples
        # if max_samples == 0.0:
        #     if n_cells < 1000:
        #         max_samples = 250
        #     else:
        #         max_samples = 0.1
        print("Using max_samples of:\t{}".format(max_samples))
        clf = IsolationForest(
            #n_estimators=500,
            max_samples=max_samples,
            warm_start=False,
            contamination=outliers_fraction,
            random_state=0,
            bootstrap=True
        )
    elif method == 'EllipticEnvelope':
        if outliers_fraction == 'auto':
            outliers_fraction = 0.1
        clf = EllipticEnvelope(
            contamination=outliers_fraction
        )
    elif method == 'OneClassSVM':
        clf = OneClassSVM(
            # nu=n_outliers,
            # kernel="rbf",
            # gamma=0.1
        )
    else:
        raise ValueError('ERROR: invalid method.')

    # Fit the model
    if method != 'LocalOutlierFactor':
        clf = clf.fit(adata.obs[metadata_columns].values)

    if method == 'LocalOutlierFactor':
        # use fit_predict to compute the predicted labels of the training
        # samples (when LOF is used for outlier detection, the estimator has
        # no predict, decision_function and score_samples methods).
        adata.obs[cell_qc_column] = clf.fit_predict(
            adata.obs[metadata_columns].values
        ) == 1
    else:
        adata.obs[cell_qc_column] = clf.predict(
            adata.obs[metadata_columns].values
        ) == 1
    adata.uns['cell_outlier_estimator'] = method

    # Update the original data to flag those cells that passed the outlier
    adata_original.obs[cell_qc_column] = False
    adata_original.obs.loc[
        adata.obs['cell_id'][adata.obs[cell_qc_column]],
        cell_qc_column
    ] = True
    adata_original.obs[['cell_id', cell_qc_column]].to_csv(
        '{}-outliers_filtered.tsv.gz'.format(options.of),
        sep='\t',
        compression=compression_opts,
        index=False,
        header=True
    )
    del adata_original.obs['cell_id']

    # Calculate cell_filtered_per_experiment
    filter_columns = [
        'experiment_id',
        'filter_type',
        'n_cells_left_in_adata'
    ]
    if options.cell_filtered_per_experiment == 'None':
        df_cell_filt_per_exp = adata.obs['experiment_id'].value_counts()
        df_cell_filt_per_exp = df_cell_filt_per_exp.rename_axis(
            'experiment_id'
        ).reset_index(name='n_cells_left_in_adata')
        df_cell_filt_per_exp['filter_type'] = 'before_filters'
        df_cell_filt_per_exp = df_cell_filt_per_exp[filter_columns]
    else:
        # Load the samples filtered per experiment file:
        df_cell_filt_per_exp = pd.read_csv(
            options.cell_filtered_per_experiment,
            sep="\t"
        )
        filt = df_cell_filt_per_exp['filter_type'] != 'after_filters'
        df_cell_filt_per_exp = df_cell_filt_per_exp.loc[filt, :]
    # Now calculate the n cells left after all filters
    adata_after_filters = adata.obs.loc[adata.obs[cell_qc_column], :]
    df_cells_filtered = adata_after_filters['experiment_id'].value_counts()
    df_cells_filtered = df_cells_filtered.rename_axis(
        'experiment_id'
    ).reset_index(name='n_cells_left_in_adata')
    df_cells_filtered['filter_type'] = '{} {} outlier_{}'.format(
        'filter__all_samples',
        'after_outlier_filter',
        method
    )
    df_cells_filtered = df_cells_filtered[filter_columns]
    df_cell_filt_per_exp = df_cell_filt_per_exp.append(
        df_cells_filtered,
        ignore_index=True
    )
    df_cells_filtered['filter_type'] = 'after_filters'
    df_cell_filt_per_exp = df_cell_filt_per_exp.append(
        df_cells_filtered,
        ignore_index=True
    )
    # Save the final dataframe
    df_cell_filt_per_exp.to_csv(
        '{}-cell_filtered_per_experiment.tsv.gz'.format(options.of),
        sep='\t',
        compression=compression_opts,
        index=False,
        header=True
    )

    # Plot the identified outliers
    metadata_columns_original = metadata_columns.copy()
    metadata_columns.append(cell_qc_column)
    # print("Making plot")

    sns_plot = sns.PairGrid(
        adata.obs[metadata_columns],
        hue=cell_qc_column,
        height=2.5,
        diag_sharey=False
    )
    sns_plot.map_upper(
        sns.scatterplot,
        marker='+',
        alpha=0.05,
        s=8,
        edgecolor=None
    )
    sns_plot.add_legend()
    for lh in sns_plot._legend.legendHandles:
        lh.set_alpha(1)
        lh._sizes = [50]
    sns_plot.map_diag(sns.kdeplot)
    sns_plot.map_lower(
        sns.kdeplot,
        levels=3
    )
    sns_plot.savefig('{}-outlier_cells.png'.format(options.of))

    # Plot the cell density
    # print("Making density plot")
    sns_plot = sns.PairGrid(
        adata.obs[metadata_columns_original],
        height=2.5,
        diag_sharey=False
    )
    sns_plot.map_upper(
        sns.kdeplot,
        cmap="viridis",
        shade=True
    )
    sns_plot.add_legend()
    for lh in sns_plot._legend.legendHandles:
        lh.set_alpha(1)
        lh._sizes = [50]
    sns_plot.map_diag(sns.kdeplot)
    sns_plot.map_lower(
        sns.kdeplot,
        cmap="viridis",
        shade=True
    )
    sns_plot.savefig('{}-cell_desity.png'.format(options.of))

    # Save the updated adata matrix
    # print("Saving data")
    adata_original.write(
        '{}.h5ad'.format(options.of),
        compression='gzip',
        compression_opts=options.anndata_compression_opts
    )


if __name__ == '__main__':
    main()
