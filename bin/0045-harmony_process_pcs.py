#!/usr/bin/env python

__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
from distutils.version import LooseVersion
import os
import random
import numpy as np
import pandas as pd
import csv
import harmonypy as hm

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
            Runs Harmony on PCs.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-pc', '--pca_file',
        action='store',
        dest='pc',
        required=True,
        help='Tab-delimited file of PCs for each cell. First column is\
            cell_barcode. Subsequent columns are PCs.'
    )

    parser.add_argument(
        '-mf', '--metadata_file',
        action='store',
        dest='mf',
        required=True,
        help='Tab-delimited metadata file, must have a column labelled\
            cell_barcode that maps to pca_file.'
    )

    parser.add_argument(
        '-mc', '--metadata_columns',
        action='store',
        dest='mc',
        required=True,
        help='Comma separated string of columns to use in metadata_file.'
    )

    parser.add_argument(
        '-t', '--theta',
        action='store',
        dest='theta',
        default='',
        help='Comma separated string of theta values (corresponding to\
            metadata_columns). If "" then sets theta to 2 for all\
            columns. Larger values of theta result in more diverse\
            clusters.\
            (default: "")'
    )

    parser.add_argument(
        '-npc', '--n_pcs',
        action='store',
        dest='npc',
        default=0,
        type=int,
        help='Number of PCs to use.\
            (default: maximum number in tsv_pcs file)'
    )

    parser.add_argument(
        '-of', '--out_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: <tsv_pcs>-harmony)'
    )

    options = parser.parse_args()

    # Fixed settings.
    verbose = True

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}-harmony'.format(
            os.path.basename(options.pc.rstrip('tsv.gz').rstrip('.'))
        )

    # Load the PCs.
    df_pca = pd.read_csv(options.pc, sep='\t', index_col='cell_barcode')

    # Check that nPCs is valid.
    n_pcs = options.npc
    if n_pcs == 0:
        n_pcs = len(df_pca.columns)
    elif n_pcs > len(df_pca.columns):
        raise Exception(
            '--number_pcs ({}) is > than n_pcs in --tsv_pcs ({}).'.format(
                n_pcs,
                len(df_pca.columns)
            )
        )
    if verbose:
        print('Using {} PCs.'.format(n_pcs))

    # Subset down to these PCs.
    df_pca = df_pca.iloc[:, range(0, n_pcs)]

    # Get the metadata_file columns that we want to adjust with Harmony.
    metadata_columns = options.mc.split(',')

    # Read in the metadata file.
    df_meta = pd.read_csv(options.mf, sep='\t', index_col='cell_barcode')
    # Ensure cell order in df_meta is the same as df_pca
    df_meta = df_meta.loc[df_pca.index, metadata_columns]
    # Also ensure that the metadata columns are categorical -- run_harmony
    # fails if not categorical
    try:
        df_meta[metadata_columns].describe().loc['unique']
    except KeyError:
        print("metadata_columns contains non-categorical attributes. Harmony does \
        not work with continuous variables. Either make attributes a string or \
        use a different column.")

    # Get the theta values for each column (if none, set to 2 for all columns).
    theta = [2] * len(metadata_columns)
    if options.theta != '':
        theta = [float(i) for i in options.theta.split(',')]

    # Run Harmony
    harmony_embeddings = hm.run_harmony(
        data_mat=df_pca.values,  # Pandas dataframe to numpy.ndarray
        meta_data=df_meta,
        vars_use=metadata_columns,
        theta=theta,
        max_iter_kmeans=500,
        verbose=verbose
    )
    # NOTE: harmony_embeddings.result() == harmony_embeddings.Z_corr
    df_harmony = pd.DataFrame(np.transpose(harmony_embeddings.Z_corr))
    harmony_cols = [
        'harmony{}'.format(i + 1) for i in range(df_harmony.shape[1])
    ]
    df_harmony.columns = harmony_cols
    df_harmony['cell_barcode'] = df_pca.index
    final_col_order = ['cell_barcode']
    final_col_order.extend(harmony_cols)
    df_harmony = df_harmony.loc[:, final_col_order]

    # Save the clustered data to a data frame.
    df_harmony.to_csv(
        '{}.tsv.gz'.format(out_file_base),
        sep='\t',
        index=False,
        quoting=csv.QUOTE_NONNUMERIC,
        # index_label='cell_barcode',
        na_rep='',
        compression=compression_opts
    )


if __name__ == '__main__':
    main()
