#!/usr/bin/env python

__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
from distutils.version import LooseVersion
import os
import random
import numpy as np
import pandas as pd
import scanpy as sc
# import csv

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
            Read anndata object. Runs BBKNN. Assumes PCS have already been
            calculated.
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
        '-bk', '--batch_key',
        action='store',
        dest='batch_key',
        default='experiment_id',
        help='Batch key for BBKNN.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-npc', '--n_pcs',
        action='store',
        dest='n_pcs',
        default=30,
        type=int,
        help='Number of PCs to use.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-ncpu', '--number_cpu',
        action='store',
        dest='ncpu',
        default=1,
        type=int,
        help='Number of CPUs to use.\
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

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='adata-normalize_pca-bbknn',
        help='Directory and basename of output files.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    # sc.set_figure_params(dpi_save = 300)

    # Set some parameters.
    n_pcs = options.n_pcs
    output_file = options.of
    batch_key = options.batch_key

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)

    # Run bbknn
    # Total number of neighbours = neighbors_within_batch x the number of
    # batches.
    sc.external.pp.bbknn(
        adata=adata,
        batch_key=batch_key,
        copy=False,
        # neighbors_within_batch=5
        n_pcs=n_pcs
    )
    adata.uns['neighbors']['params']['bbknn'] = True
    adata.uns['neighbors']['params']['n_pcs'] = n_pcs

    # Save the resulting adata.
    adata.write(
        '{}.h5ad'.format(
            # os.path.basename(options.h5.rstrip('h5ad').rstrip('.'))
            output_file
        ),
        compression='gzip'
        #compression_opts=options.anndata_compression_opts
    )

    # UMAP
    # Saved to adata.uns['umap'] and adata.obsm['X_umap']
    sc.tl.umap(
        adata,
        min_dist=0.05,  # Scanpy default = 0.05
        spread=1.0,  # Scanpy default = 1.0
        init_pos='spectral',  # Scanpy default = spectral
        n_components=n_pcs,
        # For some reason cannot access neighbors key slot, thus we
        # must keep uns['neighbors'] until we have run this.
        # neighbors_key='neighbors__{}'.format(plt__label),
        copy=False,
        random_state=0
    )
    # Add umap info to params stash
    # adata.uns['umap']['params']['tsv_reduced_dims'] = options.pc
    adata.uns['umap']['params']['n_reduced_dims_input'] = n_pcs
    adata.uns['umap']['params']['umap_init'] = 'spectral'
    # adata.uns['umap']['params']['n_neighbors'] = i__n_neighbors
    adata.uns['umap']['params']['umap_min_dist'] = 0.05
    adata.uns['umap']['params']['umap_spread'] = 1.0
    adata.uns['umap']['params']['density'] = ''

    # Save the resulting umap dimensions.
    df_reduced_dims = pd.DataFrame(
        adata.obsm['X_umap'],
        index=adata.obs_names,
        columns=[
            'UMAP{}'.format(x) for x in range(
                1, adata.obsm['X_umap'].shape[1]+1
            )
        ]
    )
    df_reduced_dims.to_csv(
        '{}-reduced_dims.tsv.gz'.format(output_file),
        sep='\t',
        index=True,
        index_label='cell_barcode',
        na_rep='',
        compression=compression_opts
    )


if __name__ == '__main__':
    main()
