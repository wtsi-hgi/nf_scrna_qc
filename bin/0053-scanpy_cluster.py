#!/usr/bin/env python

__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
import random
from distutils.version import LooseVersion
import numpy as np
import pandas as pd
import scanpy as sc
import csv
import warnings
import plotnine as plt9

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
            Read AnnData object and PCs file. Clusters the data. Saves an
            AnnData object with clusters in the clusters slot, a clusters
            file, and QC plots.
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
        '-pc', '--tsv_pcs',
        action='store',
        dest='pc',
        default='',
        help='Tab-delimited file of PCs for each cell. First column is\
            cell_barcode. Subsequent columns are PCs. If "", uses pca\
            slot in AnnData.\
            (default: "")'
    )

    parser.add_argument(
        '-cm', '--cluster_method',
        action='store',
        dest='cm',
        default='leiden',
        help='Clustering method. Valid options: [leiden|louvain].\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-npc', '--number_pcs',
        action='store',
        dest='npc',
        default=0,
        type=int,
        help='Number of PCs to use.\
            (default: maximum number in tsv_pcs file)'
    )

    parser.add_argument(
        '-r', '--resolution',
        action='store',
        dest='r',
        default=1.0,
        type=float,
        help='Resolution.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-nn', '--number_neighbors',
        action='store',
        dest='number_neighbors',
        default=25,
        type=int,
        help='Number of neighbors. If <= 0, sets to the number of unique\
            "experiment_id".\
            (default: %(default)s)'
    )

    parser.add_argument(
        '--force_recalculate_neighbors',
        action='store_true',
        dest='calculate_neighbors',
        default=False,
        help='Calculate neighbor graph even if it already exists in the\
            AnnData (which it my do so if you already ran BBKNN).\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-ncpu', '--number_cpu',
        action='store',
        dest='ncpu',
        default=4,
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
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: <h5_anndata>-<tsv_pcs>-clustered)'
    )

    options = parser.parse_args()

    # Fixed settings.
    verbose = True

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    # sc.set_figure_params(dpi_save = 300)

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}-{}-clustered'.format(
            os.path.basename(options.h5.rstrip('h5ad').rstrip('.')),
            os.path.basename(options.pc.rstrip('tsv.gz').rstrip('.'))
        )

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=options.h5)

    # Load the PCs.
    if options.pc == '':
        df_pca = pd.DataFrame(
            data=adata.obsm['X_pca'],
            index=adata.obs.index,
            columns=[
                'PC{}'.format(x) for x in
                range(1, adata.obsm['X_pca'].shape[1]+1)
            ]
        )
    else:
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

    # Add the reduced dimensions to the AnnData object.
    adata.obsm['X_pca'] = df_pca.loc[adata.obs.index, :].values.copy()

    # Check if BBKNN
    # Calculate neighbors for on the specified PCs.
    # By default saved to adata.uns['neighbors']
    #
    # First, however, check to see if adata.uns['neighbors'] already exists
    # ...and unless the user tells us not to, use that slot, not calculating
    # neighbors. This default behaviour is to accommodate the instance when
    # bbknn has been run on the data.
    if 'neighbors' not in adata.uns or options.calculate_neighbors:
        number_neighbors = options.number_neighbors
        if number_neighbors <= 0:
            number_neighbors = len(adata.obs['experiment_id'].cat.categories)
        sc.pp.neighbors(
            adata,
            use_rep='X_pca',
            n_neighbors=options.number_neighbors,
            n_pcs=n_pcs,
            copy=False,
            random_state=0
        )
    else:
        warnings.warn(
            'WARNING: found neighbors slot in adata.uns. {}'.format(
                'Not calculating neighbors (ignoring n_neighbors parameter).'
            )
        )
        # If we are using the pre-calculated neighbors drop npcs note.
        if 'n_pcs' in adata.uns['neighbors']['params']:
            n_pcs = adata.uns['neighbors']['params']['n_pcs']

    # Run the clustering, choosing either leiden or louvain algorithms
    cluster_method = options.cm
    cluster_resolution = options.r
    if cluster_method == 'leiden':
        sc.tl.leiden(
            adata,
            resolution=cluster_resolution,
            key_added=cluster_method,
            copy=False,
            random_state=0
        )
    elif cluster_method == 'louvain':
        sc.tl.louvain(
            adata,
            flavor='vtraag',
            resolution=cluster_resolution,
            key_added=cluster_method,
            copy=False,
            random_state=0
        )
    else:
        raise Exception(
            'Invalid --cluster_method: {}.'.format(
                cluster_method
            )
        )
    # Also save the clusters to the same spot so we know where they will be.
    adata.uns['cluster'] = adata.uns[cluster_method]
    adata.uns['cluster']['params']['method'] = cluster_method
    adata.obs['cluster'] = adata.obs[cluster_method]

    # Print the final number of clustered discrovered
    if verbose:
        print('{} clusters identified'.format(
            adata.obs[cluster_method].drop_duplicates().shape[0]
        ))

    # Save the clustered data to a data frame.
    cell_clustering_df = adata.obs[[cluster_method]].copy()
    cell_clustering_df.columns = ['cluster']
    cell_clustering_df['cluster_method'] = cluster_method
    cell_clustering_df['cluster_resolution'] = cluster_resolution
    cell_clustering_df.to_csv(
        '{}.tsv.gz'.format(out_file_base),
        sep='\t',
        index=True,
        quoting=csv.QUOTE_NONNUMERIC,
        index_label='cell_barcode',
        na_rep='',
        compression=compression_opts
    )

    adata.write(
        '{}.h5ad'.format(out_file_base),
        compression='gzip'
        #compression_opts=options.anndata_compression_opts
    )

    # Save dotplot of number of cells for each sample in each cluster
    df = adata.obs[['experiment_id', 'cluster']]
    df_ncells_cluster_per_sample = df.groupby(
        ['cluster', 'experiment_id']
    ).size().reset_index(name='nr_cells_cluster_smpl')
    df_ncells_cluster = df.groupby(
        ['cluster']
    ).size().reset_index(name='nr_cells_cluster')
    df_plt = pd.merge(
        df_ncells_cluster_per_sample,
        df_ncells_cluster,
        on='cluster'
    )
    df_plt['frac_cells_cluster_smpl'] = (
        df_plt['nr_cells_cluster_smpl'] / df_plt['nr_cells_cluster']
    )
    df_plt['n_cells_less_5'] = df_plt['nr_cells_cluster_smpl'] < 5
    df_plt['n_cells_zero'] = df_plt['nr_cells_cluster_smpl'] == 0

    gplt = plt9.ggplot(df_plt, plt9.aes(x='experiment_id', y='cluster'))
    gplt = gplt + plt9.theme_bw()
    gplt = gplt + plt9.theme(axis_text_x=plt9.element_text(angle=90))

    gplta = gplt + plt9.geom_point(
        plt9.aes(
            size='frac_cells_cluster_smpl',
            color='nr_cells_cluster_smpl',
            alpha='frac_cells_cluster_smpl'
        )
    )
    gplta = gplta + plt9.labs(
        title='',
        y='Cluster',
        x='',
        color='# of cells',
        size='Fraction of cells',
        alpha='Fraction of cells'
    )
    gplta.save(
        'dotplot_sample-{}.pdf'.format(out_file_base),
        width=8.5,
        height=8.5
    )

    gpltb = gplt + plt9.geom_point(
        plt9.aes(
            size='frac_cells_cluster_smpl',
            color='n_cells_less_5'
        ),
        alpha=0.75
    )
    gpltb = gpltb + plt9.scale_color_brewer(
        palette='Dark2',
        type='qual'
    )
    gpltb = gpltb + plt9.labs(
        title='',
        y='Cluster',
        x='',
        color='# of cells < 5',
        size='Fraction of cells'
    )
    gpltb.save(
        'dotplot_sample-{}_ncellsless5.pdf'.format(out_file_base),
        width=8.5,
        height=8.5
    )

    gpltb = gplt + plt9.geom_point(
        plt9.aes(
            size='frac_cells_cluster_smpl',
            color='n_cells_zero'
        ),
        alpha=0.75
    )
    gpltb = gpltb + plt9.scale_color_brewer(
        palette='Dark2',
        type='qual'
    )
    gpltb = gpltb + plt9.labs(
        title='',
        y='Cluster',
        x='',
        color='Zero cells',
        size='Fraction of cells'
    )
    gpltb.save(
        'dotplot_sample-{}_ncells0.pdf'.format(out_file_base),
        width=8.5,
        height=8.5
    )


if __name__ == '__main__':
    main()
