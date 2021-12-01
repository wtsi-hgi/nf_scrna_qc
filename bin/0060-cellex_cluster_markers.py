#!/usr/bin/env python

__author__ = 'Moritz Przybilla'
__date__ = '2020-10-27'
__version__ = '0.0.1'

import argparse
import os
import random
import numpy as np
from distutils.version import LooseVersion
import cellex
import pandas as pd
import scanpy as sc
import time
from datetime import timedelta

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)


def run_CELLEX(
    adata,
    output_file,
    top_n=100
):
    """Use CELLEX to prioritize marker genes for clusters.

    Parameters
    ----------
    adata : AnnData
        Input AnnData file.
    output_file : string
        Basename of output_file, will have -normalized_pca.h5ad appended to it.
    verbose : string
        Specify whether summary should be printed or not using True or False.

    Returns
    -------
    output_file : string
        output_file
    """
    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # create a dataframe from the count matrix
    data = pd.DataFrame(
        adata.layers['counts'].toarray(),
        index=adata.obs_names,
        columns=adata.var_names
    ).transpose()

    # establish the metadata file from the annData
    # assign the column with the cluster id
    metadata = pd.DataFrame(
        data={"cluster": adata.obs["cluster"]}, index=adata.obs_names)

    # check if the data shape matches up
    if (adata.shape[0] == len(metadata)):
        # The ESObject encapsulates the core features of CELLEX. We set
        # "verbose=True" to get some progress reports.
        # The computations may take a while depending on the data and
        # available computational power.
        eso = cellex.ESObject(data=data, annotation=metadata, verbose=True)
        eso.compute(verbose=True)
    else:
        raise Exception(
            'The number of cells provided in the CELLEX data does not \
            match the number of cells in the metadata.'
        )

    # save the CELLEX results - the specificity value and std dev
    keys = ['esmu', 'essd']
    for key in keys:
        df = eso.results[key].reset_index()
        df = df.rename(columns={'gene': 'ensembl_gene_id'})

        # Add the HUGO gene name
        df['gene_symbols'] = adata.var.loc[
            df['ensembl_gene_id'], 'gene_symbols'
        ].values
        cols = list(df.columns)
        cols = [cols[-1]] + cols[:-1]
        df = df[cols]

        # Save the full gene matrix
        df.to_csv(
            "{}-{}.tsv.gz".format(output_file, key),
            sep='\t',
            compression=compression_opts,
            index=False,
            header=True
        )

        if key == 'esmu':
            # Sort df_esmu by value per cluster and take the top n genes
            df_long = pd.melt(
                df,
                id_vars=['ensembl_gene_id', 'gene_symbols'],
                var_name='cluster'
            )
            df_long = df_long.sort_values(by=['value'], ascending=False)
            df_results = df_long.groupby(
                'cluster'
            ).head(top_n).reset_index(drop=True)
            df_results = df_results.sort_values(
                by=['cluster'],
                ascending=False
            )
            df_results.to_csv(
                "{}-{}-top_markers.tsv.gz".format(output_file, key),
                sep='\t',
                compression=compression_opts,
                index=False,
                header=True
            )

            # Take top n genes again, but filtering for genes that we excluded
            # from the highly variable gene analysis
            filt = df_long['ensembl_gene_id'].isin(
                adata.uns['df_highly_variable_gene_filter']['ensembl_gene_id']
            )
            df_long = df_long.loc[np.invert(filt), :]
            df_results = df_long.groupby(
                'cluster'
            ).head(top_n).reset_index(drop=True)
            df_results = df_results.sort_values(
                by=['cluster'],
                ascending=False
            )
            df_results.to_csv(
                "{}-{}-top_markers__user_variable_genes_exclude.tsv.gz".format(
                    output_file,
                    key
                ),
                sep='\t',
                compression=compression_opts,
                index=False,
                header=True
            )

    return eso.results['esmu']


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read anndata object. Run CELLEX and save the expression specificity
            object to file.
            """
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='sc_cellex',
        help='Basename of output files.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-vb', '--verbose',
        action='store',
        dest='vb',
        default='True',
        help='Specify if the output should be printed.'
    )

    options = parser.parse_args()

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)

    start_time = time.time()
    _ = run_CELLEX(
        adata,
        output_file=options.of
    )

    execution_summary = "Analysis execution time [{}]:\t{}".format(
        "cellex_cluster_markers.py",
        str(timedelta(seconds=time.time()-start_time))
    )
    if options.vb:
        print(execution_summary)


if __name__ == '__main__':
    main()
