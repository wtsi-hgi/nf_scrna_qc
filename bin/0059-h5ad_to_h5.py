#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-05-01'
__version__ = '0.0.1'


import argparse
import scanpy as sc
import h5py
import os
import pandas as pd


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object with cluster information. Saves normalized \
            matrix and cluster information as h5 and csv files.
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
         default='adata-log1p_cp10k',
         help='Basename of output files, assuming output in current working \
             directory.\
             (default: %(default)s)'
    )

    options = parser.parse_args()

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=options.h5)

    # Save log1p cp10k matrix
    out_file_base = options.of
    hf = h5py.File('{}.h5'.format(out_file_base), 'w')
    hf.create_dataset(
                    'X',
                    data=adata.layers['log1p_cp10k'].todense(),
                    compression='gzip'
    )
    hf.create_dataset('genes', data=adata.var.index, compression='gzip')
    hf.create_dataset('cells', data=adata.obs.index, compression='gzip')
    hf.create_dataset(
                    'cluster',
                    data=pd.DataFrame(data=adata.obs['cluster']),
                    compression='gzip')
    hf.close()


if __name__ == '__main__':
    main()
