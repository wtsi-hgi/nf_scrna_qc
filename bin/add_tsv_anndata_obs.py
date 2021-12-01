#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-04-29'
__version__ = '0.0.1'


import argparse
import pandas as pd
import scanpy as sc
import os


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Adds data from tsv file to adata obs slot.
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
        '--tsv_file',
        action='store',
        dest='tsv_file',
        help='Tab-delimited file to append to adata obs slot.'
    )

    parser.add_argument(
        '-of', '--out_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: <h5_anndata>-<tsv_file>)'
    )

    opt = parser.parse_args()

    # Get the out file base.
    out_file_base = opt.of
    if out_file_base == '':
        out_file_base = '{}-{}'.format(
            os.path.basename(opt.h5.rstrip('h5ad').rstrip('.'),
            os.path.basename(options.pc.rstrip('tsv.gz').rstript('.'))
        )

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=opt.h5)

    # Bind cluster merging information to Anndata file
    merging_progress = pd.read_table(opt.tsv_file, dtype='category')
    adata.obs = pd.concat([adata.obs, merging_progress], axis=1)

    # Save Anndata to .h5ad
    adata.write('{}.h5ad'.format(out_file_base), compression='gzip')


if __name__ == '__main__':
    main()
