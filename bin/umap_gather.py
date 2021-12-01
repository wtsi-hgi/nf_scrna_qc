#!/usr/bin/env python

__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import numpy as np
import scanpy as sc


def add_info_to_adata_file(adata_root, adata_to_add):
    """Add UMAP related info from adata_to_add to adata_root."""
    # Check that the var and obs columns match up.
    assert np.array_equal(adata_root.obs.index, adata_to_add.obs.index)
    assert np.array_equal(adata_root.var.index, adata_to_add.var.index)

    # Update the naming of the data slots based on umap params.
    update_slot_string = ''
    for key, value in adata_to_add.uns['umap']['params'].items():
        if key not in ['a', 'b']:
            # We ignore keys that are calculated by umap and not set by user
            value = str(value)
            if value == '':
                value = 'None'
            update_slot_string = '{}{}={},'.format(
                update_slot_string,
                key,
                value.replace('.', 'pt')
            )
    update_slot_string = update_slot_string.rstrip(',')

    # Optionally update density slot
    if 'density' in adata_to_add.uns['umap']['params'] and adata_to_add.uns[
        'umap'
    ]['params']['density'] != '':
        assert adata_to_add.uns['umap']['params']['density'] == 'umap_density'
        adata_root.obs[
            'umap_density__{}'.format(update_slot_string)
        ] = adata_to_add.obs.pop('umap_density')
        adata_root.uns[
            'umap_density_params__{}'.format(update_slot_string)
        ] = adata_to_add.uns.pop('umap_density_params')

    # Update the UMAP slots
    tmp_dict = adata_to_add.uns.pop('umap')
    adata_root.uns[
        'umap_params__{}'.format(update_slot_string)
    ] = tmp_dict['params']
    del tmp_dict
    adata_root.obsm[
        'X_umap__{}'.format(update_slot_string)
    ] = adata_to_add.obsm.pop('X_umap')

    # Update the neighbors slot
    adata_root.uns[
        'neighbors__{}'.format(update_slot_string)
    ] = adata_to_add.uns.pop('neighbors')

    return adata_root


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read list AnnData object and merges UMAP data slots. NOTE: all
            data is 'appended' to the first AnnData object, and only UMAP,
            connectivity, and optionally density data are added.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-h5', '--h5_anndata_list',
        action='store',
        dest='h5',
        required=True,
        help='Comma seperated list of H5 AnnData files.'
    )

    parser.add_argument(
        '-h5_root', '--h5_root',
        action='store',
        dest='h5_root',
        default='',
        required=False,
        help='Optional h5ad file to use as root.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working\
            directory.\
            (default: <first_h5_anndata_file_name>-umap_sweep.h5ad)'
    )

    options = parser.parse_args()

    # Load the AnnData files
    adata_file_list = options.h5.split(',')

    # Read in the first AnnData file.
    if options.h5_root == '':
        first_file = adata_file_list.pop(0)
        adata = sc.read_h5ad(filename=first_file)
        adata = add_info_to_adata_file(adata, adata)
    else:
        adata = sc.read_h5ad(filename=options.h5_root)

    # Loop over the remaining AnnData files and:
    # 1. Read in the parameters.
    # 2. Copy the apporpriate dataslots to the first adata object to the
    #    appropriate new slot based on the parameters recorded in the original
    #    AnnData object.
    for file in adata_file_list:
        adata_i = sc.read_h5ad(filename=file)
        adata = add_info_to_adata_file(adata, adata_i)

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}-umap'.format(
            first_file
        )

    adata.write(
        '{}.h5ad'.format(out_file_base),
        compression='gzip'
    )


if __name__ == '__main__':
    main()
