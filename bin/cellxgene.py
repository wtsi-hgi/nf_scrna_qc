#!/usr/bin/env python

__date__ = '2020-05-29'
__version__ = '0.0.1'

import argparse
import os
import numpy as np
import scipy as sp
import scanpy as sc


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object and preps for cellxgene pipeline.
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
        '--drop_extra_info',
        action='store_true',
        dest='drop_extra_info',
        default=False,
        help='Drops extra info from file to make a smaller file.\
            This also enables cellxgene to run faster and take up less memory.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: <h5_anndata>-<tsv_pcs>-cellxgene)'
    )

    options = parser.parse_args()

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}-cellxgene'.format(
            os.path.basename(options.h5.rstrip('h5ad').rstrip('.'))
        )

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=options.h5)

    # Format file as described here:
    # https://cellgeni.readthedocs.io/en/latest/visualisations.html

    # Set X to cp10k
    # warnings.warn(
    #     'WARNING: All functions in this script use adata.raw.to_adata(),',
    #     'assuming that adata.raw.to_adata stores ln(CPM+1) or ln(CP10K+1)',
    #     'normalized data filtered for genes where expression >5.'
    # )
    # # Use the raw slot which has no filters for genes expressed in
    # # >=5 cells, but is on scale log1p_cp10k.
    # # This will change X and adata.var.columns
    # adata = adata.raw.to_adata()
    # adata.X = np.expm1(adata.X)
    adata.X = np.expm1(adata.layers['log1p_cp10k'])

    # cellxgene works faster when the expression matrix is stored in CSC
    # (compressed sparse column) format instead of CSR
    # (compressed sparse row) format
    adata.X = sp.sparse.csc_matrix(adata.X)

    # For the query on cellxgene, change to gene ids rather than ensembl
    adata.var['gene_ids'] = adata.var.index
    adata.var['index_values'] = adata.var['gene_symbols'].replace(dict(zip(
        adata.var['gene_symbols'].values,
        adata.var['gene_symbols']
    )))
    adata.var = adata.var.set_index('index_values')
    adata.var_names_make_unique()

    # List of columns to drop to minimize data size for cellxgene
    if options.drop_extra_info:
        cols_drop_obs = set([])

        # Get every variable that is only one value across the whole dataset
        cols_drop_obs = cols_drop_obs.union(set(
            adata.obs.columns[adata.obs.nunique() <= 1].to_list()
        ))

        # Drop basic qc variables
        cols_drop_obs = cols_drop_obs.union(set([
            'estimated_number_of_cells',
            'mean_reads_per_cell',
            'median_genes_per_cell',
            'number_of_reads',
            'sequencing_saturation',
            'q30_bases_in_barcode',
            'q30_bases_in_rna_read',
            'q30_bases_in_sample_index',
            'q30_bases_in_umi',
            'reads_mapped_to_genome',
            'reads_mapped_confidently_to_genome',
            'reads_mapped_confidently_to_intergenic_regions',
            'reads_mapped_confidently_to_intronic_regions',
            'reads_mapped_confidently_to_exonic_regions',
            'reads_mapped_confidently_to_transcriptome',
            'reads_mapped_antisense_to_gene',
            'fraction_reads_in_cells',
            'total_genes_detected',
            'median_umi_counts_per_cell',
            'n_genes_by_counts',
            'log1p_n_genes_by_counts',
            'total_counts',
            'log1p_total_counts',
            'pct_counts_in_top_50_genes',
            'pct_counts_in_top_100_genes',
            'pct_counts_in_top_200_genes',
            'pct_counts_in_top_500_genes',
            'total_counts_gene_group__mito_transcript',
            'log1p_total_counts_gene_group__mito_transcript',
            # 'pct_counts_gene_group__mito_transcript',
            'total_counts_gene_group__mito_protein',
            'log1p_total_counts_gene_group__mito_protein',
            # 'pct_counts_gene_group__mito_protein',
            'total_counts_gene_group__ribo_protein',
            'log1p_total_counts_gene_group__ribo_protein',
            # 'pct_counts_gene_group__ribo_protein',
            'total_counts_gene_group__ribo_rna',
            'log1p_total_counts_gene_group__ribo_rna',
            # 'pct_counts_gene_group__ribo_rna',
            'normalization_factor'
        ]))
        # Drop any duplicate clusters, assuming main clusters stored in
        # 'clusters'
        if 'cluster' in adata.obs.columns:
            cols_drop_obs = cols_drop_obs.union(
                [col for col in adata.obs.columns if 'leiden' in col]
            )
            cols_drop_obs = cols_drop_obs.union(
                [col for col in adata.obs.columns if 'louvein' in col]
            )

        # Drop signatures with hvg_only
        cols_drop_obs = cols_drop_obs.union(
            [col for col in adata.obs.columns if '__hvg_only' in col]
        )

        # Drop other optional columns... don't worry this will not throw an
        # error if not present
        cols_drop_obs = cols_drop_obs.union(set([
            'batch',
            'scrublet__multiplet_scores'
        ]))

        # Drop the columns
        adata.obs = adata.obs.drop(columns=cols_drop_obs, errors='ignore')

        # Remove other data we cannot view in cellxgene
        del adata.layers
        del adata.obsp
        del adata.uns
        del adata.varm
        del adata.raw

        # NOTE: could remove extra columns from var
        cols_drop_var = set([
            'gene_group__mito_transcript',
            'gene_group__mito_protein',
            'gene_group__ribo_protein',
            'gene_group__ribo_rna',
            'n_cells_by_counts',
            'mean_counts',
            'log1p_mean_counts',
            'pct_dropout_by_counts',
            'total_counts',
            'log1p_total_counts',
            'n_cells',
            # 'highly_variable',
            'means',
            'dispersions',
            'dispersions_norm',
            'highly_variable_nbatches',
            'highly_variable_intersection',
            'mean',
            'std'
        ])
        adata.var = adata.var.drop(columns=cols_drop_var, errors='ignore')

    # Try to set default UMAP using umap_min_dist=1pt0 and umap_spread=1pt0
    umap_default = [i for i in adata.obsm.keys() if 'umap_min_dist=1pt0' in i]
    if len(umap_default) == 0:
        umap_default = [
            i for i in adata.obsm.keys() if 'umap_spread=1pt0' in i
        ]
        if len(umap_default) > 0:
            adata.obsm['X_umap'] = adata.obsm[umap_default.pop()]
    elif len(umap_default) >= 0:
        umap_default2 = [
            i for i in umap_default if 'umap_spread=1pt0' in i
        ]
        if len(umap_default2) > 0:
            adata.obsm['X_umap'] = adata.obsm[umap_default2.pop()]
        else:
            adata.obsm['X_umap'] = adata.obsm[umap_default.pop()]
    # NOTE: an alternative would be to run UMAPs here.

    # For coloring, make sure categorical variables are properly set to
    # Categorical
    # adata.obs['metadata_name'] = pd.Categorical(adata.obs['metadata_name'])
    # adata.obs['metadata_name'] = np.float32(adata.obs['metadata_name'])

    # Save the resulting anndata
    adata.write(
        '{}.h5ad'.format(out_file_base),
        compression='gzip',
        compression_opts=9  # takes ages, but we want a small file for system
    )


if __name__ == '__main__':
    main()
