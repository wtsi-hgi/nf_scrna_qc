#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-05-29'
__version__ = '0.0.1'

import argparse
import os
import re
import scanpy as sc
import pandas as pd
import numpy as np
import warnings

# If tk.Tk issues
# import matplotlib
# matplotlib.use('Agg')


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object and markers from database. Plot expression \
            expression of markers in each cluster.
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
        '--markers_database_tsv',
        action='store',
        dest='markers_database_tsv',
        default='none',
        help='Markers to plot. Must have the following columns: cell_type, \
            hgnc_symbol.'
    )

    parser.add_argument(
        '--max_n_markers',
        action='store',
        dest='max_n_markers',
        default=100,
        type=int,
        help='Max number of markers to plot (default: %(default)s).'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: None)'
    )

    options = parser.parse_args()

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    # sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    #sc.set_figure_params(dpi_save=300)

    max_n_markers = options.max_n_markers

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}'.format(
            os.path.basename(options.h5.rstrip('h5ad').rstrip('.'))
        )

    # Read in the data
    adata = sc.read_h5ad(filename=options.h5)

    data_scale = 'log1p_cp10k'
    if data_scale == 'cp10k':  # Set X to cp10k
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
    elif data_scale == 'log1p_cp10k':  # Set X to log1p_cp10k
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
        adata.X = adata.layers['log1p_cp10k']
    elif data_scale == 'counts':
        adata.X = adata.layers['counts']
        # adata = adata.raw.to_adata() # This is not counts.

    # Read in the marker database file
    df = pd.read_table(options.markers_database_tsv)

    if 'p_value_adj' in df.columns and all(pd.notnull(df['p_value_adj'])):
        df = df[['cell_type', 'hgnc_symbol', 'p_value_adj']]
        df = df.sort_values(
            by=['cell_type', 'p_value_adj'],
            ascending=[True, True]
        )
        
    elif 'auc' in df.columns and all(pd.notnull(df['auc'])):
        df = df[['cell_type', 'hgnc_symbol', 'auc']]
        df = df.sort_values(
            by=['cell_type', 'auc'],
            ascending=[True, False]
        )
    elif 'power' in df.columns and all(pd.notnull(df['power'])):
        df = df[['cell_type', 'hgnc_symbol', 'power']]
        df = df.sort_values(
            by=['cell_type', 'power'],
            ascending=[True, False]
        )

    else:
        df = df[['cell_type', 'hgnc_symbol']]

    cell_type = np.unique(df['cell_type'])
    for i in cell_type:
        marker_genes = df.loc[df['cell_type'] == i]['hgnc_symbol'].values
        marker_genes_found = marker_genes[
            np.isin(marker_genes, adata.var['gene_symbols'])
        ]
        # Clean up cell id names
        i_out = re.sub(r'\W+', '', i)  # Strip all non alphanumeric characters
        print(i, i_out)
        # Trim markers to options.max_n_markers
        if len(marker_genes_found) > max_n_markers:
            marker_genes_found = marker_genes_found[0:max_n_markers]
        elif len(marker_genes_found) == 0:
            warnings.warn("No marker genes for cluster {}".format(i_out))
            continue
        # Dotplots
        # Generate dendrogram using the marker genes...this will be used in the
        # below dotplots.
        # NOTE: With latest version of pandas, sc.tl.dendrogram throws an error
        run_dendrogram = True
        if run_dendrogram:
            sc.tl.dendrogram(
                adata,
                groupby='cluster',
                use_rep='X_pca',
                var_names=marker_genes_found,
                use_raw=False,
                cor_method='pearson',
                linkage_method='complete',
                optimal_ordering=True,
                inplace=True
            )

        _ = sc.pl.dotplot(
            adata=adata,
            var_names=marker_genes_found,
            groupby='cluster',
            gene_symbols='gene_symbols',
            dendrogram=run_dendrogram,
            show=False,
            use_raw=False,
            log=False,
            color_map='Blues',
            save='-{}-{}-{}.png'.format(
                out_file_base,
                i_out,
                data_scale
            )
        )
        _ = sc.pl.dotplot(
            adata=adata,
            var_names=marker_genes_found,
            groupby='cluster',
            gene_symbols='gene_symbols',
            dendrogram=run_dendrogram,
            show=False,
            standard_scale='var',  # Scale color between 0 and 1
            use_raw=False,
            color_map='Blues',
            save='-{}-{}-{}_standardized.png'.format(
                out_file_base,
                i_out,
                data_scale
            )
        )

        # Plot heatmaps
        # _ = sc.pl.heatmap(
        #     adata=adata_raw,
        #     var_names=marker_genes_found,
        #     groupby='cluster',
        #     dendrogram=True,
        #     gene_symbols='gene_symbols',
        #     use_raw=False,
        #     show=False,
        #     save='-{}-{}.png'.format(out_file_base, i_out)
        # )
        # _ = sc.pl.heatmap(
        #     adata=adata_raw,
        #     var_names=marker_genes_found,
        #     groupby='cluster',
        #     dendrogram=True,
        #     gene_symbols='gene_symbols',
        #     standard_scale='var',  # Scale color between 0 and 1
        #     use_raw=False,
        #     show=False,
        #     save='-{}-{}-standardized.png'.format(out_file_base, i_out)
        # )


if __name__ == '__main__':
    main()
