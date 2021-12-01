#!/usr/bin/env python


__author__ = 'Tobi Alegbe'
__date__ = '2020-06-30'
__version__ = '0.0.1'

import argparse
import os
import scanpy as sc

# Nice large palette.
COLORS_LARGE_PALLETE = [
    '#0F4A9C', '#3F84AA', '#C9EBFB', '#8DB5CE', '#C594BF', '#DFCDE4',
    '#B51D8D', '#6f347a', '#683612', '#B3793B', '#357A6F', '#989898',
    '#CE778D', '#7F6874', '#E09D37', '#FACB12', '#2B6823', '#A0CC47',
    '#77783C', '#EF4E22', '#AF1F26'
]

def save_pc_fig(
    adata,
    out_file_base,
    color_var,
    num_PCs,
    colors_quantitative=True
):
    """
    Create and save a figure containing the cells plotted in the space of \
    two principal components.
    """
    # Plot the cells in the principal component space
    pc_pairs_to_plot = ["{},{}".format(i,i+1) for i in range(1,num_PCs+1,2)]

    color_palette = 'viridis'

    sc.pl.pca(
        adata=adata,
        color=color_var,
        palette=color_palette,
        components=pc_pairs_to_plot,
        alpha=0.4,
        save='-{}-{}.png'.format(
            out_file_base,
            color_var,
        )
    )

def save_pc_genes_fig(
    adata,
    out_file_base,
    num_PCs,
    drop_legend=-1,
):
    """
    Create and save a figure containing the top genes that contribute to a \
    principal component.
    """

    # Plot the top genes contributing to each principal component
    adata_temp = adata.copy()
    adata_temp.var_names = adata_temp.var.gene_symbols
    pcs_to_plot = ','.join(map(str,list(range(1,num_PCs+1))))

    sc.pl.pca_loadings(
        adata=adata_temp,
        include_lowest=False,
        components = pcs_to_plot,
        save='-{}-n_pcs={}.png'.format(
            out_file_base,
            num_PCs,
        )
    )

def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object with principal components calculated. Plot the \
            cells in principal component coordinates and plot top genes in each \
            component ranked by their contribution.
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
        '--num_pcs',
        action='store',
        dest='num_pcs',
        default='20',
        help='Number of principal components to plot'
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

    parser.add_argument(
        '-cq', '--colors_quantitative',
        action='store',
        dest='cq',
        default='',
        help='Comma seperated list of quantitative variable names for colors.\
            (default: "")'
    )

    parser.add_argument(
        '-cc', '--colors_categorical',
        action='store',
        dest='cc',
        default='',
        help='Comma seperated list of categorical variable names for colors.\
            (default: "")'
    )

    options = parser.parse_args()

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    # sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    # sc.set_figure_params(dpi_save=300)

    # Read in the data
    adata = sc.read_h5ad(filename=options.h5)
    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}'.format(
            os.path.basename(options.h5.rstrip('h5ad').rstrip('.'))
        )

    num_PCs = int(options.num_pcs)

    # Parse the color variables.
    colors_quantitative = []
    if options.cq != '':
        colors_quantitative = options.cq.split(',')

    colors_categorical = []
    if options.cc != '':
        colors_categorical = options.cc.split(',')

    # For each color to plot, loop over the different iterations.
    for color_var in colors_quantitative:
        save_pc_fig(
            adata=adata,
            out_file_base=out_file_base,
            color_var=color_var,
            colors_quantitative=True,
            num_PCs=num_PCs
        )
    for color_var in colors_categorical:
        save_pc_fig(
            adata=adata,
            out_file_base=out_file_base,
            color_var=color_var,
            colors_quantitative=False,
            num_PCs=num_PCs
        )
    save_pc_genes_fig(
        adata=adata,
        out_file_base=out_file_base,
        num_PCs=num_PCs
        )

if __name__ == '__main__':
    main()
