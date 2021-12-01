#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-05-28'
__version__ = '0.0.1'

import argparse
import scanpy as sc
import plotnine as plt9


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object and list of phenotypes. Plot boxplots of \
            phenotypes across clusters.
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
        '--pheno_columns',
        action='store',
        dest='pheno_columns',
        default='',
        help='Pheno column to be boxplotted by cluster.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='plot_boxplot_cluster',
        help='Basename of output png file. Will have .png appended.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    adata = sc.read_h5ad(filename=options.h5)

    pheno_to_plot = options.pheno_columns.split(',')

    plt_height = 4
    plt_width = 16

    # Plot the data.
    for pheno in pheno_to_plot:
        # plt_width = adata.obs['cluster'].nunique() * 0.25

        gplt = plt9.ggplot(adata.obs)
        gplt = gplt + plt9.geom_boxplot(plt9.aes(x='cluster', y=pheno))
        gplt = gplt + plt9.theme(axis_text_x=plt9.element_text(angle=90))
        gplt.save(
            'boxplot-{}.png'.format(pheno),
            #dpi=300,
            width=plt_width,
            height=plt_height
        )

        # Add log10 transformation plot
        lab = 'log10'
        if adata.obs[pheno].min() < 0:
            adata.obs[pheno] = adata.obs[pheno] + abs(
                adata.obs[pheno].min()
            ) + 1
            lab = 'plusmin1log10'
        elif adata.obs[pheno].min() == 0:
            adata.obs[pheno] = adata.obs[pheno] + 1
            lab = 'plus1log10'
        gplt = plt9.ggplot(adata.obs)
        gplt = gplt + plt9.geom_boxplot(plt9.aes(x='cluster', y=pheno))
        gplt = gplt + plt9.theme(axis_text_x=plt9.element_text(angle=90))
        gplt = gplt + plt9.scale_y_continuous(
            trans='log10',
            # labels=comma_labels,
            minor_breaks=0
        )
        gplt.save(
            'boxplot_{}-{}.png'.format(lab, pheno),
            #dpi=300,
            width=plt_width,
            height=plt_height
        )


if __name__ == '__main__':
    main()
