#!/usr/bin/env python

__date__ = '2020-06-26'
__version__ = '0.0.1'

import argparse
import scanpy as sc
import numpy as np
import plotnine as plt9
import scipy.stats as stat

def comma_labels(x_list):
    """Change list of int to comma format."""
    result = []
    for x in x_list:
        result.append(format(int(x), ','))
    return(result)


def calculate_density(df_plot, facet_column):
    """
    Calculate the density for each category in the facet. Return an array
    that corresponds to the row of the df_plot
    """
    if facet_column == 'none':
        den_array = np.vstack([
            np.log10(df_plot['total_counts']),
            df_plot['pct_counts_gene_group__mito_transcript']
        ])
        return(stat.gaussian_kde(den_array)(den_array))

    sorted_df = df_plot.sort_values(facet_column)
    density_array = []
    for facet_key in np.unique(sorted_df[facet_column]):
        sub_df = sorted_df[sorted_df[facet_column] == facet_key]
        density_array += calculate_density(sub_df, 'none').tolist()
    sorted_df['temp_12345676'] = density_array
    return(sorted_df.loc[df_plot.index.tolist(), "temp_12345676"].values)


def plot_umi_mt_density(
    df_plot,
    output_file='plot_umi_mt_density',
    facet_column='none',
    color_var='density',
    density_contour=False
):
    """Plot plot_umi_mt_density to png.

    Parameters
    ----------
    df_plot : pandas.DataFrame
        DataFrame with the followig keys 'total_counts',
        'pct_counts_gene_group__mito_transcript'.
    output_file : string
        Basename of output file.
    facet_column : string
        Column to facet the output by.

    Returns
    -------
    NULL
    """
    if color_var == 'density':
        color_title = 'Density\n'
        # Also calculate density using a gaussian 2d kernal -- use random
        # name for plot column
        color_var = "1251234_density"
        df_plot[color_var] = calculate_density(df_plot, facet_column)
    elif color_var == 'pct_counts_gene_group__mito_transcript':
        color_title = '% MT\n'
    elif color_var == 'cell_passes_qc':
        color_title = 'Cell passed QC\n'
    else:
        color_title = color_var
    gplt = plt9.ggplot(df_plot, plt9.aes(
        x='total_counts',
        y='pct_counts_gene_group__mito_transcript',
        color=color_var
    ))
    gplt = gplt + plt9.theme_bw()
    gplt = gplt + plt9.geom_point(alpha=0.5, size=0.8)
    gplt = gplt + plt9.scale_x_continuous(
        trans='log10',
        labels=comma_labels,
        minor_breaks=0
    )
    if color_var == 'pct_counts_gene_group__mito_transcript':
        gplt = gplt + plt9.scale_color_gradient2(
            low='#3B9AB2',
            mid='#EBCC2A',
            high='#F21A00',
            midpoint=50,
            limits=[0, 100]
        )
        gplt = gplt + plt9.guides(color=plt9.guide_colorbar(ticks=False))
    elif color_var == 'cell_passes_qc':
        gplt = gplt + plt9.scale_colour_brewer(type='qual', palette='Dark2')
    elif color_var == '1251234_density':
        gplt = gplt + plt9.scale_color_cmap(cmap_name='viridis')

    if density_contour:
        gplt = gplt + plt9.geom_density_2d(alpha=0.5)
    gplt = gplt + plt9.labs(
        x='Number of molecules',
        y='Percent of molecules from MT genes',
        title='',
        color=color_title
    )
    if facet_column != 'none':
        gplt = gplt + plt9.facet_wrap(
            '~ {}'.format(facet_column),
            ncol=5
        )
        n_samples = df_plot[facet_column].nunique()
        gplt.save(
            '{}.png'.format(output_file),
            # dpi=300,
            width=4*(n_samples/2),
            height=4*(n_samples/4),
            limitsize=False
        )
    else:
        gplt.save(
            '{}.png'.format(output_file),
            # dpi=300,
            width=4,
            height=4
        )


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            QC plots.
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
        '-of', '--output_file',
        action='store',
        dest='of',
        default='plot_umi_mt_density',
        help='Basename of output png file. Will have .png appended.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '--facet_columns',
        action='store',
        dest='facet_columns',
        default='none',
        help='Column to facet plot by.\
            (default: No facet.)'
    )

    options = parser.parse_args()

    # Scanpy settings
    # sc.settings.figdir = os.getcwd()  # figure output directory 2 match base.
    # # sc.settings.n_jobs = options.ncpu  # number CPUs
    # # sc.settings.max_memory = 500  # in Gb
    # sc.set_figure_params(dpi_save=300)

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=options.h5)

    # Get a list of the facets to plot.
    facets_to_plot = options.facet_columns.split(',')
    if len(facets_to_plot) == 0:
        facets_to_plot = ['none']

    # Plot the data.
    for facet in facets_to_plot:
        plot_umi_mt_density(
            df_plot=adata.obs,
            output_file='plot_umi_mt_density.facet={}-{}'.format(
                facet,
                options.of
            ),
            facet_column=facet,
            color_var='density',
            density_contour=True
        )
        if 'cell_passes_qc' in adata.obs:
            plot_umi_mt_density(
                df_plot=adata.obs,
                output_file='plot_umi_mt_cellpassqc.facet={}-{}'.format(
                    facet,
                    options.of
                ),
                facet_column=facet,
                color_var='cell_passes_qc',
                density_contour=False
            )


if __name__ == '__main__':
    main()
