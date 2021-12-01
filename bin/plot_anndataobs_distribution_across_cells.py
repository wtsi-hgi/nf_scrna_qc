#!/usr/bin/env python

__date__ = '2020-04-08'
__version__ = '0.0.1'

import argparse
import numpy as np
import scanpy as sc
import plotnine as plt9


def comma_labels(x_list):
    """Change list of int to comma format."""
    result = []
    for x in x_list:
        result.append(format(int(x), ','))
    return(result)


def plot_histogram(
    df_plot,
    variable_column,
    output_file='plot_distribution',
    facet_column='none',
    x_log10=False
):
    """Plot plot_distribution to png.

    Parameters
    ----------
    df_plot : pandas.DataFrame
        DataFrame with <variable_column> as a column.
    variable_column : string
        String of variable_column column to plot.
    output_file : string
        Basename of output file.
    facet_column : string
        Column to facet the plot by.

    Returns
    -------
    NULL
    """
    df_plot['x'] = df_plot[variable_column]
    if x_log10:
        if np.any(df_plot['x'].values < 0):
            return 1
        elif np.any(df_plot['x'].values == 0):
            df_plot['x'] = np.log10(df_plot['x'].values + 1e-10)
            variable_column = variable_column + ' (log10)'
        else:
            df_plot['x'] = np.log10(df_plot['x'].values)
            variable_column = variable_column + ' (log10)'
    gplt = plt9.ggplot(df_plot, plt9.aes(
        x='x'
    ))
    gplt = gplt + plt9.theme_bw()
    gplt = gplt + plt9.geom_histogram(alpha=0.8)
    gplt = gplt + plt9.scale_x_continuous(
        # trans='log10',
        # labels=comma_labels,
        minor_breaks=0
    )
    gplt = gplt + plt9.scale_y_continuous(
        # trans='log10',
        # labels=comma_labels,
        minor_breaks=0
    )
    gplt = gplt + plt9.labs(
        title='',
        x=variable_column
    )
    gplt = gplt + plt9.theme(
        axis_text_x=plt9.element_text(angle=-45, hjust=0)
    )
    if facet_column != 'none':
        gplt = gplt + plt9.facet_wrap(
            '~ {}'.format(facet_column),
            ncol=5
        )
        n_facets = df_plot[facet_column].nunique()
        gplt.save(
            '{}.png'.format(output_file),
            # dpi=300,
            width=6*(n_facets/4),
            height=4*(n_facets/4),
            limitsize=False
        )
    else:
        gplt.save(
            '{}.png'.format(output_file),
            # dpi=300,
            width=4,
            height=4
        )
    return 0


def plot_ecdf(
    df_plot,
    variable_column,
    color_column='none',
    output_file='plot_distribution',
    facet_column='none',
    x_log10=False
):
    """Plot plot_distribution to png.

    Parameters
    ----------
    df_plot : pandas.DataFrame
        DataFrame with <variable_column> as a column.
    variable_column : string
        String of variable_column column to plot.
    color_column : string
        String of color column to plot.
    output_file : string
        Basename of output file.
    facet_column : string
        Column to facet the plot by.

    Returns
    -------
    NULL
    """
    n_colors = 0
    if color_column != 'none':
        gplt = plt9.ggplot(df_plot, plt9.aes(
            x=variable_column,
            color=color_column
        ))
        n_colors = df_plot[color_column].nunique()
    else:
        gplt = plt9.ggplot(df_plot, plt9.aes(
            x=variable_column
        ))
    gplt = gplt + plt9.theme_bw()
    gplt = gplt + plt9.stat_ecdf(alpha=0.8)
    if x_log10:
        gplt = gplt + plt9.scale_x_continuous(
            trans='log10',
            # labels=comma_labels,
            minor_breaks=0
        )
    else:
        gplt = gplt + plt9.scale_x_continuous(
            # trans='log10',
            # labels=comma_labels,
            minor_breaks=0
        )
    gplt = gplt + plt9.scale_y_continuous(
        # trans='log10',
        # labels=comma_labels,
        minor_breaks=0
    )
    gplt = gplt + plt9.labs(
        y='Cumulative density',
        title=''
    )
    if n_colors != 0 and n_colors > 20:
        gplt = gplt + plt9.theme(
            legend_position='none'
        )
    elif n_colors != 0 and n_colors < 9:
        gplt = gplt + plt9.scale_colour_brewer(
            palette='Dark2',
            type='qual'
        )
    if facet_column != 'none':
        gplt = gplt + plt9.facet_wrap(
            '~ {}'.format(facet_column),
            ncol=5
        )
        n_facets = df_plot[facet_column].nunique()
        gplt.save(
            '{}.png'.format(output_file),
            # dpi=300,
            width=6*(n_facets/4),
            height=4*(n_facets/4),
            limitsize=False
        )
    else:
        gplt.save(
            '{}.png'.format(output_file),
            # dpi=300,
            width=4,
            height=4
        )
    return 0


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
        '-var', '--variable_columns',
        action='store',
        dest='var',
        help='Comma seperated list of variable columns to plot.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='plot_distribution',
        help='Basename of output png file. Will have .png appended.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '--ecdf',
        action='store_true',
        dest='ecdf',
        default=False,
        help='Plot ecdf rather than histogram.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '--facet_columns',
        action='store',
        dest='facet_columns',
        default='none',
        help='Column to facet plot by. If --ecdf then used as color.\
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

    # Get a list of variables to plot.
    variables_to_plot = options.var.split(',')
    if len(variables_to_plot) == 0:
        raise Exception('No variable.')

    # Get a list of the facets to plot.
    facets_to_plot = options.facet_columns.split(',')
    if len(facets_to_plot) == 0:
        facets_to_plot = ['none']

    for facet in facets_to_plot:
        for var in variables_to_plot:
            # Ensure our facet variable is a string (e.g., if batch = [1,2])
            if facet != 'none':
                adata.obs[facet] = adata.obs[facet].apply(str)
            # Plot the data.
            if options.ecdf:
                _ = plot_ecdf(
                    df_plot=adata.obs,
                    variable_column=var,
                    output_file='{}.var={}.color={}-{}'.format(
                        'plot_ecdf',
                        var,
                        facet,
                        options.of
                    ),
                    color_column=facet,
                    x_log10=False
                )
                _ = plot_ecdf(
                    df_plot=adata.obs,
                    variable_column=var,
                    output_file='{}.var={}.color={}-{}'.format(
                        'plot_ecdf-x_log10',
                        var,
                        facet,
                        options.of
                    ),
                    color_column=facet,
                    x_log10=True
                )
            else:
                _ = plot_histogram(
                    df_plot=adata.obs,
                    variable_column=var,
                    output_file='{}.var={}.facet={}-{}'.format(
                        'plot_histogram',
                        var,
                        facet,
                        options.of
                    ),
                    facet_column=facet,
                    x_log10=False
                )
                _ = plot_histogram(
                    df_plot=adata.obs,
                    variable_column=var,
                    output_file='{}.var={}.facet={}-{}'.format(
                        'plot_histogram-x_log10',
                        var,
                        facet,
                        options.of
                    ),
                    facet_column=facet,
                    x_log10=True
                )


if __name__ == '__main__':
    main()
