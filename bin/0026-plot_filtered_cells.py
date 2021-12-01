#!/usr/bin/env python

__date__ = '2020-05-29'
__version__ = '0.0.1'

import argparse
import os
import pandas as pd
import plotnine as plt9


def comma_labels(x_list):
    """Change list of int to comma format."""
    result = []
    for x in x_list:
        result.append(format(int(x), ','))
    return(result)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Filter and merge 10x data. Save to AnnData object.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '--tsv_file',
        action='store',
        dest='tsv',
        required=True,
        help='cell_filtered_per_experiment tsv file.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output png file. Will have .png appended.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Get basename of the output file
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}'.format(
            os.path.basename(options.tsv.rstrip('tsv.gz').rstrip('\\.'))
        )

    # Load the data
    df = pd.read_csv(options.tsv, sep='\t')

    # Get the total number of input cells per sample
    df_before_filters = df[df.filter_type.isin(['before_filters'])]
    df_before_filters = df_before_filters.set_index('experiment_id')

    # Check if any difference between before and after filters.	If not,
    # return early.
    df_after_filters = df[df.filter_type.isin(['after_filters'])]
    filt = df_after_filters.n_cells_left_in_adata == df_before_filters.loc[
        df_after_filters.experiment_id,
        'n_cells_left_in_adata'
    ].values
    if all(filt):
        print("No difference detected before and after filters. No plots.")
        return()

    # Set some plotting parameters
    plt_height = 16  # 1.5 * df.experiment_id.nunique()

    # Plot the number of cells before and after all filters across experiments
    df_plt = df[df.filter_type.isin(['before_filters', 'after_filters'])]
    gplt = plt9.ggplot(df_plt, plt9.aes(
        x='experiment_id',
        y='n_cells_left_in_adata',
        # label='n_cells',
        fill='filter_type'
    ))
    gplt = gplt + plt9.theme_bw()
    gplt = gplt + plt9.geom_bar(stat='identity', position='dodge')
    # gplt = gplt + plt9.geom_text(vjust=1.6, color='white', size=3.5)
    gplt = gplt + plt9.scale_y_continuous(
        trans='log10',
        labels=comma_labels,
        minor_breaks=0
    )
    gplt = gplt + plt9.scale_fill_brewer(
        palette='Dark2',
        type='qual'
    )
    gplt = gplt + plt9.labs(
        title='',
        y='Number of cells',
        x='',
        fill=''
    )
    # NOTE: legend_position bug https://github.com/has2k1/plotnine/issues/245
    gplt = gplt + plt9.theme(
        # legend_position='bottom',
        subplots_adjust={'bottom': 0.15},
        legend_position=(.5, .05),
        legend_direction='horizontal',
        legend_title=plt9.element_blank()
    )
    gplt = gplt + plt9.coord_flip()
    gplt.save(
        '{}-n_cells_before_after.png'.format(out_file_base),
        #dpi=300,
        width=4,
        height=plt_height
    )

    # Plot the final fraction of cells filtered per experiment
    df_plt = df_after_filters.copy()
    # Invert the numbers, so instead of the number of cells that pass, get
    # the number of cells that fail at each filter.
    df_plt.n_cells_left_in_adata = df_before_filters.loc[
        df_plt.experiment_id, 'n_cells_left_in_adata'
    ].values - df_plt.n_cells_left_in_adata
    # Now calculate the fraction removed
    df_plt['fraction_cells'] = df_plt.n_cells_left_in_adata / \
        df_before_filters.loc[
            df_plt.experiment_id,
            'n_cells_left_in_adata'
        ].values
    gplt = plt9.ggplot(df_plt, plt9.aes(
        x='experiment_id',
        y='fraction_cells',
        fill='filter_type'
    ))
    gplt = gplt + plt9.theme_bw()
    gplt = gplt + plt9.geom_bar(stat='identity', position='dodge')
    if df_plt.filter_type.nunique() < 9:
        gplt = gplt + plt9.scale_fill_brewer(
            palette='Dark2',
            type='qual'
        )
    gplt = gplt + plt9.labs(
        title='',
        y='Fraction of total cells excluded',
        x='',
        fill='Filter'
    )
    # NOTE: legend_position bug https://github.com/has2k1/plotnine/issues/245
    gplt = gplt + plt9.theme(
        # legend_position='bottom',
        subplots_adjust={'bottom': 0.15},
        legend_position=(.5, .05),
        legend_direction='vertical'
    )
    gplt = gplt + plt9.coord_flip()
    gplt.save(
        '{}-fraction_before_after.png'.format(out_file_base),
        #dpi=300,
        width=4,
        height=plt_height
    )

    # Plot the number of cells falling into each filter acoss experiments.
    # NOTE: cells can fall into multiple filters.
    # Remove the rows that we do not want
    df_plt = df[~df.filter_type.isin(['before_filters', 'after_filters'])]
    df_plt = df_plt[~df_plt.filter_type.str.contains('after_filter')]
    # Invert the numbers, so instead of the number of cells that pass, get
    # the number of cells that fail at each filter.
    df_plt.n_cells_left_in_adata = df_before_filters.loc[
        df_plt.experiment_id, 'n_cells_left_in_adata'
    ].values - df_plt.n_cells_left_in_adata
    gplt = plt9.ggplot(df_plt, plt9.aes(
        x='experiment_id',
        y='n_cells_left_in_adata',
        fill='filter_type'
    ))
    gplt = gplt + plt9.theme_bw()
    gplt = gplt + plt9.geom_bar(stat='identity', position='dodge')
    if df_plt.filter_type.nunique() < 9:
        gplt = gplt + plt9.scale_fill_brewer(
            palette='Dark2',
            type='qual'
        )
    gplt = gplt + plt9.labs(
        title='',
        y='Number of cells excluded',
        x='',
        fill='Filter'
    )
    # NOTE: legend_position bug https://github.com/has2k1/plotnine/issues/245
    gplt = gplt + plt9.theme(
        # legend_position='bottom',
        subplots_adjust={'bottom': 0.15},
        legend_position=(.5, .05),
        legend_direction='vertical'
    )
    gplt = gplt + plt9.coord_flip()
    gplt.save(
        '{}-n_cells_excluded.png'.format(out_file_base),
        #dpi=300,
        width=4,
        height=plt_height
    )

    # Plot the ratio of the total number of cells removed in each filter across
    # experiments.
    # NOTE: cells can fall into multiple filters.
    df_plt['fraction_cells'] = df_plt.n_cells_left_in_adata / \
        df_before_filters.loc[
            df_plt.experiment_id,
            'n_cells_left_in_adata'
        ].values
    gplt = plt9.ggplot(df_plt, plt9.aes(
        x='experiment_id',
        y='fraction_cells',
        fill='filter_type'
    ))
    gplt = gplt + plt9.theme_bw()
    gplt = gplt + plt9.geom_bar(stat='identity', position='dodge')
    if df_plt.filter_type.nunique() < 9:
        gplt = gplt + plt9.scale_fill_brewer(
            palette='Dark2',
            type='qual'
        )
    gplt = gplt + plt9.labs(
        title='',
        y='Fraction of total cells excluded',
        x='',
        fill='Filter'
    )
    # NOTE: legend_position bug https://github.com/has2k1/plotnine/issues/245
    gplt = gplt + plt9.theme(
        # legend_position='bottom',
        subplots_adjust={'bottom': 0.15},
        legend_position=(.5, .05),
        legend_direction='vertical'
    )
    gplt = gplt + plt9.coord_flip()
    gplt.save(
        '{}-fraction_cells_excluded.png'.format(out_file_base),
        #dpi=300,
        width=4,
        height=plt_height
    )


if __name__ == '__main__':
    main()
