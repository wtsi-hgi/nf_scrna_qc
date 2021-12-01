#!/usr/bin/env python

__date__ = '2020-04-29'
__version__ = '0.0.1'

import argparse
from distutils.version import LooseVersion
import os
import pandas as pd
import csv
import plotnine as plt9
import harmonypy as hm


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Calcualte and compare LISI across a series of reduced dims and
            categorical variables.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    # parser.add_argument(
    #     '-h5', '--h5_anndata',
    #     action='store',
    #     dest='h5',
    #     required=True,
    #     help='H5 AnnData file.'
    # )

    parser.add_argument(
        '-rf', '--reduced_dims_tsv',
        action='store',
        dest='reduced_dims',
        required=True,
        help='List of tab-delimited files of reduced dimensions (e.g., PCs)\
            for each cell. First column is cell_barcode. List should be\
            split by "::" (e.g. file1.tsv.gz::file2.tsv.gz).'
    )

    parser.add_argument(
        '-lbl', '--reduced_dims_tsv_labels',
        action='store',
        dest='reduced_dims_labels',
        required=True,
        help='String of labels for each reduced_dims_tsv file. List should be\
            split by "::".'
    )

    parser.add_argument(
        '-mf', '--metadata_tsv',
        action='store',
        dest='metadata_tsv',
        required=True,
        help='Tab-delimited file of metadata for each cell. First column\
            is cell_barcode.'
    )

    parser.add_argument(
        '-mv', '--metadata_columns',
        action='store',
        dest='metadata_columns',
        default='experiment_id',
        help='Comma separated string of categorical variables to calculate\
            LISI with.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-p', '--perplexity',
        action='store',
        dest='perplexity',
        default=30.0,
        type=float,
        help='Perplexity.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: <metadata_tsv>-lisi)'
    )

    options = parser.parse_args()

    # Fixed settings.
    # verbose = True

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}-lisi'.format(
            os.path.basename(options.metadata_tsv.rstrip('tsv.gz').rstrip('.'))
        )

    # Get the columns to use
    lisi_columns = options.metadata_columns.split(',')
    # lisi_columns = ['experiment_id', 'batch']
    lisi_columns_dtype = dict(
        zip(lisi_columns, ['category']*len(lisi_columns))
    )

    # Load the metadata file
    file_meta = options.metadata_tsv
    df_meta = pd.read_csv(
        file_meta,
        sep='\t',
        index_col='cell_barcode',
        dtype=lisi_columns_dtype
    )

    # Load the reduced dims.
    files = options.reduced_dims.split('::')
    labels = options.reduced_dims_labels.split('::')
    assert len(files) == len(labels), 'ERROR: check files and labels input'

    # Make a dict of theoretical maximum LISI value for each label.
    lisi_limit = {}
    for col in lisi_columns:
        n_cat = len(df_meta[col].cat.categories)
        lisi_limit[col] = n_cat

    list_lisi = []
    for i in range(len(files)):
        df_reduced_dims = pd.read_csv(
            files[i],
            sep='\t',
            index_col='cell_barcode'
        )

        # Run lisi and save results to dataframe
        _df_lisi = pd.DataFrame(
            hm.compute_lisi(
                df_reduced_dims.loc[df_meta.index, :],
                df_meta[lisi_columns],
                lisi_columns
            ),
            columns=lisi_columns
        )
        _df_lisi['file'] = files[i]
        _df_lisi['label'] = labels[i]
        _df_lisi['cell_barcode'] = df_meta.index
        list_lisi.append(_df_lisi)

    # Make one long dataframe.
    df_lisi = pd.concat(list_lisi)
    # Make cell_barcode the first column.
    cols = list(df_lisi.columns)
    cols = [cols[-1]] + cols[:-1]

    # Save the results
    df_lisi[cols].to_csv(
        '{}.tsv.gz'.format(out_file_base),
        sep='\t',
        index=False,
        quoting=csv.QUOTE_NONNUMERIC,
        na_rep='',
        compression=compression_opts
    )

    # Compare the lisi distributions
    n_labels = len(labels)
    for lisi_column in lisi_columns:
        # Make density plot.
        gplt = plt9.ggplot(df_lisi, plt9.aes(
            fill='label',
            x='label',
            y=lisi_column,
        ))
        gplt = gplt + plt9.theme_bw(base_size=12)
        gplt = gplt + plt9.geom_violin(alpha=0.9)
        gplt = gplt + plt9.geom_boxplot(
            group='label',
            position=plt9.position_dodge(width=.9),
            width=.1,
            fill='white',
            outlier_alpha=0  # Do not know how to totally remove outliers.
        )
        # Add a line at the theoretical maximum
        gplt = gplt + plt9.geom_hline(plt9.aes(
            yintercept=lisi_limit[lisi_column]
        ))
        # gplt = gplt + plt9.facet_grid('{} ~ .'.format(label))
        gplt = gplt + plt9.labs(
            x='Reduced dimensions',
            y='LISI',
            title=''
        )
        gplt = gplt + plt9.theme(
            axis_text_x=plt9.element_text(angle=-45, hjust=0)
        )
        gplt = gplt + plt9.theme(
            legend_position='none'
        )
        if n_labels != 0 and n_labels < 9:
            gplt = gplt + plt9.scale_fill_brewer(
                palette='Dark2',
                type='qual'
            )
        gplt.save(
            '{}-{}-violin.png'.format(out_file_base, lisi_column),
            #dpi=300,
            width=4*(n_labels/4),
            height=10,
            # height=4*(n_samples/4),
            limitsize=False
        )

        # Make ecdf.
        gplt = plt9.ggplot(df_lisi, plt9.aes(
            x=lisi_column,
            color='label',
        ))
        gplt = gplt + plt9.theme_bw(base_size=12)
        gplt = gplt + plt9.stat_ecdf(alpha=0.8)
        gplt = gplt + plt9.labs(
            x='LISI',
            y='Cumulative density',
            # color='Reduction',
            title=''
        )
        if n_labels != 0 and n_labels < 9:
            gplt = gplt + plt9.scale_color_brewer(
                palette='Dark2',
                type='qual'
            )
        gplt.save(
            '{}-{}-ecdf.pdf'.format(out_file_base, lisi_column),
            #dpi=300,
            width=10,
            height=4,
            limitsize=False
        )


if __name__ == '__main__':
    main()
