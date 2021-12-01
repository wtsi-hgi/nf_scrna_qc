#!/usr/bin/env python

__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
import pandas as pd


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
        '--barcodes_list',
        action='store',
        dest='barcodes_list',
        required=True,
        help='List of barcode files in the following format: \
            <sample_id1>---barcodes.tsv.gz,<sample_id2>---barcodes.tsv.gz'
    )

    parser.add_argument(
        '--features_list',
        action='store',
        dest='features_list',
        required=True,
        help='List of features files in the following format: \
            <sample_id1>---features.tsv.gz,<sample_id2>---features.tsv.gz'
    )

    parser.add_argument(
        '--matrix_list',
        action='store',
        dest='matrix_list',
        required=True,
        help='List of matrix files in the following format: \
            <sample_id1>---matrix.mtx.gz,<sample_id2>---matrix.mtx.gz'
    )

    parser.add_argument(
        '-txd', '--tenxdata_file',
        action='store',
        dest='txd',
        required=False,
        help='File with the following headers: experiment_id\
            data_path_10x_format.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='nf_prepped__file_paths_10x.tsv',
        help='Basename of output anndata file, assuming output in current \
            working directory. Will have .h5ad appended.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Get all of the input file lists
    barcodes = options.barcodes_list.split(',')
    features = options.features_list.split(',')
    matricies = options.matrix_list.split(',')

    input_dict = {}
    for file in barcodes:
        experiment_id = file.replace('---barcodes.tsv.gz', '')
        if experiment_id not in input_dict:
            input_dict[experiment_id] = {}
        input_dict[experiment_id]['barcode'] = file
    for file in features:
        experiment_id = file.replace('---features.tsv.gz', '')
        if experiment_id not in input_dict:
            input_dict[experiment_id] = {}
        input_dict[experiment_id]['feature'] = file
    for file in matricies:
        experiment_id = file.replace('---matrix.mtx.gz', '')
        if experiment_id not in input_dict:
            input_dict[experiment_id] = {}
        input_dict[experiment_id]['matrix'] = file

    # Check to make sure we have an entry in the dict for every experiment_id
    # in the original txd file
    if options.txd != '':
        df = pd.read_csv(options.txd, sep='\t')
        for i in df['experiment_id']:
            if i not in input_dict:
                raise Exception('Missing experiment_id:\t{}'.format(
                    i
                ))

    output_dict = {}
    output_dict['experiment_id'] = []
    output_dict['data_path_10x_format'] = []
    for experiment_id, value in input_dict.items():
        os.makedirs(
            'tenx_input_data/{}'.format(experiment_id),
            exist_ok=True
        )
        output_dict['experiment_id'].append(experiment_id)
        output_dict['data_path_10x_format'].append(
            os.path.abspath('tenx_input_data/{}'.format(experiment_id))
        )
        # os.symlink(
        #     value['barcode'],
        #     'tenx_input_data/{}/barcode.tsv.gz'.format(experiment_id)
        # )
        f = 'tenx_input_data/{}/barcodes.tsv.gz'.format(experiment_id)
        if not os.path.exists(f):
            os.link(value['barcode'], f)
        f = 'tenx_input_data/{}/features.tsv.gz'.format(experiment_id)
        if not os.path.exists(f):
            os.link(value['feature'], f)
        f = 'tenx_input_data/{}/matrix.mtx.gz'.format(experiment_id)
        if not os.path.exists(f):
            os.link(value['matrix'], f)

    output_df = pd.DataFrame(output_dict)
    output_df = output_df.sort_values(
        by=['experiment_id'],
        ascending=[True]
    )
    output_df.to_csv(options.of, sep='\t')


if __name__ == '__main__':
    main()
