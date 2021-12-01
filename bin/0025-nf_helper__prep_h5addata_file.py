#!/usr/bin/env python


__author__ = 'Guillaume Noell'
__date__ = '2020-04-03'
__version__ = '0.0.1'

import argparse
import os
import pandas as pd


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Filter and merge h5ad data. Save to AnnData object.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '--h5ad_list',
        action='store',
        dest='h5ad_list',
        required=True,
        help='List of h5ad files in the following format: \
            <sample_id1>---h5ad.h5ad,<sample_id2>---h5ad.h5ad'
    )

    parser.add_argument(
        '-txd', '--h5addata_file',
        action='store',
        dest='txd',
        required=False,
        help='File with the following headers: experiment_id\
            data_path_h5ad_format.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='nf_prepped__file_paths_h5ad.tsv',
        help='Basename of output anndata file, assuming output in current \
            working directory. Will have .h5ad appended.\
            (default: %(default)s)'
    )

    options = parser.parse_args()
    print('\n0025-nf_helper__prep_h5addata_file.py options: ' + str(options) + '\n')

    # Get all of the input file lists
    h5ad = options.h5ad_list.split(',')

    input_dict = {}
    for file in h5ad:
        experiment_id = file.replace('---h5ad.h5ad', '')
        if experiment_id not in input_dict:
            input_dict[experiment_id] = {}
        input_dict[experiment_id]['h5ad'] = file

    # Check to make sure we have an entry in the dict for every experiment_id
    # in the original txd file
    if options.txd is not None:
        df = pd.read_csv(options.txd, sep='\t')
        for i in df['experiment_id']:
            if i not in input_dict:
                raise Exception('Missing experiment_id:\t{}'.format(
                    i
                ))

    output_dict = {}
    output_dict['experiment_id'] = []
    output_dict['data_path_h5ad_format'] = []
    for experiment_id, value in input_dict.items():
        os.makedirs(
            'h5ad_input_data/{}'.format(experiment_id),
            exist_ok=True
        )
        output_dict['experiment_id'].append(experiment_id)
        output_dict['data_path_h5ad_format'].append(
            os.path.abspath('h5ad_input_data/{}'.format(experiment_id))
        )
        # os.symlink(
        #     value['h5ad'],
        #     'h5ad_input_data/{}/h5ad.tsv.gz'.format(experiment_id)
        # )
        f = 'h5ad_input_data/{}/h5ad.h5ad'.format(experiment_id)
        if not os.path.exists(f):
            os.link(value['h5ad'], f)

    output_df = pd.DataFrame(output_dict)
    output_df = output_df.sort_values(
        by=['experiment_id'],
        ascending=[True]
    )
    output_df.to_csv(options.of, sep='\t')


if __name__ == '__main__':
    main()
