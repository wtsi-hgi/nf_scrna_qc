#!/usr/bin/env python

__date__ = '2020-05-26'
__version__ = '0.0.1'

import argparse
import os
import csv
import numpy as np
import pandas as pd
import scanpy as sc
from kneed import KneeLocator

import matplotlib.pyplot as plt
# from matplotlib import rcParams
# from matplotlib import gridspec

# Paper:
# https://raghavan.usc.edu//papers/kneedle-simplex11.pdf


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Calculates the knee or elbow of PC plot. Uses method described in
            https://doi.org/10.1109/ICDCSW.2011.20.
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
        '--add_n_pcs_to_elbow',
        action='store',
        dest='add_n_pcs_to_elbow',
        default=0,
        type=int,
        help='This value is added to the elbow estimate.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: <metadata_tsv>-knee)'
    )

    options = parser.parse_args()

    # Fixed settings.
    verbose = False

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}-knee'.format(
            os.path.basename(options.h5.rstrip('h5ad').rstrip('.'))
        )

    # Read in the dataframe
    adata = sc.read_h5ad(filename=options.h5)

    kneedle_dict = {}
    output_dict = {}

    # Get kneedle estimate using raw data OR fit smoothing spline. Also
    # maximize sensitivity.
    # NOTE: don't fit polynomial makes no sense for our data
    # for fit in [None, 'interp1d', 'polynomial']:
    for fit in [None, 'interp1d']:
        for i in ['variance_ratio', 'variance']:
            key = '{}-spline={}'.format(i, fit)
            sensitivity = 1000
            while True:
                if fit:
                    kneedle_dict[key] = KneeLocator(
                        np.array(list(range(0, len(adata.uns['pca'][i]))))+1,
                        adata.uns['pca'][i],
                        curve='convex',
                        direction='decreasing',
                        S=sensitivity,
                        interp_method=fit
                    )
                else:
                    kneedle_dict[key] = KneeLocator(
                        np.array(list(range(0, len(adata.uns['pca'][i]))))+1,
                        adata.uns['pca'][i],
                        curve='convex',
                        direction='decreasing',
                        S=sensitivity
                    )
                if kneedle_dict[key].knee is None:
                    sensitivity -= 5
                else:
                    if verbose:
                        print('{}\nS:\t{}\nknee:\t{}\nelbow:\t{}'.format(
                            i,
                            sensitivity,
                            round(kneedle_dict[key].knee, 3),
                            round(kneedle_dict[key].elbow, 3)
                        ))
                    output_dict[key] = {
                        'sensitivity': sensitivity,
                        'knee': kneedle_dict[key].knee,
                        'elbow': kneedle_dict[key].elbow,
                        # 'norm_knee': kneedle_dict[key].norm_knee,
                        # 'norm_elbow': kneedle_dict[key].norm_elbow,
                        'interp_method': fit,
                        'y': i
                    }
                    break

    # Save the results
    output_df = pd.DataFrame(output_dict).transpose().reset_index(drop=True)
    output_df.to_csv(
        '{}.tsv'.format(out_file_base),
        sep='\t',
        index=False,
        quoting=csv.QUOTE_NONNUMERIC,
        # compression='gzip',
        na_rep=''
    )

    # Save the results to pass if the user sets auto estimator of pc knee
    output_auto = output_df.loc[output_df['interp_method'].isna(), :]
    output_auto = output_auto.loc[output_auto['y'] == 'variance', :]
    if options.add_n_pcs_to_elbow:
        if verbose:
            print('Adding {} to elbow estimate'.format(
                options.add_n_pcs_to_elbow
            ))
        output_auto['elbow'] += options.add_n_pcs_to_elbow
    output_auto['elbow'].to_csv(
        '{}-auto_elbow_estimate.tsv'.format(out_file_base),
        sep='\t',
        index=False,
        header=False,
        quoting=csv.QUOTE_NONNUMERIC,
        na_rep=''
    )

    # Plot the results
    for key, kneedle in kneedle_dict.items():

        # Raw data and knee.
        kneedle.plot_knee()
        plt.savefig(
            '{}-{}-knee_raw.png'.format(out_file_base, key),
            #dpi=300,
            bbox_inches='tight'
        )
        plt.close()

        # Normalized data, normalized knee, and normalized distance curve.
        kneedle.plot_knee_normalized()
        plt.savefig(
            '{}-{}-knee_normalized.png'.format(out_file_base, key),
            #dpi=300,
            bbox_inches='tight'
        )
        plt.close()


if __name__ == '__main__':
    main()
