#!/usr/bin/env python3

## assess clustering using SCCAF metric
DEBUG = False

import sys
import os
import argparse
import pathlib
import pandas
import scanpy
import SCCAF
import matplotlib.pyplot as plt

OSB_COLNAM_OPTIM_PRFX = "SCCAF"
OBS_COLNAM_OPTIM_FMT = "{}_Round{:d}"

def optimize_clustering(anndat, resolution=float(4.0), n_repeats = 1, n_jobs = 1, n_cells_max = 100, plotfil_prfx = None):


def set_argument_parser():
    parser = argparse.ArgumentParser(description="Assess clustering using SCCAF")
    parser.add_argument("fnam_anndat", metavar = "anndat_fnam",
        help="Input file for AnnData object [h5ad]")
    parser.add_argument(
        "--cluster-file",
        help="Input cluster table [TSV]",
        dest="fnam_clusters")
    parser.add_argument("--output-prefix", "-p", metavar="output_prefix", dest = "out_prfx",
        default="sccaf_output",
        help="Prefix for output files (default: sccaf_output)")
    parser.add_argument("-r", "--repeats", metavar="n_repeats", dest = "n_repeats",
                        help="Number of times to iterate the assessment to build distributions of accuracies", type=int,
                        default=1)
    parser.add_argument("-n", "--nthreads", metavar = "n_threads", dest="n_threads",
                        help="Number of processors to use", type=int, default=1)
    parser.add_argument("--use-pca",
                        help="Use PCA components for assessment (needs to be available in the ann data ) (default: False)",
                        action='store_true')
    parser.add_argument(
        "-c", "--cluster-resolution",
        help="""Resolution parameter for initial Leiden clustering (only used in connection with --optimize).
            Only used in connection with --optimize and when --cluster-file is not specified.
            """,
        type = float,
        default = float(4.0)
    )

    return parser.parse_args()

if __name__ == '__main__':
    args = set_argument_parser()
    sys.stderr.write("# reading {:s}\n".format(args.fnam_anndat))
    anndat = scanpy.read(args.fnam_anndat, cache = True)
    colnam_start = OBS_COLNAM_OPTIM_FMT.format(OBS_COLNAM_OPTIM_PRFX, 0)
    if args.fnam_clusters:
        clustdat = pandas.read_table(args.fnam_clusters, usecols=[0,1], index_col=0, header=0)
        anndat.obs[colnam_start] = pandas.Index.to_series(anndat.obs.index).align(clustdat['cluster'], join = 'left')[1]
    else:
        scanpy.tl.leiden(anndat, resolution=args.cluster_resolution,
            add_key = colnam_start)

    optimize_clustering(anndat, n_repeats = args.n_repeats, n_jobs = args.n_threads, plotfil_prfx = args.out_prfx)

    sys.exit()
