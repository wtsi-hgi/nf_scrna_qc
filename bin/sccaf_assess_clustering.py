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

def assess_clustering(X_m, cluster_v, n_repeats = 1, n_jobs = 1, n_cells_max = 100, plotfil_prfx = None):
    accs =[]
    typeAcc = []
    for i in range(n_repeats):
        sys.stderr.write("# SCCAF cluster assessment iteration {:d} ...\n"
            .format(i+1))
        if DEBUG:
            sys.stderr.flush()

        y_prob, y_pred, y_test, clf, cvsm, acc = SCCAF.SCCAF_assessment(X_m, cluster_v, n_jobs=n_jobs, n=n_cells_max)
        print("type y_prob:", type(y_prob), " shape:", y_prob.shape)
        print("type y_pred: ", type(y_pred), " shape:", y_pred.shape)
        print("type y_test: ", type(y_test), " shape:", y_test.shape)
        print("type clf: ", type(clf))
        accs.append(acc)
        typeAcc.append("Test")
        accs.append(cvsm)
        typeAcc.append("CV")
        if plotfil_prfx and i < 1:
            #dfv = pandas.DataFrame(data={"y_prob": y_prob, "y_pred": y_pred, "y_test": y_test})
            #dfv.to_csv("{:s}_testdat.tsv.gz".format(plotfil_prfx), index=False)
            try:
                SCCAF.plot_roc(y_prob, y_test, clf, cvsm=cvsm, acc=acc, plot = 'roc')
                plt.savefig("{:s}_clust_roc.pdf".format(plotfil_prfx))
                SCCAF.plot_roc(y_prob, y_test, clf, cvsm=cvsm, acc=acc, plot = 'prc')
                plt.savefig("{:s}_clust_prc.pdf".format(plotfil_prfx))
                plt.close()
            except:
                sys.stderr.write("WARNING: no ROC curvers plotted.\n")

    df = pandas.DataFrame(data={"Accuracy": accs, "Type": typeAcc})
    return df

def set_argument_parser():
    parser = argparse.ArgumentParser(description="Assess clustering using SCCAF")
    parser.add_argument("fnam_anndat", metavar = "anndat_fnam",
        help="Input file for AnnData object [h5ad]")
    parser.add_argument("fnam_clusters", metavar = "clust_fnam",
        help="Input cluster table [TSV]")
    parser.add_argument("--output-prefix", "-p", metavar="output_prefix", dest = "out_prfx",
        default="sccaf_output",
        help="Prefix for output files (default: sccaf_output)")
    parser.add_argument("-r", "--repeats", metavar="n_repeats", dest = "n_repeats",
                        help="Number of times to iterate the assesment to build distributions of accuracies", type=int,
                        default=1)
    parser.add_argument("-n", "--nthreads", metavar = "n_threads", dest="n_threads",
                        help="Number of processors to use", type=int, default=1)
    parser.add_argument("--use-pca",
                        help="Use PCA components for assesment (needs to be available in the ann data ) (default: False)",
                        action='store_true')

    return parser.parse_args()

if __name__ == '__main__':
    args = set_argument_parser()
    sys.stderr.write("# reading {:s}\n".format(args.fnam_anndat))
    anndat = scanpy.read(args.fnam_anndat, cache = True)
    clustdat = pandas.read_table(args.fnam_clusters, usecols=[0,1], index_col=0, header=0)
    y = pandas.Index.to_series(anndat.obs.index).align(clustdat['cluster'], join = 'left')[1]
    X = None
    if args.use_pca:
        try:
            X = anndat.obsm['X_pca']
        except KeyError:
            X = None
            sys.stderr.write("Warning: no PCA set falling back to default.\n")
    if X is None:
        raw = getattr(anndat, 'raw')
        if raw:
            X = raw.X
        else:
            X = anndat.X
    df = assess_clustering(X, y, n_repeats = args.n_repeats, n_jobs = args.n_threads, plotfil_prfx = args.out_prfx)
    df.to_csv("{:s}_acc.txt".format(args.out_prfx), index=False)

    sys.exit()
