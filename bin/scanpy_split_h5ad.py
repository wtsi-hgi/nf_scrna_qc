#!/usr/bin/env python3

import sys
import os
from distutils.version import LooseVersion
import scipy
import scipy.io
import gzip
import pandas
import scanpy
## split h5ad file by chromium channel

# for testing use
# infnam = "/lustre/scratch123/hgi/mdt1/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/franke_data/work/2a/ebdc5079ae949777263ce3b1aca510/BF61CE54F4603C9F-adata.h5ad"

def split_h5ad_by_batch(ad, oufnprfx, colnam_batch = 'experiment_id', anndata_compression_opts=None):
    oufn_list_fnam = '{}_files.txt'.format(oufnprfx)
    oufn_list = []
    batch_labels = pandas.Categorical(ad.obs[colnam_batch].apply(lambda a: a.split('__')[0])) # <class 'pandas.core.series.Series'>
    for bl in batch_labels.categories:
        oufnam = '{0}_{1}.h5ad'.format(oufnprfx, bl)
        oufn_list.append(oufnam)
        adb = ad[batch_labels == bl,:]
        # strip unneccessary meta data - this seems to aid subsequent transformation to Seurat h5 format
        # assumes that cells failing qc already were stripped out
        # ad.obs = ad.obs[['convoluted_samplename', 'cell_passes_qc']]
        # ad.var = ad.var[['feature_types', 'genome', 'gene_symbols']]
        adb.obs = pandas.DataFrame(adb.obs.index, index = adb.obs.index, columns = ["cell_barcode"])

        vdf = adb.var[["feature_types"]]
        vdf.insert(1,"gene_ids", vdf.index)
        vdf.index = pandas.Index(adb.var['gene_symbols'].astype('str'))
        #ad.var = vdf.set_index("gene_symbols", drop = True, verify_integrity = False)
        adb.var = vdf

        del adb.uns

        if anndata_compression_opts is None:
            adb.write(oufnam)
        else:
            adb.write(
                oufnam,
                compression='gzip',
                compression_opts=anndata_compression_opts
            )
        oufn_list.append(oufnam)
    with open(oufn_list_fnam, 'w') as oufh:
        for fn in oufn_list:
            oufh.write(fn + '\n')
    return oufn_list_fnam

if __name__ == '__main__':
    nargs = len(sys.argv)
    if nargs != 3:
        sys.exit("usage: %s ./<input_file *.h5ad> <output_file_prefix>".format(sys.argv[0]))

    infnam = sys.argv[1]
    oufnprfx = sys.argv[2]

    # print(infnam)
    ad = scanpy.read(infnam)
    split_h5ad_by_batch(ad, oufnprfx, colnam_batch = 'experiment_id', anndata_compression_opts = None)
    sys.exit()
