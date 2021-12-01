#!/usr/bin/env python

__author__ = 'Guillaume Noell'
__date__ = '2021-07-08'
__version__ = '0.0.1'

import argparse
import os
from distutils.version import LooseVersion
from typing import Dict
import tables
import scipy
import scipy.io
import scipy.sparse
import gzip
import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import click
import logging

from plotnine import ggplot, geom_vline, geom_point, aes, labs, stat_smooth, facet_wrap, geom_histogram, geom_boxplot, theme, theme_bw, element_text, scale_x_continuous, scale_y_continuous, element_blank
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec

# To resolve strange TclError for interactive job
import matplotlib
matplotlib.use('Agg')  # Agg for png and pdf for pdf

@click.command()
@click.option('--samplename', required=True, type=click.STRING, help='')
@click.option('--raw_cellranger_mtx', required=True, type=click.STRING, help='') # type=click.Path(exists=False, resolve_path=True, dir_okay=True, file_okay=False), help='')
@click.option('--filtered_cellranger_mtx', required=True, type=click.STRING, help='') #,type=click.Path(exists=False, resolve_path=True, dir_okay=True, file_okay=False), help='')
@click.option('--cellbender_unfiltered_h5', required=True,type=click.STRING, help='')
@click.option('--fpr', required=True,type=click.STRING, help='')
@click.option('--n_expected_cells', required=True,type=click.INT, help='')
@click.option('--n_total_droplets_included', required=True,type=click.INT, help='')
@click.option('--out_dir', required=True,type=click.STRING, help='')

def plot_cellranger_vs_cellbender(samplename, raw_cellranger_mtx, filtered_cellranger_mtx,
                                  cellbender_unfiltered_h5, fpr, n_expected_cells,
                                  n_total_droplets_included, out_dir):
    """compare cellranger raw vs cellranger filtered vs cellbender outputs"""
    logging.info('samplename ' + str(samplename))
    logging.info('raw_cellranger_mtx ' + str(raw_cellranger_mtx))
    logging.info('filtered_cellranger_mtx ' + str(filtered_cellranger_mtx))
    logging.info('cellbender_unfiltered_h5 ' + str(cellbender_unfiltered_h5))
    logging.info('fpr ' + str(fpr))
    logging.info('n_expected_cells ' + str(n_expected_cells))
    logging.info('n_total_droplets_included ' + str(n_total_droplets_included))
    logging.info('out_dir ' + str(out_dir))

    # Make the output directory if it does not exist.
    if out_dir == '':
        out_dir = os.getcwd()
    else:
        os.makedirs(out_dir, exist_ok=True)
        out_dir = out_dir + '/fpr_' + fpr
        os.makedirs(out_dir, exist_ok=True)
        os.makedirs(out_dir + '/' + samplename, exist_ok=True)
        logging.info(out_dir)
    # logging.info(df.head())
    
    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)
       
    # read cellranger raw
    adata_cellranger_raw = sc.read_10x_mtx(
        raw_cellranger_mtx, var_names='gene_symbols', make_unique=True,
        cache=False, cache_compression=compression_opts)

    # First filter out any cells that have 0 total counts
    zero_count_cells_cellranger_raw = adata_cellranger_raw.obs_names[np.where(adata_cellranger_raw.X.sum(axis=1) == 0)[0]]
    # sc.pp.filter_cells(adata, min_counts=1, inplace=True) # Minimum number of counts required for a cell to pass filtering.
    logging.info("_cellranger_raw: Filtering {}/{} cells with 0 counts.".format(
        len(zero_count_cells_cellranger_raw),
        adata_cellranger_raw.n_obs))
    adata_cellranger_raw = adata_cellranger_raw[adata_cellranger_raw.obs_names.difference(zero_count_cells_cellranger_raw, sort=False)]
    
    sc.pp.calculate_qc_metrics(adata_cellranger_raw, inplace=True)

    logging.info('cellranger raw n barcodes(.obs) x cells(.var) .X.shape:'); logging.info(adata_cellranger_raw.X.shape)
    logging.info('cellranger raw .obs:'); logging.info(adata_cellranger_raw.obs)
    logging.info('cellranger raw .var:'); logging.info(adata_cellranger_raw.var)
    
    df_total_counts = pd.DataFrame(data= adata_cellranger_raw.obs.sort_values(by=['total_counts'], ascending=False).total_counts)
    df_total_counts['barcode_row_number'] = df_total_counts.reset_index().index + 1 
    df_total_counts['barcodes'] = df_total_counts.index
    df_total_counts_cellranger_raw = df_total_counts
    df_total_counts_cellranger_raw['dataset']='Cellranger Raw'
    
    logging.info(df_total_counts)
    # read cellranger filtered
    adata_cellranger_filtered = sc.read_10x_mtx(
        filtered_cellranger_mtx, var_names='gene_symbols', make_unique=True,
        cache=False, cache_compression=compression_opts)

    # First filter out any cells that have 0 total counts
    zero_count_cells_cellranger_filtered = adata_cellranger_filtered.obs_names[np.where(adata_cellranger_filtered.X.sum(axis=1) == 0)[0]]
    # sc.pp.filter_cells(adata, min_counts=1, inplace=True) # Minimum number of counts required for a cell to pass filtering.
    logging.info("_cellranger_filtered: Filtering {}/{} cells with 0 counts.".format(
        len(zero_count_cells_cellranger_filtered),
        adata_cellranger_filtered.n_obs))
    adata_cellranger_filtered = adata_cellranger_filtered[adata_cellranger_filtered.obs_names.difference(zero_count_cells_cellranger_filtered, sort=False)]
    
    sc.pp.calculate_qc_metrics(adata_cellranger_filtered, inplace=True)
    
    logging.info('cellranger filtered n barcodes(.obs) x cells(.var) .X.shape:'); logging.info(adata_cellranger_filtered.X.shape)
    logging.info('cellranger filtered .obs:'); logging.info(adata_cellranger_filtered.obs.columns)
    logging.info(adata_cellranger_filtered.obs)
    logging.info('cellranger filtered .var:'); logging.info(adata_cellranger_filtered.var)

    df_total_counts = pd.DataFrame(data= adata_cellranger_filtered.obs.sort_values(by=['total_counts'], ascending=False).total_counts)
    df_total_counts['barcodes'] = df_total_counts.index
    df_total_counts['barcode_row_number'] = df_total_counts.reset_index().index + 1 
    df_total_counts_cellranger_filtered = df_total_counts
    df_total_counts_cellranger_filtered['dataset'] = 'Cellranger Filtered'

    logging.info(df_total_counts)
    # read cellbender output
    adata_cellbender = anndata_from_h5(cellbender_unfiltered_h5,
                                          analyzed_barcodes_only=True)

    # First filter out any cells that have 0 total counts
    zero_count_cells_cellbender_filtered = adata_cellbender.obs_names[np.where(adata_cellbender.X.sum(axis=1) == 0)[0]]
    # sc.pp.filter_cells(adata, min_counts=1, inplace=True) # Minimum number of counts required for a cell to pass filtering.
    logging.info("_cellbender_filtered: Filtering {}/{} cells with 0 counts.".format(
        len(zero_count_cells_cellbender_filtered),
        adata_cellbender.n_obs))
    adata_cellbender = adata_cellbender[adata_cellbender.obs_names.difference(zero_count_cells_cellbender_filtered, sort=False)]
    
    sc.pp.calculate_qc_metrics(adata_cellbender, inplace=True)

    logging.info('cellbender cellbender.n barcodes(.obs) x cells(.var) .X.shape:'); logging.info(adata_cellbender.X.shape)
    logging.info('cellbender cellbender.obs:'); logging.info(adata_cellbender.obs)
    logging.info('cellbender cellbender.var:'); logging.info(adata_cellbender.var)

    df_total_counts = pd.DataFrame(data= adata_cellbender.obs.sort_values(by=['total_counts'], ascending=False).total_counts)
    df_total_counts['barcodes'] = df_total_counts.index
    df_total_counts['barcode_row_number'] = df_total_counts.reset_index().index + 1 
    df_total_counts_cellbender = df_total_counts
    df_total_counts_cellbender['dataset'] = 'Cellbender'
    
    logging.info(df_total_counts)
    
    # df_total_counts_cellranger_filtered.rename(columns={"total_counts": "cellranger_filtered_total_counts"})
    df_cellranger_cellbender = pd.merge(df_total_counts_cellranger_filtered, df_total_counts_cellbender, 
                                        how='outer', left_index=True, right_index=True,
                                        suffixes=('_cellranger', '_cellbender')).sort_values(by=['total_counts_cellbender'],
                                                                                             ascending=False)
    logging.info(df_cellranger_cellbender)
    df_cellranger_cellbender[['cellranger', 'cellbender']] = np.where(df_cellranger_cellbender[['total_counts_cellranger', 'total_counts_cellbender']].isnull(), 0, 1)

    #df_cellranger_cellbender.to_csv('df_cellranger_cellbender.csv', index=True, index_label='barcode')

    grouped = df_cellranger_cellbender[['cellranger', 'cellbender']].groupby(["cellranger", "cellbender"]).size().reset_index(name='counts')
    logging.info(grouped.columns)
    #grouped.to_csv('cellranger_cellbender.csv', index=False)
    
    df_cellranger_cellbender['barcode_row_number'] = df_cellranger_cellbender.reset_index().index + 1 


    ### plot UMI counts descending order 
    df_merged = pd.concat([df_total_counts_cellranger_raw, df_total_counts_cellranger_filtered, df_total_counts_cellbender])
    #df_merged.to_csv('df_merged.csv', index=True, index_label='barcode') 

    df_vline =  pd.DataFrame(data={'x': [int(n_expected_cells), int(n_total_droplets_included)],
                                   'color': ['expected-cells', 'total-droplets-included']})
    
    gplt = ggplot(df_merged, aes(x='barcode_row_number', y='total_counts')) \
        + geom_point() \
        + geom_vline(df_vline, aes(xintercept='x', color='color')) \
        + theme_bw() + facet_wrap('dataset') \
        + labs(x='Barcodes (ordered by descending cell total couts)',color='Cellbender input',
               y='Cell total counts', title='Cells filtered out by Cellranger or Cellbender') \
        + scale_y_continuous(trans='log10',minor_breaks=0) + scale_x_continuous(trans='log10',minor_breaks=0)
    gplt.save(out_dir + '/' + samplename + '/barcode_vs_total_counts.png', width=12,height=5, dpi=300) # dpi=300,

    df_cellranger_cellbender_count = grouped # pd.read_csv('cellranger_cellbender.csv')
    
    df = pd.merge(df_merged, df_cellranger_cellbender[['cellranger','cellbender']], how='left', left_index=True, right_index=True)
    df = pd.merge(df, df_cellranger_cellbender_count, how='left', left_on=['cellranger','cellbender'], right_on=['cellranger','cellbender'])
    df["counts"].fillna(df['counts'].isnull().sum(), inplace=True)
    df["counts"] = df["counts"].astype(int)
    # df.replace({"counts": {""}  }, inplace=True)
    
    df["filtered"] = df["cellranger"].astype(str) + '-' +  df["cellbender"].astype(str)
    df.replace({"filtered": {"nan-nan":'Cellranger Raw only', "1.0-1.0":"Cellranger Filtered + Cellbender", "1.0-0.0":"Cellranger Filtered only",  "0.0-1.0":"Cellbender only",  "0.0-0.0":"0.0-0.0"}  }, inplace=True)
    df["filtered"] = df["filtered"] + ', n=' + df["counts"].astype(str)
    df['filtered'].value_counts()
    df.replace({"dataset": {"cellbender":"Cellbender output","cellranger_raw":"Cellranger Raw output","cellranger_filtered":"Cellranger Filtered output"}  }, inplace=True)
    
    
    gplt = ggplot(df, aes(x='filtered', y='total_counts', color='filtered')) \
        + geom_boxplot() \
        + theme_bw() \
        + facet_wrap('dataset') \
        + theme(axis_text_x=element_blank()) \
        + scale_y_continuous(trans='log10',minor_breaks=0) \
        + labs(color='n cells in intersection of datasets', x='', y='Cell total counts', title='Total cell counts compared across datasets (facets)')
    gplt.save(out_dir + '/' + samplename + '/boxplots_cellranger_vs_cellbender.png', width=12,height=5, dpi=300) # dpi=300,


    # plot difference cellbender filtered vs cellranger filtered for common cells between the 2 datasets
    df_cellranger_cellbender = df_cellranger_cellbender[df_cellranger_cellbender['cellranger'] == 1]
    df_cellranger_cellbender = df_cellranger_cellbender[df_cellranger_cellbender['cellbender'] == 1]

    # Subset the datasets to the relevant barcodes.
    adata_cellbender_common = adata_cellbender[df_cellranger_cellbender.index.values]
    adata_cellranger_filtered_common = adata_cellranger_filtered[df_cellranger_cellbender.index.values]
    # Put count matrices into 'layers' in anndata for clarity.
    adata = adata_cellbender_common
    adata.layers['counts_cellbender'] = adata_cellbender_common.X.copy()
    adata.layers['counts_raw'] = adata_cellranger_filtered_common.X.copy()
    # Get the differences in counts per cell
    X_raw_minus_cb = adata.layers['counts_raw'] - adata.layers['counts_cellbender']
    X_raw_minus_cb = abs(X_raw_minus_cb)
    # Get the top most different genes
    df_diff_genes = pd.DataFrame(data=adata.var.gene_symbols.values)
    df_diff_genes['ensembl_id'] = adata.var.index
    df_diff_genes['gene_symbol'] = adata.var.gene_symbols.values
    df_diff_genes['dif_across_cells'] = np.asarray(X_raw_minus_cb.sum(axis=0)).reshape(-1)
    df_diff_genes = df_diff_genes.sort_values('dif_across_cells',ascending=False).head(n=100)
    #df_diff_genes.to_csv('df_diff_genes.csv', index=True)
    top_genes = df_diff_genes['ensembl_id']
    top_genes_symbols = df_diff_genes['gene_symbol']
    logging.info('top_genes:')
    logging.info(top_genes)

    logging.info(adata_cellbender_common.var.index)
    adata_cellbender_common = adata_cellbender[df_cellranger_cellbender.index.values, top_genes].to_df()
    adata_cellbender_common['barcode'] = adata_cellbender_common.index
    adata_cellbender_common = pd.melt(adata_cellbender_common, ignore_index = True, id_vars=['barcode'],
                                      var_name='ensembl_id', value_name='count')
    adata_cellbender_common = pd.merge(adata_cellbender_common, df_diff_genes[['ensembl_id','gene_symbol']],
                                        how='left', left_on='ensembl_id', right_on='ensembl_id')
    adata_cellbender_common = adata_cellbender_common.sort_values(by=['barcode', 'ensembl_id'], ascending=False)
    adata_cellbender_common['dataset'] = 'Cellbender'
    #adata_cellbender_common.to_csv('adata_cellbender_common.csv', index=True)
    
    logging.info(adata_cellranger_filtered.var.index)
    adata_cellranger_filtered_common = adata_cellranger_filtered[df_cellranger_cellbender.index.values,
                                                                 top_genes_symbols].to_df()
    adata_cellranger_filtered_common['barcode'] = adata_cellranger_filtered_common.index
    adata_cellranger_filtered_common = pd.melt(adata_cellranger_filtered_common, ignore_index = True,
                                               id_vars=['barcode'], var_name='gene_symbol', value_name='count')
    adata_cellranger_filtered_common = pd.merge(adata_cellranger_filtered_common, df_diff_genes[['ensembl_id','gene_symbol']],
                                        how='left', left_on='gene_symbol', right_on='gene_symbol')
    adata_cellranger_filtered_common['dataset'] = 'Cellranger Filtered'
    adata_cellranger_filtered_common = adata_cellranger_filtered_common.sort_values(by=['barcode', 'ensembl_id'], ascending=False)
    adata_cellranger_filtered_common = adata_cellranger_filtered_common[adata_cellbender_common.columns]
    #adata_cellranger_filtered_common.to_csv('adata_cellranger_filtered_common.csv', index=True)
    
    logging.info(adata_cellranger_raw.var.index)
    adata_cellranger_raw_common = adata_cellranger_raw[df_cellranger_cellbender.index.values,
                                                                 top_genes_symbols].to_df()
    adata_cellranger_raw_common['barcode'] = adata_cellranger_raw_common.index
    adata_cellranger_raw_common = pd.melt(adata_cellranger_raw_common, ignore_index = True,
                                               id_vars=['barcode'], var_name='gene_symbol', value_name='count')
    adata_cellranger_raw_common = pd.merge(adata_cellranger_raw_common, df_diff_genes[['ensembl_id','gene_symbol']],
                                        how='left', left_on='gene_symbol', right_on='gene_symbol')
    adata_cellranger_raw_common['dataset'] = 'Cellranger Raw'
    adata_cellranger_raw_common = adata_cellranger_raw_common.sort_values(by=['barcode', 'ensembl_id'], ascending=False)
    adata_cellranger_raw_common = adata_cellranger_raw_common[adata_cellbender_common.columns]
    #adata_cellranger_raw_common.to_csv('adata_cellranger_raw_common.csv', index=True)

    logging.info(adata_cellranger_raw_common['gene_symbol']== adata_cellbender_common['gene_symbol'])
    logging.info(adata_cellranger_raw_common['ensembl_id']== adata_cellbender_common['ensembl_id'])


    adata_filtered_cellbender_diff = adata_cellbender_common.copy()
    adata_filtered_cellbender_diff['count'] = adata_cellranger_filtered_common['count'] - adata_cellbender_common['count']
    adata_filtered_cellbender_diff['dataset'] = 'Cellranger Filtered - Cellbender'

    adata_raw_cellbender_diff = adata_cellbender_common.copy()
    adata_raw_cellbender_diff['count'] = adata_cellranger_raw_common['count'] - adata_cellbender_common['count']
    adata_raw_cellbender_diff['dataset'] = 'Cellranger Raw - Cellbender'
    
    df_merged = pd.concat([adata_cellbender_common, adata_cellranger_filtered_common,
                           adata_cellranger_raw_common,
                           adata_filtered_cellbender_diff, adata_raw_cellbender_diff], ignore_index=True)
    
    gplt = ggplot(df_merged, aes(x='gene_symbol',y='count')) \
        + geom_boxplot() \
        + theme_bw() \
        + theme(axis_text_x = element_text(angle = 90, hjust = 1, size= 6)) \
        + facet_wrap('dataset', scales = 'free', ncol = 1) \
        + labs(x='Genes (top 100 Genes most different between Cellranger Filtered counts and Cellbender filtered counts)', y='Cell total counts', title='Total cell counts compared across most different genes (x-axis) and datasets (facets)')
    gplt.save(out_dir + '/' + samplename + '/boxplot_topgenes_cellranger_vs_cellbender.png', width=10,height=20, dpi=300) # dpi=300,
    logging.info('script done.')



    
def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary."""
    d = {}
    with tables.open_file(file) as f:
        # read in everything
        for array in f.walk_nodes("/", "Array"):
            d[array.name] = array.read()
    return d

def anndata_from_h5(
    file: str,
    analyzed_barcodes_only: bool = True
 ) -> 'anndata.AnnData':
    """Load an output h5 file into an AnnData object for downstream work.

    Args:
        file: The h5 file
        analyzed_barcodes_only: False to load all barcodes, so that the size of
            the AnnData object will match the size of the input raw count
            matrix. True to load a limited set of barcodes: only those
            analyzed by the algorithm. This allows relevant latent
            variables to be loaded properly into adata.obs and adata.obsm,
            rather than adata.uns.

    Returns:
        adata: The anndata object, populated with inferred latent variables
            and metadata.
    """
    d = dict_from_h5(file)
    X = scipy.sparse.csc_matrix(
        (d.pop('data'), d.pop('indices'), d.pop('indptr')),
        shape=d.pop('shape')
    ).transpose().tocsr()

    if analyzed_barcodes_only:
        if 'barcodes_analyzed_inds' in d.keys():
            X = X[d['barcodes_analyzed_inds'], :]
            d['barcodes'] = d['barcodes'][d['barcodes_analyzed_inds']]
        elif 'barcode_indices_for_latents' in d.keys():
            X = X[d['barcode_indices_for_latents'], :]
            d['barcodes'] = d['barcodes'][d['barcode_indices_for_latents']]
        else:
            print(
                'Warning: analyzed_barcodes_only=True, but the key ',
                '"barcodes_analyzed_inds" or "barcode_indices_for_latents" ',
                'is missing from the h5 file. ',
                'Will output all barcodes, and proceed as if ',
                'analyzed_barcodes_only=False'
            )

    print(d.keys())

    # Construct the count matrix.
    if 'gene_names' in d.keys():
        gene_symbols = d.pop('gene_names').astype(str)
    else:
        gene_symbols = d.pop('name').astype(str)
    adata = anndata.AnnData(
        X=X,
        obs={'barcode': d.pop('barcodes').astype(str)},
        var={
            'gene_ids': d.pop('id').astype(str),
            'gene_symbols': gene_symbols
        }
    )
    adata.obs.set_index('barcode', inplace=True)
    adata.var.set_index('gene_ids', inplace=True)

    # Add other information to the adata object in the appropriate slot.
    for key, value in d.items():
        try:
            value = np.asarray(value)
            if len(value.shape) == 0:
                adata.uns[key] = value
            elif value.shape[0] == X.shape[0]:
                if (len(value.shape) < 2) or (value.shape[1] < 2):
                    adata.obs[key] = value
                else:
                    adata.obsm[key] = value
            elif value.shape[0] == X.shape[1]:
                if value.dtype.name.startswith('bytes'):
                    adata.var[key] = value.astype(str)
                else:
                    adata.var[key] = value
            else:
                adata.uns[key] = value
        except Exception:
            print(
                'Unable to load data into AnnData: ', key, value, type(value)
            )

    if analyzed_barcodes_only:
        cols = adata.obs.columns[
            adata.obs.columns.str.startswith('barcodes_analyzed')
            | adata.obs.columns.str.startswith('barcode_indices')
        ]
        for col in cols:
            try:
                del adata.obs[col]
            except Exception:
                pass

    return adata

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO) # set logging level
    plot_cellranger_vs_cellbender()
