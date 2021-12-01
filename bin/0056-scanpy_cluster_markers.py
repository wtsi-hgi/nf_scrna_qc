#!/usr/bin/env python

__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
from distutils.version import LooseVersion
import csv
import random
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
# os.environ['PYTHONHASHSEED']=str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)


def nd(arr):
    """Make nd array."""
    return np.asarray(arr).reshape(-1)


def dotplot(adata, markers):
    """Customize dotplot function."""
    adata_raw = adata.raw.to_adata()
    cluter_ids = np.unique(
        adata_raw.obs.cluster.values.astype(int)
    ).astype(str)

    features = adata_raw.var.gene_symbols.values
    midx = [np.where(i == features)[0][0] for i in markers]

    # for each cluster for each gene get two things
    # 1 percent of cells in the cluster expressing that gene
    # 2 average expression of that gene (for cells that are expressing it)
    per = np.zeros((len(cluter_ids), len(markers)))
    avg = np.zeros((len(cluter_ids), len(markers)))
    mtx = adata_raw.X  # .todense()

    for cn, c in enumerate(cluter_ids):
        tmp_mtx = mtx[adata_raw.obs.cluster.values == c]
        sub_mtx = tmp_mtx[:, midx]
        avg[cn] = nd(sub_mtx.mean(axis=0))
        per[cn] = (sub_mtx > 0).sum(axis=0) / sub_mtx.shape[0]

    fig, ax = plt.subplots(
        figsize=(len(markers) * 0.35, len(cluter_ids) * 0.3 + 1)
    )
    xidx = np.arange(len(markers))
    yidx = np.arange(len(cluter_ids))

    xlabels = markers
    ylabels = cluter_ids

    X, Y = np.meshgrid(xidx, yidx)

    for dn, d in enumerate(per):
        # a = ax.scatter(X[dn], Y[dn], s=d*500+10, c=avg[dn], cmap='Blues')
        a = ax.scatter(X[dn], Y[dn], s=d*100, c=avg[dn], cmap='Blues')

    ax.set_xticks(xidx)
    ax.set_yticks(yidx)

    ax.set_xticklabels(xlabels, rotation=90, ha='center')
    ax.set_yticklabels(ylabels)

    ax.set_xlabel('Gene')
    ax.set_ylabel('Cluster')
    ax.figure.colorbar(a, ax=ax, label='$ln(CPM+1)$')

    handles = [
        mpl.lines.Line2D([
            0], [0], marker='o', color='w', label='  0%',
            markerfacecolor='black', markersize=7
        ),
        mpl.lines.Line2D(
            [0], [0], marker='o', color='w', label=' 25%',
            markerfacecolor='black', markersize=10
        ),
        mpl.lines.Line2D(
            [0], [0], marker='o', color='w', label=' 50%',
            markerfacecolor='black', markersize=12
        ),
        mpl.lines.Line2D(
            [0], [0], marker='o', color='w', label=' 75%',
            markerfacecolor='black', markersize=13.5
        ),
        mpl.lines.Line2D(
            [0], [0], marker='o', color='w', label='100%',
            markerfacecolor='black', markersize=17
        )
    ]
    ax.legend(
        handles=handles,
        loc='center left',
        bbox_to_anchor=(1.15, 0.5),
        title='Percent of cells'
    )

    # plt.show()

    return(fig, ax)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Gets cluster marker genes.
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
        help='H5 AnnData file where clusters have been saved to cluster slot.'
    )

    parser.add_argument(
        '-rgm', '--rank_genes_method',
        action='store',
        dest='rgm',
        default='wilcoxon',
        help='Method used to rank marker genes. Valid options:\
            [wilcoxon|logreg].\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-ncpu', '--number_cpu',
        action='store',
        dest='ncpu',
        default=4,
        type=int,
        help='Number of CPUs to use.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: <h5_anndata>-<tsv_pcs>-clustered)'
    )

    options = parser.parse_args()

    # Fixed settings.
    # verbose = True

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    # sc.set_figure_params(dpi_save = 300)

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}'.format(
            os.path.basename(options.h5.rstrip('h5ad').rstrip('.'))
        )

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=options.h5)

    # See https://github.com/theislab/scanpy/issues/967
    warnings.warn(
        ('WARNING: All functions in this script set use_raw=True, '
            'assuming that adata.raw.to_adata stores ln(CPM+1) or ln(CP10K+1) '
            'normalized data.')
    )

    # NOTE: You should be using the ln(CPM+1) data here. By default these
    #       functions use the .raw attribute of AnnData if present which is
    #       assumed to be ln(CPM+1).

    # Identify cell type makers.
    if options.rgm == 'wilcoxon':
        # If considering simple Wilcoxon vs t-test, choose Wilcoxon as
        # recommended in this paper:
        # https://www.nature.com/articles/nmeth.4612.
        # NOTE: adata.uns['rank_genes_groups'].keys() = ['params',
        #       'scores', 'names', 'logfoldchanges', 'pvals', 'pvals_adj']
        sc.tl.rank_genes_groups(
            adata,
            groupby='cluster',
            groups='all',
            reference='rest',
            use_raw=True,
            method='wilcoxon',
            n_genes=1000,
            corr_method='bonferroni'
        )
    elif options.rgm == 'logreg':
        # Rank genes by logistic regression (multi-variate appraoch),
        # suggested by: https://doi.org/10.1101/258566
        # NOTE: adata.uns['rank_genes_groups'].keys() = ['params',
        #       'scores', 'names']
        sc.tl.rank_genes_groups(
            adata,
            groupby='cluster',
            groups='all',
            reference='rest',
            use_raw=True,
            method='logreg',
            n_genes=1000,
            max_iter=5000  # passed to sklearn.linear_model.LogisticRegression
        )
    else:
        # Other options for cell type marker identification in scanpy:
        # MAST, limma, DESeq2, diffxpy, logreg
        raise Exception('Method not implemented')

    # Save the ranks.
    results_dict = dict()
    for cluster_i in adata.uns['rank_genes_groups']['names'].dtype.names:
        # print(cluster_i)
        # Get keys that we want from the dataframe.
        data_keys = list(
            set(['names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']) &
            set(adata.uns['rank_genes_groups'].keys())
        )
        # Build a table using these keys.
        key_i = data_keys.pop()
        results_dict[cluster_i] = pd.DataFrame(
            row[cluster_i] for row in adata.uns['rank_genes_groups'][key_i]
        )
        results_dict[cluster_i].columns = [key_i]
        for key_i in data_keys:
            results_dict[cluster_i][key_i] = [
                row[cluster_i] for row in adata.uns['rank_genes_groups'][key_i]
            ]
        results_dict[cluster_i]['cluster'] = cluster_i
    marker_df = pd.concat(results_dict, ignore_index=True)

    # Clean up naming.
    marker_df = marker_df.rename(
        columns={'names': 'ensembl_gene_id'},
        inplace=False
    )

    # Add gene_symbols.
    # NOTE: Because rank_genes_groups was run on the ln(CPM+1) data,
    #       we must be sure to use the gene symbols from that data since the
    #       data in adata.X may contain fewer genes, for instance if the
    #       matrix was filtered before scaling.
    marker_df = marker_df.set_index('ensembl_gene_id', inplace=False)
    marker_df = marker_df.join(
        adata.raw.var[['gene_symbols']],
        how='left'
    )
    if (np.invert(marker_df.gene_symbols.notnull()).sum() > 0):
        filt = np.invert(marker_df.gene_symbols.notnull())
        print(marker_df.loc[filt, :])
        raise Exception('Missing gene_symbols in marker_df.')
    elif (np.invert(marker_df.gene_symbols.notna()).sum() > 0):
        filt = np.invert(marker_df.gene_symbols.notna())
        print(marker_df.loc[filt, :])
        raise Exception('Missing gene_symbols in marker_df.')

    marker_df = marker_df.reset_index(drop=False)
    marker_df = marker_df.rename(
        columns={'index': 'ensembl_gene_id'},
        inplace=False
    )
    if 'pvals_adj' in marker_df.columns:  # if logreg, no pvals_adj in marker
        marker_df = marker_df.sort_values(
            by=['cluster', 'scores', 'pvals_adj'],
            ascending=[True, False, True]
        )
    else:
        marker_df = marker_df.sort_values(
            by=['cluster', 'scores'],
            ascending=[True, False]
        )

    # Save the marker dataframe.
    marker_df.to_csv(
        '{}-cluster_markers.tsv.gz'.format(out_file_base),
        sep='\t',
        index=False,
        quoting=csv.QUOTE_NONNUMERIC,
        na_rep='',
        compression=compression_opts
    )
    # Save a filtered marker dataframe for pvals_adj
    marker_df_tmp = marker_df
    additional_out_tag = 'filter'
    if 'pvals_adj' in marker_df_tmp.columns:
        additional_out_tag = '{}__fdr0pt05'.format(additional_out_tag)
        marker_df_tmp = marker_df_tmp.loc[
            (marker_df_tmp['pvals_adj'] < 0.05), :
        ]
        marker_df_tmp.to_csv(
            '{}-cluster_markers-{}.tsv.gz'.format(
                out_file_base,
                additional_out_tag
            ),
            sep='\t',
            index=False,
            quoting=csv.QUOTE_NONNUMERIC,
            na_rep='',
            compression=compression_opts
        )
    # Drop any genes that are in the list of genes to for highly variable genes
    additional_out_tag = '{}__user_variable_genes_exclude'.format(
        additional_out_tag
    )
    filt = marker_df_tmp['ensembl_gene_id'].isin(
        adata.uns['df_highly_variable_gene_filter']['ensembl_gene_id']
    )
    marker_df_tmp = marker_df_tmp.loc[np.invert(filt), :]
    marker_df_tmp.to_csv(
        '{}-cluster_markers-{}.tsv.gz'.format(
            out_file_base,
            additional_out_tag
        ),
        sep='\t',
        index=False,
        quoting=csv.QUOTE_NONNUMERIC,
        na_rep='',
        compression=compression_opts
    )

    # Plot cell type makers.
    # Annoyingly, prefix hardcoded as rank_genes_groups_<cluster_id>.
    _ = sc.pl.rank_genes_groups(
        adata,
        gene_symbols='gene_symbols',
        n_genes=25,
        sharey=False,
        show=False,
        save='-{}.pdf'.format(out_file_base)
    )

    # Plot cell type markers in dotplot.
    # Annoyingly, prefix hardcoded as dotplot.
    # sc.pl.rank_genes_groups_dotplot(
    #     adata,
    #     n_genes=25,
    #     sharey=False,
    #     show=False,
    #     save='-{}.pdf'.format(out_file_base)
    # )

    # Sort by scores: same order as p-values except most methods return scores.
    marker_df = marker_df.sort_values(by=['scores'], ascending=False)
    # Make dataframe of the top 3 markers per cluster
    marker_df_plt = marker_df.groupby('cluster').head(3)
    if (np.invert(marker_df_plt.gene_symbols.notnull()).sum() > 0):
        filt = np.invert(marker_df_plt.gene_symbols.notna())
        print(marker_df_plt.loc[filt, :])
        raise Exception('Missing gene_symbols in marker_df_plt.')
    elif (np.invert(marker_df_plt.gene_symbols.notna()).sum() > 0):
        filt = np.invert(marker_df_plt.gene_symbols.notna())
        print(marker_df_plt.loc[filt, :])
        raise Exception('Missing gene_symbols in marker_df_plt.')
    # Drop markers that are not good... the markers from above are just the
    # top n ranked markers.
    # NOTE: Not sure how to do this when we only have scores, for instance
    #       with logreg marker discovery.
    if 'pvals_adj' in marker_df_plt.columns:
        marker_df_plt.loc[(marker_df_plt['pvals_adj'] < 0.05), :]

    # Dict for plotting both ensembl_gene_id and gene_symbols on the same plot.
    marker_dict_plt = {}
    for idx, row in marker_df_plt.iterrows():
        marker_dict_plt[row['gene_symbols']] = row['ensembl_gene_id']

    # NOTE: You should be using the ln(CPM+1) data here. "$ln(CPM+1)$"
    # Plot cell type markers in dotplot.
    fig, ax = dotplot(
        adata,
        markers=marker_df_plt['gene_symbols'].values
    )
    fig.savefig(
        'dotplot-{}-nodendrogram.pdf'.format(out_file_base),
        #dpi=300,
        bbox_inches='tight'
    )
    # Annoyingly, the prefix is hardcoded as dotplot. Should save to axes.
    # NOTE: We use adata_raw in order to properly get the gene symbols on the
    #       x axis. In theory, the below command should work.
    #           sc.pl.dotplot(adata,
    #               var_names=marker_df_plt['ensembl_gene_id'].to_list(),
    #               use_raw=True, gene_symbols='gene_symbols')
    #       ...however it does not. It seems like properly adding gene_symbols
    #       functionality is an ongoing issue in scanpy.
    #       See: https://github.com/theislab/scanpy/issues/455
    adata_raw = adata.raw.to_adata()
    # Generate dendrogram using the marker genes... this will be used in the
    # below dotplots.
    # NOTE: With latest version of pandas, sc.tl.dendrogram throws an error.
    run_dendrogram = True
    if run_dendrogram:
        sc.tl.dendrogram(
            adata_raw,
            groupby='cluster',
            use_rep='X_pca',
            var_names=marker_dict_plt,
            use_raw=False,
            cor_method='pearson',
            linkage_method='complete',
            optimal_ordering=True,
            inplace=True
        )
    _ = sc.pl.dotplot(
        adata_raw,
        var_names=marker_dict_plt,
        groupby='cluster',
        dendrogram=run_dendrogram,
        use_raw=False,
        show=False,
        color_map='Blues',
        save='_ensembl-{}.pdf'.format(out_file_base)
    )
    _ = sc.pl.dotplot(
        adata_raw,
        var_names=marker_df_plt['gene_symbols'].to_list(),
        groupby='cluster',
        dendrogram=run_dendrogram,
        gene_symbols='gene_symbols',
        use_raw=False,
        show=False,
        color_map='Blues',
        save='-{}.pdf'.format(out_file_base)
    )
    _ = sc.pl.dotplot(
        adata_raw,
        marker_df_plt['gene_symbols'].to_list(),
        groupby='cluster',
        dendrogram=run_dendrogram,
        gene_symbols='gene_symbols',
        standard_scale='var',  # Scale color between 0 and 1
        use_raw=False,
        show=False,
        color_map='Blues',
        save='-{}-standardized.pdf'.format(out_file_base)
    )
    _ = sc.pl.heatmap(
        adata_raw,
        marker_df_plt['gene_symbols'].to_list(),
        groupby='cluster',
        dendrogram=run_dendrogram,
        gene_symbols='gene_symbols',
        use_raw=False,
        show=False,
        save='-{}.pdf'.format(out_file_base)
    )
    _ = sc.pl.heatmap(
        adata_raw,
        marker_df_plt['gene_symbols'].to_list(),
        groupby='cluster',
        dendrogram=run_dendrogram,
        gene_symbols='gene_symbols',
        standard_scale='var',  # Scale color between 0 and 1
        use_raw=False,
        show=False,
        save='-{}-standardized.pdf'.format(out_file_base)
    )


if __name__ == '__main__':
    main()
