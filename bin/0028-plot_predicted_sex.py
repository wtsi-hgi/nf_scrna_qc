#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-05-20'
__version__ = '0.0.1'

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import plotnine as plt9
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (10, 10)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read anndata object. Read chrX and chrY genes, plot scatterplot
            of mean expression of those signature across samples. Anndata
            should have experiment_id and sex columns.
            """
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '-Y', '--chrY_genes',
        default='',
        dest='Y',
        required=False,
        help='TSV file of Y genes. If none, uses all genes on Y chr.'
    )

    parser.add_argument(
        '-X', '--chrX_genes',
        default='',
        dest='X',
        required=False,
        help='TSV file of X genes. If none, uses XIST.'
    )

    parser.add_argument(
        '-o', '--output_file',
        default='scatterplot-sex_sample_swap_check',
        dest='o',
        help='Basename for output files.'
    )

    options = parser.parse_args()

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)

    # If we have a flag for cells that pass QC then filter down to them
    if 'cell_passes_qc' in adata.obs:
        adata = adata[adata.obs['cell_passes_qc'], :]
        del adata.obs['cell_passes_qc']

    # Read Chr X and Chr Y genes
    if options.X != '':
        X = pd.read_csv(options.X, sep="\t")
        X = X['ensembl_gene_id']
        X_lab = "Mean X chr gene expression (counts)"
    else:
        X = ['ENSG00000229807']
        # X = ['XIST']
        X_lab = "Mean XIST gene expression (counts)"
    if options.Y != '':
        Y = pd.read_csv(options.Y, sep="\t")
        Y = Y['ensembl_gene_id']
    else:
        Y = [
            "ENSG00000184895",
            "ENSG00000129824",
            "ENSG00000067646",
            "ENSG00000176679",
            "ENSG00000099715",
            "ENSG00000168757",
            "ENSG00000099721",
            "ENSG00000092377",
            "ENSG00000099725",
            "ENSG00000233803",
            "ENSG00000229549",
            "ENSG00000228927",
            "ENSG00000258992",
            "ENSG00000238074",
            "ENSG00000236424",
            "ENSG00000114374",
            "ENSG00000067048",
            "ENSG00000183878",
            "ENSG00000154620",
            "ENSG00000129864",
            "ENSG00000129862",
            "ENSG00000165246",
            "ENSG00000129873",
            "ENSG00000182415",
            "ENSG00000172468",
            "ENSG00000169953",
            "ENSG00000286265",
            "ENSG00000012817",
            "ENSG00000198692",
            "ENSG00000280969",
            "ENSG00000242875",
            "ENSG00000234414",
            "ENSG00000244395",
            "ENSG00000242389",
            "ENSG00000169807",
            "ENSG00000169800",
            "ENSG00000226941",
            "ENSG00000169789",
            "ENSG00000183753",
            "ENSG00000188120",
            "ENSG00000205944",
            "ENSG00000169763",
            "ENSG00000172352",
            "ENSG00000183795",
            "ENSG00000187191",
            "ENSG00000205916",
            "ENSG00000185894",
            "ENSG00000172288"
        ]
        # Same as above, but hugo names.
        # Y = [
        #     "SRY",
        #     "RPS4Y1",
        #     "ZFY",
        #     "TGIF2LY",
        #     "PCDH11Y",
        #     "TSPY2",
        #     "AMELY",
        #     "TBL1Y",
        #     "PRKY",
        #     "TSPY4",
        #     "TSPY8",
        #     "TSPY3",
        #     "TSPY1",
        #     "TSPY9P",
        #     "TSPY10",
        #     "USP9Y",
        #     "DDX3Y",
        #     "UTY",
        #     "TMSB4Y",
        #     "VCY",
        #     "VCY1B",
        #     "NLGN4Y",
        #     "CDY2B",
        #     "CDY2A",
        #     "HSFY1",
        #     "HSFY2",
        #     "AC007244.1",
        #     "KDM5D",
        #     "EIF1AY",
        #     "RPS4Y2",
        #     "RBMY1B",
        #     "RBMY1A1",
        #     "RBMY1D",
        #     "RBMY1E",
        #     "PRY2",
        #     "RBMY1F",
        #     "RBMY1J",
        #     "PRY",
        #     "BPY2",
        #     "DAZ1",
        #     "DAZ2",
        #     "PRYP3",
        #     "CDY1B",
        #     "BPY2B",
        #     "DAZ3",
        #     "DAZ4",
        #     "BPY2C",
        #     "CDY1"
        # ]

    # Make the plot
    adata.var['X_chr-gene'] = np.in1d(adata.var.index, X)
    adata.var['Y_chr-gene'] = np.in1d(adata.var.index, Y)
    adata.obs['X_chr-sum'] = adata[:, adata.var['X_chr-gene']].X.todense(
        ).sum(axis=1)
    adata.obs['Y_chr-sum'] = adata[:, adata.var['Y_chr-gene']].X.todense(
        ).sum(axis=1)
    if 'sex' not in adata.obs.columns:
        adata.obs['sex'] = 'not reported'
    df = adata.obs[['experiment_id', 'sex', 'Y_chr-sum', 'X_chr-sum']]
    df = df.groupby(['experiment_id', 'sex']).mean().dropna().reset_index()

    # Save scatterplot with mean expression per sample
    plt = plt9.ggplot(df) + plt9.aes(x='X_chr-sum', y='Y_chr-sum', color='sex')
    plt = plt + plt9.theme_bw()
    plt = plt + plt9.scale_colour_brewer(type='qual', palette='Dark2')
    plt = plt + plt9.geom_point(alpha=0.45)
    plt = plt + plt9.ylab("Mean Y chr gene expression (counts)")
    plt = plt + plt9.xlab(X_lab)
    plt.save(
        '{}.png'.format(options.o),
        #dpi=300,
        width=4,
        height=4
    )


if __name__ == '__main__':
    main()
