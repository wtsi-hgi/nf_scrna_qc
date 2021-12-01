# Output Files and Directories

## example directory structure
```
results
├── azimuth
├── handover
│   └── minimal_dataset
├── multiplet.method=scrublet
│   └── plots
├── normalize=total_count.vars_to_regress=pct_counts_gene_group__mito_transcript.hvg_exclude=genes_remove_hvg_v001.scores=genes_score_v001
│   ├── plots
│   └── reduced_dims-vars_to_regress=pct_counts_gene_group__mito_transcript-pca.n_pcs=18
│       └── plots
└── plots
```

## output files
```
results_nf_qc_cluster
├── adata-cell_filtered_per_experiment.tsv.gz
├── adata.h5ad
├── adata-mads.tsv
├── azimuth/
├── file_cellmetadata.tsv
├── handover/
│   └── minimal_dataset/
│       ├── files.tsv
├── multiplet.method=scrublet/
│   └── plots/
├── normalize=total_count.vars_to_regress=pct_counts_gene_group__mito_transcript.hvg_exclude=genes_remove_hvg_v001.scores=genes_score_v001/
│   ├── adata-metadata.tsv.gz
│   ├── adata-normalized_pca-counts.h5ad
│   ├── adata-normalized_pca.h5ad
│   ├── adata-normalized_pca-knee.tsv
│   ├── adata-pcs.tsv.gz
│   ├── plots/
│   │   ├── adata-normalized_pca-knee-variance_ratio-spline=interp1d-knee_normalized.png
│   │   ├── adata-normalized_pca-knee-variance_ratio-spline=interp1d-knee_raw.png
│   │   ├── adata-normalized_pca-knee-variance_ratio-spline=None-knee_normalized.png
│   │   ├── adata-normalized_pca-knee-variance_ratio-spline=None-knee_raw.png
│   │   ├── adata-normalized_pca-knee-variance-spline=interp1d-knee_normalized.png
│   │   ├── adata-normalized_pca-knee-variance-spline=interp1d-knee_raw.png
│   │   ├── adata-normalized_pca-knee-variance-spline=None-knee_normalized.png
│   │   ├── adata-normalized_pca-knee-variance-spline=None-knee_raw.png
│   │   ├── filter_genes_dispersion-adata.pdf
│   │   ├── highest_expr_genes-adata.pdf
│   │   ├── pca_variance_ratio-adata-log.pdf
│   │   └── pca_variance_ratio-adata.pdf
│   └── reduced_dims-vars_to_regress=pct_counts_gene_group__mito_transcript-pca.n_pcs=18/
│       ├── plots/
│       │   ├── pca_loadings-pca-n_pcs=18.png
│       │   ├── pca-pca-experiment_id.png
│       │   ├── pca-pca-n_cells.png
│       │   ├── pca-pca-pct_counts_gene_group__mito_transcript.png
│       │   └── pca-pca-total_counts.png
│       └── reduced_dims.tsv.gz
└── plots/
    ├── adata-cell_desity.png
    ├── adata-cell_filtered_per_experiment-fraction_before_after.png
    ├── adata-cell_filtered_per_experiment-fraction_cells_excluded.png
    ├── adata-cell_filtered_per_experiment-n_cells_before_after.png
    ├── adata-cell_filtered_per_experiment-n_cells_excluded.png
    ├── adata-mads-n_genes_by_counts.png
    ├── adata-mads-pct_counts_gene_group__mito_protein.png
    ├── adata-mads-pct_counts_gene_group__mito_transcript.png
    ├── adata-mads-pct_counts_gene_group__ribo_protein.png
    ├── adata-mads-pct_counts_gene_group__ribo_rna.png
    ├── adata-mads-total_counts.png
    ├── adata-outlier_cells.png
    ├── plot_ecdf.var=pct_counts_gene_group__mito_transcript.color=experiment_id-adata.png
    ├── plot_ecdf.var=total_counts.color=experiment_id-adata.png
    ├── plot_ecdf-x_log10.var=pct_counts_gene_group__mito_transcript.color=experiment_id-adata.png
    ├── plot_ecdf-x_log10.var=total_counts.color=experiment_id-adata.png
    ├── plot_histogram.var=pct_counts_gene_group__mito_transcript.facet=experiment_id-adata.png
    ├── plot_histogram.var=total_counts.facet=experiment_id-adata.png
    ├── plot_histogram-x_log10.var=pct_counts_gene_group__mito_transcript.facet=experiment_id-adata.png
    ├── plot_histogram-x_log10.var=total_counts.facet=experiment_id-adata.png
    ├── plot_nfeature_mt_cellpassqc.facet=experiment_id-adata.png
    ├── plot_nfeature_mt_density.facet=experiment_id-adata.png
    ├── plot_umi_mt_cellpassqc.facet=experiment_id-adata.png
    ├── plot_umi_mt_density.facet=experiment_id-adata.png
    ├── plot_umi_ngene_cellpassqc.facet=experiment_id-adata.png
    ├── plot_umi_ngene_mt_density.facet=experiment_id-adata.png
    ├── plot_umi_ngene_mt.facet=experiment_id-adata.png
    └── scatterplot-sex_sample_swap_check.png 
```
