#!/usr/bin/env nextflow

import org.yaml.snakeyaml.Yaml
nextflow.enable.dsl = 2

VERSION = "0.0.1" // Do not edit, controlled by bumpversion.

def findDeep(Map m, String key) {
    if (m.containsKey(key)) return m[key]
    m.findResult { k, v -> v instanceof Map ? findDeep(v, key) : null }
}

def make_abspath_if_relpath(dirpath, refpath, filnam) {
  if (dirpath =~ /^\//) {abspath = dirpath} else {abspath = refpath.plus('/' + dirpath)}
  if (filnam.length() > 0) {abspath = abspath.plus('/' + filnam)}
  abspath
}

// Modules to include.
include {
    wf__multiplet;
} from "./modules/multiplet.nf"

include {
    prep_merge_samples;
    merge_samples;
    prep_merge_samples_from_h5ad;
    merge_samples_from_h5ad;
    filter_outlier_cells;
    plot_filtered_cells;
    plot_predicted_sex;
    plot_qc;
    plot_distributions;
    estimate_pca_elbow;
    get_estimate_pca_elbow;
    normalize_and_pca;
    subset_pcs;
    plot_pcs;
    harmony;
    bbknn;
    lisi;
} from "./modules/qc.nf"


include {
    convert_seurat;
    wf__cluster;
    wf__cluster as wf__cluster_harmony;
    wf__cluster as wf__cluster_bbknn;
} from "./modules/cluster.nf"

include {
    wf__umap;
    wf__umap as wf__umap_harmony;
    wf__umap as wf__umap_bbknn;
    // umap_calculate_and_plot;
    // umap_calculate_and_plot as umap_calculate_and_plot__harmony;
} from "./modules/umap.nf"

include {
  wf__celltype_assignment;
} from "./modules/celltype.nf"

include {
  //gather_handover_dataset;
  prep_handover_raw_input;
  wf__handover_dataset
} from "./modules/handover.nf"

// Set default parameters.
params.output_dir           = "${launchDir}/nf-qc_cluster_results"
params.help                 = false
params.run_multiplet        = false
params.mode                 = "conventional"
params.layer                = "none"
params.file_sample_qc       = "no_file__file_sample_qc"
params.file_cellmetadata    = "no_file__file_cellmetadata"
params.file_anndata         = "no_file__file_anndata"
params.genes_exclude_hvg    = "no_file__genes_exclude_hvg"
params.genes_score          = "no_file__genes_score"
params.anndata_compression_opts = 9
// NOTE: The default parameters below were chosen to show the flexiblity of
//       this pipeline. They were not chosen because these are the values one
//       should use for final analysis.
//
// Default key to add in metadata
params.metadata_key_column = [
    value: "experiment_id"
]
// Default parameters for qc plots.
params.plots_qc = [
    facet_columns: [value: ["experiment_id"]],
    variable_columns_distribution_plots: [value: [
        "total_counts,pct_counts_gene_group__mito_transcript"
    ]]
]
// Default parameters for reduced dimension calculations.
// run_downstream_analysis: If false don't run clustering or umaps
params.reduced_dims = [
    run_downstream_analysis: false,
    vars_to_regress: [value: [
        "",
        "total_counts,pct_counts_gene_group__mito_transcript"
    ]],
    n_dims: [
        auto_estimate: false,
        add_n_to_estimate: 0,
        value: [15, 30]
    ]
]
// Default parameters for harmony.
params.harmony = [
    run_process: false,
    variables_and_thetas: [value: [
        [variable: "experiment_id", theta: "1.0"],
        [variable: "experiment_id,phase", theta: "1.0,0.2"]
    ]]
]
// Default parameters for bbknn
params.bbknn = [
    run_process: false,
    batch_variable: [value: ["experiment_id"]]
]
// Default parameters for lisi
params.lisi = [
    run_process: false,
    variables: [value: ["experiment_id,phase"]]
]
// Default parameters for cluster calculations.
params.cluster = [
    number_neighbors: [value: [15, 20]],
    methods: [value: ["leiden"]],
    resolutions: [value: [1.0, 3.0]],
    variables_boxplot: [value: [
        "total_reads,pct_counts_gene_group__mito_transcript"
    ]],
    known_markers: [run_process: false],
    convert_seurat: false
]
// Default parameters for cluster resolution validation using logistic
// regression implemented in keras 
params.cluster_validate_resolution = [
    run_process: true,
    sparsity: [value: [0.0001]],
    train_size_cells: [value: [-1]]
]
// Default parameters for cluster marker gene calculations.
params.cluster_marker = [
    methods: [value: ["wilcoxon"]]
]
// Default parameters for umap calculations.
params.umap = [
    run_process: true,
    colors_quantitative: [value: "total_counts"],
    colors_categorical: [value: "experiment_id,phase"],
    n_neighbors: [value: [15, 30]],
    umap_init: [value: ["X_pca", "spectral"]],
    umap_min_dist: [value: [0.5]],
    umap_spread: [value: [1.0, 5.0]]
]
// Default parameters for sccaff cluster optimization/assessment
params.sccaf = [
  run_assessment: true,
  run_optimization: false,
  min_accuracy: 0.92f,
  default_leiden_res: 4.0
]
// Default parameters for Azimuth
params.azimuth = [
    run_process: true,
]
//Default parameters for dataset handover
params.data_handover = [
    run_process: false,
    cellranger_raw_files_table: '',
    cellbender_files_table:'',
    deconvolution_dir:''
]

// Define the help messsage.
def help_message() {
    log.info """
    ============================================================================
     single cell qc and clustering ~ v${VERSION}
    ============================================================================

    Runs basic single cell qc and clustering

    Usage:
    nextflow run main.nf -profile <local|lsf> -params-file params.yaml [options]

    Mandatory arguments:
        --mode              Specify whether you want to run "conventional"
                            scRNA-seq analysis workflow starting with raw
                            10X data or a "subclustering" analysis of an
                            already existing annData object, or
                            "conventional_h5ad" which runs as "conventional" but takes
                            .h5ad inputs files (annData objects) rather than 10x (cf. --file_paths_h5ad).

        --file_paths_10x    Tab-delimited file containing experiment_id and
                            path_data_10xformat columns.
                            Not required if --mode is set to "conventional_h5ad" (use --file_paths_h5ad instead).

        --file_paths_h5ad   Required only and only if --mode is set to "conventional_h5ad".
                            Tab-delimited file containing experiment_id and
                            path_data_h5adformat columns (file paths tos annData .h5ad objects).

        --file_metadata     Tab-delimited file containing sample metadata.

        --file_anndata      AnnData object which contains a raw count matrix
                            in the aData.layers slot. This is used for sub-
                            clustering of a specific set of cells as for
                            instance T or B cells. Should be provided as .h5ad file.
                            If this file is included, --file_paths_10x and
                            --file_metadata will be ignored.

        --layer             Specify a layer in your annData object that should be used
                            for the subclustering analysis at the PCA step.

    Optional arguments:
        --file_sample_qc    YAML file containing sample quality control
                            filters.

        --file_cellmetadata Tab-delimited file containing experiment_id and
                            data_path_cellmetadata columns.

        --genes_exclude_hvg
                            Tab-delimited file with genes to exclude from
                            highly variable gene list. Must contain
                            ensembl_gene_id column. If no filter, then pass an
                            empty file.

        --genes_score
                            Tab-delimited file with genes to use to score
                            cells. Must contain ensembl_gene_id and score_id
                            columns.  If one score_id == "cell_cycle", then
                            requires a grouping_id column with "G2/M" and "S".
                            If no filter, then pass an empty file.

        --output_dir        Directory name to save results to. (Defaults to
                            'nf-qc_cluster')

        -params-file        YAML file containing analysis parameters. See
                            example in test/params.yml.

    Profiles:
        local               local execution
        lsf                 lsf cluster execution
    """.stripIndent()
}


// Boot message - either help message or the parameters supplied.
if (params.help){
    help_message()
    exit 0
} else {
    log.info """
    ============================================================================
     single cell qc and clustering ~ v${VERSION}
    ============================================================================
    mode                          : ${params.mode}
    file_paths_10x                : ${params.file_paths_10x}
    file_paths_h5ad                : ${params.file_paths_h5ad}
    file_metadata                 : ${params.file_metadata}
    file_anndata                  : ${params.file_anndata}
    file_sample_qc                : ${params.file_sample_qc}
    file_cellmetadata             : ${params.file_cellmetadata}
    layer                         : ${params.layer}
    genes_exclude_hvg             : ${params.genes_exclude_hvg}
    genes_score                   : ${params.genes_score}
    output_dir (output folder)    : ${params.output_dir}
    """.stripIndent()
    // A dictionary way to accomplish the text above.
    // def summary = [:]
    // summary['file_paths_10x'] = params.file_paths_10x
    // log.info summary.collect { k,v -> "${k.padRight(20)} : $v" }.join("\n")
}

log.info "$params"

// Initalize Channels.
// Channel: example init
// Channel
//     .fromPath( params.file_paths_10x )
//     .println()
// Channel: required files
// Channel
//     .fromPath(params.file_paths_10x)
//     .splitCsv(header: true, sep: "\t", by: 1)
//     .map{row -> tuple(row.experiment_id, file(row.data_path_10x_format))}
//     .view()
// Here, this channel is conntected to the 10x input files. They are considered individually
// which is why .splitCsv is there. Then, the respective files for the classic 10x format is added.


if (params.mode == "conventional") {
    refdir_10x = file(params.file_paths_10x).getParent()
    channel__file_paths_10x = Channel
	.fromPath(params.file_paths_10x)
	.splitCsv(header: true, sep: "\t", by: 1)
	.map{row -> tuple(
        row.experiment_id,
        file(make_abspath_if_relpath(row.data_path_10x_format, refdir_10x, "barcodes.tsv.gz")),
        file(make_abspath_if_relpath(row.data_path_10x_format, refdir_10x, "features.tsv.gz")),
        file(make_abspath_if_relpath(row.data_path_10x_format, refdir_10x, "matrix.mtx.gz"))
    )}
    //n_experiments = file(params.file_paths_10x).countLines()
}
if (params.mode == "conventional_h5ad") {
  refdir_h5ad = file(params.file_paths_h5ad).getParent()
    channel__file_paths_h5ad = Channel
	.fromPath(params.file_paths_h5ad)
	.splitCsv(header: true, sep: "\t", by: 1)
	.map{row -> tuple(
    row.experiment_id,
    file(make_abspath_if_relpath(row.h5ad_filepath, refdir_h5ad, ""))
    )}
}

// Initialize known markers channel
// cluster__known_markers is a list of tsv files, first serialize
// the array then run plot_known_markers
// This channel takes a list of genes in, which shoul be processed in the differential gene expression visualization.
if (params.cluster.known_markers.run_process) {
    channel__cluster__known_markers = Channel
        .fromList(params.cluster.known_markers.value)
        .map{row -> tuple(row.file_id, file(row.file))}
} else {
    channel__cluster__known_markers = tuple('', '')
}

// If a sample QC file is supplied, then load it. This is to enable
// filter_outlier_cells.
def params_sample_qc = [:]
if (params.file_sample_qc != "no_file__file_sample_qc") {
    parameter_yaml = new FileInputStream(new File(params.file_sample_qc))
    new Yaml().load(parameter_yaml).each {
        k, v -> if (k == "sample_qc") { params_sample_qc[k] = v}
    }
}
// If the flag to tell us not to filter outliers is not defined then define it
if (findDeep(params_sample_qc, "filter_outliers") != null) {
    if (!params_sample_qc.sample_qc.cell_filters.filter_outliers.containsKey("run_process")) {
        params_sample_qc.sample_qc.cell_filters.filter_outliers.run_process = false
    }
} else {
    params_sample_qc.sample_qc = [
        cell_filters: [
            filter_outliers: [
                run_process: false,
                method: 'IsolationForest',
                metadata_columns: 'pct_counts_gene_group__mito_transcript,log1p_total_counts,log1p_n_genes_by_counts'
            ]
        ]
    ]
}

if (params.data_handover.run_process) {
  channel__file_paths_h5ad_deconv = Channel
	.fromPath(params.data_handover.deconvolution_files_table)
	.splitCsv(header: true, sep: "\t", by: 1)
	.map{row -> file(row.h5ad_filepath)}
  //row.experiment_id contains mangled file name {expid}__[donor{donor_no}|unsassigned|doublet]}

  channel__file_paths_raw = Channel
    .fromPath(params.data_handover.cellranger_raw_files_table)
    .splitCsv(header: true, sep: "\t", by: 1)
    .map{row -> tuple(
        row.experiment_id,
        file(row.data_path_raw_h5))}
}


// Run the workflow.
workflow {
    main:
    log.info "Running main workflow."
    log.info """\n
      -------------------------------
      -- parameters for multiplets --
      -------------------------------
      params.mode = '${params.mode}'
      params_sample_qc.sample_qc.cell_filters.filter_multiplets.run_process = ${params_sample_qc.sample_qc.cell_filters.filter_multiplets.run_process}
      params.file_cellmetadata = '${params.file_cellmetadata}'
    """.stripIndent()

    // Optionally run multiplet filters.
    // set file_cellmetadata
    if (params.mode == "conventional" &&
	params_sample_qc.sample_qc.cell_filters.filter_multiplets.run_process &
        params.file_cellmetadata == "no_file__file_cellmetadata"
    ) {
	log.info "Running multiplet filters."
        wf__multiplet(
            params.output_dir,
            channel__file_paths_10x,
            params_sample_qc.sample_qc.cell_filters.filter_multiplets.expected_multiplet_rate,
            params_sample_qc.sample_qc.cell_filters.filter_multiplets.n_simulated_multiplet,
            params_sample_qc.sample_qc.cell_filters.filter_multiplets.multiplet_threshold_method,
            params_sample_qc.sample_qc.cell_filters.filter_multiplets.scale_log10
        )
        // NOTE: file__cellmetadata is already defined as path, so no need
        // to call file or path.
        file_cellmetadata = wf__multiplet.out.file__cellmetadata
        multiplet_calls = wf__multiplet.out.multiplet_calls
    } else {
        // For some reason cannot use path here, must use file.
        file_cellmetadata = file(params.file_cellmetadata)
        multiplet_calls = null
    }

    // Check whether we conduct the conventional or sub-clustering workflow
    if (params.mode == "conventional" | params.mode == "conventional_h5ad") {
        // This is the conventional workflow
        log.info """\n
        ----------------------------------------------------------------------------
        Running conventional Workflow
        ----------------------------------------------------------------------------
         """.stripIndent()

        if (
            params.mode == "conventional_h5ad" &
		    params.file_anndata == "no_file__file_anndata"
        ) {
	        log.info "Running conventional_h5ad inputs merge."
	        // channel__file_paths_h5ad.view()
            file_sample_qc = file(params.file_sample_qc)
            prep_merge_samples_from_h5ad(channel__file_paths_h5ad)
            merge_samples_from_h5ad(
                    params.output_dir,
                    params.file_metadata,
                    file_sample_qc,
                    file_cellmetadata,
                    params.metadata_key_column.value,
                    prep_merge_samples_from_h5ad.out.h5ad.collect(),
                            params.anndata_compression_opts
            )
            file__anndata_merged = merge_samples_from_h5ad.out.anndata
            file__cells_filtered = merge_samples_from_h5ad.out.cells_filtered
        }
        else if (
            params.mode == "conventional" &
		    params.file_anndata == "no_file__file_anndata"
        ) {
	        log.info "Running conventional inputs merge."
            // Merge the samples, perform cell + gene filtering, add metadata.
            file_sample_qc = file(params.file_sample_qc)
            prep_merge_samples(channel__file_paths_10x)
            merge_samples(
                params.output_dir,
                params.file_paths_10x,
                params.file_metadata,
                file_sample_qc,
                file_cellmetadata,
                params.metadata_key_column.value,
                prep_merge_samples.out.barcodes.collect(),
                prep_merge_samples.out.features.collect(),
                prep_merge_samples.out.matrix.collect(),
                params.anndata_compression_opts
            )
            file__anndata_merged = merge_samples.out.anndata
            file__cells_filtered = merge_samples.out.cells_filtered
            // If add an additional layer of automatic outlier removal (in
            // addition to what was already done in merge_samples, then
            // do that).
	    }




        if (
            params_sample_qc.sample_qc.cell_filters.filter_outliers.run_process
        ) {
            log.info """Running automatic outlier cell filtering."""
            filter_outlier_cells(
                params.output_dir,
                file__anndata_merged,
                file__cells_filtered,
                params_sample_qc.sample_qc.cell_filters.filter_outliers.metadata_columns,
                params_sample_qc.sample_qc.cell_filters.filter_outliers.method,
                params_sample_qc.sample_qc.cell_filters.filter_outliers.outliers_fraction,
                params_sample_qc.sample_qc.cell_filters.filter_outliers.max_samples,
                params.anndata_compression_opts
            )
            file__anndata_merged = filter_outlier_cells.out.anndata
            file__cells_filtered = filter_outlier_cells.out.cells_filtered
        }

        // Plot the filtered cells per sample.
        plot_filtered_cells(
            params.output_dir,
            file__cells_filtered
        )

        // Predict sex from gene expression and check against phenotypes.
        plot_predicted_sex(
            params.output_dir,
            file__anndata_merged
        )


        // Make QC plots of the merged data.
        plot_qc(
            params.output_dir,
            file__anndata_merged,
            params.plots_qc.facet_columns.value
        )


        plot_distributions(
            params.output_dir,
            file__anndata_merged,
            params.plots_qc.facet_columns.value,
            params.plots_qc.variable_columns_distribution_plots.value
        )

        // Azimuth cell-type assignment branches off here
        if ( params.azimuth.run_process ) {
            wf__celltype_assignment(
                params.output_dir,
                file__anndata_merged
                )
            if (multiplet_calls) {
                multiplet_calls.collect().set{ch_multiplet_calls}
                // ch_multiplet_calls.subscribe onNext: { println 'multiplet calls', it }, onComplete: { println 'Done' }
            } else {
                ch_multiplet_calls = null
            }
        }

        if ( params.data_handover.run_process ) {
          prep_handover_raw_input(channel__file_paths_raw)
          wf__celltype_assignment.out.predicted_celltypes
              .concat(channel__file_paths_h5ad_deconv, prep_handover_raw_input.out.raw_feature_bc_mx)
              .collect()
              .set { ch_gather_trigger }
          // ch_gather_trigger.subscribe onNext: { println 'trigger handover', it }, onComplete: { println 'Done' }

          wf__handover_dataset (
              ch_gather_trigger,
              params.output_dir,
              file__anndata_merged,
              params.data_handover.cellranger_raw_files_table,
              params.data_handover.cellbender_files_table,
              params.data_handover.deconvolution_files_table,
              ch_multiplet_calls,
              params.data_handover.Deconvolution_path,
              params.output_dir
          )
        }

        // Normalize, regress (optional), scale, and calculate PCs.
        genes_exclude_hvg = file(params.genes_exclude_hvg)
        genes_score = file(params.genes_score)
        analysis_mode = params.mode
        layer = params.layer
        normalize_and_pca(
            params.output_dir,
            file__anndata_merged,
            analysis_mode,
            layer,
            genes_exclude_hvg,
            genes_score,
            params.reduced_dims.vars_to_regress.value
        )
    }
    else if (params.mode == "subclustering") {
        // This is the subclustering workflow
        // Normalize, regress (optional), scale, and calculate PCs.
        log.info """\n
        ----------------------------------------------------------------------------
        Running subclustering Workflow"
        ----------------------------------------------------------------------------
         """.stripIndent()

        genes_exclude_hvg = file(params.genes_exclude_hvg)
        genes_score = file(params.genes_score)
        file_anndata = file(params.file_anndata)
        analysis_mode = params.mode
        layer = params.layer
        normalize_and_pca(
            params.output_dir,
            file_anndata,
            analysis_mode,
            layer,
            genes_exclude_hvg,
            genes_score,
            params.reduced_dims.vars_to_regress.value
        )
    }
    // Make Seurat dataframes of the normalized anndata
    // convert_seurat(
    //     normalize_and_pca.out.outdir,
    //     normalize_and_pca.out.anndata
    // )
    // Estimate number of PCs to use using eblow from variance explained

    estimate_pca_elbow(
        normalize_and_pca.out.outdir,
        normalize_and_pca.out.anndata,
        params.reduced_dims.n_dims.add_n_to_estimate
    )
    // If the auto estimate PCS, then apply auto esimation.
    if (params.reduced_dims.n_dims.auto_estimate) {
        log.info "n_pcs = automatically estimated."
        n_pcs = estimate_pca_elbow.out.auto_elbow
    } else {
        log.info "n_pcs = Channel.from(params.reduced_dims.n_dims.value)"
        n_pcs = Channel.from(params.reduced_dims.n_dims.value)
    }
    // An alternative version to the one above where the automatically
    // estimated number of PCs is appended to the user supplied list
    // if (params.reduced_dims.n_dims.auto_estimate) {
    //     n_pcs_auto = estimate_pca_elbow.out.auto_elbow
    //     n_pcs = n_pcs_auto.concat(
    //         Channel.from(params.reduced_dims.n_dims.value)
    //     ).unique()
    // } else {
    //     n_pcs = Channel.from(params.reduced_dims.n_dims.value)
    // }
    // Subset PCs to those for anlaysis
        subset_pcs(
            normalize_and_pca.out.outdir,
            normalize_and_pca.out.anndata,
            normalize_and_pca.out.metadata,
            normalize_and_pca.out.pcs,
            normalize_and_pca.out.param_details,
            n_pcs
        )

        // Plot the PCs and gene's contributing to each PC
        plot_pcs(
            subset_pcs.out.outdir,
            subset_pcs.out.anndata,
            n_pcs,
            params.umap.colors_quantitative.value,
            params.umap.colors_categorical.value
        )

        // "Correct" PCs using Harmony
        if (params.harmony.run_process) {
            harmony(
                normalize_and_pca.out.outdir,
                normalize_and_pca.out.anndata,
                normalize_and_pca.out.metadata,
                normalize_and_pca.out.pcs,
                normalize_and_pca.out.param_details,
                n_pcs,
                params.harmony.variables_and_thetas.value
            )
        }

        // Run BBKNN
        // if (params.bbknn.run_process) {
        //     bbknn(
        //         normalize_and_pca.out.outdir,
        //         normalize_and_pca.out.anndata,
        //         normalize_and_pca.out.metadata,
        //         normalize_and_pca.out.pcs,
        //         normalize_and_pca.out.param_details,
        //         n_pcs,
        //         params.bbknn.batch_variable.value
        //     )
        // }
        // TODO: There is a bug below where lisi will be called for each
        // normalize_and_pca call. It just means there will be some duplicate
        // output files in each normalize_and_pca dir and a bit of wasted CPU.
        if (params.lisi.run_process) {
            lisi_input = subset_pcs.out.reduced_dims_params.collect()
            if (params.harmony.run_process) {
                lisi_input = lisi_input.mix(
                    harmony.out.reduced_dims_params.collect()
                )
            }
            if (params.bbknn.run_process) {
                lisi_input = lisi_input.mix(
                    bbknn.out.reduced_dims_params.collect()
                )
            }
            lisi(
                normalize_and_pca.out.outdir,
                normalize_and_pca.out.metadata,
                params.lisi.variables.value,
                lisi_input.collect()
            )
        }
        // Scatter-gather UMAP plots
        if (
            params.reduced_dims.run_downstream_analysis &
            params.umap.run_process
        ) {
            wf__umap(
                subset_pcs.out.outdir,
                subset_pcs.out.anndata,
                subset_pcs.out.metadata,
                subset_pcs.out.pcs,
                subset_pcs.out.reduced_dims,
                "False",
                params.umap.n_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.umap.colors_quantitative.value,
                params.umap.colors_categorical.value
            )
            // Use the data with all of the umaps calculated for downstream
            // clustering, so that we have all of the umap dims in adata.
            cluster_subset_pcs__outdir = wf__umap.out.outdir
            cluster_subset_pcs__anndata = wf__umap.out.anndata
            cluster_subset_pcs__metadata = wf__umap.out.metadata
            cluster_subset_pcs__pcs = wf__umap.out.pcs
            cluster_subset_pcs__reduced_dims = wf__umap.out.reduced_dims
        } else if (params.reduced_dims.run_downstream_analysis) {
            // If running downstream analysis and no umaps, set input for
            // downstream analysis
            cluster_subset_pcs__outdir = subset_pcs.out.outdir
            cluster_subset_pcs__anndata = subset_pcs.out.anndata
            cluster_subset_pcs__metadata = subset_pcs.out.metadata
            cluster_subset_pcs__pcs = subset_pcs.out.pcs
            cluster_subset_pcs__reduced_dims = subset_pcs.out.reduced_dims
        }


        if (params.harmony.run_process & params.umap.run_process) {

            wf__umap_harmony(
                harmony.out.outdir,
                harmony.out.anndata,
                harmony.out.metadata,
                harmony.out.pcs,
                harmony.out.reduced_dims,
                "False",
                params.umap.n_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.umap.colors_quantitative.value,
                params.umap.colors_categorical.value
            )
            // Use the data with all of the umaps calculated for downstream
            // clustering, so that we have all of the umap dims in adata.
            cluster_harmony__outdir = wf__umap_harmony.out.outdir
            cluster_harmony__anndata = wf__umap_harmony.out.anndata
            cluster_harmony__metadata = wf__umap_harmony.out.metadata
            cluster_harmony__pcs = wf__umap_harmony.out.pcs
            cluster_harmony__reduced_dims = wf__umap_harmony.out.reduced_dims
        } else if (params.harmony.run_process) {
            cluster_harmony__outdir = harmony.out.outdir
            cluster_harmony__anndata = harmony.out.anndata
            cluster_harmony__metadata = harmony.out.metadata
            cluster_harmony__pcs = harmony.out.pcs
            cluster_harmony__reduced_dims = harmony.out.reduced_dims
        }
        // NOTE: for BBKNN, we specifically pass the PCs to the reduced dims
        ///      slot not the UMAPS.
        // NOTE: for BBKNN n_neighbors is not needed since already calculated
        if (params.bbknn.run_process & params.umap.run_process) {
            wf__umap_bbknn(
                bbknn.out.outdir,
                bbknn.out.anndata,
                bbknn.out.metadata,
                bbknn.out.pcs,
                bbknn.out.reduced_dims,
                "True",  // Don't look at the reduced_dims parameter
                ["-1"],  // params.cluster.number_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.umap.colors_quantitative.value,
                params.umap.colors_categorical.value
            )
            // Use the data with all of the umaps calculated for downstream
            // clustering, so that we have all of the umap dims in adata.
            cluster_bbknn__outdir = wf__umap_bbknn.out.outdir
            cluster_bbknn__anndata = wf__umap_bbknn.out.anndata
            cluster_bbknn__metadata = wf__umap_bbknn.out.metadata
            cluster_bbknn__pcs = wf__umap_bbknn.out.pcs
            cluster_bbknn__reduced_dims = wf__umap_bbknn.out.reduced_dims
        } else if (params.bbknn.run_process) {
            cluster_bbknn__outdir = bbknn.out.outdir
            cluster_bbknn__anndata = bbknn.out.anndata
            cluster_bbknn__metadata = bbknn.out.metadata
            cluster_bbknn__pcs = bbknn.out.pcs
            cluster_bbknn__reduced_dims = bbknn.out.reduced_dims
        }
        // START LEGACY CODE ----------------------------------------------
        // NOTE: Legacy code due to it being hard to compare
        // Make UMAPs of the reduced dimensions - no scatter gather.
        // umap_calculate_and_plot(
        //     subset_pcs.out.outdir,
        //     subset_pcs.out.anndata,
        //     subset_pcs.out.reduced_dims,
        //     params.umap.colors_quantitative.value,
        //     params.umap.colors_categorical.value,
        //     params.umap.n_neighbors.value,
        //     params.umap.umap_init.value,
        //     params.umap.umap_min_dist.value,
        //     params.umap.umap_spread.value
        // )
        // umap_calculate_and_plot__harmony(
        //     harmony.out.outdir,
        //     harmony.out.anndata,
        //     harmony.out.reduced_dims,
        //     params.umap.colors_quantitative.value,
        //     params.umap.colors_categorical.value,
        //     params.umap.n_neighbors.value,
        //     params.umap.umap_init.value,
        //     params.umap.umap_min_dist.value,
        //     params.umap.umap_spread.value
        // )
        // END LEGACY CODE ------------------------------------------------
        // Cluster the results, varying the resolution.
        // Also, generate UMAPs of the results.i
        if (params.reduced_dims.run_downstream_analysis) {
            wf__cluster(
                cluster_subset_pcs__outdir,
                cluster_subset_pcs__anndata,
                cluster_subset_pcs__metadata,
                cluster_subset_pcs__pcs,
                cluster_subset_pcs__reduced_dims,
                "False",  // use_pcs_as_reduced_dims
                params.cluster.number_neighbors.value,
                params.cluster.methods.value,
                params.cluster.resolutions.value,
                params.cluster.variables_boxplot.value,
                channel__cluster__known_markers,
                params.cluster_validate_resolution.sparsity.value,
                params.cluster_validate_resolution.train_size_cells.value,
                params.cluster_marker.methods.value,
                params.umap.n_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.sccaf.min_accuracy
            )

        }

        if (params.harmony.run_process) {
            wf__cluster_harmony(
                cluster_harmony__outdir,
                cluster_harmony__anndata,
                cluster_harmony__metadata,
                cluster_harmony__pcs,
                cluster_harmony__reduced_dims,
                "False",  // use_pcs_as_reduced_dims
                params.cluster.number_neighbors.value,
                params.cluster.methods.value,
                params.cluster.resolutions.value,
                params.cluster.variables_boxplot.value,
                channel__cluster__known_markers,
                params.cluster_validate_resolution.sparsity.value,
                params.cluster_validate_resolution.train_size_cells.value,
                params.cluster_marker.methods.value,
                params.umap.n_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.sccaf.min_accuracy
            )
        }

        if (params.bbknn.run_process) {
            wf__cluster_bbknn(
                cluster_bbknn__outdir,
                cluster_bbknn__anndata,
                cluster_bbknn__metadata,
                cluster_bbknn__pcs,
                cluster_bbknn__reduced_dims,
                "True",  // use_pcs_as_reduced_dims
                ["-1"],  // params.cluster.number_neighbors.value,
                params.cluster.methods.value,
                params.cluster.resolutions.value,
                params.cluster.variables_boxplot.value,
                channel__cluster__known_markers,
                params.cluster_validate_resolution.sparsity.value,
                params.cluster_validate_resolution.train_size_cells.value,
                params.cluster_marker.methods.value,
                ["-1"],  // params.umap.n_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.sccaf.min_accuracy
            )
        }
}



workflow.onComplete {
    log.info """\n On complete actions \n"""

    try {
         if (params.transfer_to_web.run_process){
            println(["python3", "$workflow.projectDir/bin/transfer_data.py", "--pipeline", "QC", "--name", "${params.transfer_to_web.tranche.tranche_run_name}", "--directory", "${params.output_dir}","--out_directory", "${params.output_dir}"].execute().text)
            println(["python3", "$workflow.projectDir/bin/transfer_data.py", "--pipeline", "Deconvolution", "--name", "${params.transfer_to_web.tranche.tranche_run_name}", "--directory", "${params.transfer_to_web.tranche.Deconvolution_path}","--out_directory", "${params.output_dir}"].execute().text)
            println(["python3", "$workflow.projectDir/bin/transfer_data.py", "--pipeline", "Fetch", "--name", "${params.transfer_to_web.tranche.tranche_run_name}", "--directory", "${params.transfer_to_web.tranche.Fetch_path}","--out_directory", "${params.output_dir}"].execute().text)
            println(["python3", "$workflow.projectDir/bin/transfer_data.py", "--pipeline", "Cellbender", "--name", "${params.transfer_to_web.tranche.tranche_run_name}", "--directory", "${params.transfer_to_web.tranche.Cellbender_path}","--out_directory", "${params.output_dir}"].execute().text)
            println(["rsync", "-vr", "${params.output_dir}/handover/${params.transfer_to_web.tranche.tranche_run_name}", "${params.transfer_to_web.destination}"].execute().text)
            log.info "Transfer made to the website- check the results on https://apps.hgi.sanger.ac.uk/scrna/"
         }
         else{
             log.info "Not transfering data"
         }
      }catch(Exception ex) {
         log.info('Transfer is not activated - the results will not be reflected on https://apps.hgi.sanger.ac.uk/scrna/');
         log.info """\n
         to activate transfer please provide folowing parameters in params.yml file:
         ----------------------------------------------------------------------------
         transfer_to_web:
            run_process: 'true'
            destination: 'ubuntu@prod_swarm:/volume/scRNA_test_app/scrna_static_and_media_files/bin'
            tranche:
                tranche_run_name: 'Test 4' #its important that the name is unique, if not this will overwrite the existing run
                Fetch_path: 'false'
                Cellbender_path: 'path/to/pipeline/run/directory'
                Deconvolution_path: 'path/to/pipeline/run/directory'
        ----------------------------------------------------------------------------
         """.stripIndent()
      }

    // executed after workflow finishes
    // ------------------------------------------------------------------------
    log.info """\n
    ----------------------------------------------------------------------------
     pipeline execution summary
    ----------------------------------------------------------------------------
    Completed         : ${workflow.complete}
    Duration          : ${workflow.duration}
    Success           : ${workflow.success}
    Work directory    : ${workflow.workDir}
    Exit status       : ${workflow.exitStatus}
    """.stripIndent()
}
