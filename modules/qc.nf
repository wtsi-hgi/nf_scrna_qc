#!/usr/bin/env nextflow

def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}


process prep_merge_samples {
    input:
        tuple(
            val(experiment_id),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix)
        )

    output:
        path("${experiment_id}---barcodes.tsv.gz", emit: barcodes)
        path("${experiment_id}---features.tsv.gz", emit: features)
        path("${experiment_id}---matrix.mtx.gz", emit: matrix)

    script:
        """
        ln --physical ${file_10x_barcodes} ${experiment_id}---barcodes.tsv.gz
        ln --physical ${file_10x_features} ${experiment_id}---features.tsv.gz
        ln --physical ${file_10x_matrix} ${experiment_id}---matrix.mtx.gz
        """
}

process prep_merge_samples_from_h5ad {
    input:
        tuple(
            val(experiment_id),
            path(file_h5ad)
        )

    output:
        path("${experiment_id}---h5ad.h5ad", emit: h5ad)

    script:
        """
        ln --physical ${file_h5ad} ${experiment_id}---h5ad.h5ad
        """
}

process merge_samples {
    // Takes a list of raw 10x files and merges them into one anndata object.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache true        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode         // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file_paths_10x)
        path(file_metadata)
        path(file_params)
        path(file_cellmetadata)
        val(metadata_key)
        file(file_10x_barcodes)
        file(file_10x_features)
        file(file_10x_matrix)
        val(anndata_compression_opts)

    // NOTE: use path here and not file see:
    //       https://github.com/nextflow-io/nextflow/issues/1414
    output:
        path("${runid}-adata.h5ad", emit: anndata)
        path(
            "${runid}-adata-cell_filtered_per_experiment.tsv.gz",
            emit: cells_filtered
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // String filename = './parameters.yml'
        // yaml.dump(file_params , new FileWriter(filename))
        // Customize command for optional files.
        cmd__params = ""
        if (file_params.name != "no_file__file_sample_qc") {
            cmd__params = "--params_yaml ${file_params}"
        }
        cmd__cellmetadata = ""
        if (file_cellmetadata.name != "no_file__file_cellmetadata") {
            cmd__cellmetadata = "--cell_metadata_file ${file_cellmetadata}"
        }
        files__barcodes = file_10x_barcodes.join(',')
        files__features = file_10x_features.join(',')
        files__matrix = file_10x_matrix.join(',')
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "merge_samples: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        0025-nf_helper__prep_tenxdata_file.py \
            --barcodes_list ${files__barcodes} \
            --features_list ${files__features} \
            --matrix_list ${files__matrix} \
            --tenxdata_file ${file_paths_10x} \
            --output_file nf_prepped__file_paths_10x.tsv
        0025-scanpy_merge.py \
            --tenxdata_file nf_prepped__file_paths_10x.tsv \
            --sample_metadata_file ${file_metadata} \
            --sample_metadata_columns_delete "sample_status,study,study_id" \
            --metadata_key ${metadata_key} \
            --number_cpu ${task.cpus} \
            --output_file ${runid}-adata \
            --anndata_compression_opts ${anndata_compression_opts} \
            ${cmd__params} \
            ${cmd__cellmetadata}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

process merge_samples_from_h5ad {
    // Takes a list of h5ad files and merges them into one anndata object.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache true        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode         // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file_metadata)
        path(file_params)
        path(file_cellmetadata)
        val(metadata_key)
        file(file_h5ad)
        val(anndata_compression_opts)

    // NOTE: use path here and not file see:
    //       https://github.com/nextflow-io/nextflow/issues/1414
    output:
        path("${runid}-adata.h5ad", emit: anndata)
        path(
            "${runid}-adata-cell_filtered_per_experiment.tsv.gz",
            emit: cells_filtered
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // String filename = './parameters.yml'
        // yaml.dump(file_params , new FileWriter(filename))
        // Customize command for optional files.
        cmd__params = ""
        if (file_params.name != "no_file__file_sample_qc") {
            cmd__params = "--params_yaml ${file_params}"
        }
        cmd__cellmetadata = ""
        if (file_cellmetadata.name != "no_file__file_cellmetadata") {
            cmd__cellmetadata = "--cell_metadata_file ${file_cellmetadata}"
        }
        files__h5ad = file_h5ad.join(',')
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "merge_samples: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        0025-nf_helper__prep_h5addata_file.py \
            --h5ad_list ${files__h5ad} \
            --output_file nf_prepped__file_paths_h5ad.tsv
        0025-scanpy_merge_from_h5ad.py \
            --h5addata_file nf_prepped__file_paths_h5ad.tsv \
            --sample_metadata_file ${file_metadata} \
            --sample_metadata_columns_delete "sample_status,study,study_id" \
            --metadata_key ${metadata_key} \
            --number_cpu ${task.cpus} \
            --output_file ${runid}-adata \
            --anndata_compression_opts ${anndata_compression_opts} \
            ${cmd__params} \
            ${cmd__cellmetadata}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

process filter_outlier_cells {
    // Takes annData object, plots predicted outlier cells
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__cells_filtered)
        val(metadata_columns)
        val(method)
        val(outliers_fraction)
        val(max_samples)
        val(anndata_compression_opts)

    output:
        path("${runid}-adata.h5ad", emit: anndata)
        path(
            "${runid}-adata-cell_filtered_per_experiment.tsv.gz",
            emit: cells_filtered
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // Append run_id to output file.
        outfile = "${runid}-adata"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "filter_outlier_cells: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        0026-filter_outlier_cells.py \
            --h5_anndata ${file__anndata} \
            --cell_filtered_per_experiment_file ${file__cells_filtered} \
            --outliers_fraction 0 \
            --metadata_columns ${metadata_columns} \
            --cell_qc_column cell_passes_qc \
            --method ${method} \
            --outliers_fraction ${outliers_fraction} \
            --max_samples ${max_samples} \
            --output_file ${outfile} \
            --anndata_compression_opts ${anndata_compression_opts}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process plot_filtered_cells {
    // Takes annData object, plots filtered cells
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__filtered_cells)

    output:
        val(outdir, emit: outdir)
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_filtered_cells: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        0026-plot_filtered_cells.py \
            --tsv_file ${file__filtered_cells} \
            --output_file ${runid}-adata-cell_filtered_per_experiment
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process plot_predicted_sex {
    // Takes annData object, plots the predicted sex fron gene expression
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)

    output:
        val(outdir, emit: outdir)
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // Append run_id to output file.
        outfile = "${runid}-scatterplot-sex_sample_swap_check"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_predicted_sex: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        0028-plot_predicted_sex.py \
            --h5_anndata ${file__anndata} \
            --output_file ${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process plot_distributions {
    // Takes annData object, generates basic qc plots
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        each facet_columns
        each variable_columns_distribution_plots

    output:
        val(outdir, emit: outdir)
        path("plots/*.png")
        path("plots/*.pdf") optional true
        path("*.tsv") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        // Append run_id to output file.
        outfile = "${runid}-${outfile}"
        // Figure out if we are facetting the plot and update accordingly.
        cmd__facet_columns = ""
        if (facet_columns != "") {
            cmd__facet_columns = "--facet_columns ${facet_columns}"
        }
        // Run distribution across cells if a value is specified
        cmd__anndataobs = ""
        cmd__anndataobs_ecdf = ""
        if (variable_columns_distribution_plots != "") {
            cmd__anndataobs = "plot_anndataobs_distribution_across_cells.py"
            cmd__anndataobs = "${cmd__anndataobs} --h5_anndata ${file__anndata}"
            cmd__anndataobs = "${cmd__anndataobs} --output_file ${outfile}"
            cmd__anndataobs = "${cmd__anndataobs} --variable_columns ${variable_columns_distribution_plots}"
            cmd__anndataobs = "${cmd__anndataobs} ${cmd__facet_columns}"
            cmd__anndataobs_ecdf = "${cmd__anndataobs} --ecdf"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_distributions: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        ${cmd__anndataobs}
        ${cmd__anndataobs_ecdf}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process plot_qc {
    // Takes annData object, generates basic qc plots
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        each facet_columns

    output:
        val(outdir, emit: outdir)
        path("plots/*.png")
        path("plots/*.pdf") optional true
        path("*.tsv") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        // Append run_id to output file.
        outfile = "${runid}-${outfile}"
        // Figure out if we are facetting the plot and update accordingly.
        cmd__facet_columns = ""
        if (facet_columns != "") {
            cmd__facet_columns = "--facet_columns ${facet_columns}"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_qc: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        plot_qc_umi_nfeature_mt.py \
            --h5_anndata ${file__anndata} \
            --output_file ${outfile} \
            ${cmd__facet_columns}
        plot_qc_umi_mt_density.py \
            --h5_anndata ${file__anndata} \
            --output_file ${outfile} \
            ${cmd__facet_columns}
        plot_qc_nfeature_mt_density.py \
            --h5_anndata ${file__anndata} \
            --output_file ${outfile} \
            ${cmd__facet_columns}
        0027-calculate_mads.py \
            --h5_anndata ${file__anndata} \
            --qc_key 'pct_counts_gene_group__mito_transcript,pct_counts_gene_group__mito_protein,pct_counts_gene_group__ribo_protein,pct_counts_gene_group__ribo_rna,total_counts,n_genes_by_counts' \
            --output_file ${outfile}-mads
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process normalize_and_pca {
    // Takes annData object, nomalizes across samples, calculates PCs.
    // NOTE: Once normalization is set, it would be faster to normalize per
    //       sample and then merge.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        val(analysis_mode)
        val(layer)
        path(file__genes_exclude_hvg)
        path(file__genes_score)
        each vars_to_regress

    output:

        val(outdir, emit: outdir)

        val("${outdir}", emit: outdir3)
        path("${runid}-adata-normalized_pca.h5ad", emit: anndata)
        path("${runid}-adata-metadata.tsv.gz", emit: metadata)
        path("${runid}-adata-pcs.tsv.gz", emit: pcs)

        path(
            "${runid}-adata-normalized_pca-counts.h5ad",
            emit: anndata_filtered_counts
        )
        val("${param_details}", emit: param_details)
        path("plots/*.pdf")
        path("plots/*.png") optional true
        // tuple(
        //     val(outdir),
        //     path("${runid}-adata-normalized_pca.h5ad"),
        //     path("${runid}-adata-metadata.tsv.gz"),
        //     path("${runid}-adata-pcs.tsv.gz"),
        //     emit: results
        // )

    script:
        runid = random_hex(16)
        analysis_mode = "${analysis_mode}"
        if (analysis_mode == "subclustering"){
            layer = "${layer}"
        }
        // Add any variables we are regressing to the output dir.
        param_details="vars_to_regress=none"
        if (vars_to_regress == "") {
            cmd__vars_to_regress = ""
        } else {
            param_details = "vars_to_regress=${vars_to_regress}"
            cmd__vars_to_regress = "--vars_to_regress ${vars_to_regress}"
        }

        // todo - mo11 - these paths are confusing

        outdir = "${outdir_prev}/normalize=total_count.${param_details}"
        // Add details on the genes we are exlcuding from hgv list.
        file_vge = "${file__genes_exclude_hvg.getSimpleName()}"
        outdir = "${outdir}.hvg_exclude=${file_vge}"
        // Add details on the scores we are using.
        file_score = "${file__genes_score.getSimpleName()}"
        outdir = "${outdir}.scores=${file_score}"


        // this is where the subfolder 1 is determined

        // Customize command for optional files.
        cmd__genes_exclude_hvg = ""
        if (file__genes_exclude_hvg.name != "no_file__genes_exclude_hvg") {
            cmd__genes_exclude_hvg = "--variable_genes_exclude ${file__genes_exclude_hvg}"
        }
        cmd__genes_score = ""
        if (file__genes_score.name != "no_file__genes_score") {
            cmd__genes_score = "--score_genes ${file__genes_score}"
        }
        // Basic details on the run.
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"

        """
        echo "normalize_pca: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        0035-scanpy_normalize_pca.py \
            --h5_anndata ${file__anndata} \
            --overwrite_x_with_layer ${layer} \
            --output_file ${runid}-adata \
            --number_cpu ${task.cpus} \
            ${cmd__vars_to_regress} \
            ${cmd__genes_exclude_hvg} \
            ${cmd__genes_score}
        mkdir plots

        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
        // Old version with bash evaluation of optional commands
        //
        // echo "normalize_pca: ${process_info}"
        // # If there are entries in the variable_genes_exclude file, add it to
        // # the call.
        // cmd__vg_exclude="--variable_genes_exclude ${file__genes_exclude_hvg}"
        // val=\$(cat ${file__genes_exclude_hvg} | wc -l)
        // if [ \$val -eq 0 ]; then cmd__vg_exclude=""; fi
        // # If there are entries in the score_genes file, add it to the call.
        // cmd__score_genes="--score_genes ${file__genes_score}"
        // val=\$(cat ${file__genes_score} | wc -l)
        // if [ \$val -eq 0 ]; then cmd__score_genes=""; fi
        // 0035-scanpy_normalize_pca.py \
        //     --h5_anndata ${file__anndata} \
        //     --output_file ${runid}-adata \
        //     --number_cpu ${task.cpus} \
        //     ${cmd__vars_to_regress} \
        //     \${cmd__vg_exclude} \
        //     \${cmd__score_genes}
        // mkdir plots
        // mv *pdf plots/ 2>/dev/null || true
        // mv *png plots/ 2>/dev/null || true
}


process estimate_pca_elbow {
    // Takes annData object, estiamtes the elbow in PC var explained.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        val(add_n_to_estimate)

    output:
        val(outdir, emit: outdir)
        path("${runid}-${outfile}.tsv", emit: pca_elbow_estimate)
        env(AUTO_ELBOW, emit: auto_elbow)
        path("plots/*.png")
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        log.info("""outdir = ${outdir}""")
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        outfile = "${outfile}-knee"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "estimate_pca_elbow: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        0030-estimate_pca_elbow.py \
            --h5_anndata ${file__anndata} \
            --add_n_pcs_to_elbow ${add_n_to_estimate} \
            --output_file ${runid}-${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        AUTO_ELBOW=\$(cat ${runid}-${outfile}-auto_elbow_estimate.tsv)
        """
}


process get_estimate_pca_elbow {
    //
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    input:
        path(file__estimate)

    output:
        stdout emit: estimate

    script:
        """
        cat ${file__estimate}
        """
}


process subset_pcs {
    // Takes PCs (rows = cell barcodes) and subsets down to a specified number.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    //saveAs: {filename -> filename.replaceAll("${runid}-", "")},
    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("normalized_pca.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("${param_details}.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        val(params__pcs)
        each n_pcs

    output:
        val(outdir, emit: outdir)
        path(file__anndata, emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path("${runid}-reduced_dims.tsv.gz", emit: reduced_dims)
        // NOTE: passing the param details as an unpublished file is messy,
        // but I could not get collect of ${param_details} and file to work.
        path(
            "reduced_dims-${param_details}.tsv.gz",
            emit: reduced_dims_params
        )
        // val(n_pcs, emit: n_pcs)
        // tuple(
        //     val(outdir),
        //     path("${runid}-reduced_dims.tsv.gz"),
        //     emit: results
        // )

    script:
        runid = random_hex(16)
        param_details = "${params__pcs}-pca"
        param_details = "${param_details}.n_pcs=${n_pcs}"
        outdir = "${outdir_prev}/reduced_dims-${param_details}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "subset_pcs: ${process_info}"
        echo "publish_directory: ${outdir}"
        0045-subset_pca_file.py \
            --tsv_pcs ${file__pcs} \
            --number_pcs ${n_pcs} \
            --output_file ${runid}-reduced_dims
        cp ${runid}-reduced_dims.tsv.gz \
            reduced_dims-${param_details}.tsv.gz
        """
}

process plot_pcs {
    // Takes annData object with PCs and returns plots
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        each n_pcs
        val(colors_quantitative)
        val(colors_categorical)

    output:
        val(outdir, emit: outdir)
        path("plots/*.png")
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        // outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        outfile = "pca"
        cmd__colors_quant = ""
        if (colors_quantitative != "") {
            cmd__colors_quant = "--colors_quantitative ${colors_quantitative}"
        }
        cmd__colors_cat = ""
        if (colors_categorical != "") {
            cmd__colors_cat = "--colors_categorical ${colors_categorical}"
        }
        // drop_legend_n = "-1"
        // if (cmd__colors_cat.contains("experiment_id")) {
        //     drop_legend_n = "8"
        // }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "pca_plot: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        pca_plot.py \
            --h5_anndata ${file__anndata} \
            --num_pcs ${n_pcs} \
            ${cmd__colors_quant} \
            ${cmd__colors_cat} \
            --output_file ${runid}-${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process harmony {
    // Takes PCs (rows = cell barcodes) and metadata (rows = cell barcodes),
    // runs Harmony
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("normalized_pca.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("${param_details}.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        val(params__pcs)
        each n_pcs
        each variables_and_thetas

    output:
        val(outdir, emit: outdir)
        path(file__anndata, emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path("${runid}-reduced_dims.tsv.gz", emit: reduced_dims)
        // NOTE: passing the param details as an unpublished file is messy,
        // but I could not get collect of ${param_details} and file to work.
        path(
            "reduced_dims-${param_details}.tsv.gz",
            emit: reduced_dims_params
        )
        // val(n_pcs, emit: n_pcs)
        // tuple(
        //     val(outdir),
        //     path("${runid}-reduced_dims.tsv.gz"),
        //     val(n_pcs),
        //     emit: results
        // )

    script:
        runid = random_hex(16)
        param_details = "${params__pcs}-harmony"
        param_details = "${param_details}.n_pcs=${n_pcs}"
        param_details = "${param_details}.variables=${variables_and_thetas.variable}"
        outdir = "${outdir_prev}/reduced_dims-${param_details}"
        outdir = "${outdir}.thetas=${variables_and_thetas.theta}"
        theta_str = "${variables_and_thetas.theta}".replaceAll("\\.", "pt")
        param_details = "${param_details}.thetas=${theta_str}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "harmony: ${process_info}"
        echo "publish_directory: ${outdir}"
        0045-harmony_process_pcs.py \
            --pca_file ${file__pcs} \
            --metadata_file ${file__metadata} \
            --metadata_columns ${variables_and_thetas.variable} \
            --theta ${variables_and_thetas.theta} \
            --n_pcs ${n_pcs} \
            --out_file ${runid}-reduced_dims
        cp ${runid}-reduced_dims.tsv.gz \
            reduced_dims-${param_details}.tsv.gz
        """
        // NOTE: below code for harmony in R
        // 0045-harmony_process_pcs.R \
        //     --pca_file ${file__pcs} \
        //     --metadata_file ${file__metadata} \
        //     --metadata_columns ${variables_and_thetas.variable} \
        //     --theta ${variables_and_thetas.theta} \
        //     --n_pcs ${n_pcs} \
        //     --out_file ${runid}-reduced_dims
}


process bbknn {
    // Calulates bbknn neighbors and saves UMAPS of these
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("normalized_pca.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("${param_details}.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        val(params__pcs)
        each n_pcs
        each batch_var

    output:
        val(outdir, emit: outdir)
        path("${runid}-${outfile}-bbknn.h5ad", emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path("${runid}-reduced_dims.tsv.gz", emit: reduced_dims)
        // NOTE: passing the param details as an unpublished file is messy,
        // but I could not get collect of ${param_details} and file to work.
        path(
            "reduced_dims-${param_details}.tsv.gz",
            emit: reduced_dims_params
        )

    script:
        runid = random_hex(16)
        param_details = "${params__pcs}-bbknn"
        param_details = "${param_details}.batch=${batch_var}"
        param_details = "${param_details}.n_pcs=${n_pcs}"
        outdir = "${outdir_prev}/reduced_dims-${param_details}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "bbknn: ${process_info}"
        echo "publish_directory: ${outdir}"
        0045-bbknn.py \
            --h5_anndata ${file__anndata} \
            --batch_key ${batch_var} \
            --n_pcs ${n_pcs} \
            --output_file ${runid}-${outfile}-bbknn
        cp ${runid}-${outfile}-bbknn-reduced_dims.tsv.gz \
            reduced_dims-${param_details}.tsv.gz
        mv ${runid}-${outfile}-bbknn-reduced_dims.tsv.gz \
            ${runid}-reduced_dims.tsv.gz
        """
}


process lisi {
    // Takes a list of reduced_dims and calculates lisi
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("normalized_pca.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__metadata)
        val(variables)
        file(file__reduced_dims)
        //tuple(val(label__reduced_dims), file(file__reduced_dims))

    output:
        val(outdir, emit: outdir)
        path(file__metadata, emit: metadata)
        path("${runid}-${outfile}-lisi.tsv.gz", emit: clusters)
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__metadata job.
        outfile = "${file__metadata}".minus(".tsv.gz")
            .split("-").drop(1).join("-")
        file__reduced_dims = file__reduced_dims.join("::")
        label__reduced_dims = file__reduced_dims
            .replaceAll("reduced_dims-", "")
            .replaceAll(".tsv.gz", "")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "lisi: ${process_info}"
        echo "publish_directory: ${outdir}"
        sleep 5m
        rm -fr plots
        0047-lisi.py \
            --reduced_dims_tsv ${file__reduced_dims} \
            --reduced_dims_tsv_labels ${label__reduced_dims} \
            --metadata_tsv ${file__metadata} \
            --metadata_columns ${variables} \
            --perplexity 30 \
            --output_file ${runid}-${outfile}-lisi
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}
