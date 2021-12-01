#!/usr/bin/env nextflow


include {
    umap_calculate_and_plot as umap_calculate_and_plot__cluster;
} from "./umap.nf"

include {
    sccaf_assess_clustering
    sccaf_optimize_clustering
} from "./sccaf.nf"

def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}


process cluster {
    // Clusters results.
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
                    } else if(filename.endsWith("reduced_dims.tsv.gz")) {
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
        path(file__reduced_dims)
        each number_neighbors
        each method
        each resolution
        // tuple(val(outdir_prev), path(file__reduced_dims))

    output:
        val(outdir, emit: outdir)
        path("${runid}-${outfile}-clustered.h5ad", emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path(file__reduced_dims, emit: reduced_dims)
        path("${runid}-${outfile}-clustered.tsv.gz", emit: clusters)
        val(outdir_prev, emit: outdir__reduced_dims)
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true

    script:
        runid = random_hex(16)
        resolution_str = "${resolution}" //.replaceAll("\\.", "pt")
        outdir = "${outdir_prev}/cluster"
        outdir = "${outdir}.number_neighbors=${number_neighbors}"
        outdir = "${outdir}.method=${method}"
        outdir = "${outdir}.resolution=${resolution_str}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cluster: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        0053-scanpy_cluster.py \
            --h5_anndata ${file__anndata} \
            --tsv_pcs ${file__reduced_dims} \
            --number_neighbors ${number_neighbors} \
            --cluster_method ${method} \
            --resolution ${resolution} \
            --number_cpu ${task.cpus} \
            --output_file ${runid}-${outfile}-clustered
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process plot_phenotype_across_clusters {
    // Takes annData object, plots distribution of obs value across clusters
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
        each variables

    output:
        val(outdir, emit: outdir)
        path("plots/*.png")
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        // Append run_id to output file.
        outfile = "${runid}-${outfile}-cluster_boxplot"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_phenotype_across_clusters: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        0055-plot_anndataobs_across_clusters.py \
            --h5_anndata ${file__anndata} \
            --pheno_columns ${variables} \
            --output_file ${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process serialize_known_markers {
    // Serializes known markers for analysis
    // ------------------------------------------------------------------------
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    input:
        tuple(
            val(file_id),
            file(marker_file)
        )

    output:
        path("${outfile}", emit: marker_file)

    script:
        runid = random_hex(16)
        outfile = "${file_id}"
        cmd__run = "touch NO_FILE.tsv"
        if (file_id != "") {
            outfile = "${file_id}.tsv"
            cmd__run = "ln --physical ${marker_file} ${outfile}"
        } else {
            outfile = "NO_FILE"
            cmd__run = "touch ${outfile}"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "serialize_known_markers: ${process_info}"
        ${cmd__run}
        """
}


process plot_known_markers {
    // Plots markers from previous studies as dotplots
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    label 'long_job'     // use GPU
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    //saveAs: {filename -> filename.replaceAll("${runid}-", "")},
    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        each file(marker_file)

    output:
        val(outdir, emit: outdir)
        path("plots_known_markers/*.pdf") optional true
        path("plots_known_markers/*.png") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        // outfile = "${file__anndata}".minus(".h5ad")
        //     .split("-").drop(1).join("-")
        outfile = "${marker_file.name}".minus(".tsv")
        // Only run this script if there is a value
        cmd__run = ""
        if (outfile != "NO_FILE") {
            cmd__run = "0055-plot_known_markers.py"
            cmd__run = "${cmd__run} --h5_anndata ${file__anndata}"
            cmd__run = "${cmd__run} --markers_database ${marker_file}"
            cmd__run = "${cmd__run} --output_file ${outfile}"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_marker_dotplot: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots_known_markers
        ${cmd__run}
        mkdir plots_known_markers
        mv *pdf plots_known_markers/ 2>/dev/null || true
        mv *png plots_known_markers/ 2>/dev/null || true
        """
}


process cluster_validate_resolution_keras {
    // Validate the resolution for clusters.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    //maxForks 2         // hard to control memory usage. limit to 3 concurrent
    label 'gpu'        // use GPU
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    //saveAs: {filename -> filename.replaceAll("${runid}-", "")},
    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("clustered.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("reduced_dims.tsv.gz")) {
                        null
                    } else if(filename.endsWith("clustered.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "${task.publish_mode}",
                overwrite: "true"

    when:
      params.cluster_validate_resolution.run_process

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        path(file__reduced_dims)
        path(file__clusters)
        each sparsity
        each train_size_cells
        val(outdir__reduced_dims)

    output:
        val(outdir, emit: outdir)
        path(file__anndata, emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path(file__reduced_dims, emit: reduced_dims)
        path(file__clusters, emit: clusters)
        path("${runid}-${outfile}.h5", emit: model)
        path("${runid}-${outfile}.yml", emit: model_yaml)
        path("${runid}-${outfile}-weights.h5", emit: model_weights)
        path("${runid}-${outfile}-model_report.tsv.gz", emit: model_report)
        path(
            "${runid}-${outfile}-test_result.tsv.gz",
            emit: model_test_result
        )
        path(
            "${runid}-${outfile}-weights.tsv.gz",
            emit: model_weights_tsv
        )
        tuple(
            val("${outdir__reduced_dims}"),
            file("${runid}-${outfile}-model_report.tsv.gz",),
            file("${runid}-${outfile}-test_result.tsv.gz",),
            emit: plot_input
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/validate_resolution"
        // outdir = "${outdir}.method=${method}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        // Add sparsity information.
        sparsity_txt = "${sparsity}".replaceAll("\\.", "pt")
        outfile = "${outfile}-sparsity_l1=${sparsity_txt}"
        // Add info on number of cells for training
        // cmd__train_cells = ""
        // if (train_size_cells > 0) {
        //     train_size_cells_txt = "${train_size_cells}"
        //     cmd__train_cells = "--train_size_cells ${train_size_cells}"
        // }
        outfile = "${outfile}-train_size_cells=${train_size_cells}"
        // Job info
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        tf_memory = "${task.memory}".replaceAll(" GB", "")
        """
        echo "cluster_validate_resolution: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        0057-scanpy_cluster_validate_resolution-keras.py \
            --h5_anndata ${file__anndata} \
            --sparsity_l1 ${sparsity} \
            --number_epoch 25 \
            --batch_size 32 \
            --train_size_cells ${train_size_cells} \
            --memory_limit ${tf_memory} \
            --output_file ${runid}-${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process plot_resolution_validate {
    // Plot the AUC from validation models across resolutions
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    //saveAs: {filename -> filename.replaceAll("${runid}-", "")},
    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("clustered.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("reduced_dims.tsv.gz")) {
                        null
                    } else if(filename.endsWith("clustered.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "${task.publish_mode}",
                overwrite: "true"

    when:
      params.cluster_validate_resolution.run_process

    input:
        tuple(
            val(outdir_prev),
            path(files__model_report),
            path(files__y_prob_df)
        )
        // val(outdir_prev)
        // path(files__model_report)
        // path(files__y_prob_df)

    output:
        val(outdir, emit: outdir)
        path(
            "${runid}-${outfile}-merged_model_report.tsv.gz",
            emit: merged_model_report
        )
        path(
            "${runid}-${outfile}-merged_test_result.tsv.gz",
            emit: merged_test_result
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        files__model_report = files__model_report.join('::')
        files__y_prob_df = files__y_prob_df.join('::')
        outfile = "resolution_tuning"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_resolution: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        0058-plot_resolution_boxplot.py \
            --model_reports ${files__model_report} \
            --h_line 0.75 \
            --output_file ${runid}-${outfile}-hline0pt8
        0058-plot_resolution_boxplot.py \
            --model_reports ${files__model_report} \
            --output_file ${runid}-${outfile}
        0058-plot_resolution_curve.py \
             --y_prob_dfs ${files__y_prob_df} \
             --output_file ${runid}-${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process cluster_markers {
    // Find markers for clusters.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    //saveAs: {filename -> filename.replaceAll("${runid}-", "")},
    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("clustered.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("reduced_dims.tsv.gz")) {
                        null
                    } else if(filename.endsWith("clustered.tsv.gz")) {
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
        path(file__reduced_dims)
        path(file__clusters)
        each method

    output:
        val(outdir, emit: outdir)
        path(file__anndata, emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path(file__reduced_dims, emit: reduced_dims)
        path(file__clusters, emit: clusters)
        path(
            "${runid}-${outfile}-cluster_markers.tsv.gz",
            emit: cluster_markers
        )
        path(
            "${runid}-${outfile}-cluster_markers-filter__*.tsv.gz"
        ) optional true
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/cluster_markers"
        outdir = "${outdir}.method=${method}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        // Only run this script if there is a value
        cmd__enrich = ""
        // if (method == "wilcoxon") {
        //     f = "${runid}-${outfile}-cluster_markers"
        //     cmd__enrich = "0056-marker_enrichmet.R"
        //     cmd__enrich = "${cmd__enrich} --markers_table ${f}.tsv.gz;"
        //     cmd__enrich = "${cmd__enrich} 0056-marker_enrichmet.R"
        //     cmd__enrich = "${cmd__enrich} --markers_table ${f}-filter__fdr0pt05.tsv.gz;"
        //     cmd__enrich = "${cmd__enrich} 0056-marker_enrichmet.R"
        //     cmd__enrich = "${cmd__enrich} --markers_table ${f}-filter__fdr0pt05__user_variable_genes_exclude.tsv.gz;"
        // }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cluster: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        0056-scanpy_cluster_markers.py \
            --h5_anndata ${file__anndata} \
            --rank_genes_method ${method} \
            --number_cpu ${task.cpus} \
            --output_file ${runid}-${outfile}
        ${cmd__enrich}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

process cellex_cluster_markers {
    // Find markers for clusters using CELLEX.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    label 'long_job'     // use GPU
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    //saveAs: {filename -> filename.replaceAll("${runid}-", "")},
    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)

    output:
        path(
            "${runid}-${outfile}-essd.tsv.gz"
        )
        path(
            "${runid}-${outfile}-esmu*.tsv.gz"
        )

    script:
    // The script generates two files, adata-normalized_pca-bbknn-umap-clustered.esmu.csv.gz
    // and adata-normalized_pca-bbknn-umap-clustered.essd.csv.gz.
    // These are two matrices in the format genes x clusters. The first file, esmu.csv.gz,
    // contains the expression specificity scores for each gene within each cluster (scale from 0 to 1),
    // while the second file essd.csv.gz comprises the respective standard deviation for each gene and cluster.
        runid = random_hex(16)
        outdir = "${outdir_prev}/cluster_markers.method=cellex"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        // Only run this script if there is a value
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "cluster: ${process_info}"
        echo "publish_directory: ${outdir}"
        0060-cellex_cluster_markers.py \
            --h5_anndata ${file__anndata} \
            --output_file ${runid}-${outfile} \
            --verbose False
        """
}

process merge_clusters {
    // Merges clusters.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    //saveAs: {filename -> filename.replaceAll("${runid}-", "")},
    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        each maximum_de
        each auc_difference

    output:
        val(outdir, emit: outdir)
        path("${runid}-${outfile}-merged_clusters.h5ad", emit: anndata)
        path("${runid}-${outfile}-merged_clusters.tsv.gz", emit: clusters)
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true

    script:
        runid = random_hex(16)
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        outdir = "${outdir_prev}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "merge_clusters: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        0059-h5ad_to_h5.py \
            --h5_anndata ${file__anndata} \
            --output_file ${runid}-temp
        0059-seurat_cluster_merge.R \
            --input_file ${runid}-temp.h5 \
            --output_file_basename ${runid}-${outfile}-merged_clusters \
            --maximum_de ${maximum_de} \
            --auc_difference ${auc_difference}
        add_tsv_anndata_obs.py \
            --h5_anndata ${file__anndata} \
            --tsv_file ${runid}-${outfile}-merged_clusters.tsv.gz \
            --out_file ${runid}-${outfile}-merged_clusters
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process prep_cellxgene {
    // Preps adata file for cellxgene
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false           // cache results from run
    scratch false           // use tmp directory
    echo echo_mode          // echo output from script

    //saveAs: {filename -> filename.replaceAll("${runid}-", "")},
    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("clustered.h5ad")) {
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

    output:
        val(outdir, emit: outdir)
        path(
            "${runid}-${outfile}-cellxgene.h5ad",
            emit: cluster_markers
        )

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "prep_cellxgene: ${process_info}"
        echo "publish_directory: ${outdir}"
        cellxgene.py \
            --h5_anndata ${file__anndata} \
            --drop_extra_info \
            --output_file ${runid}-${outfile}-cellxgene
        """
}


process convert_seurat {
    // Converts anndata h5 file to a Seurat data object.
    // TODO: automatically add reduced_dims to Seurat data object.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    //saveAs: {filename -> filename.replaceAll("${runid}-", "")},
    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    when:
      params.cluster.convert_seurat

    input:
        val(outdir_prev)
        path(file__anndata)

    output:
        val(outdir, emit: outdir)
        path("${runid}-${outfile}.rds.gz", emit: seurat_data)
        path("${outdir_relative}/*", emit: matrix_data)
        // tuple(path("${outdir_relative}/*"), emit: matrix_data)

    script:
        runid = random_hex(16)
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        outdir_relative = "${runid}-matrices-${outfile}"
        outdir = "${outdir_prev}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "convert_seurat: ${process_info}"
        echo "publish_directory: ${outdir}"
        outdir_relative_full_path=\$(pwd)"/${outdir_relative}"
        convert-anndata_10x.py \
            --h5_anndata ${file__anndata} \
            --output_dir ${outdir_relative}
        ln -s \${outdir_relative_full_path}/matrix-X.mtx.gz ${outdir_relative}/matrix.mtx.gz
        convert-10x_seurat.R \
            --in_dir \${outdir_relative_full_path} \
            --metadata_file \${outdir_relative_full_path}/metadata-barcodes.tsv.gz \
            --count_matrix_file \${outdir_relative_full_path}/matrix-counts.mtx.gz \
            --out_file ${runid}-${outfile}
        unlink ${outdir_relative}/matrix.mtx.gz
        """
}


workflow wf__cluster {
    take:
        outdir
        anndata
        metadata
        pcs
        reduced_dims
        use_pcs_as_reduced_dims
        cluster__number_neighbors
        cluster__methods
        cluster__resolutions
        cluster__boxplot_variables
        cluster__known_markers
        cluster_validate_resolution__sparsity
        cluster_validate_resolution__train_size_cells
        cluster_marker__methods
        n_neighbors
        umap_init
        umap_min_dist
        umap_spread
        sccaf_minacc
    main:
        // Cluster the results, varying the resolution ------------------------
        cluster(
            outdir,
            anndata,
            metadata,
            pcs,
            reduced_dims,
            cluster__number_neighbors,
            cluster__methods,
            cluster__resolutions
        )
        // Boxplot of phenotype across clusters -------------------------------
        plot_phenotype_across_clusters(
            cluster.out.outdir,
            cluster.out.anndata,
            cluster__boxplot_variables
        )
        // Dotplot of marker genes across clusters ----------------------------
        // Serialize marker genes
        serialize_known_markers(
            cluster__known_markers
        )
        plot_known_markers(
            cluster.out.outdir,
            cluster.out.anndata,
            serialize_known_markers.out.marker_file
        )
        // Validate the resolution
            cluster_validate_resolution_keras(  //#mo11
             cluster.out.outdir,
             cluster.out.anndata,
             cluster.out.metadata,
             cluster.out.pcs,
             cluster.out.reduced_dims,
             cluster.out.clusters,
             cluster_validate_resolution__sparsity, // "0.0001",
             cluster_validate_resolution__train_size_cells,
             cluster.out.outdir__reduced_dims
         )

        sccaf_assess_clustering(
          cluster.out.outdir,
          cluster.out.anndata,
          cluster.out.clusters
        )

        sccaf_optimize_clustering(
           cluster.out.outdir,
           cluster.out.anndata,
           cluster.out.clusters,
           sccaf_minacc
         ) // - mo11


        // Plot the AUC across the resolutions
        // NOTE: cannot just run a collect() in output, because there might
        // not be a unique call - e.g. harmony with multiple theta
        plot_resolution_validate(
            cluster_validate_resolution_keras.out.plot_input.groupTuple()
         ) //#mo11
        // Make Seurat dataframes of the clustered anndata
        convert_seurat(
            cluster.out.outdir,
            cluster.out.anndata
        )
        // Generate UMAPs of the results.
        umap_calculate_and_plot__cluster(
            cluster.out.outdir,
            cluster.out.anndata,
            cluster.out.pcs,
            cluster.out.reduced_dims,
            use_pcs_as_reduced_dims,
            "",
            "cluster",
            n_neighbors,
            umap_init,
            umap_min_dist,
            umap_spread
        )
        // Find marker genes for clusters
        cluster_markers(
            cluster.out.outdir,
            cluster.out.anndata,
            cluster.out.metadata,
            cluster.out.pcs,
            cluster.out.reduced_dims,
            cluster.out.clusters,
            cluster_marker__methods
        )
        // Find marker genes for clusters using CELLEX
        cellex_cluster_markers(
            cluster.out.outdir,
            cluster.out.anndata
        )
        // // Merge clusters
        // merge_clusters(
        //     cluster.out.outdir,
        //     cluster.out.anndata,
        //     ['5'],
        //     ['0.1']
        // )
        // Prep adata file for cellxgene website
        prep_cellxgene(
            cluster.out.outdir,
            cluster.out.anndata
        )
}
