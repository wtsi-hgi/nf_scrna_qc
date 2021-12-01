#!/usr/bin/env nextflow


def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}


process run_scrublet {
    // Runs scrublet for each sample.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("multiplet_calls_published.txt")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "${task.publish_mode}",
                overwrite: "true"

    //each smplid_and_datapath
    input:
        val(outdir_prev)
        tuple(
            val(experiment_id),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix)
        )
        val(expected_multiplet_rate)
        val(n_simulated_multiplet)
        val(multiplet_threshold_method)
        val(scale_log10)

    output:
        val(outdir, emit: outdir)
        val(experiment_id, emit: experiment_id)
        path("${runid}-${outfile}-scrublet.tsv.gz", emit: multiplet_calls)
        path(
            "${runid}-${outfile}-multiplet_calls_published.txt",
            emit: multiplet_calls_published
        )
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/multiplet"
        outdir = "${outdir}.method=scrublet"
        outfile = "${experiment_id}"

        // Check to see if we should use use log10 of the doublet simulations
        // to derive the threshold
        cmd__scale_log10 = ""
        if (scale_log10 == "True") {
            cmd__scale_log10 = "--scale_log10"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "run_scrublet: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        TMP_DIR=\$(mktemp -d -p \$(pwd))
        ln --physical ${file_10x_barcodes} \$TMP_DIR
        ln --physical ${file_10x_features} \$TMP_DIR
        ln --physical ${file_10x_matrix} \$TMP_DIR
        0015-run_scrublet.py \
            --tenxdata_dir \$TMP_DIR \
            --expected_multiplet_rate ${expected_multiplet_rate} \
            --n_simulated_multiplet ${n_simulated_multiplet} \
            --multiplet_threshold_method ${multiplet_threshold_method} \
            ${cmd__scale_log10} \
            --output_file ${runid}-${outfile}
        echo -e "${experiment_id}\t${outdir}/${outfile}-scrublet.tsv.gz" > \
            ${runid}-${outfile}-multiplet_calls_published.txt
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process make_cellmetadata_pipeline_input {
    // Makes a input tsv file for the main pipeline.
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
        path("*multiplet_calls_published.txt")

    output:
        val(outdir, emit: outdir)
        path('file_cellmetadata.tsv', emit: file__cellmetadata)

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "make_pipeline_input_file: ${process_info}"
        echo "publish_directory: ${outdir}"
        # Note: the default paste delim is tab
        cat *multiplet_calls_published.txt \
            | awk 'BEGIN{print "experiment_id\tdata_path_cellmetadata"}1' \
            > file_cellmetadata.tsv
        """
}


workflow wf__multiplet {
    take:
        output_dir
        channel__file_paths_10x
        expected_multiplet_rate
        n_simulated_multiplet
        multiplet_threshold_method
        scale_log10
    main:
        // Identify multiplets using scrublet.
        run_scrublet(
            output_dir,
            channel__file_paths_10x,
            expected_multiplet_rate,
            n_simulated_multiplet,
            multiplet_threshold_method,
            scale_log10
        )
        // Generate input file for merge based in multiplets
        make_cellmetadata_pipeline_input(
            output_dir,
            run_scrublet.out.multiplet_calls_published.collect()
        )
    emit:
        // Return merged input data file.
        outdir = make_cellmetadata_pipeline_input.out.outdir
        file__cellmetadata = make_cellmetadata_pipeline_input.out.file__cellmetadata
        multiplet_calls = run_scrublet.out.multiplet_calls
}
