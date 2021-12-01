#!/usr/bin/env nextflow

process prep_handover_raw_input {
    input:
        tuple(
            val(experiment_id),
            path(file__feature_bc_matrix_h5),
        )

    output:
        path("${experiment_id}__raw_feature_bc_matrix.h5", emit: raw_feature_bc_mx)

    script:
        """
        ln --physical ${file__feature_bc_matrix_h5} ${experiment_id}__raw_feature_bc_matrix.h5
        #cp ${file__feature_bc_matrix_h5} ${experiment_id}__raw_feature_bc_matrix.h5
        """
}

process gather_handover_dataset {

  publishDir  path: "${outdir}",
              mode: "${task.publish_mode}",
              overwrite: "true"

  when:
    params.data_handover.run_process

  input:
    path(execution_trigger)
    val(outdir_prev)
    path(file__anndata_merged)
    path(file__cellranger_raw_files_table_tsv)
    path(file__cellbender_files_table_tsv)
    path(file__deconv_files_table_tsv)
    path(multiplet_calls)
    path(deconvolution_path)
    path(qc_output_dir)

  output:
    path("${subdir}/*", emit:outfiles_dataset)
    path("${subdir2}/*", emit:outfiles_dataset2)
    val(outdir, emit: outdir_dataset)

  script:
    outdir = "${outdir_prev}/handover"
    subdir = "minimal_dataset"
    subdir2='minimal_dataset_summary'
    if (multiplet_calls) {
      argstr = " --scrublet-output-dir=${qc_output_dir}/multiplet.method=scrublet"
    } else {
      argstr = ""
    }
    """
      echo "execution_trigger: ${execution_trigger}"
      echo "gather_minimal_dataset.py --cellranger-rawfiles-table=${file__cellranger_raw_files_table_tsv}"
      echo "multiplet_calls: ${multiplet_calls}"
      echo "outdir = ${outdir}"
      echo "azimuth-output-dir=${qc_output_dir}/azimuth"
      echo "deconvolution-output-dir=${deconvolution_path}/results/split_donor_h5ad"
      

      gather_minimal_dataset.py \
        --output-dir=${subdir} \
        --cellranger-rawfiles-table=${file__cellranger_raw_files_table_tsv} \
        --cellbender-files-table=${file__cellbender_files_table_tsv} \
        --deconvolution-files-table=${file__deconv_files_table_tsv} \
        --deconvolution-output-dir=${deconvolution_path}/results/split_donor_h5ad \
        --azimuth-output-dir=${qc_output_dir}/azimuth \
        --qc-merged-h5ad=${file__anndata_merged} ${argstr} \
    """
}



workflow wf__handover_dataset {
  take:
    predicted_celltypes
    outdir
    file__anndata_merged
    file__cellranger_raw_files_table_tsv
    file__cellbender_files_table_tsv
    file__deconv_files_table_tsv
    multiplet_calls
    deconvolution_path
    qc_output_dir

  main:
    gather_handover_dataset(
      predicted_celltypes,
      outdir,
      file__anndata_merged,
      file__cellranger_raw_files_table_tsv,
      file__cellbender_files_table_tsv,
      file__deconv_files_table_tsv,
      multiplet_calls,
      deconvolution_path,
      qc_output_dir
    )
}
