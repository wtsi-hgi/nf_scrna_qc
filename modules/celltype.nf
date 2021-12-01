#!/usr/bin/env nextflow

process split_h5ad_by_batch {
  input:
    path(file__anndata) // anndata h5ad file

  output:
    path("${outfil_prfx}_*.h5ad", emit:files_anndata_batch)
    path(outfile, emit: file_batch_list)

  script:
    process_info = "${task.cpus} (cpus), ${task.memory} (memory)"
    outfil_prfx = "${file__anndata}".minus(".h5ad")
    outfile = "${outfil_prfx}".plus("_files.txt")
    """
    echo "split_h5ad_by_batch: ${process_info}"
    echo "outfile = ${outfile}"
    scanpy_split_h5ad.py ${file__anndata} ${outfil_prfx}
    """
}

process assign_cell_types_azimuth {

  publishDir  path: "${outdir}",
              saveAs: {filename -> "${outfil_prfx}_" + filename},
              mode: "${task.publish_mode}",
              overwrite: "true"

  when:
    params.azimuth.run_process

  stageInMode 'copy'
  // stageInMode 'copy' is needed because SeuratDisk:::Convert()
  // generates the output file apparently from the absolute path of input file.
  // Symbolic links have the output file written to the link target directory
  // where it cannot be found by the azimuth.R script.

  input:
    val outdir_prev
    path file_h5ad_batch

  output:
    path(celltype_table, emit:predicted_celltypes)
    path "ncells_by_type_barplot.pdf"
    path "query_umap.pdf"
    path "prediction_score_umap.pdf"
    path "prediction_score_vln.pdf"
    path "mapping_score_umap.pdf"
    path "mapping_score_vln.pdf"

  script:
  process_info = "${task.cpus} (cpus), ${task.memory} (memory)"
  outdir = "${outdir_prev}/azimuth"
  // output file prefix: strip random hex number form beginning of file name
  outfil_prfx = "${file_h5ad_batch}".minus(".h5ad").split("-").drop(1).join("-")
  //outfil_prfx = "${file_h5ad_batch}".minus(".h5ad")
  celltype_table = "${outfil_prfx}_predicted_celltype_l2.tsv.gz"
  """
    echo "assign_cell_types_azimuth: ${process_info}"
    echo "input: ${file_h5ad_batch}"
    echo "outfil_prfx: ${outfil_prfx}"
    azimuth.R ./${file_h5ad_batch}
    gzip -c predicted_celltype_l2.tsv > ${celltype_table}
  """
}


workflow wf__celltype_assignment {
  take:
    outdir
    file__anndata

  main:
      split_h5ad_by_batch(
        file__anndata
      )
      split_h5ad_by_batch.out.files_anndata_batch
        .flatMap()
        .set{ch_batch_files}
      assign_cell_types_azimuth(
        outdir,
        ch_batch_files
      )

  emit:
    predicted_celltypes = assign_cell_types_azimuth.out.predicted_celltypes
}
