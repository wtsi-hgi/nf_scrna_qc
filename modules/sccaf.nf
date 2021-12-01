#!/usr/bin/env nextflow

process sccaf_assess_clustering {
  publishDir path: "${outdir}",
             mode: "${task.publish_mode}",
             overwrite: "true"

  when:
    params.sccaf.run_assessment

  input:
    val(outdir_prev)
    path(file__anndata) // anndata h5ad file
    path(file__external_clustering_tsv)

  output:
    path("${outfile_roc_pdf}")
    path("${outfile_prc_pdf}")
    path("${outfile_acc_txt}")

  script:
    outdir = "${outdir_prev}/sccaf/clustering_assessment"
    process_info = "${task.cpus} (cpus), ${task.memory} (memory)"
    outfil_prfx = "${file__external_clustering_tsv}".minus(".tsv.gz").plus("_assess")
    outfile_roc_pdf = outfil_prfx.plus("_clust_roc.pdf")
    outfile_prc_pdf = outfil_prfx.plus("_clust_prc.pdf")
    outfile_acc_txt = outfil_prfx.plus("_acc.txt")
    """
    echo "sccaf_assess_clustering: ${process_info}"
    echo "publish_directory: ${outdir}"
    echo "params.sccaf.run_assessment = ${params.sccaf.run_assessment}"

    sccaf_assess_clustering.py \\
       --nthreads ${task.cpus} \\
       --use-pca \\
       --output-prefix ${outfil_prfx} \\
       ${file__anndata} \\
       ${file__external_clustering_tsv}
    """
}

process sccaf_optimize_clustering {
  publishDir path: "${outdir}",
             saveAs: {filename -> outfil_prfx.plus(filename)},
             mode: "${task.publish_mode}",
             overwrite: "true"
             //saveAs: {filename -> outfil_prfx.plus(filename)},

  when:
    params.sccaf.run_optimization

  input:
    val(outdir_prev)
    path(file__anndata)
    path(file__external_clustering_tsv)
    val(min_accuracy)

    output:
      path("${outfile}")
      path("roc-curve.png")
      path("optim.png")
      path("rounds.txt")

  script:
    outdir = "${outdir_prev}/sccaf/clustering_optimization"
    outdir_rel = "optim_out"
    process_info = "${task.cpus} (cpus), ${task.memory} (memory)"
    outfil_prfx = "${file__external_clustering_tsv}".minus(".tsv.gz").plus("_optim_")
    outfile = outfil_prfx.plus("clustopt.h5ad")
    """
    echo "sccaf_optimize_clustering: ${process_info}"
    echo "publish_directory: ${outdir}"
    echo "params.sccaf.run_optimization = ${params.sccaf.run_optimization}"
    echo "min_accuracy = ${min_accuracy}"

    sccaf \
          --output-file ${outfile} \
          --cores ${task.cpus} \
          --input-file ${file__anndata} \
          --external-clustering-tsv ${file__external_clustering_tsv} \
          --use-pca \
          --optimise \
          --min-accuracy ${min_accuracy} \
          --produce-rounds-summary
    #      --optimisation-plots-output ${outdir_rel}
    """
}

workflow wf__optimize_clustering {
  take:
    outdir
    anndata
    external_clustering
    min_accuracy

  main:
    sccaf_optimize_clustering(outdir, anndata, exernal_clustering, min_accuracy)
}
