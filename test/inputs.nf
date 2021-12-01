params {
    // --mode Specify whether you want to run "conventional" scRNA-seq analysis workflow starting with raw 10X data or
    // a "subclustering" analysis of an already existing annData object.
    // "conventional_h5ad" mode runs like "conventional" but takes list of h5ad scanpy objects rather than raw cellranger 10x data.

    mode = "conventional"
    // file_paths_10x is required also for mode == "conventional_h5ad" & run_multiplet = true
    file_paths_10x = "${projectDir}/test/data/table_10x_file_paths.tsv"
    // Tab-delimited file containing sample metadata.
    file_metadata = "${projectDir}/test/data/metadata.tsv"

    //mode = "conventional_h5ad"
    // tsv table with columns "experiment_id" "donor" and "h5ad_filepath"
    file_paths_h5ad = "${projectDir}/test/data/deconv_h5ad_files.tsv"
    // Tab-delimited file containing sample metadata.
    //file_metadata = "${projectDir}/test/data/metadata_deconv.tsv"

    //    --genes_score Tab-delimited file with genes to use to score cells.
    // Must contain ensembl_gene_id and score_id columns. If one score_id == "cell_cycle",
    // then requires a grouping_id column with "G2/M" and "S". If no filter, then pass an empty file.
    genes_score = "${projectDir}/example_input_files/genes_score_v001.tsv"

    //    --genes_exclude_hvg Tab-delimited file with genes to exclude from highly variable gene list.
    // Must contain ensembl_gene_id column. If no filter, then pass an empty file.
    genes_exclude_hvg = "${projectDir}/example_input_files/genes_remove_hvg_v001.tsv"

}
