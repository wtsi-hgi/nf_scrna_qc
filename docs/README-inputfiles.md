## Input files

**Note:** input files should be specified as absolute file paths.

*required input files:*

1. **--file_paths_10x**:  Tab-delimited file containing experiment_id and data_path_10x_format columns (i.e., list of input samples). *Reqired*.
2. **--file_metadata**:  Tab-delimited file containing sample metadata. This will automatically be subset down to the sample list from 1. *Reqired*.

*optional input files:*

3. **--file_sample_qc**:  YAML file containing sample qc and filtering parameters. Optional. NOTE: in the example config file, this is part of the YAML file for `-params-file`.
4. **--genes_exclude_hvg**:  Tab-delimited file with genes to exclude from
highly variable gene list. Must contain ensembl_gene_id column. Optional.
5. **--genes_score**:  Tab-delimited file with genes to use to score cells. Must contain ensembl_gene_id and score_idvcolumns. If one score_id == "cell_cycle", then requires a grouping_id column with "G2/M" and "S" (see example file in `example_input_files`). Optional.
6. **-params-file**:  YAML file containing analysis parameters. Optional.
7. **--run_multiplet**:  Flag to run multiplet analysis. Optional.
8. **--file_cellmetadata**:  Tab-delimited file containing experiment_id and data_path_cellmetadata columns. For instance this file can be used to pass per cell doublet annotations. Optional.
