#!/usr/bin/env Rscript

SCRIPT_NAME <- "harmony_process_pcs.R"

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("harmony"))


#' Command line interface wrapper
#'
#' @importFrom optparse make_option
#' @importFrom optparse OptionParser
#' @importFrom optparse parse_args
command_line_interface <- function() {
    optionList <- list(
        optparse::make_option(c("-f", "--pca_file"),
            type = "character",
            help = paste0(
                "Tab-delimited input file of PCs. Columns: 'cell_barcode'",
                " followed by PC1,PC2,etc."
            )
        ),

        optparse::make_option(c("--metadata_file"),
            type = "character",
            help = paste0(
              "Tab-delimited metadata file, must have a column labelled",
              " cell_barcode that maps to pca_file."
            )
        ),

        optparse::make_option(c("--metadata_columns"),
            type = "character",
            help = paste0(
              "Comma separated string of columns to use in metadata_file."
            )
        ),

        optparse::make_option(c("--theta"),
            type = "character",
            default = "",
            help = paste0(
                "Comma separated string of theta values (corresponding to",
                " metadata_columns). If '' then sets theta to 2 for all",
                " columns. Larger values of theta result in more diverse",
                " clusters.",
                " [default: %default]"
            )
        ),

        optparse::make_option(c("--n_pcs"),
            type = "numeric",
            default = 0,
            help = paste0(
                "Number of PCs to use from pca_file. If 0 then use all.",
                " [default: %default]"
            )
        ),


        optparse::make_option(c("--out_file"),
            type = "character",
            default = "",
            help = paste0(
                "Name (and possibly path) of output file. Will have tsv.gz",
                " appended to it. If '' then add '-harmony.tsv.gz' to",
                " pca_file.",
                " [default: %default]"
            )
        )

        # optparse::make_option(c("--verbose"),
        #     type = "logical",
        #     action = "store_true",
        #     default = FALSE,
        #     help = paste0(
        #         "Verbose mode (write extra info to std.err).",
        #         " [default: %default]"
        #     )
        # )
    )

    parser <- optparse::OptionParser(
        usage = "%prog",
        option_list = optionList,
        description = paste0(
            "Runs Harmony using a PC and metadata file from scanpy."
        )
    )

    # a hack to fix a bug in optparse that won"t let you use positional args
    # if you also have non-boolean optional args:
    getOptionStrings <- function(parserObj) {
        optionStrings <- character()
        for (item in parserObj@options) {
            optionStrings <- append(optionStrings,
                                    c(item@short_flag, item@long_flag))
        }
        optionStrings
    }
    optStrings <- getOptionStrings(parser)
    arguments <- optparse::parse_args(parser, positional_arguments = TRUE)

    # read in the parameters
    param <- list()
    for (i in names(arguments$options)) {
        param[[i]] <- arguments$options[[i]]
    }

    # Read in the PCA file.
    f_pca <- param[["pca_file"]]
    mtx_pca <- data.table::fread(
        cmd = paste("gunzip -c", f_pca),
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE
    )
    if (!("cell_barcode" %in% colnames(mtx_pca))) {
        stop("Missing cell_barcode column.")
    }
    mtx_pca <- as.matrix(mtx_pca, rownames = "cell_barcode")

    # Check the input pcs.
    n_pcs <- param[["n_pcs"]]
    if (n_pcs > length(colnames(mtx_pca))) {
        cat(
            "n_pcs =", n_pcs,
            " length(colnames(mtx_pca)) =", length(colnames(mtx_pca)), "\n"
        )
        stop("Invalid n_pcs value: n_pcs > length(colnames(mtx_pca)).")
    } else if (n_pcs < 0) {
        stop("Invalid n_pcs vaue: n_pcs < 0.")
    } else if (n_pcs == 0) {
        n_pcs <- length(colnames(mtx_pca))
    } # else n_pcs is a valid value


    # Get the metadata_file columns that we want to adjust with Harmony.
    metadata_columns <- unlist(strsplit(param[["metadata_columns"]], ","))


    # Read in the metadata file.
    f_meta <- param[["metadata_file"]]
    df_meta <- data.frame(data.table::fread(
        cmd = paste("gunzip -c", f_meta),
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE
    ))
    if (!("cell_barcode" %in% colnames(df_meta))) {
        stop("Missing cell_barcode column.")
    }
    rownames(df_meta) <- df_meta[["cell_barcode"]]
    df_meta <- df_meta[rownames(mtx_pca), metadata_columns]


    # Get the theta values for each column (if none, set to 2 for all columns).
    theta <- rep(2, length(metadata_columns))
    if (param[["theta"]] != '') {
        theta <- as.numeric(unlist(strsplit(param[["theta"]], ",")))
    }


    harmony_embeddings <- harmony::HarmonyMatrix(
        data_mat = mtx_pca[, seq(1, n_pcs)],
        meta_data = df_meta,
        do_pca = FALSE,
        # npcs = n_pcs, # only relevant if do_pca == TRUE
        verbose = TRUE,
        # epsilon.harmony = -Inf, # Set to -Inf to never stop early.
        vars_use = metadata_columns,
        theta = theta
    )
    colnames(harmony_embeddings) <- gsub(
        "PC",
        "harmony",
        colnames(harmony_embeddings)
    )
    out_col_order <- colnames(harmony_embeddings)
    out_col_order <- c("cell_barcode", out_col_order)
    harmony_embeddings <- as.data.frame(harmony_embeddings)
    harmony_embeddings[["cell_barcode"]] <- rownames(harmony_embeddings)


    base <- param[["out_file"]]
    if (base == "") {
        base <- paste0(
            gsub(".tsv.gz", "", param[["pca_file"]]),
            "-harmony"
        )
    }
    gzfh <- gzfile(paste(base, ".tsv.gz", sep = ""), "w", compression = 9)
    write.table(
        harmony_embeddings[out_col_order],
        gzfh,
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE,
        sep = "\t",
        na = ""
    )
    close(gzfh)

    return(0)
}


main <- function() {
    # Run analysis
    run_time <- system.time(df_results <- command_line_interface())
    execution_summary <- paste0(
        "Analysis execution time", " [", SCRIPT_NAME, "]:\t",
        run_time[["elapsed"]]/3600, # proc.time sec to hours
        " hours."
    )
    print(execution_summary)
    #return(execution_summary)
    #return(0) # For nextflow, more meaningful to return execution_summary
}


# code for development
dev <- function() {
    f_pca <- "adata-pcs.tsv.gz"
    mtx_pca <- data.table::fread(
        cmd = paste("gunzip -c", f_pca),
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE
    )
    mtx_pca <- as.matrix(mtx_pca, rownames = "cell_barcode")

    f_meta <- "adata-metadata.tsv.gz"
    df_meta <- data.frame(data.table::fread(
        cmd = paste("gunzip -c", f_meta),
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE
    ))
    rownames(df_meta) <- df_meta[["cell_barcode"]]
    df_meta <- df_meta[rownames(mtx_pca), c("experiment_id", "bead_version")]

    n_pcs <- 10
    harmony_embeddings <- harmony::HarmonyMatrix(
        data_mat = mtx_pca[, seq(1, n_pcs)],
        meta_data = df_meta,
        do_pca = FALSE,
        #npcs = n_pcs,
        verbose = TRUE,
        # epsilon.harmony = -Inf, # Set to -Inf to never stop early.
        # vars_use = c("experiment_id"),
        # theta = c(1)
        vars_use = c("experiment_id", "bead_version"),
        theta = c(1, 0.2)
    )
    colnames(harmony_embeddings) <- gsub(
        "PC",
        "harmony",
        colnames(harmony_embeddings)
    )
    out_col_order <- colnames(harmony_embeddings)
    out_col_order <- c("cell_barcode", out_col_order)
    harmony_embeddings <- as.data.frame(harmony_embeddings)
    harmony_embeddings[["cell_barcode"]] <- rownames(harmony_embeddings)

    base <- "harmony_dev"
    gzfh = gzfile(paste(base, ".tsv.gz", sep = ""), "w", compression = 9)
    write.table(
        harmony_embeddings[out_col_order],
        gzfh,
        row.names = FALSE,
        col.names = TRUE,
        quote = TRUE,
        sep = "\t",
        na = ""
    )
    close(gzfh)
}


# like python if __name__ == '__main__'
if (sys.nframe() == 0) {
    #dev()
    main()
}
