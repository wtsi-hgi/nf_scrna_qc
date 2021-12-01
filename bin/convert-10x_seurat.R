#!/usr/bin/env Rscript

# disable strings as factors, but re-enable upon exit
old <- options(stringsAsFactors = FALSE)
on.exit(options(old), add = TRUE)

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("Seurat"))


#' Command line interface wrapper
#'
#' @importFrom optparse make_option
#' @importFrom optparse OptionParser
#' @importFrom optparse parse_args
#' @importFrom data.table fread
#' @importFrom Seurat Read10X
#' @importFrom Seurat CreateSeuratObject
command_line_interface <- function() {
    optionList <- list(
        optparse::make_option(c("--in_dir"),
            type = "character",
            help = paste0(
                "Input directory containing barcodes.tsv.gz, features.tsv.gz,",
                " and matrix.mtx.gz."
            )
        ),

        optparse::make_option(c("--metadata_file"),
            type = "character",
            default = "",
            help = paste0(
              "Tab-delimited metadata file. Should be a dataframe where the",
              " rows are cell names and the columns are additional metadata",
              " fields. Should contain a cell_barcode column. NOTE: should",
              " be full path to file.",
              " [default: None]"
            )
        ),

        optparse::make_option(c("--count_matrix_file"),
            type = "character",
            default = "",
            help = paste0(
              "If matrix.mtx.gz is not the count matrix and count_matrix_file",
              " is provided, counts slot will be updated to properly reflect",
              " the count matrix. Otherwise the count matrix in the count",
              " slot will be wrong. Note barcodes.tsv.gz and features.tsv.gz",
              " from the main input must also match up to this file.",
              " NOTE: should be full path to file.",
              " [default: None]"
            )
        ),

        optparse::make_option(c("--out_file"),
            type = "character",
            default = "sc_df",
            help = paste0(
                "Name (and possibly path) of output file. Will have .rds.gz",
                " appended to it.",
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
            "Makes Seurat data object with metadata on samples and cells."
        )
    )

    # A hack to fix a bug in optparse that won't let you use positional args
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

    message(paste0(
        "WARNING: Unless ",
        param[["in_dir"]],
        "/matrix.mtx.gz is a counts matrix or --count_matrix_file is",
        " specified, then 'counts' slot in the Seurat data object",
        " will not actually be counts."
    ))

    message(paste0(
        "NOTE: If a \"scan() expected 'an integer', got <int> is thrown\", it",
        " is likely because you have run out of memory."
    ))

    # Read in the data.
    in_dir <- param[["in_dir"]]
    df_tmp <- Seurat::Read10X(
        in_dir,
        gene.column = 2  # 1 = ensembl_ids, 2 = gene_symbols
    )

    metadata_file <- param[["metadata_file"]]
    if (metadata_file != "") {
        df_meta <- data.frame(data.table::fread(
            cmd = paste("gunzip -c", metadata_file),
            sep = "\t",
            header = T,
            stringsAsFactors = F
        ))
        rownames(df_meta) <- df_meta[["cell_barcode"]]

        # Initalize the Seurat data object.
        df <- Seurat::CreateSeuratObject(
            df_tmp,
            min.cells = 0,
            min.features = 0,
            assay = "RNA",
            meta.data = df_meta[unlist(df_tmp@Dimnames[2]),]
        )
    } else {
        # Initalize the Seurat data object.
        df <- Seurat::CreateSeuratObject(
            df_tmp,
            min.cells = 0,
            min.features = 0,
            assay = "RNA"
        )
    }

    # Optionally add in the counts data.
    count_matrix_file <- param[["count_matrix_file"]]
    if (count_matrix_file != "") {
        # Make a temp directory.
        dir_tmp <- tempdir(check = TRUE)

        # Symlink all the files.
        for (i in c("barcodes.tsv.gz", "features.tsv.gz")) {
            unlink(paste0(dir_tmp, "/", i))
            file.symlink(
                from = paste0(in_dir, "/", i),
                to = paste0(dir_tmp, "/", i)
            )
        }
        unlink(paste0(dir_tmp, "/matrix.mtx.gz"))
        file.symlink(
            from = count_matrix_file,
            to = paste0(dir_tmp, "/matrix.mtx.gz")
        )
        # print(dir_tmp)
        # print(list.files(dir_tmp))

        # Read in the data.
        df_tmp2 <- Seurat::Read10X(
            dir_tmp,
            gene.column = 2  # 1 = ensembl_ids, 2 = gene_symbols
        )

        # Initalize the Seurat data object.
        df2 <- Seurat::CreateSeuratObject(
            df_tmp2,
            min.cells = 0,
            min.features = 0,
            assay = "RNA"
        )

        # Update the counts slot of the main matrix.
        # mtx_new_counts <- df2$RNA@counts
        # mtx_new_counts[
        #     which(rownames(mtx_new_counts) %in% rownames(df$RNA@counts)),
        #     which(colnames(mtx_new_counts) %in% colnames(df$RNA@counts))
        # ]
        mtx_new_counts <- df2$RNA@counts[
            rownames(df$RNA@counts),
            colnames(df$RNA@counts)
        ]
        df$RNA@counts <- mtx_new_counts
    }

    print(str(df))

    # Save the Seurat data object.
    saveRDS(
        df,
        file = paste0(param[["out_file"]], ".rds.gz"),
        compress = TRUE
    )

    return(0)
}


main <- function() {
    SCRIPT_NAME <- "convert-10x_seurat.R"

    # run analysis
    run_time <- system.time(df_results <- command_line_interface())
    message(paste0(
        "Analysis execution time", " [", SCRIPT_NAME, "]:\t",
        run_time[["elapsed"]]/3600, # proc.time sec to hours
        " hours."
    ))
    #return(0)
}


# like python if __name__ == '__main__'
if (sys.nframe() == 0) {
    main()
}
