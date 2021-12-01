#!/usr/bin/env Rscript

## adapted from https://azimuth.hubmapconsortium.org

## expects input *.h5 file as argument
## writes plots to a single PDF file Rplots.pdf

library(Seurat)
library(SeuratDisk)
library(Matrix)
library(hdf5r)
library(ggplot2)
library(tools)

# reference files are expected in the following directory
REFERENCE_DIR <- "/opt/PBMC_reference/"

#' Load file input into a \code{Seurat} object
#'
#' Take a file and load it into a \code{\link[Seurat]{Seurat}} object. Supports
#' a variety of file types and always returns a \code{Seurat} object
#'
#' @details \code{LoadFileInput} supports several file types to be read in as
#' \code{Seurat} objects. File type is determined by extension, matched in a
#' case-insensitive manner See sections below for details about supported
#' filtypes, required extension, and specifics for how data is loaded
#'
#' @param path Path to input data
#'
#' @return A \code{\link[Seurat]{Seurat}} object
#'
#' @importFrom tools file_ext
#' @importFrom SeuratDisk Connect
#' @importFrom SeuratObject CreateSeuratObject Assays GetAssayData
#' DefaultAssay<-
#' @importFrom Seurat Read10X_h5 as.sparse Assays DietSeurat as.Seurat
#'
#' @section 10X H5 File (extension \code{h5}):
#' 10X HDF5 files are supported for all versions of Cell Ranger; data is read
#' in using \code{\link[Seurat]{Read10X_h5}}. \strong{Note}: for multi-modal
#' 10X HDF5 files, only the \emph{first} matrix is read in
#'
#' @section Rds File (extension \code{rds}):
#' Rds files are supported as long as they contain one of the following data
#' types:
#' \itemize{
#'  \item A \code{\link[Seurat]{Seurat}} V3 object
#'  \item An S4 \code{\link[Matrix]{Matrix}} object
#'  \item An S3 \code{\link[base]{matrix}} object
#'  \item A \code{\link[base]{data.frame}} object
#' }
#' For S4 \code{Matrix}, S3 \code{matrix}, and \code{data.frame} objects, a
#' \code{Seurat} object will be made with
#' \code{\link[Seurat]{CreateSeuratObject}} using the default arguments
#'
#' @section h5Seurat File (extension \code{h5seurat}):
#' h5Seurat files and all of their features are fully supported. They are read
#' in via \code{\link[SeuratDisk]{LoadH5Seurat}}. \strong{Note}: only the
#' \dQuote{counts} matrices are read in and only the default assay is kept
#'
#' @inheritSection LoadH5AD AnnData H5AD File (extension \code{h5ad})
#' @export
#'
LoadFileInput <- function(path) {
  # TODO: add support for loom files
  on.exit(expr = gc(verbose = FALSE))
  type <- tolower(x = tools::file_ext(x = path))
  return(switch(
    EXPR = type,
    'h5' = {
      mat <- Read10X_h5(filename = path)
      if (is.list(x = mat)) {
        mat <- mat[[1]]
      }
      CreateSeuratObject(counts = mat, min.cells = 1, min.features = 1)
    },
    'rds' = {
      object <- readRDS(file = path)
      if (inherits(x = object, what = c('Matrix', 'matrix', 'data.frame'))) {
        object <- CreateSeuratObject(counts = as.sparse(x = object), min.cells = 1, min.features = 1)
      } else if (inherits(x = object, what = 'Seurat')) {
        if (!'RNA' %in% Assays(object = object)) {
          stop("No RNA assay provided", call. = FALSE)
        } else if (Seurat:::IsMatrixEmpty(x = GetAssayData(object = object, slot = 'counts', assay = 'RNA'))) {
          stop("No RNA counts matrix present", call. = FALSE)
        }
        object <- CreateSeuratObject(
          counts = GetAssayData(object = object[["RNA"]], slot = "counts"),
          min.cells = 1,
          min.features = 1,
          meta.data = object[[]]
        )
      } else {
        stop("The RDS file must be a Seurat object", call. = FALSE)
      }
      object
    },
    'h5ad' = LoadH5AD(path = path),
    'h5seurat' = {
      if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
        stop("Loading h5Seurat files requires SeuratDisk", call. = FALSE)
      }
      hfile <- suppressWarnings(expr = SeuratDisk::Connect(filename = path))
      on.exit(expr = hfile$close_all())
      if (!'RNA' %in% names(x = hfile[['assays']])) {
        stop("Cannot find the RNA assay in this h5Seurat file", call. = FALSE)
      } else if (!'counts' %in% names(x = hfile[['assays/RNA']])) {
        stop("No RNA counts matrix provided", call. = FALSE)
      }
      object <- as.Seurat(
        x = hfile,
        assays = list('RNA' = 'counts'),
        reductions = FALSE,
        graphs = FALSE,
        images = FALSE
      )
      object <- CreateSeuratObject(
        counts = GetAssayData(object = object[["RNA"]], slot = "counts"),
        min.cells = 1,
        min.features = 1,
        meta.data = object[[]]
      )
    },
    stop("Unknown file type: ", type, call. = FALSE)
  ))
}

#' Load a diet H5AD file
#'
#' Read in only the counts matrix and (if present) metadata of an H5AD file and
#' return a \code{Seurat} object
#'
#' @inheritParams LoadFileInput
#'
#' @return A \code{Seurat} object
#'
#' @importFrom hdf5r H5File h5attr
#' @importFrom SeuratObject AddMetaData CreateSeuratObject
#'
#' @keywords internal
#'
#' @section AnnData H5AD File (extension \code{h5ad}):
#' Only H5AD files from AnnData v0.7 or higher are supported. Data is read from
#' the H5AD file in the following manner
#' \itemize{
#'  \item The counts matrix is read from \dQuote{/raw/X}; if \dQuote{/raw/X} is
#'  not present, the matrix is read from \dQuote{/X}
#'  \item Feature names are read from feature-level metadata. Feature level
#'  metadata must be an HDF5 group, HDF5 compound datasets are \strong{not}
#'  supported. If counts are read from \code{/raw/X}, features names are looked
#'  for in \dQuote{/raw/var}; if counts are read from \dQuote{/X}, features
#'  names are looked for in \dQuote{/var}. In both cases, feature names are read
#'  from the dataset specified by the \dQuote{_index} attribute, \dQuote{_index}
#'  dataset, or \dQuote{index} dataset, in that order
#'  \item Cell names are read from cell-level metadata. Cell-level metadata must
#'  be an HDF5 group, HDF5 compound datasets are \strong{not} supported.
#'  Cell-level metadata is read from \dQuote{/obs}. Cell names are read from the
#'  dataset specified by the \dQuote{_index} attribute, \dQuote{_index} dataset,
#'  or \dQuote{index} dataset, in that order
#'  \item Cell-level metadata is read from the \dQuote{/obs} dataset. Columns
#'  will be returned in the same order as in the \dQuote{column-order}, if
#'  present, or in alphabetical order. If a dataset named \dQuote{__categories}
#'  is present, then all datasets in \dQuote{__categories} will serve as factor
#'  levels for datasets present in \dQuote{/obs} with the same name (eg. a
#'  dataset named \dQuote{/obs/__categories/leiden} will serve as the levels for
#'  \dQuote{/obs/leiden}). Row names will be set as cell names as described
#'  above. All datasets in \dQuote{/obs} will be loaded except for
#'  \dQuote{__categories} and the cell names dataset
#' }
#'
LoadH5AD <- function(path) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Loading H5AD files requires hdf5r", call. = FALSE)
  }
  adata <- hdf5r::H5File$new(filename = path, mode = 'r')
  on.exit(expr = adata$close_all())
  Exists <- function(name) {
    name <- unlist(x = strsplit(x = name[1], split = '/', fixed = TRUE))
    hpath <- character(length = 1L)
    exists <- TRUE
    for (i in seq_along(along.with = name)) {
      hpath <- paste(hpath, name[i], sep = '/')
      exists <- adata$exists(name = hpath)
      if (isFALSE(x = exists)) {
        break
      }
    }
    return(exists)
  }
  IsMetaData <- function(md) {
    return(Exists(name = md) && inherits(x = adata[[md]], what = 'H5Group'))
  }
  GetIndex <- function(md) {
    return(
      if (adata[[md]]$attr_exists(attr_name = '_index')) {
        hdf5r::h5attr(x = adata[[md]], which = '_index')
      } else if (adata[[md]]$exists(name = '_index')) {
        '_index'
      } else if (adata[[md]]$exists(name = 'index')) {
        'index'
      } else {
        stop("Cannot find the rownames for '", md, "'", call. = FALSE)
      }
    )
  }
  GetRowNames <- function(md) {
    return(adata[[md]][[GetIndex(md = md)]][])
  }
  LoadMetadata <- function(md) {
    factor.cols <- if (adata[[md]]$exists(name = '__categories')) {
      names(x = adata[[md]][['__categories']])
    } else {
      NULL
    }
    index <- GetIndex(md = md)
    col.names <- names(x = adata[[md]])
    if (adata[[md]]$attr_exists(attr_name = 'column-order')) {
      tryCatch(
        expr = {
          col.order <- hdf5r::h5attr(x = adata[[md]], which = 'column-order')
          col.names <- c(
            intersect(x = col.order, y = col.names),
            setdiff(x = col.names, y = col.order)
          )
        },
        error = function(...) {
          return(invisible(x = NULL))
        }
      )
    }
    col.names <- col.names[!col.names %in% c('__categories', index)]
    df <- sapply(
      X = col.names,
      FUN = function(i) {
        x <- adata[[md]][[i]][]
        if (i %in% factor.cols) {
          x <- factor(x = x, levels = adata[[md]][['__categories']][[i]][])
        }
        return(x)
      },
      simplify = FALSE,
      USE.NAMES = TRUE
    )
    return(as.data.frame(x = df, row.names = GetRowNames(md = md)))
  }
  if (Exists(name = 'raw/X')) {
    md <- 'raw/var'
    counts <- as.matrix(x = adata[['raw/X']])
  } else if (Exists(name = 'X')) {
    md <- 'var'
    counts <- as.matrix(x = adata[['X']])
  } else {
    stop("Cannot find counts matrix")
  }
  if (!IsMetaData(md = md)) {
    stop("Cannot find feature-level metadata", call. = FALSE)
  } else if (!IsMetaData(md = 'obs')) {
    stop("Cannot find cell-level metadata", call. = FALSE)
  }
  cat("loading meta data ...\n")
  metadata <- LoadMetadata(md = 'obs')
  cat("loading row names ...\n")
  rownames(x = counts) <- GetRowNames(md = md)
  colnames(x = counts) <- rownames(x = metadata)
  object <- CreateSeuratObject(counts = counts)
  if (ncol(x = metadata)) {
    object <- AddMetaData(object = object, metadata = metadata)
  }
  return(object)
}

#' Load the reference RDS files
#'
#' Read in a reference \code{\link[Seurat]{Seurat}} object and annoy index. This
#' function can read either from URLs or a file path. In order to read properly,
#' there must be the following files:
#' \itemize{
#'  \item \dQuote{ref.Rds} for the downsampled reference \code{Seurat}
#'  object (for mapping)
#'  \item \dQuote{idx.annoy} for the nearest-neighbor index object
#' }
#'
#' @param path Path or URL to the two RDS files
#' @param seconds Timeout to check for URLs in seconds
#'
#' @return A list with two entries:
#' \describe{
#'  \item{\code{map}}{
#'   The downsampled reference \code{\link[Seurat]{Seurat}}
#'   object (for mapping)
#'  }
#'  \item{\code{plot}}{The reference \code{Seurat} object (for plotting)}
#' }
#'
#' @importFrom SeuratObject Idents<-
#' @importFrom Seurat LoadAnnoyIndex
#' @importFrom httr build_url parse_url status_code GET timeout
#' @importFrom utils download.file
#' @importFrom Matrix sparseMatrix
#'
#' @export
#' @examples
#' \dontrun{
#' # Load from a URL
#' ref <- LoadReference("https://seurat.nygenome.org/references/pbmc")
#' # Load from a directory
#' ref2 <- LoadReference("/var/www/html")
#' }
#'
LoadReference <- function(path, seconds = 10L) {
  ref.names <- list(
    map = 'ref.Rds',
    ann = 'idx.annoy'
  )
  if (substr(x = path, start = nchar(x = path), stop = nchar(x = path)) == '/') {
    path <- substr(x = path, start = 1, stop = nchar(x = path) - 1)
  }
  uri <- httr::build_url(url = httr::parse_url(url = path))
  if (grepl(pattern = '^://', x = uri)) {
    if (!dir.exists(paths = path)) {
      stop("Cannot find directory ", path, call. = FALSE)
    }
    mapref <- file.path(path, ref.names$map)
    annref <- file.path(path, ref.names$ann)
    exists <- file.exists(c(mapref, annref))
    if (!all(exists)) {
      stop(
        "Missing the following files from the directory provided: ",
        Oxford(unlist(x = ref.names)[!exists], join = 'and')
      )
    }
  } else {
    ref.uris <- paste(uri, ref.names, sep = '/')
    names(x = ref.uris) <- names(x = ref.names)
    online <- vapply(
      X = ref.uris,
      FUN = Online,
      FUN.VALUE = logical(length = 1L),
      USE.NAMES = FALSE
    )
    if (!all(online)) {
      stop(
        "Cannot find the following files at the site given: ",
        Oxford(unlist(x = ref.names)[!online], join = 'and')
      )
    }
    mapref <- url(description = ref.uris[['map']])
    annref <- tempfile()
    download.file(url = ref.uris[['ann']], destfile = annref, quiet = TRUE)
    on.exit(expr = {
      close(con = mapref)
      unlink(x = annref)
    })
  }
  # Load the map reference
  map <- readRDS(file = mapref)
  # Load the annoy index into the Neighbor object in the neighbors slot
  map[["refdr.annoy.neighbors"]] <- LoadAnnoyIndex(
    object = map[["refdr.annoy.neighbors"]],
    file = annref
  )
  # Validate that reference contains required dims
  if (ncol(x = map[["refDR"]]) < getOption(x = "Azimuth.map.ndims")) {
    stop("Provided reference doesn't contain at least ",
         getOption(x = "Azimuth.map.ndims"), " dimensions. Please either
         regenerate reference with requested dimensionality or adjust ",
         "the Azimuth.map.ndims option.")
  }
  # Create plotref
  ad <- Tool(object = map, slot = "AzimuthReference")
  # plotref.dr <- GetPlotRef(object = ad) # UseMethod(generic = 'GetPlotRef', object = ad)
  plotref.dr <- slot(object = ad, name = "plotref")
  cm <- sparseMatrix(
    i = 1, j = 1, x = 0, dims = c(1, nrow(x = plotref.dr)),
    dimnames = list("placeholder", Cells(x = plotref.dr))
  )
  plot <- CreateSeuratObject(
    counts = cm
  )
  plot[["refUMAP"]] <- plotref.dr
  plot <- AddMetaData(object = plot, metadata = Misc(object = plotref.dr, slot = "plot.metadata"))
  gc(verbose = FALSE)
  return(list(
    map = map,
    plot = plot
  ))
}

#' Transform an NN index
#'
#' @param object Seurat object
#' @param meta.data Metadata
#' @param neighbor.slot Name of Neighbor slot
#' @param key Column of metadata to use
#'
#' @return \code{object} with transfomed neighbor.slot
#'
#' @importFrom SeuratObject Indices
#'
#' @keywords internal
#'
NNTransform <- function(
  object,
  meta.data,
  neighbor.slot = "query_ref.nn",
  key = 'ori.index'
) {
  on.exit(expr = gc(verbose = FALSE))
  ind <- Indices(object[[neighbor.slot]])
  ori.index <- t(x = sapply(
    X = 1:nrow(x = ind),
    FUN = function(i) {
      return(meta.data[ind[i, ], key])
    }
  ))
  rownames(x = ori.index) <- rownames(x = ind)
  slot(object = object[[neighbor.slot]], name = "nn.idx") <- ori.index
  return(object)
}

#######################
######## main #########
#######################

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("input h5 file required as input.\n", call=FALSE)
}

inputfile.h5ad = args[1]
cat("inputfile.h5ad = ", inputfile.h5ad, "\n")
# input_dir = args[1]

options("Azimuth.map.ndims" = 50) # used in helpers.R::LoadReference()

# Ensure Seurat v4.0 or higher is installed
if (packageVersion(pkg = "Seurat") < package_version(x = "4.0.0")) {
  stop("Mapping datasets requires Seurat v4 or higher.", call. = FALSE)
}

# Ensure glmGamPoi is installed
if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    BiocManager::install("glmGamPoi")
  }
}

# copy input file and convert to h5seurat
Convert(inputfile.h5ad, dest="h5seurat", overwrite = TRUE)
inputfile.h5seurat <- paste0(file_path_sans_ext(inputfile.h5ad), ".h5seurat")
cat("inputfile.h5seurat = ", inputfile.h5seurat, "\n")

# Download the Azimuth reference and extract the archive

# Load the reference
# Change the file path based on where the reference is located on your system.
## reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")
reference <- LoadReference(REFERENCE_DIR)

# Load the query object for mapping
# Change the file path based on where the query file is located on your system.
# cat("Loading file", input_file, "\n")
#query <- LoadFileInput(path = input_file)
#query <- Seurat::Read10X(
#  input_dir
#  gene.column = 2 # 1: ensembl_ids, 2: gene_symbols
#  )
cat("Loading file", inputfile.h5seurat, "\n")
query <- LoadH5Seurat(inputfile.h5seurat)
cat("query file loaded.\n")
# Calculate nCount_RNA and nFeature_RNA if the query does not
# contain them already
if (!all(c("nCount_RNA", "nFeature_RNA") %in% c(colnames(x = query[[]])))) {
    calcn <- as.data.frame(x = Seurat:::CalcN(object = query))
    colnames(x = calcn) <- paste(
      colnames(x = calcn),
      "RNA",
      sep = '_'
    )
    query <- AddMetaData(
      object = query,
      metadata = calcn
    )
    rm(calcn)
}

# Calculate percent mitochondrial genes if the query contains genes
# matching the regular expression "^MT-"
if (any(grepl(pattern = '^MT-', x = rownames(x = query)))) {
  query <- PercentageFeatureSet(
    object = query,
    pattern = '^MT-',
    col.name = 'percent.mt',
    assay = "RNA"
  )
}

## cells must have been filtered by this stage
# Filter cells based on the thresholds for nCount_RNA and nFeature_RNA
# you set in the app
# cells.use <- query[["nCount_RNA", drop = TRUE]] <= 37340 &
#  query[["nCount_RNA", drop = TRUE]] >= 500 &
#  query[["nFeature_RNA", drop = TRUE]] <= 5321 &
#  query[["nFeature_RNA", drop = TRUE]] >= 50

# If the query contains mitochondrial genes, filter cells based on the
# thresholds for percent.mt you set in the app
#if ("percent.mt" %in% c(colnames(x = query[[]]))) {
#  cells.use <- cells.use & (query[["percent.mt", drop = TRUE]] <= 80 &
#    query[["percent.mt", drop = TRUE]] >= 0)
#}

# Remove filtered cells from the query
#query <- query[, cells.use]

# Preprocess with SCTransform
query <- SCTransform(
  object = query,
  assay = "RNA",
  new.assay.name = "refAssay",
  residual.features = rownames(x = reference$map),
  reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
  method = 'glmGamPoi',
  ncells = 2000,
  n_genes = 2000,
  do.correct.umi = FALSE,
  do.scale = FALSE,
  do.center = TRUE
)

# Find anchors between query and reference
anchors <- FindTransferAnchors(
  reference = reference$map,
  query = query,
  k.filter = NA,
  reference.neighbors = "refdr.annoy.neighbors",
  reference.assay = "refAssay",
  query.assay = "refAssay",
  reference.reduction = "refDR",
  normalization.method = "SCT",
  features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
  dims = 1:50,
  n.trees = 20,
  mapping.score.k = 100
)

# Transfer cell type labels and impute protein expression
#
# Transferred labels are in metadata columns named "predicted.*"
# The maximum prediction score is in a metadata column named "predicted.*.score"
# The prediction scores for each class are in an assay named "prediction.score.*"
# The imputed assay is named "impADT" if computed

refdata <- lapply(X = "celltype.l2", function(x) {
  reference$map[[x, drop = TRUE]]
})
names(x = refdata) <- "celltype.l2"
if (FALSE) {
  refdata[["impADT"]] <- GetAssayData(
    object = reference$map[['ADT']],
    slot = 'data'
  )
}
query <- TransferData(
  reference = reference$map,
  query = query,
  dims = 1:50,
  anchorset = anchors,
  refdata = refdata,
  n.trees = 20,
  store.weights = TRUE
)

# Calculate the embeddings of the query data on the reference SPCA
query <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = reference$map,
  query = query,
  reductions = "pcaproject",
  reuse.weights.matrix = TRUE
)

# Calculate the query neighbors in the reference
# with respect to the integrated embeddings
query[["query_ref.nn"]] <- FindNeighbors(
  object = Embeddings(reference$map[["refDR"]]),
  query = Embeddings(query[["integrated_dr"]]),
  return.neighbor = TRUE,
  l2.norm = TRUE
)

# The reference used in the app is downsampled compared to the reference on which
# the UMAP model was computed. This step, using the helper function NNTransform,
# corrects the Neighbors to account for the downsampling.
query <- NNTransform(
  object = query,
  meta.data = reference$map[[]]
)

# Project the query to the reference UMAP.
query[["proj.umap"]] <- RunUMAP(
  object = query[["query_ref.nn"]],
  reduction.model = reference$map[["refUMAP"]],
  reduction.key = 'UMAP_'
)


# Calculate mapping score and add to metadata
query <- AddMetaData(
  object = query,
  metadata = MappingScore(anchors = anchors),
  col.name = "mapping.score"
)

# save mapped data set
# save(query, file = "azimuth.bin")
# load("azimuth.bin" )


# VISUALIZATIONS
cat("Generating Visualizations ...\n")
pdf(NULL) # switch off automatic generation of Rplots.pdf

# First predicted metadata field, change to visualize other predicted metadata
id <- "celltype.l2"[1]
predicted.id <- paste0("predicted.", id)
predicted.id.score <- paste0(predicted.id, ".score")

# write a table of cell-type assignments, prediction and mapping scores:
fnam.table <- paste0(gsub(".", "_", predicted.id, fixed = TRUE),".tsv")
data <- FetchData(object = query, vars = c(predicted.id, predicted.id.score, "mapping.score"), slot = "data")
write.table(data, fnam.table, quote = FALSE, sep="\t")
#gzip(fnam.table, overwrite = TRUE)

# make a barplot for the number of cells in each celltype category
ctyp.counts <- function(vec, label) {
  tb <- table(as.factor(vec))
  v <- vector(mode = "character", length = length(tb))
  v[] <- label
  tabf <- cbind(v, as.data.frame(tb))
  names(tabf) <- c("threshold", "cell_type", "count")
  tabf
}
tdf <- ctyp.counts(data[,1], label = "all")
for (thresh in c(.5, .8)) {
  tag <- paste0("score>",format(thresh, digits=1))
  tdf <- rbind(
    tdf,
    ctyp.counts(data[data[,2] > thresh, 1],
      label = paste0("score>",format(thresh, digits=1))
    )
  )
}
p <- ggplot(tdf, aes(x=cell_type, y=count, fill = threshold))
p <- p + geom_col(position = "dodge")
p <- p + theme(axis.text.x = element_text(angle = 90))
ggsave("ncells_by_type_barplot.pdf")

# DimPlot of the reference
#ref.plt <- DimPlot(object = reference$plot, reduction = "refUMAP", group.by = id, label = TRUE) + NoLegend()
#ggsave(file = "ref_umap.pdf", plot = ref.plt)
#ggsave("ref_umap.pdf")

# DimPlot of the query, colored by predicted cell type
DimPlot(object = query, reduction = "proj.umap", group.by = predicted.id, label = TRUE) + NoLegend()
ggsave("query_umap.pdf")

# Plot the score for the predicted cell type of the query
FeaturePlot(object = query, features = paste0(predicted.id, ".score"), reduction = "proj.umap")
ggsave("prediction_score_umap.pdf")
VlnPlot(object = query, features = paste0(predicted.id, ".score"), group.by = predicted.id) + NoLegend()
ggsave("prediction_score_vln.pdf")

# Plot the mapping score
FeaturePlot(object = query, features = "mapping.score", reduction = "proj.umap")
ggsave("mapping_score_umap.pdf")
VlnPlot(object = query, features = "mapping.score", group.by = predicted.id) + NoLegend()
ggsave("mapping_score_vln.pdf")

# Plot the prediction score for the class CD16 Mono
#FeaturePlot(object = query, features = "CD16 Mono", reduction = "proj.umap")
#ggsave("prediction_score_CD16Mono_umap.pdf")
#VlnPlot(object = query, features = "CD16 Mono", group.by = predicted.id) + NoLegend()
#ggsave("prediction_score_CD16Mono_vln.pdf")

# Plot an RNA feature
# FeaturePlot(object = query, features = "GNLY", reduction = "proj.umap")
# VlnPlot(object = query, features = "GNLY", group.by = predicted.id) + NoLegend()
