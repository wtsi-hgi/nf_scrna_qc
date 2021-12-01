#!/usr/bin/env Rscript

library(Seurat)
library(hdf5r)
library(optparse)

unfactor <- function(f){
  out <- as.numeric(levels(f))[f]
}


main <- function() {
    optionList <- list(
        optparse::make_option(c("-i", "--input_file"),
            type = "character",
            help = paste0(
                "Input h5 file"
            )
        ),

        optparse::make_option(c("--maximum_de"),
            type = "numeric",
            default = 5,
            help = paste0(
                "Maximum differentially expressed genes.",
                " [default: %default]"
            )
        ),

        optparse::make_option(c("--auc_difference"),
            type = "numeric",
            default = 0.1,
            help = paste0(
                "Difference for AUC, truncated de genes.",
                " [default: %default]"
            )
        ),

        optparse::make_option(c("--n_cpu"),
            type = "numeric",
            default = 1,
            help = paste0(
                "Number of CPUs to use.",
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

    input.file=arguments$options$input_file
    output.file.basename=arguments$options$out_file
    max_de=arguments$options$maximum_de
    auc_diff=arguments$options$auc_difference
    #with: AUC < 0 + auc_diff and AUC > 1 - auc_diff

    # Increase the number of CPUs for Seurat
    # need to increase the RAM
    # see https://satijalab.org/seurat/v3.0/future_vignette.html
    future::plan("multiprocess", workers = arguments$options$n_cpu)
    options(future.globals.maxSize = 100000 * 1024^2)

    # Data
    print("reading h5")
    X <- H5File$new(input.file, mode="r")
    M=X[['X']][1:X[['X']]$dims[1],1:X[['X']]$dims[2]]
    cells=X[['cells']][1:X[['cells']]$dims]
    genes=X[['genes']][1:X[['genes']]$dims]
    df=as.data.frame(M)
    colnames(df) <- cells
    rownames(df) <- genes
    print("done h5")

    # Metadata
    obs=data.frame("cluster"=X[['cluster']][1, 1:X[['cluster']]$dims[2]])
    rownames(obs) <- cells
    print("metadata")

    # Seurat object
    seur <- CreateSeuratObject(counts=df, assay = "RNA", meta.data = obs)
    Idents(object=seur) <- factor(seur$cluster)
    print("made suerat")

    # Initialize
    clusters_active=unfactor(seur@active.ident)
    clusters_unique=sort(unique(unfactor(seur@active.ident)))
    print(paste0("clusters_unique = ", clusters_unique))
    update_min=0
    k=1

    merging_progress=list()
    merging_matrix=list()
    plt=list()

    merging_progress[[k]] <- clusters_active


    while(update_min<max_de) { 
      print(paste0("k=",k))
      mmm=matrix(data=NA, nrow=length(clusters_unique),
                 ncol=length(clusters_unique))

      for(j in 1:length(clusters_unique)){
        for(i in 1:length(clusters_unique)){
          print(paste0("i=", i, "j=", j))

          if(clusters_unique[j]==clusters_unique[i]){

            mmm[j,i]=NA

          } else {

            Y = FindMarkers(seur,
                            ident.1 = clusters_unique[j],
                            ident.2 = clusters_unique[i],
                            min.pct = 0.25,
                            random.seed = 1,
                            test.use = "roc")

            # matrix with DE genes between each pair of clusters
            mmm[j,i]=sum(Y$myAUC<auc_diff, Y$myAUC>(1-auc_diff))
          }

        }

      }

      df=data.frame(mmm)
      rownames(df) <- clusters_unique
      colnames(df) <- clusters_unique

      merging_matrix[[k]]=df
      plt[[k]] <- ggplot(
          data.frame("de_genes"= as.vector(t(merging_matrix[[k]]))),
            aes(x=de_genes)) + geom_histogram(binwidth=0.9,
            colour="black", fill="white")
      plt[[k]] <- plt[[k]] + ggtitle(paste0("Merging step: ", k))

      update_min=min(mmm, na.rm = TRUE) #find minimum de genes
      cond=which(mmm == update_min, arr.ind = TRUE) #find clusters with min de gene

      if(length(cond) == 2){
        clusters_active=replace(clusters_active,
                                clusters_active==clusters_unique[cond[1]],
                                clusters_unique[cond[2]]) # merge clusters

      } else if (length(cond) > 2){
        cond_sampled=sample(1:nrow(cond), 1)
        clusters_active=replace(clusters_active,
            clusters_active==clusters_unique[cond[cond_sampled,][1]],
            clusters_unique[cond[cond_sampled,][2]]
        )
      }

      Idents(object=seur)=factor(clusters_active)
      clusters_unique.update=sort(unique(unfactor(seur@active.ident)))

      k=k+1
      clusters_unique=clusters_unique.update
      merging_progress[[k]]=clusters_active
    }


    # Save results
    print("saving results")
    mp=as.data.frame(do.call(cbind,merging_progress))
    colnames(mp) <- paste0("merge_step_",0:(k-1))
    mp = cbind("cell_barcode" = colnames(seur), mp)
    gzfh <- gzfile(paste0(base, ".tsv.gz"), "w", compression = 9)
    write.table(
        mp,
        gzfh,
        row.names = FALSE,
        col.names = TRUE,
        quote = TRUE,
        sep = "\t",
        na = ""
    )
    close(gzfh)

    # Save plots
    pdf("hist_de_genes.pdf", onefile = TRUE)
    for (i in seq(length(plt))) {
      print(plt[[i]])
    }
    dev.off()

}


# like python if __name__ == '__main__'
if (sys.nframe() == 0) {
    #dev()
    main()
}
