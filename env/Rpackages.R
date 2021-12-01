pkgs <- c("igraph", "optparse","data.table", "ggplot2", "devtools", "BiocManager")
github.pkgs <- c("immunogenomics/harmony", "YuLab-SMU/enrichplot")
bioc.pkgs <- c("clusterProfiler", "org.Hs.eg.db")

cran.mirrors <- c("https://cran.ma.imperial.ac.uk/", "https://www.stats.bris.ac.uk/R/")

#user.lib <- c("/home/container_user/R/lib")
#.libPaths(c(.libPaths(), user.lib))
user.lib <- NULL

install.packages(pkgs, lib = user.lib, repos=cran.mirrors)
devtools::install_github(c("immunogenomics/harmony", "YuLab-SMU/enrichplot"), lib = user.lib)
BiocManager::install(bioc.pkgs, version = "3.12", update = TRUE, ask = FALSE, lib = user.lib)
