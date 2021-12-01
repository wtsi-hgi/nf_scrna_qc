
# Notes

## built dockerfile and push to dockerhub

```
nohup sudo docker build . -t wtsihgi/nf_qc_cluster:2.1 &
tail -f nohup.out
docker login
docker push wtsihgi/nf_qc_cluster:2.1
```

## General

We use conda for version control for command line programs and python pacakges.

## R

CURRENTLY THIS PIPELINE REQUIRES NO R PACKAGES. THE BELOW CODE SHOULD NOT BE RUN.

Conda, however, does not work well with R packages. One could use conda to install a base version of R. Alternatively, one can pull from from a rocker image (e.g., `rocker/r-ver:3.6.2` instead of `buildpack-deps:curl` or a minimal install like `alpine` ).

For R package management we use `renv`.

```R
# initialize renv
renv::init()

# update renv after intalling any packages
renv::settings$snapshot.type("simple")
renv::snapshot()

# add renv.lock, .Rprofile, and renv/activate.R to git
```

Note: for the enrichment script (currently not run), the following packages must be installed
```R
- bioconductor-dose=3.14.0
- bioconductor-clusterprofiler=3.16.0
- bioconductor-org.hs.eg.db=3.11.1
```
