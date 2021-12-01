# Workarounds for problems
- if the automatic build of Singularty containers fails:
  - pull the required images from Dockerhub individually with [Singularity](https://sylabs.io/singularity/) and place them in the directory pointed to by the environment variable NXF_SINGULARITY_CACHEDIR before executing [nextflow run](https://www.nextflow.io):
  ```bash
  export NXF_SINGULARITY_CACHEDIR=$(realpath ./singularity_images)
  mkdir -p $NXF_SINGULARITY_CACHEDIR
  singularity pull --dir ${NXF_SINGULARITY_CACHEDIR} docker://wtsihgi//nf_qc_cluster-2.4
  singularity pull --dir ${NXF_SINGULARITY_CACHEDIR} docker://wtsihgi/seurat_azimuth_pbmc:1.1

  ```
