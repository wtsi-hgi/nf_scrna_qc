# Workarounds for problems
- if the automatic build of Singularty containers fails:
  - pull the required images from Dockerhub individually with [Singularity](https://sylabs.io/singularity/) and place them in the directory pointed to by the environment variable NXF_SINGULARITY_CACHEDIR before executing [nextflow run](https://www.nextflow.io):
  ```bash
  export NXF_SINGULARITY_CACHEDIR=$(realpath ./singularity_images)
  mkdir -p $NXF_SINGULARITY_CACHEDIR
  singularity pull --dir ${NXF_SINGULARITY_CACHEDIR} docker://wtsihgi/nf_qc_cluster:2.4
  singularity pull --dir ${NXF_SINGULARITY_CACHEDIR} docker://wtsihgi/seurat_azimuth_pbmc:1.1
  singularity pull --dir ${NXF_SINGULARITY_CACHEDIR} docker://wtsihgi/nf_qc_cluster:sccaf_1.5
  ```

will create the .sif files in $NXF_SINGULARITY_CACHEDIR dir:
```
nf_qc_cluster_2.4.sif
seurat_azimuth_pbmc_1.1.sif
nf_qc_cluster_sccaf_1.5.sif
```

Copy these images in Nextflow singularity cache dir, defined in configuration fines as singularity.cacheDir, and renames files to change extension from .sif to .img 
