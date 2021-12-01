# running pipeline tests

## test dataset
The test dataset contains UMI counts of pooled PBMC samples from 10 donors. The data where down-sampled to 200 cells per donor and 324 cells for which the donor could not be assigned.

## requirements for processing the test dataset
For the processing of the test dataset a linux machine is required with
- 4 GB of memory
- 4 cores
- 5 GB of disk storage

and the following software installed:
- [Singularity](https://sylabs.io/singularity/) container engine
- [Nextflow](https://www.nextflow.io)

## how to run the pipeline on the test dataset
To process the test dataset, execute the following commands:
```bash
  git clone https://github.com/wtsi-hgi/nf_qc_cluster.git
  nextflow run ./nf_qc_cluster/main.nf \
    --file_sample_qc ./nf_qc_cluster/test/params.yml \
    -params-file ./nf_qc_cluster/test/params.yml \
    -c ./nf_qc_cluster/nextflow.config \
    -c ./nf_qc_cluster/test/inputs.nf \
    -profile test
```

This will create the following directories in the calling directory
- ./test_output
  - contains the [output files and directories](README-outputfiles.md) of the pipeline
- ./test_work
  - contains the nexflow work directories for each process
- ./test_reports
  - tracing ouput
- ./singularity_images
  - cache for [Singularity](https://sylabs.io/singularity/) image files.

If the automatic build of Singularity containers fails, see [problems & workarounds](README-workarounds.md).
