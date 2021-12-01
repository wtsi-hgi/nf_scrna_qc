#!/usr/bin/env bash

test_dir="test_dir_$(date +"%Y_%m_%d_%I_%M_%p")"
mkdir $test_dir && cd $test_dir
echo testing in $PWD

git clone https://github.com/wtsi-hgi/nf_qc_cluster.git

# check nextflow is in PATH
nextflow -version

# check singularity is in PATH
singularity --version

# create and use local temporary directory
mkdir -p $PWD/tmp
export TMPDIR=$PWD/tmp
export SINGULARITY_TMPDIR=$PWD/tmp
export SINGULARITY_CACHEDIR=$PWD/tmp
export NXF_SINGULARITY_CACHEDIR=$PWD/tmp
export SINGULARITY_LOCALCACHEDIR=$PWD/tmp

# run test 
#   (the 'test' profile runs all NF tasks locally, and using singularity containers pulled from dockerhub).
nextflow run ./nf_qc_cluster/main.nf \
	--file_sample_qc \
	./nf_qc_cluster/test/params.yml \
	-params-file ./nf_qc_cluster/test/params.yml \
	-c ./nf_qc_cluster/nextflow.config \
	-c ./nf_qc_cluster/test/inputs.nf \
	-profile test
