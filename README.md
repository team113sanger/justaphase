# justaphase
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.04.5-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

justaphase is a bioinformatics pipeline written in [Nextflow](http://www.nextflow.io) that can be used for species-agnostic assembly and annotaton of MNVS from Caveman SNV calls. 

## Pipeline summary

In brief, the pipeline takes a set of bam files and corresponding Caveman VCFs and finds those variants which can be phased to produce MNVs. It utilises CASM-Smartphase and Whatshap for phasing and reconstructs haplotypes with fur_phaser_py, an internally developed tool for merging SNVs with the same phase group.

## Inputs 
- `vcf_files`: Path to a set of Caveman VCF files to be phased
- `bam_files`: Path to a set of corresponding bam files to use in phasing
- `genome`: The reference genome file `.fa` file and corresponding index files (`fa.fai`, `fa.dict`)
- `baitset`: A bed file descibing the regoions that were used for targeted capture 
- `vep_cache`: Path to a local ensembl vep cache to use in variant annotation
-`custom_files`: Path to custom files used in generating vep annotations


## Usage 

The recommended way to launch this pipeline on LSF is using a wrapper script (e.g. `bsub < my_wrapper.sh`) that submits nextflow as a job and records the version (**e.g.** `-r 0.1.1`)  and the `.json` parameter file supplied for a run.

An example wrapper script:
```
#!/bin/bash
#BSUB -q normal
#BSUB -R "select[mem>8000] rusage[mem=8000] span[hosts=1]"
#BSUB -M 8000
#BSUB -oo nf_out.o
#BSUB -eo nf_out.e

PARAMS_FILE="<YOUR_PARAMS_FILE>"

# Load module dependencies
module load nextflow-23.10.0
module load singularity/3.11.4

# Create a nextflow job that will spawn other jobs

nextflow run 'https://github.com/team113sanger/justaphase' \
-r 0.2.0 \
-params-file $PARAMS_FILE \
-c nextflow.config \
-profile farm22 
```

The pipeline can configured to run on either Sanger OpenStack secure-lustre instances or farm22 by changing the profile speicified:
`-profile secure_lustre` or `-profile farm22`. 



