#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process FOO {
  container "quay.io/biocontainers/samtools:1.9--h91753b0_8"
  input:
  tuple val(metadata), path(bamfile)
  
  output:
  tuple val(metadata), path("*_output_bamfile.bam"), emit: example

  script:
  """
  samtools --version
  echo $bamfile > ${metadata.sample_id}_output_bamfile.bam
  """
}

workflow {

    bamfiles = Channel.fromPath(params.bamfiles, checkIfExists: true) \
    | map { file -> tuple([sample_id: file.baseName], file)}
    bamfiles.view()
    
    
    // tn_pairs = Channel.fromPath(params.manifest, checkIfExists: true) \
    // | splitCsv(sep:"\t", header:['tumor', 'normal','t_bam','n_bam']) 
    // | map{ row -> tuple([pair_id: row.normal+ "_" + row.tumor],row.t_bam,row.n_bam) }
    // tn_pairs.view()
    
    FOO(bamfiles)
    FOO.out.example.view()


}