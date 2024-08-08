#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
process FOO {
  input:
  tuple val(metadata), path(bamfile)
  
  output:
  
  tuple val(metadata), path("*_output_bamfile.bam"), emit: example

  script:
  """
  echo $bamfile > ${metadata.sample_id}_output_bamfile.bam
  """
}

workflow {

    // Cohort files 
    bamfiles = Channel.fromPath(params.bamfiles, checkIfExists: true) \
    | map { file -> tuple([sample_id: file.baseName], file)}
    bamfiles.view()
    FOO(bamfiles)
    FOO.out.example.view()

}