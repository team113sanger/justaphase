process RUN_WHATSHAP {
    publishDir "${params.outdir}", mode: 'copy'
    container "quay.io/biocontainers/whatshap:2.3--py39h1f90b4d_0"
    
    input:
    tuple val(meta), path(vcf_file), path(bam_file), path(bai)
    path genome_fa
    path genome_fai

    output: 
    tuple val(meta), path("*.phased.vcf.gz"), emit: phased_vcf
    
    script:
    """
    whatshap phase \
    -o "${meta.contrast}.phased.vcf.gz" \
    --reference=${genome_fa} \
    --ignore-read-groups \
    --sample TUMOUR \
    ${vcf_file} \
    ${bam_file}
    """

}