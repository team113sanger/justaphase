process RUN_WHATSHAP {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    tuple val(meta), path(vcf_file), path(bam_file), path(bai)
    tuple path(reference), path(index)

    output: 
    tuple val(meta), path("*.phased.vcf.gz"), emit: phased_vcf
    
    script:
    def refbase = reference[0].baseName
    """
    whatshap phase \
    -o "${meta.contrast}.phased.vcf.gz" \
    --reference=$refbase \
    --ignore-read-groups \
    --sample TUMOUR \
    ${vcf_file} \
    ${bam_file}
    """

}