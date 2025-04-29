
process INDEX_PHASED_VARS {
    publishDir "${params.outdir}", mode: 'copy'
    container "quay.io/biocontainers/bcftools:1.19--h8b25389_1"
    label "phase"

    input: 
    tuple val(meta), path(vcf_file)

    output: 
    tuple val(meta), path(vcf_file), path("*.csi"), emit: indexed_vcf

    script:
    """
    bcftools index $vcf_file
    """

}



process FIND_ADJACENT_VARIANTS {
    publishDir "${params.outdir}", mode: 'copy'
    label "casm_smartphase"
    container "quay.io/wtsicgp/casm-smart-phase:0.1.8"
    
    input:
    tuple val(meta), path(vcf), path(index)

    output:
    tuple val(meta), path("*.adjacent_snvs.bed"), path(vcf), path(index), emit: vcf_bed_pair
    
    script:
    """
    casmsmartphase generate-bed \
    -f $vcf \
    -o "${meta.contrast}.adjacent_snvs.bed" \
    --markhz
    """

}

process FIND_MNV_CANDIDATES {
    container "quay.io/biocontainers/bcftools:1.19--h8b25389_1"
    label "phase"
    
    input:
    tuple val(meta), val(bed_values), path(vcf_file), path(index)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.csi"), emit: lines
    tuple val(meta), path("*_mnv_lines.tsv"), emit: queries

    script:
    """
    bcftools view \
    -r "${bed_values.chr}:${bed_values.start}-${bed_values.stop}" $vcf_file \
    -O z -o "${bed_values.chr}:${bed_values.start}-${bed_values.stop}.vcf.gz"
    bcftools index "${bed_values.chr}:${bed_values.start}-${bed_values.stop}.vcf.gz"
    bcftools view -H "${bed_values.chr}:${bed_values.start}-${bed_values.stop}.vcf.gz" > "${meta.sample_id}_${bed_values.chr}:${bed_values.start}-${bed_values.stop}_mnv_lines.tsv"
    """

}


process MERGE_SORT_AND_UNHEAD {
    label "phase"
    container "quay.io/biocontainers/bcftools:1.19--h8b25389_1"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(meta), path(region), path(indices)

    output:
    path("*mnv_lines.tsv"), emit: lines

    script:
    """
    bcftools view -H $region > "${meta.sample_id}_mnv_lines.tsv"
    """

}

process EXTRACT_BAITSET_VARIANTS {
    publishDir "${params.outdir}", mode: 'copy'
    container "quay.io/biocontainers/bedtools:2.31.1--h13024bc_3"
    input: 
    tuple val(meta), path(vcf_file) 
    path(baitset)
    
    output:
    tuple val(meta), path("*.baitset_only.vcf.gz"), emit: bait_variants
    
    script: 
    def outfname = "${vcf_file}".replace(".vcf.gz", "") + ".baitset_only.vcf.gz"
    """
    bedtools intersect -header -a $vcf_file -b $baitset > $outfname
    """
}

process COMPOSE_MNV_VARIANTS {
    publishDir "${params.outdir}", mode: 'copy'
    container "quay.io/team113sanger/fur_phaser_py:0.9.2"

    input:
    tuple val(meta), path(subset), path(vcf_file), path(vcf_index)

    output: 
    tuple val(meta), path("*.vcf.gz"), emit: mnv_file

    script:
    """
    python3 -m fur_phaser_py \
    -c $subset -p $vcf_file \
    -o "${meta.contrast}.vcf.gz"
    """
}