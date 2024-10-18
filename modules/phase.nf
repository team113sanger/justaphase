
process INDEX_PHASED_VARS {
    publishDir "${params.outdir}", mode: 'copy'
    label "phase"
    module "bcftools-1.19/python-3.11.6"

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
    module "/software/CASM/modules/modulefiles/casm-smart-phase/0.1.8"
    
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
    module "bcftools-1.19/python-3.11.6"
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
    module "bcftools-1.19/python-3.11.6"
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

process EXTACT_BAITSET_VARIANTS {
    module "bedtools2-2.31.1/python-3.10.10"
    input: 
    tuple val(meta), path(vcf_file) 
    path(baitset)
    
    output:
    tuple val(meta), path(".baitset_only.vcf.gz"), emit: bait_variants
    
    script: 
    """
    bedtools intersect -header -a $vcf_file -b $baitset > "${meta.sample_id}.baitset_only.vcf.gz"
    """
}

process COMPOSE_MNV_VARIANTS {
    container 'docker://gitlab-registry.internal.sanger.ac.uk/dermatlas/fur_phaser_py/feature/build_env:1d53ec06'

    input:
    tuple val(meta), path(subset)
    tuple val(meta), path(vcf_file), path(vcf_index)

    output: 
    path("*.vcf.gz"), emit: mnv_file

    script:
    """
    python3 /opt/repo/src/fur_phaser_py/phaser.py -c $subset -p $vcf_file \
    -o "${meta.contrast}.vcf.gz"
    """
}