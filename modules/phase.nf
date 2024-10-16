
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

process FIND_MNV_VARIANTS {
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
    bcftools view -H $region > "${meta.sample_id}_mnv_lines.tsv"
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

process INTERSECT {
    module "bedtools2-2.31.1/python-3.10.10"
    input: 
    tuple val(meta), path(vcf_file) 
    path(baitset)
    
    output:
    tuple val(meta), path(".baitset_only.vcf.gz")
    
    script: 
    """
    bedtools intersect -header -a $vcf_file -b $baitset > "${meta.sample_id}.baitset_only.vcf.gz"
    """
}

process COMPOSE_MNV_VARIANTS {
    container 'docker://gitlab-registry.internal.sanger.ac.uk/dermatlas/fur_phaser_py/feature/build_env'
    input:
    tuple val(meta), path(subset), path(index)
    tuple val(meta), path(vcf_file), path(vcf_index)

    output: 
    path("*.vcf")

    script:
    """
    python3 /opt/repo/src/fur_phaser_py/phaser.py -c $subset -v $vcf_file
    """
}