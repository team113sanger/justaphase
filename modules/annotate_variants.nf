

process ANNOTATE_VARIANTS {
    publishDir "${params.outdir}", mode: "copy"
    container "ensemblorg/ensembl-vep:release_104.1"

    input:
    tuple val(meta), path(vcf_file)
    path(vep_cache)
    path(custom_files)
    val(custom_args)
    tuple path(reference)//, path(index)
    val(species)
    val(assembly_name)
    val(db_version)
    
    output: 
    tuple val(meta), path("*vep.vcf.gz"),path("*vep.vcf.gz.tbi"), emit: vep_annotation
    script: 
    def outfname = "${vcf_file}".replace(".vcf.gz", "") + ".vep.vcf.gz"
    // def refbase = reference[0].baseName
    if (species == "homo_sapiens"){
    """
    vep -i ${vcf_file} \
    --cache \
    -t SO \
    --dir_cache ${vep_cache} \
    --db_version ${db_version} \
    --output_file ${outfname} \
    --buffer_size 20000 \
    --format vcf \
    --species ${species} \
    --offline \
    --symbol \
    --biotype \
    --vcf \
    --sift s \
    --no_stats \
    --assembly ${assembly_name} \
    --flag_pick_allele \
    --canonical \
    --hgvs \
    --shift_hgvs 1 \
    --fasta ${reference} \
    $custom_args \
    --compress_output bgzip \
    --mane \
    --protein \
    --numbers \
    --polyphen p \
    --domains \
    --transcript_version \
    --show_ref_allele \
    --fork 4 
    tabix -p vcf ${outfname}
    """
    }
    else {

    """
    vep -i ${vcf_file} \
    --output_file ${outfname} \
    --cache \
    --dir ${vep_cache} \
    --fasta ${reference} \
    --db_version ${db_version} \
    --species ${species} \
    --assembly ${assembly_name} \
    --offline \
    $custom_args \
    -t SO \
    --format vcf \
    --buffer_size 20000 \
    --offline \
    --symbol \
    --biotype \
    --vcf \
    --sift s \
    --no_stats \
    --flag_pick_allele_gene \
    --canonical \
    --hgvs \
    --shift_hgvs 1 \
    --compress_output bgzip \
    --mane \
    --protein \
    --numbers \
    --fork 4 \
    --domains
    tabix -p vcf ${outfname}
    """

    }

    stub: 
    """
    echo stub > item.vep.vcf.gz
    echo stub > item.vep.vcf.gz.tbi
    """
    }
