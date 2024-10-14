
process RUN_VEP {
    module "ensembl-vep/104"
    input: 
    tuple val(meta), path(INPUT_VCF)
    val(species)
    val(assembly)
    path(cache_dir)
    path(custom_libs)


    output: 
    path("*.vcf")
    
    script:
    """
    vep \
    --cache \
    --dir $cache_dir \
    -i "$INPUT_VCF" \
    --db_version 104 \
    -t SO \
    --format vcf \
    --force_overwrite \
    -o "$OUTPUT_VCF" \
    --buffer_size 20000 \
    --species $species \
    --offline \
    --symbol \
    --biotype \
    --vcf \
    --sift s \
    --no_stats \
    --assembly $build \
    --flag_pick_allele_gene \
    --canonical \
    --hgvs \
    --shift_hgvs 1 \
    --fasta $genome \
    --custom "${custom_libs},${args}" \
    --compress_output bgzip \
    --mane \
    --numbers \
    --fork 4 \
    --domains
    """
}
