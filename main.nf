include { RUN_WHATSHAP } from './modules/whatshap.nf'
include { INDEX_PHASED_VARS; FIND_ADJACENT_VARIANTS; FIND_MNV_CANDIDATES; MERGE_SORT_AND_UNHEAD; COMPOSE_MNV_VARIANTS; EXTRACT_BAITSET_VARIANTS } from './modules/phase.nf'
include { ANNOTATE_VARIANTS } from './modules/annotate_variants.nf'



workflow  {
    genome = file(params.genome_files, checkIfExists: true)
    baitset = file(params.baitset, checkIfExists: true)
    vep_cache = file(params.vep_cache, checkIfExists: true)
    custom_files = Channel.of(params.custom_files.split(';'))
    .map(it -> file(it, checkIfExists: true))
    .collect()
    custom_args = Channel.of(params.custom_args.split(';'))
    .collect()
    .map { '--custom ' + it.join(' --custom ') }

    
    bamfile_ch  = Channel.fromPath(params.bam_files,checkIfExists: true) \
    | map { file -> 
    tuple([sample_id: file.baseName.replace(".sample.dupmarked", "")], file)}
    
    vcf_file_ch = Channel.fromPath(params.vcf_files)
    .splitCsv()
    .map(it -> file(it[0], checkIfExists: true))
    .map { file -> tuple([sample_id: file.baseName.split('_vs_')[0]], file)}
    
    input_ch = vcf_file_ch.join(bamfile_ch) \
    | map{ meta, vcf, bam -> 
    tuple(meta + [contrast: vcf.baseName.replace(".vcf", "")], vcf, bam)}
    | map { meta, vcf, bam  -> 
            index = bam + ".bai"
            tuple(meta, vcf, bam, index)}


    RUN_WHATSHAP(input_ch, genome)
    INDEX_PHASED_VARS(RUN_WHATSHAP.out.phased_vcf)
    FIND_ADJACENT_VARIANTS(INDEX_PHASED_VARS.out.indexed_vcf)
    FIND_ADJACENT_VARIANTS.out.vcf_bed_pair \
    | splitCsv(elem: 1, header: ['chr', 'start', 'stop'], sep: '\t')
    | set { intervals }
    FIND_MNV_CANDIDATES(intervals) 
    mnv_ch = FIND_MNV_CANDIDATES.out.queries.groupTuple()
    variant_sets = mnv_ch.join(INDEX_PHASED_VARS.out.indexed_vcf)
    variant_sets.view()
    COMPOSE_MNV_VARIANTS(variant_sets)
    EXTRACT_BAITSET_VARIANTS(COMPOSE_MNV_VARIANTS.out.mnv_file,baitset)
    ANNOTATE_VARIANTS(EXTRACT_BAITSET_VARIANTS.out.bait_variants, 
                    vep_cache, 
                    custom_files,
                    custom_args,
                    genome,
                    params.species,
                    params.assembly,
                    params.db_version)





}