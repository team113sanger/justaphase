include { RUN_WHATSHAP } from './modules/whatshap.nf'
include { INDEX_PHASED_VARS; FIND_ADJACENT_VARIANTS; FIND_MNV_VARIANTS; MERGE_SORT_AND_UNHEAD; COMPOSE_MNV_VARIANTS } from './modules/phase.nf'



workflow  {
    genome = file(params.genome_files, checkIfExists: true)
    
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
    FIND_MNV_VARIANTS(intervals) 
    mnv_ch = FIND_MNV_VARIANTS.out.queries.groupTuple()
    COMPOSE_MNV_VARIANTS(mnv_ch,INDEX_PHASED_VARS.out.indexed_vcf)
    // groups_ch.collectFile(name: "test.txt",
    //                       storeDir:"/lustre/scratch125/casm/team113da/users/bf14/variant_caller_benchmarking/whatshap/fur_whatshap", 
    //                       newLine: false)
    
}