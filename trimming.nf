#!/usr/bin/env nextflow

params.reads = "$baseDir/data/*0.01_{1,2}.fq.gz"
params.outdir = "$baseDir/results"
params.threads = 2

println """\
         NEXTFLOW COURSE
================================
Reads            : $params.reads
Output-folder    : $params.outdir
"""

read_pairs_ch = Channel.fromFilePairs(params.reads)
                      

process trimmomatic {
    publishDir "$params.outdir/trimming", mode: 'copy'
    
    input:
    tuple val(sample), file(reads) from read_pairs_ch

    output:
    tuple val(sample), file("${sample}_1P.fq"), file("${sample}_2P.fq") into (paired_fq)
    tuple val(sample), file("${sample}_1U.fq"), file("${sample}_2U.fq") into (unpaired_fq)
    file("${sample}.trimmomatic.stats.log") into (trimmomatic_stats)

    script:
    """
    mkdir -p $params.outdir/trimmed-reads/
    trimmomatic PE -threads $params.threads ${reads[0]} ${reads[1]} -baseout ${sample} SLIDINGWINDOW:4:15    AVGQUAL:30 > ${sample}.trimmomatic.stats.log
    
    mv ${sample}_1P ${sample}_1P.fq
    mv ${sample}_2P ${sample}_2P.fq
    mv ${sample}_1U ${sample}_1U.fq
    mv ${sample}_2U ${sample}_2U.fq
    """
}

