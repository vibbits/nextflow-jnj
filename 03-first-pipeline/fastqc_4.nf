#!/usr/bin/env nextflow

params.reads = "$baseDir/data/*0.01_{1,2}.fq.gz"
params.env = "$params.abspath/environment.yml"

println "reads: $params.reads"

/**
 * Quality control fastq
 */


read_pairs_ch = Channel .fromFilePairs(params.reads)


process fastqc_raw_reads {
    conda "$params.env"

    input:
    set sample_id, file(reads) from read_pairs1_ch
    tuple val(x), file('latin.txt') from values
    script:
    """
    mkdir -p $params.outdir/quality-control
    fastqc ${reads}
    """
}