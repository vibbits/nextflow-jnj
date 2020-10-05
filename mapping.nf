#!/usr/bin/env nextflow

params.reads = "$baseDir/results/trimmed-reads/*_{1,2}P.fq"
params.outdir = "$baseDir/results"
params.genome = "$baseDir/reference/Homo_sapiens_assembly38_chr22.fa"
params.threads = 2

println """\
        NEXTFLOW COURSE
================================
Reads            : $params.reads
Genome           : $params.genome
Output-folder    : $params.outdir
"""

read_pairs_ch = Channel.fromFilePairs(params.reads)
genome_file = file(params.genome)
                      

process mapping {
    publishDir "$params.outdir/mapped-reads", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample), file(reads) from read_pairs_ch

    output:
    file "${sample}.mappings.bam" into (bam_ch)

    script:
    // CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    // readGroup = "@RG\\tID:${idRun}\\t${CN}PU:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
    
    """
    mkdir -p $params.outdir/mapped-reads/
    
    bwa mem -M ${genome_file} ${reads[0]} ${reads[1]} | samtools view -b - -o ${sample}.mappings.bam  
    """
}