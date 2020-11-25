#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

params.outdir = "$launchDir/results"

// Process trimmomatic
process trimmomatic {
    publishDir "$params.outdir/trimmed-reads", mode: 'copy', overwrite: true
    label 'low'

    // Same input as fastqc on raw reads, comes from the same channel. 
    input:
    tuple val(sample), path(reads) 

    output:
    tuple val("${sample}"), path("${sample}{1,2}_P.fq")

    script:
    """
    mkdir -p $params.outdir/trimmed-reads/
    trimmomatic PE -threads $params.threads ${reads[0]} ${reads[1]} ${sample}1_P.fq ${sample}1_U.fq ${sample}2_P.fq ${sample}2_U.fq $params.slidingwindow $params.avgqual 
    """
}

workflow TRIM {
    take: 
    input
    
    main:
		trim_fq = trimmomatic(input)

    emit:
    trim_fq
}