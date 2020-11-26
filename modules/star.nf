#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

params.outdir = "$launchDir/results"
params.genomeDir = "$launchDir/data/index/"

process star_idx {
    publishDir "$params.genomeDir", mode: 'copy', overwrite: true
    label 'high'
    container "quay.io/biocontainers/star:2.6.1d--0"

    input:
    path genome
    path gtf

    script:
    """
    mkdir -p $params.genomeDir

    STAR --runThreadN $task.cpus \\
      --runMode genomeGenerate \\
      --genomeDir $params.genomeDir \\
      --genomeFastaFiles $genome \\
      --genomeSAindexNbases $params.genomeSAindexNbases \\
      --sjdbGTFfile $gtf
    """
}

process star_alignment {
    publishDir "$params.outdir/mapped-reads", mode: 'copy', overwrite: true
    label 'low'
    container "quay.io/biocontainers/star:2.6.1d--0"

    input:
    tuple val(sample), path(reads) 
    path (index)
    path (gtf)

    //output:
    //tuple val(sample), path("*Aligned.out.bam") , emit: bam
    //tuple val(meta), path("*Log.final.out")   , emit: log_final
    //tuple val(meta), path("*Log.out")         , emit: log_out
    //tuple val(meta), path("*Log.progress.out"), emit: log_progress
    //path  "*.version.txt"                     , emit: version

    //tuple val(sample), path("*sortedByCoord.out.bam")  , emit: bam_sorted
    //tuple val(meta), path("*toTranscriptome.out.bam"), optional:true, emit: bam_transcript
    //tuple val(meta), path("*fastq.gz")               , optional:true, emit: fastq
    //tuple val(meta), path("*.tab")                   , optional:true, emit: tab

    script:
    """
    mkdir -p $params.outdir/mapped-reads/

    STAR  \\
        --genomeDir $index/ \\  
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --runThreadN $task.cpus \\
        --outSAMtype BAM SortedByCoordinate \\
        --sjdbGTFfile ${gtf} 
    """
}


// Running a workflow with the defined processes here.  
workflow IDX {
  take:
  genome
  gtf

  main: 
  star_idx(genome, gtf)
}

workflow MAP {
  take:
  reads
  index
  gtf

  main: 
  star_alignment(reads,index, gtf)

  //emit:
  //bam_sorted_output = star_alignment.out.bam_sorted
}
