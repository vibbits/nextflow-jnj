#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

params.outdir = "$launchDir/results"
params.genomeDir = "$launchDir/data/index/"

process star_index {
    publishDir "$params.genomeDir", mode: 'copy', overwrite: true
    label 'low'

    input:
    path genome_ch
    path gtf_ch

    output:
    path "star", emit: index

    script:
    """
    mkdir -p $params.genomeDir

    STAR --runThreadN $task.cpus \\
      --runMode genomeGenerate \\
      --genomeDir $params.genomeDir \\
      --genomeFastaFiles $genome_ch \\
      --genomeSAindexNbases $params.genomeSAindexNbases \\
      --sjdbGTFfile $gtf_ch 
    """
}


process star_alignment {
    publishDir "$params.outdir/mapped-reads", mode: 'copy', overwrite: true
    label 'high'

    input:
    tuple val(sample), path(reads) 
    path (index)
    path (gtf)

    output:
    //tuple val(sample), path("*Aligned.out.bam") , emit: bam
    //tuple val(meta), path("*Log.final.out")   , emit: log_final
    //tuple val(meta), path("*Log.out")         , emit: log_out
    //tuple val(meta), path("*Log.progress.out"), emit: log_progress
    //path  "*.version.txt"                     , emit: version

    tuple val(sample), path("*sortedByCoord.out.bam")  , emit: bam_sorted
    //tuple val(meta), path("*toTranscriptome.out.bam"), optional:true, emit: bam_transcript
    //tuple val(meta), path("*fastq.gz")               , optional:true, emit: fastq
    //tuple val(meta), path("*.tab")                   , optional:true, emit: tab

    script:
    """
    mkdir -p $params.outdir/mapped-reads/
    STAR --runMode alignReads \\
        --genomeDir $index \\  
        --runThreadN $task.cpus \\
        --outFileNamePrefix ${sample} \\
        --outSAMtype BAM SortedByCoordinate \\
        --sjdbGTFfile ${gtf} \\
        --readFilesIn ${reads[0]} ${reads[1]} 
    """
}


// Running a workflow with the defined processes here.  
workflow IDX {
  take:
  genome
  gtf

  main: 
  index = star_index(genome, gtf)

  emit:
  index
}

workflow MAP {
  take:
  reads
  index
  gtf

  main: 
  aligned_reads = star_alignment(reads,index, gtf)

  emit:
  aligned_reads //= star_alignment.out.samples_bam
}
