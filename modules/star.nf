#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

params.outdir = "$launchDir/results"
params.genomeDir = "$launchDir/data/index/"

process star_idx {
    label 'high'
    container "quay.io/biocontainers/star:2.6.1d--0"

    input:
    path genome
    path gtf
    
    output:
    path "genome_dir"

    script:
    """
    mkdir genome_dir
    
    STAR --runThreadN $task.cpus \\
      --runMode genomeGenerate \\
      --genomeDir genome_dir \\
      --genomeFastaFiles $genome \\
      --genomeSAindexNbases $params.genomeSAindexNbases \\
      --sjdbGTFfile $gtf
    """
}

process star_alignment {
    publishDir "$params.outdir/mapped-reads", mode: 'copy', overwrite: true
    label 'high'
    container "quay.io/biocontainers/star:2.6.1d--0"

    input:
    tuple val(sample), path(reads) 
    path genomeDir
    path genome
    path gtf

    script:
    """
    mkdir -p $params.outdir/mapped-reads/

    STAR  \\
        --genomeDir $genomeDir \\  
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
  
  emit:
  star_idx.out
}

workflow MAP {
  take:
  reads
  genomeDir
  genome
  gtf

  main: 
  star_alignment(reads, genomeDir, genome, gtf)

  //emit:
  //bam_sorted_output = star_alignment.out.bam_sorted
}
