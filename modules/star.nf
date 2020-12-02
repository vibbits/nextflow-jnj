#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

params.outdir = "$launchDir/results"

process star_idx {
    label 'high'
    container "quay.io/biocontainers/star:2.6.1d--0"

    input:
    path genome
    path gtf
    
    output:
    path "index_dir/", emit: index

    script:
    """
    mkdir index_dir
    
    STAR --runThreadN $task.cpus \\
      --runMode genomeGenerate \\
      --genomeDir index_dir/ \\
      --genomeFastaFiles $genome \\
      --genomeSAindexNbases $params.genomeSAindexNbases \\
      --sjdbGTFfile $gtf
    """
}

process star_alignment {
    publishDir "$params.outdir/mapped-reads/", mode: 'copy', overwrite: true, pattern: "*.bam"  
    label 'high'
    container "quay.io/biocontainers/star:2.6.1d--0"

    input:
    // (trim_fq, IDX.out, gtf)
    tuple val(sample), path(reads) 
    path indexDir
    path gtf

    output:
    tuple val(sample), path("*.bam") ,  emit: align_bam

    script:
    """
    mkdir -p $params.outdir/mapped-reads/


    STAR  \\
        --readFilesIn ${reads} \\
        --runThreadN $task.cpus \\
        --outSAMtype BAM SortedByCoordinate \\
        --sjdbGTFfile ${gtf} \\
        --outFileNamePrefix $sample. \\
        --genomeDir ${indexDir}
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
  index = star_idx.out.index
}

workflow MAP {
  take:
  reads
  indexDir
  gtf

  main: 
  star_alignment(reads, indexDir, gtf)

  emit:
  align_bam = star_alignment.out.align_bam
}
