#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

process star_index {
    publishDir "$dirgenome/idx/", mode: 'copy', overwrite: true
    
    input:
    path dirgenome
    file genome
    file gtf

    script:
    """
    mkdir -p $dirgenome/idx/

    STAR --runThreadN $params.threads \\
      --runMode genomeGenerate \\
      --genomeDir $dirgenome/idx/ \\
      --genomeFastaFiles $genome \\
      --sjdbGTFfile $gtf \\
      --sjdbOverhang $params.lengthreads \\
      --genomeSAindexNbases 10
    """
}


process star_alignment {
    publishDir "$params.outdir/mapped-reads", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample), file(reads) 
    path index
    file gtf

    output:
    file "${sample}.mappings.bam", emit: samples_bam

    script:
    """
    mkdir -p $params.outdir/mapped-reads/
    STAR --runMode alignReads 
        --genomeDir ${index} \\  
        --runThreadN ${params.threads} \\
        --outFileNamePrefix $prefix_star
        --outSAMtype BAM SortedByCoordinate \\
        --sjdbGTFfile $gtf \\
        --sjdbOverhang $params.lengthreads
        --readFilesIn ${reads}
    """
}

// Running a workflow with the defined processes here.  
workflow IDX {
  take:
  dirgenome, genome, gtf

  main: 
  star_index(dirgenome, genome, gtf)
}

workflow MAP {
  take:
  reads, index, gtf

  main: 
  star_alignment(reads, index, gtf)

  output:
  samples_bam = star_alignment.out.samples_bam
}
