#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

process star_index {
    publishDir "$params.dirgenome", mode: 'copy', overwrite: true
    
    input:
    file genome
    file gtf

    //output: 
    //path "$params.dirgenome/idx", emit: dir_index

    script:
    """
    mkdir -p $params.dirgenome 

    STAR --runThreadN $params.threads \\
      --runMode genomeGenerate \\
      --genomeDir $params.dirgenome \\
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
    file gtf

    output:
    //file "${sample}.mappings.bam", emit: samples_bam

        //  output:
        //  set file("*Log.final.out"), file ('*.bam') into star_aligned
        //  file "*.out" into alignment_logs
        //  file "*SJ.out.tab"
        //  file "*Log.out" into star_log
        //  file "where_are_my_files.txt"
        //  file "*Unmapped*" optional true
        //  file "${prefix}Aligned.sortedByCoord.out.bam.bai" into bam_index_rseqc, bam_index_genebody

    script:
    """
    mkdir -p $params.outdir/mapped-reads/
    STAR --runMode alignReads 
        --genomeDir $params.dirgenome \\  
        --runThreadN ${params.threads} \\
        --outFileNamePrefix ${sample}
        --outSAMtype BAM SortedByCoordinate \\
        --sjdbGTFfile $gtf \\
        --sjdbOverhang $params.lengthreads
        --readFilesIn ${reads[0]} ${reads[1]} 
    """
}

// Running a workflow with the defined processes here.  
workflow IDX {
  take:
  genome
  gtf

  main: 
  star_index(genome, gtf)
}

workflow MAP {
  take:
  reads
  gtf

  main: 
  star_alignment(reads, gtf)

//  emit:
//  samples_bam = star_alignment.out.samples_bam
}
