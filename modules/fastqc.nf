#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

params.outdir = "$launchDir/results"

process fastqc {
  publishDir "$params.outdir/quality-control-$sample/", mode: 'copy', overwrite: true
    
  input:
  tuple val(sample), file(reads)

  output:
  path("*_fastqc.{zip,html}") 

  script:
  """
  fastqc ${reads}
  """
}

workflow QC {
    take: 
    input
    
    main:
		out =  fastqc(input)
    
    emit:
    out
}