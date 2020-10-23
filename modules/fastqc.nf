#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process fastqc {

  input:
  tuple val(sample), file(reads) 

  script:
  """
  mkdir -p $params.outdir/quality-control-dsl2
  fastqc --outdir $params.outdir/quality-control ${reads}
  """
}
