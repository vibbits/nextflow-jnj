#!/usr/bin/env nextflow

nextflow.preview.dsl=2

params.reads = "$baseDir/data/*0.01_{1,2}.fq.gz"
params.outdir = "$baseDir/results"


//include {fastqc} from './modules/fastqc'

process fastqc {

  input:
  tuple val(sample), file(reads)

  script:
  """
  mkdir -p $params.outdir/quality-control-dsl2
  fastqc --outdir $params.outdir/quality-control-dsl2 ${reads}
  """
}


workflow {
	read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )

	fastqc( read_pairs_ch) 
}
