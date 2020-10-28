#!/usr/bin/env nextflow

params.reads = "$baseDir/data/*0.01_{1,2}.fq.gz"
// params.outdir = "$baseDir/nextflow"
params.env = "$params.abspath/environment.yml"


 Channel
         .fromFilePairs(params.reads) // standard option for paired-end reads
         .ifEmpty { 'No read pairs found' } // if no reads: raise error
         .into { read_pairs1_ch; read_pairs2_ch } // set two channels: fastqc and mapping

genome_file = file(params.genome)

/**
 * Quality control fastq
 */
process fastqc_raw_reads {
  conda "$params.env"
  // publishDir "$params.outdir/quality-control", mode: 'copy', overwrite: true

  input:
  set sample_id, file(reads) from read_pairs1_ch

  script:
  """
  mkdir -p $params.outdir/quality-control
  fastqc --outdir $params.outdir/quality-control ${reads}
  """
}