#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

// Similar to DSL1, the input data is defined in the beginning.
params.reads = "$launchDir/data/*{1,2}.fq.gz"
params.outdir = "$launchDir/results"
params.threads = 2
params.slidingwindow = "SLIDINGWINDOW:4:15"
params.avgqual = "AVGQUAL:30"


println """\
      LIST OF PARAMETERS
================================
Reads            : $params.reads
Output-folder    : $params.outdir/
Threads          : $params.threads
"""

// Also channels are being created. 
read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists:true)

// Definition of a process, notice the absence of the 'from channel'.
// A process being defined, does not mean it's invoked (see workflow)
process fastqc_raw {
  publishDir "$params.outdir/quality-control-$sample/", mode: 'copy', overwrite: true
    
  input:
  tuple val(sample), file(reads)

  script:
  """
  mkdir -p $params.outdir/quality-control-$sample
  fastqc --outdir $params.outdir/quality-control-$sample ${reads}
  """
}


// Process trimmomatic
process trimmomatic {
    publishDir "$params.outdir/trimmed-reads", mode: 'copy'

    // Same input as fastqc on raw reads, comes from the same channel. 
    input:
    tuple val(sample), file(reads) 

    output:
    tuple val(sample), file("${sample}_1P.fq"), file("${sample}_2P.fq"), emit: paired_fq
    tuple val(sample), file("${sample}_1U.fq"), file("${sample}_2U.fq"), emit: unpaired_fq
    path "${sample}.trimmomatic.stats.log", emit: trimmomatic_stats

    script:
    """
    mkdir -p $params.outdir/trimmed-reads/
    trimmomatic PE -threads $params.threads ${reads[0]} ${reads[1]} -baseout ${sample} $params.slidingwindow $params.avgqual > ${sample}.trimmomatic.stats.log
    
    mv ${sample}_1P ${sample}_1P.fq
    mv ${sample}_2P ${sample}_2P.fq
    mv ${sample}_1U ${sample}_1U.fq
    mv ${sample}_2U ${sample}_2U.fq
    """
}

process fastqc_trim {
  publishDir "$params.outdir/quality-control-$sample/", mode: 'copy', overwrite: true
    
  input:
  tuple val(sample), file(reads)

  script:
  """
  mkdir -p $params.outdir/quality-control-$sample
  fastqc --outdir $params.outdir/quality-control-$sample ${reads}
  """
}

// Running a workflow with the defined processes here.  
// TO DO :  fix double fastqc by invoking fastqc once with collect and mix.
workflow {
	read_pairs_ch.view()
	fastqc_raw(read_pairs_ch) 
  paired_fq = trimmomatic(read_pairs_ch)
  fastqc_trim(paired_fq.mix())
}
