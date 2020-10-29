# Creating our first pipeline
In this chapter we will build a basic RNA-seq pipeline consisting of quality controls, trimming of reads and mapping to a reference genome (and counting). We will build the pipeline step by step, starting from quality control with FastQC. We will very shortly introduce DSL1 in the beginning. After exploring Nextflow's flexibility on the quality control process, we will extend our pipeline script with the other processes while switching to the newer DSL2. 

## Quality control in DSL1

The following script can be found and run in `03-first-pipeline/fastqc_1.nf`:
```
#!/usr/bin/env nextflow

params.reads = "$launchDir/data/*.fq.gz"

/**
 * Quality control fastq
 */

reads_ch = Channel
    .fromPath( params.reads )
    .view()
    
process fastqc_raw_reads {

    input:
    file read from reads_ch 
    
    script:
    """
    fastqc ${read}
    """
}
```

The first line of our script is always a shebang line, declaring the environment where the OS can find the software (i.e. Nextflow). Generally, the input files are first assigned into parameters which allows flexibility flexibility in the pipeline. Input files are then assigned to channels and they serve as input for the process. 

Note:  
- `$launchDir`: The directory where the main script is located (version >20, replaces `$baseDir`). 
- Flexibility in the language (writing of spaces, enters with channels, asigning channel values to a variable or using `.set{}`, etc.) 

---
QUESTION:  
Run script with the following line: `nextflow run 03-first-pipeline/fastqc_1.nf -bg > log`.  
What does the `-bg > log` mean?  
FastQC generates an html- and zip-file for each read. Where are the output files?  

ANSWER: the hash at the beginning of each process reveals where you can find the result of each process. 

---


--- 
Exercise: 
- Adapt file for handling read pairs. Result: [`fastqc_2.nf`]
- Print parameters using `println` & check if the files exist when creating the channels [`checkIfExists`](https://www.nextflow.io/docs/latest/channel.html?highlight=fromfilepairs). Invoke the checkIfExists-error by running the nextflow script with wrong reads:
`nextflow run 03-first-pipeline/fastqc_3.nf --reads wrongfilename`. 


- Create a directory where the files can be stored (hint: [`publishDir`](https://www.nextflow.io/docs/latest/process.html?highlight=publishdir#publishdir))

**Warning**: Files are copied into the specified directory in an asynchronous manner, thus they may not be immediately available in the published directory at the end of the process execution. For this reason files published by a process must not be accessed by other downstream processes.

## Moving towards DSL2
Nextflow recently went through a big make-over. The premise of the next version, using DSL2, is to make the pipelines more modular and simplify the writing of complex data analysis pipelines. As an example, have a look at [this complex pipelines](https://github.com/nf-core/rnaseq/blob/master/main.nf) as part of the [nf-core](https://nf-co.re/pipelines) pipelines. 

The nf-script wil start with the following line:
```
nextflow.enable.dsl=2
```

When using DSL1 each channel could only be consumed once, this is ommited in DSL2. Once created, a channel can be consumed indefinitely. 

A new term is introduced now: `workflow`. In the workflow, the processes are called as functions with input arguments being the channels. 

Regarding the processes, the new DSL separates the definition of a process from its invocation. This means that in DSL1 the process was defined and also run when it the script was invoked, however in DSL2, the definition of a process does not necessarily mean that it will be run. Moreover, within processes there are no more references to channels (i.e. `from` and `into`). The channels are passed as inputs to the processes which are defined and invoked in the `workflow`. 

Let's have a look (`nextflow run 03-first-pipeline/dsl2-fastqc.nf`):

```
#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

// Similar to DSL1, the input data is defined in the beginning.
params.reads = "$launchDir/data/*{1,2}.fq.gz"
params.outdir = "$launchDir/results"

println """\
      LIST OF PARAMETERS
================================
Reads            : $params.reads
Output-folder    : $params.outdir/
"""

// Also channels are being created. 
read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists:true)

// Definition of a process, notice the absence of the 'from channel'.
// A process being defined, does not mean it's invoked (see workflow)
process fastqc {
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
workflow {
	read_pairs_ch.view()
	fastqc(read_pairs_ch) 
}

```


https://www.nextflow.io/docs/latest/docker.html