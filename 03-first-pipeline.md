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
Nextflow recently went through a big make-over. The premise of the next version, using DSL2, is to make the pipelines more modular and simplify the writing of complex data analysis pipelines. 

The nf-script wil start with the following line:
```
nextflow.enable.dsl=2
```

When using DSL1 each channel could only be consumed once, this is ommited in DSL2. Once created, a channel can be consumed indefinitely. 

A new term is introduced now: `workflow`. In the workflow, the processes are called as functions with input arguments being the channels. 

Regarding the processes, the new DSL separates the definition of a process from its invocation. This means that in DSL1 the process was defined and also run when it the script was invoked, however in DSL2, the definition of a process does not necessarily mean that it will be run. Moreover, within processes there are no more references to channels (i.e. `from` and `into`). The channels are passed as inputs to the processes which are defined and invoked in the `workflow`. 

Let's have a look (`nextflow run 03-first-pipeline/dsl2-fastqc.nf`):

```nextflow
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

Now we will add the next step in our pipeline, which is **trimming and filtering the low quality reads**. For this process, we will use the tool `trimmomatic`. The solution is available in `03-first-pipeline/dsl2-trimming.nf`. Here we're introducing a new option: `emit`. Defining a process output with `emit` allows us to use it as a channe lin the external scope.  


At this point we're interested in the result of the `trimmomatic` process. Hence, we want to verify the quality of the reads with another `fastqc` process. Rerun `fastqc` on the filtered read sequences by adding it in the workflow of `03-first-pipeline/dsl2-trimming.nf`. Use the parameter `-resume`. 
  
- Hmm, error? `Process fastqc has been already used -- If you need to reuse the same component include it with a different name or include in a different workflow context`. It means that processes can only be used once in a workflow. This means that we need to come up with a smarter solution. 

## Subworkflows and modules
The workflow keyword allows the definition of **sub-workflow** components that enclose the invocation of one or more processes and operators. It also allows you to use this workflow from within another workflow. The workflow that does not cary any name is considered to be the main workflow and will be executed implicitly. 

However, if we want to be truly modular, we can write a library of modules. A module can contain the definition of a function, process and workflow definitions. Navigate to the modules folder and find a script called `fastqc.nf`. This script consists of a process and a workflow. This module can be imported into our pipeline script (main workflow) like this:
```
include {QC} from './modules/fastqc.nf'
```
This line is quite specific. The workflow is defined within the curly brackets, the origin of the module defined by a relative path must start with `./`. 

When including a module component it’s possible to specify a name alias. This allows the inclusion and the invocation of the same component multiple times in your script using different names. For example:
``` 
include { QC as fastqc_raw; QC as fastqc_trim } from "${launchDir}/modules/fastqc" //addParams(OUTPUT: fastqcOutputFolder)
```
Now we're ready to use a process, defined in a module, multiple times in a workflow. 

- Hmm, error? `Workflow 'QC' declares 1 input channels but 2 were specified`. Notice that it won't work because twice same process, so do collect/mix. 


## Configuration files, executors and portability
Pipeline configuration properties are defined in a file named `nextflow.config` in the pipeline execution directory. This file can be used to define which executor to use, the processes' environment variables, pipeline parameters etc. In the example below we start with defining the processes' allowed memory, cpu-usage and execution time. 

```
process {
     memory='1G'
     cpus='1'
     time='6h'
}
```

Imagine that you want to separate analysis parameters in a separate file, this is possible by creating a `params.config` file and include it in the `nextflow.config` file as such: 
```
includeConfig "/path/to/params.config"
```

While a *process* defines *what* command or script has to be executed, the *executor* determines *how* that script is actually run on the target system. In the Nextflow framework architecture, the executor is the component that determines the system where a pipeline process is run and it supervises its execution.

If not otherwise specified, processes are executed on the local computer. The local executor is very useful for pipeline development and testing purposes, but for real world computational pipelines an HPC or cloud platform is often required.

In other words you can write your pipeline script once and have it running on your computer, a cluster resource manager or the cloud by simply changing the executor definition in the Nextflow configuration file.

![executors](img/executors-schedulers.PNG)


As discussed before, Nextflow is especially useful thanks to its portability, i.e. the native support for containers and environment managers. There are several options for attaching containers to your pipeline. Either you define a dedicated container for each process individually, or you define one container for all processes together. 

```
process.container = 'biocontainer/example:latest'
singularity.cacheDir = "/path/to/singularity"
```

https://www.nextflow.io/docs/latest/docker.html


## Modules
Put mapping into modules. 

## Publishing final results
- DAG
- Report
