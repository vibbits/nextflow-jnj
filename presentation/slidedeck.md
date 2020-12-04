---
marp: true

theme: white

---

![nextflow](../img/nextflow-logo.png)
# Session Johnson & Johnson - 04 December 2020
Tuur Muyldermans - tuur.muyldermans@vib.be
Alexander Botzki - alexander.botzki@vib.be

<!--

# Objective
- Understand Nextflow syntax
- Understand workflow pipelines
- Write simple pipelines yourself!
-->

---

# Overview:
- Introduction
- Basic concepts: processes, channels and operators
- Creating our first Nextflow script(s)
- Managing configurations: parameters, portability, execution
- Creating reports



---

# 1. Introduction

<!-- 
In the first chapter we will elaborate on how Nextflow is designed, its advantages and disadvantages, the basic components, etc.
-->

---

## Bash scripts

```
#!/bin/bash

blastp -query sample.fasta -outfmt 6 \
	| head -n 10 \
	| cut -f 2 \
	| blastdbcmd -entry - > sequences.txt
```


<!--
Writing workflows to automate processes is not something new. In the data/ folder we've written a bash script that downloads the data that we will use throughout this tutorial. These bash scripts are probably one of the oldest forms of workflows. 


Starting with a shebang line, the `blastp` command is piped through multiple times to eventually result in an output file `sequences.txt`. The downside of this very basic and intuitive pipeline is that it has a sequential flow. In response to that, pipeline tools were built which are aimed to deal with more complex situations.
-->


---

<!--
Nextflow is designed around the idea that Linux has many simple but powerful command-line and scripting tools that, when chained together, facilitate complex data manipulations. First published in 2017 by CRG in Barcelona by Paolo di Tommaso. For obvious reasons it caught a lot of attention, a community grew and Nextflow developed to what it is now. 
-->

> **Nextflow** is a reactive workflow framework and a programming Domain Specific Language that eases the writing of data-intensive computational pipelines.

![language subset](../img/java-groovy-nextflow.png)


<!--
Nextflow scripting is an extension of the Groovy programming language, which in turn is a super-set of the Java programming language. Groovy can be considered as Python for Java in a way that simplifies the writing of code and is more approachable. 
-->

---
## Why (not)? (1/2)

+ Parallelization: processes are automatically scheduled based on available resources 
+ Scalability: simple scaling from local to HPC-cluster usage
+ Portability: run across different platforms
+ Reproducible: native support for containers, conda environments, and interaction with Git.
+ Continuous checkpoints for resuming / expanding pipelines (which is usually the case for workflow pipelines)
+ Re-usability: DSL2 and its modules will allow re-using of other scripts
+ Community: nf-core, Gitter, etc. 



---

## Why (not)? (2/2)

Alternatives: [link](https://github.com/pditommaso/awesome-pipeline/)
- Syntax of the Groovy language, yet another language
- Flexibility also comes with cost of complexity
- Nitpicking details in failure of scripts

<!--
Some thoughts or disadvantages from my personal point of view, it takes some time to get used to the syntax of the Groovy language. As flexible as it is, as complex it gets. Often it's difficult to trace down the exact problem of a failure of a pipeline script, especially in the beginning. It's probably not the first thing you should be concerned of if you're doing a one-time analysis. 
-->

---

# 2. Basic concepts

--- 

![width:900 ](../img/nextflow-conceptually.png)


<!--
Nextflow consists of three main components: channels, operators and processes. 
- *Channels*: connect processes/operators with each other.  
- *Operators*: transform the content of channels by applying functions or transformations. Usually operators are applied on channels to get the input of a process in the right format.  
- *Processes*: define the piece of script that is actually being run (e.g. an alignment process with STAR)  

Since the introduction of the new DSL2, *workflows* can be added to this list. This will be discussed in the next chapter.
-->

---
Example: inspect `02-basic-concepts/firstscript.nf`

---
## 2.1 Channels
- Input of the analysis is stored in a channel (files, strings, numbers, etc.) 
- unidirectional async queues that allows the processes to communicate with each other
- Channels can be used by operators or serve as an input for the processes


```
# Channel consisting of strings
strings_ch = Channel.from('This', 'is', 'a', 'channel')

# Channel consisting of a single file
file_ch = Channel.fromPath('data/sequencefile.fastq')

# Channel consisting of multiple files by using a wildcard *
multfiles_ch = Channel.fromPath('data/*.fastq')

# Channel consisting of multiple paired-end files by using wildcard * and options {x,y}
paired_ch = Channel.fromFilePairs('data/*{1,2}.fastq')
```
Further reading: [Nextflow's documentation](https://www.nextflow.io/docs/latest/channel.html#).

<!--
The input of the analysis is stored in a channel, these are generally files like sequencing, reference fasta, annotation files, etc. however the input can be of any kind like numbers, strings, lists, etc. To have a complete overview, we refer to the official documentation[[4](https://www.nextflow.io/docs/latest/channel.html#)]. Here are some examples of how a channel is being created:
-->

---

## 2.2 Operators
- Transform content of channels
- A plethora of operators exists, only a handful used extensively
- Examples: `.view()`, `.ifEmpty()`, `.splitFasta()`, `.print()`, etc. etc. etc.


---

- `collect`: e.g. when using a channel consisting of multiple independent files (e.g. fastq-files) and need to be assembled for a next process. 
```
Channel
    .from( 1, 2, 3, 4 )
    .collect()
    .view()

# outputs
[1,2,3,4]
```
Further reading: [Nextflow's documentation](https://www.nextflow.io/docs/latest/operator.html?highlight=view#)

---
- `mix`: e.g. when assembling items from multiple channels into one channel for a next process (e.g. multiqc)

```
c1 = Channel.from( 1,2,3 )
c2 = Channel.from( 'a','b' )
c3 = Channel.from( 'z' )

c1 .mix(c2,c3)

# outputs
1
2
3
'a'
'b'
'z'
```
Further reading: [Nextflow's documentation](https://www.nextflow.io/docs/latest/operator.html?highlight=view#)

---
## 2.3 Processes
<!--
Processes are the backbone of the pipeline. They represent each individual subpart of the analysis and contain hence one of the (many) processes in a pipeline: think about fastqc, trimmomatic, star or hisat alignment, counting. In the code-snippet below, you can see that it consists of a couple of blocks: directives, input, output, when clause and the script. 
-->

```
process < name > {

   [ directives ]

   input:
    < process inputs >

   output:
    < process outputs >

   when:
    < condition >

   [script|shell|exec]:
   < user script to be executed >

}
```


---
## Processes

- Executed independently
- Isolated from any other process
- FIFO queues

<!--
Each process is executed independently and isolated from any other process. They communicate via asynchronous FIFO queues, i.e. one process will wait for the output of another and then runs reactively when the channel has contents. 
-->

--- 

## 2.4 Running our first pipeline:

```
nextflow run firstscript.nf
```
<!--
In our case, we will replace `example.nf` with `02-basic-consepts/firstscript.nf`. First, inspect the script `02-basic-consepts/firstscript.nf` and notice how the channels are being created, passed on to the process' inputs, processed by the script section and then given to the output. 
--> 

Output:

<!--When we run this script, the result file will not be present in our folder structure. Question: look at the output... Can you guess where to find the result? 
-->

```
N E X T F L O W  ~  version 20.07.1
Launching `02-basic-concepts/firstscript.nf` [elegant_curie] - revision: 9f886cc00a
executor >  local (2)
executor >  local (2)
[5e/195314] process > valuesToFile (2) [100%] 2 of 2 ✔
results file: /path/to/work/51/7023ee62af2cb4fdd9ef654265506a/result.txt
results file: /path/to/work/5e/195314955591a705e5af3c3ed0bd5a/result.txt
```
<!--
The output consists of:
- Version of nextflow 
- Information regarding the script that has ran with an identifier name
- Hash with process ID, progress and caching information
- Optional output printed to the screen as defined in the script (if present)
The results are stored in the results file as described in the two last lines. By default the results of a process are stored in the `work/` directory in subfolders with names defined by the hashes.
-->

---

Output-files stored in the work-directory. 

Besides the output, also a bunch of hidden `.command.*` files are present:
```
-... user group    0 Nov 26 15:20 .command.begin*
-... user group 1797 Nov 26 15:20 .command.err*
-... user group 1826 Nov 26 15:20 .command.log*
-... user group    0 Nov 26 15:20 .command.out*
-... user group 3187 Nov 26 15:20 .command.run*
-... user group   53 Nov 26 15:20 .command.sh*
-... user group    3 Nov 26 15:20 .exitcode*
```

<!--
- .exitcode, contains 0 if everything is ok, another value if there was a problem.
- .command.log, contains the log of the command execution. Often is identical to .command.out
- .command.out, contains the standard output of the command execution
- .command.err, contains the standard error of the command execution
- .command.begin, contains what has to be executed before .command.sh
- .command.sh, contains the block of code indicated in the process
- .command.run, contains the code made by nextflow for the execution of .command.sh and contains environmental variables, eventual invocations of linux containers etc
-->
---
## FIFO-principle 
```
nextflow run 02-basic-consepts/fifo.nf
```
Output:
```
N E X T F L O W  ~  version 20.07.1
Launching `02-basic-concepts/fifo.nf` [nauseous_mahavira] - revision: a71d904cf6
[-        ] process > whosfirst [  0%] 0 of 2
This is job number 6
This is job number 3
This is job number 7
This is job number 8
This is job number 5
This is job number 4
This is job number 1
This is job number 2
This is job number 9
executor >  local (10)
[4b/aff57f] process > whosfirst (10) [100%] 10 of 10
```
<!--
Earlier, we described that Nextflow uses an asynchronous FIFO principle. Let's exemplify this by running the script `02-basic-consepts/fifo.nf` and inspect the order that the channels are being processed. 

Note also the implicit parallelisation *.fastq in a channel (one channel) will spit it out over multiple processes simultaneously. No need of making a fors–loop.
-->

---
## Self-written scripts
<!--
A script, as part of the process, can be written in any language (bash, Python, Perl, Ruby, etc.). This allows to add self-written scripts in the pipeline. 
They can be defined in the process script itself as given in the example here, or they could also exist as a script in another folder and be run here like python script.py
--> 
- Any language (bash, Python, Perl, Ruby, etc.)
- Defined in the process or command to run the script 
```
#!/usr/bin/env nextflow
 
process python {
    
    """
    #!/usr/bin/python3

    firstWord = 'hello'
    secondWord = 'folks'
    print(f'{firstWord} {secondWord}')
    """
}
```

---
# 3. Creating our first pipeline
<!--
In this chapter we will build a basic RNA-seq pipeline consisting of quality controls, trimming of reads and mapping to a reference genome (excl. counting). We will build the pipeline step by step, starting from quality control with FastQC. We will start with DSL1 for this process. After exploring Nextflow's flexibility on the quality control process, we will switch to the newer DSL2 language and extend our pipeline script with the other processes. 
-->

---

![RNAseq.PNG](../img/RNAseq.PNG)


---
## 3.1 Using DSL1 (example: FastQC)

Inspect: `03-first-pipeline/fastqc_1.nf`

- Shebang line
- Assign input into `params`
- Comments are always useful
- Create a channel for input files - serve as input for process
- Operator on a channel `.view()`
<!--
The first line of our script is always a shebang line, declaring the environment where the OS can find the software (i.e. Nextflow). Generally, the input files and parameters of the processes are first assigned into *parameters* which allows flexibility in the pipeline. Input files are then assigned to channels and they serve as input for the process. 
-->

<!--
Note:  
- `$launchDir`: The directory where the main script is located (replaces `$baseDir` in version >20). 
- Flexibility in the language (writing of spaces, enters with channels, asigning channel values to a variable or using `.set{}`, etc.) 
-->

---

**Implicit parallellisation**:  
Run the following and keep an eye on the workload distribution (`htop`). 
```
nextflow run 03-first-pipeline/fastqc_1.nf
```

---
Let's add some new features:
1. `nextflow run 03-first-pipeline/fastqc_1.nf -bg > log`
<!--
ANSWER: run in the background and push output of nextflow to the log file. No need of explicitly using nohup, screen or tmux.
ANSWER: the hash at the beginning of each process reveals where you can find the result of each process.
-->

---

Let's add some new features:
1. `nextflow run 03-first-pipeline/fastqc_1.nf -bg > log`
2. Adapt file for handling read pairs. Which parameter would you use on runtime to overwrite the inputfiles?

<!--
ANSWER: nextflow run --reads data/WT_lib1_R1.fq.gz 03-first-pipeline/fastqc_1.nf .-->

---
Let's add some new features:
1. `nextflow run 03-first-pipeline/fastqc_1.nf -bg > log`
2. Adapt file for handling read pairs. Which parameter would you use on runtime to overwrite the inputfiles?
3. Print parameters using `println` & check if the files exist when creating the channels.
<!--
ANSWER: println section to describe your workflow
ANSWER: checkIfExists
ANSWER: invoke the checkIfExists-error by running the nextflow script with wrong reads: `nextflow run 03-first-pipeline/fastqc_3.nf --reads wrongfilename`.
-->

---
Let's add some new features:
1. `nextflow run 03-first-pipeline/fastqc_1.nf -bg > log`. Where are the output files?
2. Adapt file for handling read pairs. Which parameter would you use on runtime to overwrite the inputfiles?
3. Print parameters using `println` & check if the files exist when creating the channels.
4. Create a directory where the files can be stored with `publishDir`.
<!--
ANSWER: If the output is to be used by another process, and the files are being moved, they won't be accessible for the next process and hence you're pipeline will fail complaining about files not being present. Files published by a process must not be accessed by other downstream processes

DIFFERENT MODES:
 - Without any additional arguments, a hyperlink will be created to the files stored in the `work/` directory, with mode set to copy (`mode: 'copy'`) the files will be made available in the defined directory. 
  - Can you guess what might happen if we set the mode to move? (`mode: 'move'`)  
  - Can you guess what might happen if we set `pattern` to `'*.zip'`?

ANSWER: Only zip files will be created to the output. Patterns will specify what will be outputted. ie zip and html. So we can control which files are being outputted in the publishdir.
-->

---
## 3.2 Moving towards DSL2
- Make the pipelines more modular 
- Simplify the writing of complex data analysis pipelines

---

**Overview of the changes (DSL2 vs DSL1)**:
- `nextflow.enable.dsl=2`

- In DSL2: 
  - Channel can be used indefinitely (vs only once in DSL1) 
  - Process can (still) only be used once
  - DSL2 separates the definition of a process from its invocation
  - Within processes no more references to channels (i.e. `from` and `into`)
- Introduction of `workflow` & modules

<!--
Here is a list of the major changes: 
- Following the shebang line, the nf-script wil start with the following line: `nextflow.enable.dsl=2` (not to be mistaken with *preview*).
- When using DSL1 each channel could only be consumed once, this is ommited in DSL2. Once created, a channel can be consumed indefinitely. 
- A process on the other hand can still only be used once in DSL2 
- A new term is introduced: `workflow`. In the workflow, the processes are called as functions with input arguments being the channels. 
- Regarding the processes, the new DSL separates the definition of a process from its invocation. This means that in DSL1 the process was defined and also run when the script was invoked, however in DSL2, the definition of a process does not necessarily mean that it will be run. 
- Moreover, within processes there are no more references to channels (i.e. `from` and `into`). The channels are passed as inputs to the processes which are defined and invoked in the `workflow`. 
-->
---

### Quality control with `FastQC` (DSL2)
Inspect: `03-first-pipeline/dsl2-fastqc.nf`

```
// Running a workflow with the defined processes here.  
workflow {
	read_pairs_ch.view()
	fastqc(read_pairs_ch) 
}
```

---

### Trimming with `trimmomatic` (DSL2)
Inspect: `03-first-pipeline/dsl2-trimming.nf`

- Introducing: output and `emit`. 
- Accessing an output of a process with `<processname>.out.<emitname>`
<!--Defining a process output with `emit` allows us to use it as a channel in the external scope. -->

```
nextflow run 03-first-pipeline/dsl2-trimming.nf
```

---

### Quality control on trimmed reads with `FastQC` (DSL2)
Rerun and uncomment the last fastqc-process in the workflow:
```
nextflow run 03-first-pipeline/dsl2-trimming.nf
```
Output:
```
Error: Process fastqc has been already used -- 
If you need to reuse the same component include it with a 
different name or include in a different workflow context
```

<!--
At this point we're interested in the result of the `trimmomatic` process. Hence, we want to verify the quality of the reads with another `fastqc` process. Re-run `fastqc` on the filtered read sequences by adding it in the workflow of `03-first-pipeline/dsl2-trimming.nf`. Use the parameter `-resume` to restart the pipeline from where it stopped the last time. 

Hmm, error? `Process fastqc has been already used -- If you need to reuse the same component include it with a different name or include in a different workflow context`. It means that processes can only be used once in a workflow. This means that we need to come up with a smarter solution (see below). 
--> 

---

## 3.3 Sub-workflows 

- Sub-workflows with workflow keyword definition
- Allows use of sub-workflows within workflow
- Main workflow does not carry a name and is executed implicitly



<!--
The workflow keyword (which is basically a name for a workflow)  allows the definition of **sub-workflow** components that enclose the invocation of one or more processes and operators. It also allows you to use this workflow from within another workflow. The workflow that does not cary any name is considered to be the main workflow and will be executed implicitly. An alternative workflow entry can be specified using the -entry command line option. 
-->

---
```
workflow star{
  take:
  arg1
  arg2
  arg3

  main:
  star_index(arg1, arg2)
  star_alignment(arg1, arg2, arg3)
}

workflow hisat2{
  take:
  arg1
  arg2

  main:
  hisat_index(arg1)
  hisat_alignment(arg1, arg2)
}

workflow {
  star(arg1, arg2, arg3)
  hisat2(arg1, arg2)
}

```
<!--
- An alternative workflow entry can be specified using `-entry`
-->

---

## 3.4 Modules

<!--
However, if we want to be truly modular, we can write a library of modules and import a specific component from that library. A module can contain the definition of a function, process and workflow definitions. Navigate to the modules folder and find a script called `fastqc.nf`. This script consists of a process and a workflow. This module can be imported into our pipeline script (main workflow) like this:
-->
- Write components in a module (component = process, workflow, functions)
- Import a specific component from that module: 
```
include {QC} from './modules/fastqc.nf'
```
<!--
This line is quite specific. The workflow is defined within the curly brackets, the origin of the module defined by a relative path must start with `./`. -->


Example: `modules/fastqc.nf`. 

```
include { QC as fastqc_raw; QC as fastqc_trim } from "${launchDir}/modules/fastqc")
```

<!--
When including a module component it’s possible to specify a name alias. This allows the inclusion and the invocation of the same component multiple times in your script using different names. Now we're ready to use a process, defined in a module, multiple times in a workflow. 

-->

---
Inspect `03-first-pipeline/dsl2-subworkflow.nf`:
- Fastqc process has been removed from this script
- Fastqc is imported from `modules/fastqc.nf`
- Trimmomatic process is still internally described
```
nextflow run 03-first-pipeline/dsl2-subworkflow.nf
```

<!--
Run `03-first-pipeline/dsl2-subworkflow.nf` which contains the `trimmomatic` process internally and imports the `fastqc` process from the modules library `module/fastqc.nf`. -->

---
`03-first-pipeline/dsl2-subworkflow.nf`:  

![dsl2-subworkflow](../img/dsl2-subworkflow.PNG)

---

# RNAseq workflow example
![RNAseq.PNG](../img/RNAseq.PNG)

<!--
Similarly as described above, we can extend this pipeline and map our trimmed reads on a reference genome. First, we'll have to index our reads and afterwards we can map our reads. In the folder `modules/` find the script `star.nf` which contains two processes: `star_index` and `star_alignment`. 
These modules are called from the main script `03-first-pipeline/dsl2-RNAseq.nf`. 

This pipeline is still subject to optimizations which will be elaborated in the next section. 
-->

---

# 4. Managing configurations

---

## Configuration files  
- Pipeline configuration properties are defined in `nextflow.config`
- Technical parameters (profiles):
  - Executor, CPUs, memory, etc. 
- Pipeline specific parameters:
  - Input files, process related parameters (e.g. trimmomatic or STAR)
- Separate these variables from pipeline - more modular  

```
params.reads = "$launchDir/data/*{1,2}.fq.gz"

process {
    memory='1G'
    cpus='1'
}
```

<!--

Pipeline configuration properties are defined in a file named `nextflow.config` in the pipeline execution directory. This file can be used to define technical and project depeding parameters, e.g. which executor to use, the processes' environment variables, pipeline parameters etc. Hence, the configuration file allows to separate these variables from the script and makes the scripts more flexible depending on the project you're using it for and the infrastructure you're running it on.  

In the example below we start with defining the processes' allowed memory- and cpu-usage. This list can be extended with parameters that are more relevant for HPC usage: time, queue, executor, etc. 
-->

--- 

- Add labels that allow different resources for a particular process

```
// Define technical resources below:
process {
    withLabel: 'low' {
        memory='1G'
        cpus='1'
        time='6h'
    }
    withLabel: 'med' {
        memory='2G'
        cpus='2'
    }
    withLabel: 'high' {
        memory = '8G'
        cpus='8'
    }
}
```
<!--
It's also possible to create labels that can be chosen and used for each process separately. In the example below we can use the label `high` as a directive in a process and hence allow more resources for that particular process. 

-->

--- 

- Separate analysis parameteres in a separate file
```
includeConfig "/path/to/params.config"
```
<!--
Imagine that you want to separate analysis parameters in a separate file, this is possible by creating a `params.config` file and include it in the `nextflow.config` file as such: 
-->

---
## Portability & reproducibility 
- Support for Conda, Docker & Singularity 
<!--
As discussed before, Nextflow is especially useful thanks to its portability and reproducibility, i.e. the native support for containers and environment managers. There are two options for attaching containers to your pipeline. Either you define a dedicated container image for each process individually, or you define one container for all processes together in the configurations file. 
-->


---
**Using Docker**:
<!--
In the former case, simply define the container image name in the process directives. In the snippet below, we defined a container that already exists in [DockerHub](https://hub.docker.com/r/biocontainers/fastqc). Dockerhub is also the default location where Nextflow will search for the existence of this container if it doesn't exist locally. -->
1. In process directives:
```
process quality-control {
    container 'biocontainers/fastqc:v0.11.9_cv7'

    """
    fastqc ...
    """
}
```
2. In `nextflow.config` file:
```
process.container = 'vibbioinfocore/analysispipeline:latest'
```
<!--
We're referring to a Docker container image that exists on Dockerhub. Notice however that all the tools and dependencies necessary during your pipeline, need to be present in this image. To run the pipeline script with this Docker container image, use the following command: `nextflow run example.nf -with-docker`. 
-->

---

- Run pipeline with Docker container: `nextflow run example.nf -with-docker`
- Or add the following to `nextflow.config`-file: 
`docker.enabled = true`


Note: to set the correct user- and group-settings: `docker.runOptions = '-u \$(id -u):\$(id -g)'`

<!--
Ultimately, the parameter `-with-docker` does not need to be defined and it should use the Docker container in the background at all times, for this purpose also use the `docker.enabled = true` option in the config file. Another interesting parameter to consider addin to the configuration file is the `docker.runOptions = '-u \$(id -u):\$(id -g)'`. This allows us to create files with permissions on user-level instead of the default root-level files. -->

--- 

**Using Singularity**:

```
nextflow run example.nf -with-singularity [singularity-image-file]
```
<!--
Similarly with a singularity image for which you do not have to adapt the pipeline script. You can run with Singularity container using the following command-line parameter: `-with-singularity [singularity-image-file]`, where the image is downloaded from Dockerhub as well, built on runtime and then stored in a folder `singularity/`. Re-using a singularity image is possible with:
-->

Or extend `nextflow.config`-file with: 
```
singularity.cacheDir = "/path/to/singularity" // centralised caching directory
process.container = 'singularity.img' // define the image
singularity.enabled = true  // enable running with singularity 
```
<!--
If you want to avoid entering the Singularity image as a command line parameter, you can define it in the Nextflow configuration file. For example you can add the following lines in the `nextflow.config` file:
-->

---
## Executors
![executors](../img/executors-schedulers.png)

<!--
While a *process* defines *what* command or script has to be executed, the *executor* determines *how* that script is actually run on the target system. In the Nextflow framework architecture, the executor is the component that determines the system where a pipeline process is run and it supervises its execution.

If not otherwise specified, processes are executed on the local computer. The local executor is very useful for pipeline development and testing purposes, but for real world computational pipelines an HPC or cloud platform is often required.

Hence, you can write your pipeline script once and have it running on your computer, a cluster resource manager or the cloud by simply changing the executor definition in the Nextflow configuration file. As these configurations are often a one-time effort, managed by a local IT/admin person, we refer to the [official documentation](https://www.nextflow.io/docs/latest/executor.html). 
-->

---
# 5. Creating reports
<!-- 
Nextflow has an embedded function for reporting a number of informations about the resources needed by each job and the timing. Just by adding a parameter on run-time, different kind of reports can be created. 
-->

---
1. Workflow report (html)
<!--
After running the nextflow pipeline script with the option `-with-report`, find the html report in the folder from where you launched the pipeline. 
-->
```
nextflow run example.nf -with-docker -with-report
```
<!--
This report describes the usage of resources and job durations and give an indication of bottlenecks in the pipeline. 
-->

2. DAG: visualization of the pipeline (dependency: graphviz)
<!--Use the option `-with-dag` to create a visualization of the workflow. By default it will create a `.dot`-file that contains a description of the workflow, however use e.g. `rnaseq.PNG` as an argument to create a figure. This visualization is a nice overview of the workflow processes and how they are linked together and can be especially useful as a starting point to unravel more complex pipelines.
-->
```
nextflow run example.nf -with-dag <filename.PNG>
```
---

# Questions

---

# Further reading & references:

- Nextflow's official documentation ([link](https://www.nextflow.io/docs/latest/index.html))
- Reach out to the community on Gitter ([link](https://gitter.im/nextflow-io/nextflow))
- Curated collection of patterns ([link](https://github.com/nextflow-io/patterns))
- Workshop focused on DSL2 developed by CRG Bioinformatics Core ([link](https://github.com/biocorecrg/ELIXIR_containers_nextflow))
- Tutorial exercises (DSL1) developed by Seqera ([link](https://github.com/seqeralabs/nextflow-tutorial))
- Curated ready-to-use analysis pipelines by NF-core ([link](https://nf-co.re/))
- Model example pipeline on Variant Calling Analysis with NGS RNA-Seq data developed by CRG ([link](https://github.com/CRG-CNAG/CalliNGS-NF))