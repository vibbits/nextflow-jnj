# Creating our first pipeline
In this chapter we will build a basic RNA-seq pipeline consisting of quality controls, trimming of reads and mapping to a reference genome (and counting). We will build the pipeline step by step, starting from quality control with FastQC. We will very shortly introduce DSL1 in the beginning. After exploring Nextflow's flexibility on the quality control process, we will extend our pipeline script with the other processes while switching to the newer DSL2. 

The following script can be found and run in `03-first-pipeline/fastqc_1.nf`:
```
#!/usr/bin/env nextflow

params.reads = "$baseDir/../data/*.fq.gz"

/**
 * Quality control fastq
 */

read_ch = Channel
    .fromPath( params.reads )
    .view()

process fastqc_raw_reads {

    input:
    file read from read_ch
   
    script:
    """
    fastqc $read
    """
}
```

The first line of our script is always a shebang line, declaring the environment where the OS can find the software (i.e. Nextflow). Generally, the input files are first assigned into parameters which allows flexibility on running the nextflow pipeline. Input files are then assigned to channels and they serve as input for the process. 

- `$baseDir`: the folder from where the pipeline is run. 
- Flexibility in the language (writing of spaces, enters with channels, etc.) 

QUESTION: where are the output files (html- & zip-files)? 