## Creating reports 
Nextflow has an embedded function for reporting a number of informations about the resources needed by each job and the timing. Just adding a parameter will give you a nice html report.


1. Report
After running the nextflow pipeline script with the option `-with-report`, find the html report in the folder from where you launched the pipeline. 
```
nextflow run 05-reports/dsl2-RNAseq.nf -with-docker -bg -with-report -resume > log
```

2. DAG 
Use the option `-with-dag` to create a visualization of the workflow. By default it will create a `.dot`-file that contains a description of the workflow, however use e.g. `rnaseq.PNG` as an argument to create a figure. This visualization is a nice overview of the workflow processes and how they are linked together and can be especially useful as a starting point to unravel more complex pipelines.
```
nextflow run 05-reports/dsl2-RNAseq.nf -with-docker -bg -with-dag -resume > log
```
(Library graphviz needs to be installed for this purpose.)