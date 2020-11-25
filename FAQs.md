# 6. FAQs

Q: How to kill a pipeline? 
Look into the file `.nextflow.pid` and extract the number of your process. 
Use it to kill your process with `kill <number>`. 

Q: Add a help section for the pipeline (additional information)?
Use the following parameter on runtime: `nextflow run <pipeline.nf> --help`. In the main.nf script add the following snippet and extend with any additional information regarding the pipeline:
```
if (params.help) {
    log.info 'This is some information'
    log.info 'about the pipeline.'
    exit 1
}
```

