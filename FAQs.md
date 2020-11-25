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

Q: Why use `.flatten()` on a channel?  
Flattening a channel will divide all files in that channel in different channels and distribute it over the infrastructure and hence create implicit parallelization. 

Q: Is it possible to pipe processes in a workflow? 
The pipe operator has been enabled in the DSL2. The concept is the same as in Linux pipes and allows to write a piped workflow in the workflow-process. 

Q: Is there a possibility to have an overview of all pipelines run in Nextflow?  
Use `nextflow log` to have an overview of all previously run pipelines. 

Q: How to disable the caching process and overwrite the results?  
Nextflow's cache can be disabled for a specific process by setting the `cache` parameter to any of the following values:
- `true`: By default, cache keys are created indexing input files meta-data information (name, size and last update timestamp attributes). Set to `false` if you want to disable caching. 
- `deep`: Cache keys are created indexing input files content.
- `lenient`: (Best in HPC and shared file systems) Cache keys are created indexing input files path and size attributes

Q: What is the difference between `script:` and `shell:`?   
