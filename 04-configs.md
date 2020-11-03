# Managing configurations

## Configuration files  
Pipeline configuration properties are defined in a file named `nextflow.config` in the pipeline execution directory. This file can be used to define which executor to use, the processes' environment variables, pipeline parameters etc. In the example below we start with defining the processes' allowed memory- and cpu-usage. This list can be extended with parameters that are more relevant for HPC usage: time, queue, executor, etc. 

```
process {
    memory='1G'
    cpus='1'
}
```
It's also possible to create labels that can be chosen and used for each process separately. In the example below we can use the label `high` as a directive in a process and hence allow more resources for that particular process. 

```
process {
    memory='1G'
    cpus='1'
    withLabel: 'high'	
   	{ 
		memory='2G'
   	 	cpus='2'
	} 
}
```

Imagine that you want to separate analysis parameters in a separate file, this is possible by creating a `params.config` file and include it in the `nextflow.config` file as such: 
```
includeConfig "/path/to/params.config"
```
## Portability

As discussed before, Nextflow is especially useful thanks to its portability and reproducibility, i.e. the native support for containers and environment managers. There are two options for attaching containers to your pipeline. Either you define a dedicated container image for each process individually, or you define one container for all processes together in the configurations file. 

In the former case, simply define the container image name in the process directives. Here we defined a container that already exists in [DockerHub](https://hub.docker.com/r/biocontainers/fastqc):
```
process quality-control {
    container 'biocontainers/fastqc:latest'

    """
    fastqc ...
    """
}
```

In the latter case, write the following line in the `nextflow.config` file:
```
process.container = 'rnaseq:latest'
```
We're referring to a Docker container image that already exists locally on our computer, nl. `rnaseq:latest`. It was created by building it from the Dockerfile: `docker build -t rnaseq:latest .`. Notice however that all the tools and dependencies necessary during your pipeline, need to be present in this image. To run your pipeline script with this Docker container image, use the following command: `nextflow run example.nf -with-docker`. 

Similarly with a singularity image:
```
singularity.cacheDir = "/path/to/singularity"
```

## Executors
While a *process* defines *what* command or script has to be executed, the *executor* determines *how* that script is actually run on the target system. In the Nextflow framework architecture, the executor is the component that determines the system where a pipeline process is run and it supervises its execution.

If not otherwise specified, processes are executed on the local computer. The local executor is very useful for pipeline development and testing purposes, but for real world computational pipelines an HPC or cloud platform is often required.

Hence, you can write your pipeline script once and have it running on your computer, a cluster resource manager or the cloud by simply changing the executor definition in the Nextflow configuration file. As these configurations are often a one-time effort, managed by a local IT/admin person, we refer to the [official documentation](https://www.nextflow.io/docs/latest/executor.html). 

![executors](img/executors-schedulers.png)

