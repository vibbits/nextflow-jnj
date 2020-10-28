# Nextflow 
This tutorial aims to get you familiarized with Nextflow. After this course you should be able to understand workflow pipelines that are written in Nextflow and write simple pipelines yourself! Here's an overview of the materials that we will cover:

- General introduction to Nextflow 
- Basic concepts: processes, channels, operators
- Creating our first Nextflow script
- Executors (local, queue vs cloud) (memory, time and CPU)
- Managing environment: Docker/Singularity/Conda

[comment]: <> (- Git integration: e.g. running nextflow pipelines from github directly)


After a brief introduction and example in DSL1, we will immediately change to the newer version of Nextflow with DSL2. It's important to have an idea of the differences between the two versions as most of the pipelines today are written in DSL1, but being transformed or enriched with new pipelines written in DSL2. 

## Prerequisites
This course requires familiarity with the command-line. You should feel confident interacting with the command-line and launch commands and tools from there. 

## Installations
Nextflow should be installed on your machine. Please consult the official documentation [here](https://www.nextflow.io/docs/latest/getstarted.html#installation).