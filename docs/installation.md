# Installation and configuration 

## Compute infrastructure

This pipeline is designed to run on a distributed compute system, such as a traditional HPC cluster system. 
We have tested the pipeline on two Slurm clusters, with node configurations of 16cores/128GB Ram and 20cores/256GB Ram, respectively. 

While smaller nodes will probably work, it may require some tweaking on your end. Most importantly, if you plan on using the transcriptome 
assembly or the whole genome alignment branches of the pipeline, available memory may become limiting. For transcriptome assembly, 128GB should be fine in most cases. 
Genome alignments will easily consume up to 500GB of RAM. If you cluster does not support those, it's best to not turn on this part. 

## Installing Nextflow

The pipeline is built using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool to run tasks across 
distributed compute systems in a very portable manner. Therefore the first thing to do is to install Nextflow. 

Nextflow runs on most systems (Linux, Mac OSX etc). It can be installed by running the following commands:

```
### Make sure that Java v8+ is installed:
java -version

```

### Install Nextflow

Get Netxtflow from [here](https://github.com/nextflow-io/nextflow/releases)

Just put the downloaded file into your $PATH and rename it to 'nextflow'

## Running the code directly from github

In order to run the pipeline from github directly, you can do the following:

```bash
nextflow -c nextflow.config ikmb/esga -params-file config.yaml
```

Please see the section on how to create a custom config file for your cluster below. 

## Creating a config file for your cluster

This pipeline uses at minimum two configuration files. The file [conf/base.config](../conf/base.config) contains information about the resource requirements 
of the individual stages of the pipeline. Normally, you do not have to touch this.

In addition, you will need a config file that specifies which resource manager your cluster uses. An example for a Slurm cluster which uses 
Singularity (see below) is included as [conf/ccga_medcluster.config](../conf/ccga_medcluster.config). Detailed instructions about resource managers and 
available options can be found [here](https://www.nextflow.io/docs/latest/executor.html).

Create your own config file and pass it to the pipeline during start-up using the `-c` flag

```bash
nextflow -c nextflow.config ikmb/esga -params-file config.yaml`
```

An example of a simple SLURM config file:`

```bash
#!/usr/bin/env nextflow

/*
 * -------------------------------------------------------
 *  Your local config file
 * -------------------------------------------------------
 */

executor {
  name="slurm"
  queueSize=150
}

process {

  executor = 'slurm'
  queue = 'all'
  memory = { 8.GB * task.attempt }
  cpus = { 1 * task.attempt }
  time = { 2.h * task.attempt }
  errorStrategy = { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries = 3
  maxErrors = '-1'

}

singularity {
        enabled = true
}

params {
  max_memory = 120.GB
  max_cpus = 16
  max_time = 120.h
}

```

## Installing all other software 

This pipeline uses a lot of different bioinformatics software - you can find a full list with versions in the included 
file [environment.yml](../environment.yml). You won't have to install any of these tools - assuming your cluster offers singularity or docker support:

### Singularity

Make sure you have [Singularity](https://github.com/sylabs/singularity). If Singularity is not available on 
your cluster, please ask your admins to install it - you will likely need version 3.0 or later. 

To enable use of singularity, simply add the following to your custom config file (see below):

```bash
singularity {
	enabled = true
}
```

Depending on your cluster and configuration of singularity, you may also have to provide some additional run options. 
A typical example would be that your data is stored on a network-mounted drive, which is not automatically detected by singularity. In this case, you can do:

```bash
singularity {
	enabled=true
	runOptions="-B /path/to/network/drive"
	cacheDir="/path/to/cache"
}
```

Note, that the cacheDir option will make sure that the container is only downloaded once and can be re-used for future pipeline runs. Otherwise, nextflow will re-download it every time you start an annotation project. 
