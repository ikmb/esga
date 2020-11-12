# Installation and configuration 

## Compute infrastructure

This pipeline is designed to run on a distributed compute system, such as a traditional HPC cluster system. 
We have tested the pipeline on two Slurm clusters, with node configurations of 16cores/128GB Ram and 20cores/256GB Ram, respectively. 

While smaller nodes will probably work, it may require some tweaking on your end. Most importantly, if you plan on using the transcriptome 
assembly branch of the pipeline, available memory may become limiting (however, 128GB Ram should be fine for typical datasets; 256GB are perhaps 
necessary if you plan on using a larger sample size). 

## Installing Nextflow 

The pipeline is built using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool to run tasks across 
distributed compute systems in a very portable manner. Therefore the first thing to do is to install Nextflow. 

Nextflow runs on most systems (Linux, Mac OSX etc). It can be installed by running the following commands:

```
# Make sure that Java v8+ is installed:
java -version

```

# Install Nextflow

You will currently need Nextflow version 20.01 to run this pipeline as certain features were removed in later releases and have not been updated in our code.

[Nextflow 20.01](https://github.com/nextflow-io/nextflow/releases/download/v20.01.0/nextflow-20.01.0-all)

Just put the downloaded file into your $PATH and rename it to 'nextflow'

## Obtaining the code 

You can check out the code to a location on your system (i.e. $HOME/git/). This is recommended as you will have to make some minor changes 
so the pieline knows how to run on your system (see below). 

``` 
cd $HOME/git/ 

git clone git@github.com:ikmb-denbi/genome-annotation.git
``` 
 
Launching the pipeline then works as follows:

```bash
nextflow run $HOME/git/genome-annotation/main.nf <run options>
```

## Creating a config file for your cluster

This pipeline uses at minimum two configuration files. The file [conf/base.config](../conf/base.config) contains information about the resource requirements 
of the individual stages of the pipeline. Normally, you do not have to touch this.

In addition, you will need a config file that specifies which resource manager your cluster uses. An example for a Slurm cluster which uses 
Singularity (see below) is included as [conf/slurm_ikmba.config](../conf/slurm_ikmba.config). Detailed instructions about resource managers and 
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
your cluster, please ask your admins to install it. 

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
