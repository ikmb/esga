## What happens in this pipeline?

This pipeline accepts a genome assembly and a range of sequence inputs from which so-called gene builds are computed - often also referred to as 'annotations'.

These correspond to the structure of different (protein-coding) genes in your assembly and can inform a wide range of scientific questions - from gene evolution to disease studies. 

Please note that any computational gene build is prone to a wide range of errors, stemming from noise in the input data as well as limitations in the prediction algorithms. The more obscure the organism, the harder it will
likely be to achieve a highly sensitive annotation in just "one go". 

In any case, the product of this (or any other) annotation pipeline is a best-guess interpretation and will require manual curation to achieve the best possible result suitable for highly detailed genetic/genomic studies.

The steps in this pipeline can be broken down based on the types of input they use. 

* Proteins - are aligned against the genome assembly using SPALN. The pipeline distinguished between "targeted" proteins (stemming from your organism) and "generic" proteins (from other species)
* Transcripts - are aligned against the genome assembly using Minimap2. These can stem from a range of technologies, including IsoSeq, ESTs and assembled RNA-seq data
* Splice junctions - are extracted from aligned RNA-seq data
* Repeats - are annotated using RepeatMasker. If no repeats are available, they will be modelled "de-novo" from the assembly (only really works for assemblies using long-read technologies...).
* Gene synteny - is computed from whole genome alignments and mapped gene models from closely related high-quality reference genomes.

The resulting data will be processed by annotation "engines" (AUGUSTUS, Pasa and EvidenceModeler) and produce distinct gene builds for each of these tools.

All of the data types and tools are optional, with the expception of AUGUSTUS. More information on this is available from the [usage](usage.md) instructions

![](../images/Pipeline_dag.svg) 
