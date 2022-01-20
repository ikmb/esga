# How to run the pipeline 

The typical command for running the pipeline is as follows:

```
nextflow run ikmb/esga --genome 'genome.fasta' --proteins 'proteins.fasta' --reads 'data/*_R{1,2}.fastq' --transcripts 'transcripts.fa' --rm_lib repeats.fa
```

This will run all the steps in the pipeline (Proteins, ESTs/Transcriptome, RNAseq). The types of evidences you provide determine which parts of the pipeline are actually run. 


### Parameters file

In the next section, you will find a list of all user-configurable pipeline options. 
You can of course provide each option as a command line parameter. But this can get a bit tedious. As an alterantive, you can provide a configuration file using the YAML format. An example is included under [../assets/config.yaml](../assist/config.yaml). 
To provide a config file as an option, use `-params-file my_config.yaml`. The revised command will then read:

`nextflow -c nextflow.config run ikmb/esga -params-file config.yaml`

The default (reduced) [YAML options file](../assets/config.yaml):

```yaml
genome: ""
proteins: false
proteins_targeted: false
transcripts: false
references: false
rm_lib: false
rm_species: false
reads: false
aug_species: "human"
ncrna: false
protein_aligner: spaln
rnaseq_aligner: star
aug_options: "--singlestrand=true --alternatives-from-evidence=on --minexonintronprob=0.08 --minmeanexonintronprob=0.4 --maxtracks=3"
spaln_taxon: false
spaln_q: 5
utr: false
trinity: false
pasa: false
evm: false
npart_size: 200000000
max_intron_size: 50000
min_contig_size: 5000
singleEnd: false
rnaseq_stranded: false
```
  
An explanation of these options follows below.

### 1. Mandatory arguments 

#### `--genome`
Location of the genome you want to annotate. This file should be in FASTA format. Additionally, we recommend to make sure that you clean the fasta headers in a way that they do not contain any special characters, unnecessary spaces or other "meta" data. 
Please also be aware that some public databases do not allow the submission of assemblies that have leading or trailing 'N's in any of its scaffolds.

### 2. Evidences. At least one of:

#### `--reads` 
Location of your input FastQ files. For a set of PE libraries named using typical Illumina naming convention:

```bash
--reads 'path/to/data/*_R{1,2}_001.fastq.gz'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

#### `--transcripts` 
Location of a single FASTA file with all EST sequences or assembled transcriptome(s) (short reads or IsoSeq) from the species of interest. If you have multiple files, concatenate 
them into a single file first and make sure that the sequence names are not duplicated (this happens when you try to merge two Trinity assemblies, for example). 

#### `--proteins_targeted`
Location of a single FASTA file with exactly one proteome from your species of interest or a very closely related one. Make sure that this file contains more than just a few sequences as it expects to have a near complete proteome here. 

#### `--proteins` 
Location of a single FASTA file with protein sequences from related taxa. If you have multiple files, concatenate them into a single file first. 
The [included](../assets/Eumetazoa_UniProt_reviewed_evidence.fa) set of curated eumetazoan proteins set would be appropriate for this. 

#### `--references` 
This option points to one or more genome sequences in FASTA format, accompanied by a gene annotation in GTF format. The genome and the annotation need to share a base name for this to work.

For example, if you point to the human genome in FASTA format, named homo_sapiens.fa , this option expects there to be an annotation file called homo_sapiens.gtf right next to it. 

```bash

nextflow run ikmb/esga --references /path/to/ref_genome/my_genome.fa [...]

```

```bash

nextflow run ikmb/esga --refefences '/path/to/ref_genomes/*.fa' [...]

```

Each genome sequence will be aligned to the target assembly using [Satsuma2](https://github.com/bioinfologics/satsuma2) to produce a pairwise chain file. This chain file is then used by [Kraken](https://github.com/GrabherrGroup/kraken) to lift the original 
annotation onto the target genome. The resulting mapped models will not be corrected for splice junctions, so they are not fully valid annotations. However, the mapping of CDS and exon features can be used to inform the subsequent gene finding process.

Please make sure to use closely related genomes for this (to annotate a primate, any other primate or even mammal, would be fine). Also note that this process can consume a large amount of memory depending on the genome size(s) 
(>>250GB for vertebrates). We have had this part crash on nodes with 250GB Ram when aligning a chunk of a mammalian genome to the whole genome of another mammal. So you are likely going to need >= 250GB Ram nodes in your cluster 
to successfully use this part of ESGA when working on larger vertebrates. 

Suitable sources for this step are, for example, genomes and annotations from the [EnsEMBL](https://www.ensembl.org/index.html) database.

### 3. Programs to run 
By default, the pipeline will run all parts for which the required types of input are provided. However, some parts need to specifically "switched on" as they require longer run times and may not be strictly necessary. For example,
yes you can run the Trinity transcriptome assembly part of the pipeline (see below), but if you already have a set of assembled transcripts (--transcripts), this may not be necessary in combination with the RNA-seq hints
that are always generated when RNA-seq data is available. Likewise, you can select to only run AUGUSTUS-based gene predictions and skip the evidence modeler stage, and so on. 

The only non-optional part is the AUGUSTUS stage as this is the core around which ESGA was originally built. 

#### `--pasa` [ true | false (default) ]
Run the PASA pipeline to build gene models from aligned transcripts (requires --transcripts and/or --reads & --trinity). Please note that due to our need to use container technologies, PASA runs off a SQLite database. Consequently, this step is quite slow for larger genome/transcriptome data sets (several days!). 

#### `--evm` [ true | false (default) ]
Run the evidence-modeler gene build pipeline, combining all the various outputs produced by this workflow. 

#### `--trinity` [ true | false (default) ]
Run the de-novo transcriptome assembler Trinity to produce transcript information. Requires --reads. Will be switched on by default if --pasa is requested and --reads are available. Note that Trinity may consume large amounts of RAM, depending
on the size of the input data. However, we suspect that 120GB RAM will be fine in most instances - a typical minimum spec for a modern cluster node. 

#### `--ncrna`[ true | false (default) ]
Independently predict non-coding RNAs using RFam version 14. The resulting models will not be merged into the main gene build but can be used for manual curation in e.g. WebApollo. Please note that the node executing the nextflow process
must have access to the internet to download the RFam files on-the-fly (some clusters do not allow internet connections from compute or even login nodes).

### 4. Tool options

#### `--rnaseq_aligner` [ hisat | star (default) ]
Alignment software to use for RNAseq reads. Default is STAR.

#### `--protein_aligner` [ gth | spaln (default) ]
Alignment software to use for protein sequences. Default is `spaln` since it is much faster than `gth`. This speed advantage may be offset by some questionable alignments. In that case, re-run with gth to check if the alignments improve. 

### 5. Program parameters

#### `--spaln_q <int>` [ 5 (default) ]
Algorithm to be used for SPALN alignment. See Spaln [documentation](https://github.com/ogotoh/spaln#Exec) for details. Options 5 to 7 are typically what you want when aligning proteins against an entire genome.

#### `--spaln_taxon` [ default = "Tetrapod" ]
Name of the taxonomic group to choose the internal SPALN parameters. See column 2 of [this](https://github.com/ogotoh/spaln/blob/master/table/gnm2tab) list. The default, 'Tetrapod', should work for all tetrapods.

#### `--max_intron_size <int>` [ 20000 (default) ]
The default value is set to 20000 - for something like a nematode, this would be too long; for some lower vertebrates it would probably be fine, although 
a few introns may be much longer. Genes containing such extraordinarily large introns will then probably be mis-annotated. Information on plausible intron sizes can be obtained from the literature for many taxonomic groups. 

#### `--rm_species`
Use this taxonomic group or species to identify and mask repeats. Valid names can, in most cases, be guessed, and follow the nomenclature provided through the [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy). This option draws from available data included in DFam 3.2, which contains HMM profiles for over 273.0000 repeat families from [347 species](https://www.dfam.org/releases/Dfam_3.2/relnotes.txt).

As can be seen, this dataset is not necessarily comprehensive in its taxonomic coverage. If your species of interest falls well outside of the many organisms included with DFam, consider using a FASTA file of repeats instead (--rm_lib). If no such file exists, leave the options `--rm_lib` and `--rm_species` empty to trigger the de-novo repeat-modelling stage in ESGA using RepeatModeler. 

On de-novo repeat-modelling:
Detecting repeats de-novo  requires a sufficient number of related repeats to be present in your assembly. If your genome was assembled from short reads, this strategy may not return anything -
short read assemblies tend to collapse repeats. In this case, the pipeline will fall back to the built-in DFam library that ships with RepeatMasker.

#### `--rm_lib`[ fasta file | false ]
ESGA can run RepeatMasker using repeat sequences in FASTA format. Possible sources include self-computed repeats (using RepeatModeler) or curated repeat libraries from 
GRINST (www.grinst.org, commercial). 

If you have a copy of the complete Repeatmasker library (and an installation of RM), you can extract the repeat annotation from a species like so: 

```
# Get the Tree of all available species: 
perl /[...]/RepeatMasker/4.0.8/util/queryRepeatDatabase.pl -tree

# Select the name of the species (e.g. "Ostreoida") from the output tree and do:
perl /[...]/RepeatMasker/4.0.8/util/queryRepeatDatabase.pl -species Ostreoida > RMdb_Ostreoida.fa
``` 
You will need to delete the first line of the resulting FASTA file since RepeatMasker prints some basic information into it for no clear reason - but this otherwise breaks any FASTA parser. 

Then run the pipeline with the option "--rm_lib RMdb_Ostreoida.fa". 

#### `--aug_species` [ default = 'human' ]
Species model for AUGUSTUS. A list of valid identifiers can be found [here](https://github.com/Gaius-Augustus/Augustus/blob/master/docs/RUNNING-AUGUSTUS.md).

#### `--aug_config` [ default = 'assets/augustus/augustus_default.cfg' ]
Location of AUGUSTUS extrinsic hint configuration file. By default, this pipeline uses a [config file](../assets/augustus/augustus_default.cfg) that we found to work well for predicting gene models in mammalian genomes using the kinds of extrinsic hints constructed by this pipeline. However, you can pass your own option. However, you should start from our template to make sure that the types of hints match the ones seen by Augustus. 

#### `--aug_options` [ default = "" ]
AUGUSTUS has numerous options, not all of which are exposed through our pipeline. If you have good reason to use a specific command line flag that is not configurable through ESGA, you can use this option to set it manually. 
For example, to allow Augustus to predict overlapping genes (default: no), you could specifiy `--aug_options '--singleStrand=true'`

Our template contains some of the non-standard options that we found to help improve the results produced by AUGUSTUS. 

#### `--aug_training` [ true | false (default) ]

Enable training of an AUGUSTUS profile. This option requires either a species-specific proteome (--proteins_targeted) or a transcriptome (--pasa). ESGA will prefer a targeted proteome - however, this proteome should contain well over 1000 proteins, else the training will fail. 

This option also interacts with `--aug_species` . If the species name already exists, ESGA will perform a re-training of that model. If the species name does not exist, a new profile will be created. 

Please not that this routine is somewhat basic and reliant on the quality of your input data. Resulting models should work ok, but are likely not going to be as good as those shipping with AUGUSTUS. Especially sparse RNAseq data may cause issues when using the PASA route for producing training models. 

#### `--utr` [ true | false (default ]
Enabling prediction of UTRs during AUGUSTUS ab-initio gene finding can help produce more acurate gene models. However, this option is best used with available RNA-seq data and should only ever be switched on if the AUGUSTUS profile
was trained to predict UTRs - else the pipeline will fatally fail (just remove --utr in that case and -resume).
    
#### `--evm_weights` [ default = 'assets/evm/evm_weights.txt' ]
A file specifying the weights given to individual inputs when running EvidenceModeler. By default a [pre-configured file](../assets/evm/evm_weights.txt) is used. However, you can pass your own version - ideally starting from our template to make sure that all types of evidence are considered. 

#### `--pri_prot_targeted <int>` [ 5 (default ]
Priority of protein-based hints for Augustus gene predictions from the closest reference proteome. Higher priority hints are considered first and override lower-priority hints.

#### `--pri_prot <int>` [ 4 (default ]
Priority of protein-based hints for Augustus gene predictions. Higher priority hints are considered first and override lower-priority hints. 

#### `--pri_est <int>` [ 3 (default) ]
Priority of transcript-based hints for Augustus gene predictions. Higher priority hints are considered first and override lower-priority hints.

#### `--pri_rnaseq <int>` [ 4 (default) ]
Priority of RNAseq-based hints for Augustus gene predictions. Higher priority hints are considered first and override lower-priority hints.

#### `--pri_wiggle <int>` [ 2 (default) ]
Read coverage from RNA-seq experiments may be used to help AUGUSTUS in particular predict acurate UTRs and/or isoforms.  Higher priority hints are considered first and override lower-priority hints.

#### `--pri_trans <int>` [ 4 (default) ]
Priority for trans-mapped annotations (i.e. lift-over gene models, converted to CDS and exon hints).

### 6. How to tune the speed of the pipeline - data splitting

TL;DR: This is probably not necessary, the default options should work fine. 

One of the advantages of using Nextflow is that it allows you to speed up a pipeline by splitting some of the input files into smaller chunks before 
running specific programs. Then that program can be run on each smaller chunk in parallel in a compute cluster. 
When all instances of the program are finished, Nextflow can correctly put together all the results in a single output for that program. Depending on the size and contiguity of your target genome and the size of the evidence data, you may want 
to tweak one or several of the parameters below. If unsure, leave at the defaults.

#### `--aug_chunk_length`[default = 3000000]
The length of genomic sequence to run individual AUGUSTUS jobs on. Chunks are overlapped by 1/6*aug_chunk_length. The aug_chunk_length should be chosen so that the longest expected gene will fit in it. The default of 3MB should be sufficient for most metazoans. 
Decreasing this number will (likely) create more parallel jobs and may miss genes entirely if they are longer that the specified length. 

#### `--nproteins` [ default = 200 ]
Number of sequences in each protein alignment job. Larger values will usually create longer run times, but decrease the number of parallel jobs and load on the file system. 

#### `--npart_size` [ default = 200000000 ]
Size in bp of the pieces into which the genome is split for parallelization optimization. By default, this is set to `200000000`, i.e. 200Mb. This function will *not* break scaffolds, but simply tries to distribute the
assembly into chunks of this size - some chunks may be bigger, some smaller. This function uses [fasta-splitter.pl](http://kirill-kryukov.com/study/tools/fasta-splitter/).

Setting this to larger values will create fewer parallel jobs, so the run time is likely going to increase. 

#### `--min_contig_size` [ default = 5000 ]
Small contigs generally will not contribute anything useful to the annotation, but can increase runtime dramatically. Contigs smaller than this size are removed from the assembly prior to annotation. 

### 7. Other options 

#### `--email` [ you@somewhere.com | false (default)]
If you specify an Email address, the pipeline will send a notification upon completion. However, for this to work, the node running the nextflow process must have a configured Email server. 

#### `--singleEnd` [ true | false (default) ]
By default, the pipeline expects paired-end RNA-seq data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq.gz'
```

It is not possible to run a mixture of single-end and paired-end files in one run. Please see the [troubleshooting](troubleshooting.md) document for more information on this. 

#### `--rnaseq_stranded` [ true | false (default) ]
Whether your RNAseq library was sequenced with a fw strand-specific protocol. Our assumption is that this would be dUTP as is typical for Illumina applications. 

#### `--outdir` [ default = 'results' ]
The output directory where the results will be saved. 

#### `--run_name` [default = a random name ]
Give this pipeline run a specific name. Otherwise, nextflow will generate a random one. 

### 7. Nextflow parameters (indicate with single dash "-")

#### `-profile`
Use this parameter to choose a configuration profile. Each profile is designed for a different combination of compute environment and installation estrategy (see [Installation instructions](../docs/installation.md)).

#### `-params-file config.yaml`
All the above options can be passed from either the command line or through a configuration file. A suitable template is included under assets/config.yaml.

#### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

#### `-c`
Specify the path to a specific config file (this is a core NextFlow command).

### Job Resources
#### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

#### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples. 

