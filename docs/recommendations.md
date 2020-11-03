# Recommendations

Evidence data for annotation can be obtained from a number of sources. Below follow some general suggestions that we found to work well.

## Targeted proteins

Sources for this could be an existing annotation from your species of interest or a very closely related one. During this stage the proteins will essentially be laminated onto the assembly to define gene loci. So if this set is evolutionarily too distant, 
this step might produce sub-par results. 

You could either check public databases like Uniprot for curated proteins of your species, or turn to another annotation project such as [EnsEMBL](ftp://ftp.ensembl.org/pub/current_fasta/). However, please keep in mind
that annotating a genome guided by another annotation risks perpetuation annotation errors. That said, from our experience EnsEMBL gene builds are quite robust and scientifically usable.

## Proteins (other)

Other proteins for annotation should cover a wider taxonomic range to help the pipeline find genes that were perhaps missing from the targeted stage. Since ESGA is aimed for vertebrates specifically, and metazoans at most, we have included a set of reviewed
metazoan proteins that could be used here. Otherwise, you have some additional options:

[Uniprot](https://www.uniprot.org/) is generally a good place to start. The search interface allows you to limit results by taxnomic group and level of experimental support for 
each sequence. As a rule of thumb, try to focus on sequences that are flagged as "full length" and that have experimental support from either
protein or transcriptome sequencing. 

[RefSeq](https://www.ncbi.nlm.nih.gov/protein/) is another useful database. However, sequences here tend to include a lot more computational predictions,
which runs the risk of perpetuating annotation errors. 

## Transcripts

Transcripts provide information of those genomic region that, under a given condition, are actively transcribed. Depending on the library/sequencing strategy, these sequences will preferentially stem from protein-coding loci or instead represent the entire
transcriptional landscape. For the purpose of annotation, poly-A (coding) transcrips are to be preferred. Typical sources of this data are either (traditional) EST libraries or, more commonly, RNA-sequencing (see below). 

EST data can be downloaded from [GenBank](https://www.ncbi.nlm.nih.gov/nucleotide) or [ENA](https://www.ebi.ac.uk/ena). 
Make sure that the sequences are from your species of interest; on the nucleotide level even relatively small evolutionary timescales could
greatly dimminish the alignment quality/rate. 

An alternative to traditional ESTs are sequence data from de-novo assembled transcriptomes. For shrot reads, we recommend genome-guided [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Genome-Guided-Trinity-Transcriptome-Assembly). 

For long reads, you might want to check out IsoSeq (PacBio).

## RNA-seq data

This is perhaps the most important type of annotation evidence, as it stems specifically from your organism of interest and should provide
a good resolution of splice junctions. 

We recommend poly-A selected, strand-specific paired-end RNA-seq libraries for the purpose of annotation (any length is fine, but longer reads tend to work best as they will map more uniquely). Ideally, these libraries cover different develpomental stages and body parts.
For mammals, a typical sequencing depth per library is around 25 million PE reads (this is more related to the size of the transcriptome than 
of the genome!)

While more is often better, the pipeline is designed to process "hundreds or millions" of PE reads at most. So aim for some representative tissues, one library each (no replicates). Especially the generation of exon hints from aligned
reads takes a very long time if too much data is provided, typically without much noticable gain to the annotation quality. 

## Repeats

By default, Repeatmasker will run with the built-in DFam hmm profile for (mostly) primates. It is thus generally advisable to instead provide 
repeat annotations in FASTA format. Possible sources include self-computed repeats (using RepeatModeler - which ESGA will run for you if no other repeat informatino is provided) or curated repeat libraries from 
GRINST (www.grinst.org, commercial). 

If you have a copy of the complete Repeatmasker library (and an installation of RM), you can extract the repeat annotation from a species like so: 

```
# Get the Tree of all available species: 
perl /[...]/RepeatMasker/4.0.8/util/queryRepeatDatabase.pl -tree

# Select the name of the species (e.g. "Ostreoida") from the output tree and do:
perl /[...]/RepeatMasker/4.0.8/util/queryRepeatDatabase.pl -species Ostreoida > RMdb_Ostreoida.fa
``` 

However, if you choose to let RepeatModeler identify repeat sequences in your assembly, all you need to do is...nothing. If neither rm_lib nor rm_species are defined, repeat modelling will be performed automatically. 
This will produce a fasta file of candidate repeats that are fed to RepeatMasker. Be advised that the success of this strategy depends on the method of genome sequencing and assembly. Especially short-read based assemblies are 
prone to collapsing complex repeat structures so that these will be invisible to RepeatModeler and potentially adversely affect the annotation quality down the line. 

## Training a new Augustus model

Augustus ships with a number of high-quality prediction models from a range of taxonomic groups. Usually, the easiest approach is thus to
use a model that is taxonomically somewhat close to the species you are trying to annotate. For example, if annotating a bird, the built-in model for
chicken should work just fine. Remember, Augustus uses a range of hints produced by this pipeline to inform its gene finding anyway. As long as the basics of
the model are appropriate for you organism, this is a good approach.

However, if you are not getting satisfying results or find that no pre-existing model is likely appropriate for your genome of choice, you can enable an
automatic training routine. This will use available transcriptome data to either refine an existing model (--model exists) or built one from scratch 
(--model does not yet exist). This is somewhat experimental and depends a lot on the quality of the input data. It also takes a pretty long time (several
days for larger genomes) as it needs to first re-construct gene models from the aligned transcriptome data, select all the models that are probably
full length and finally use these gene structures to train Augustus. 

## Improving an initial gene build

Chances are the first run of the pipeline will not produce satisfactory results. Apart from trying to add additional data (e.g. more RNAseq, or ESTs, etc), an obvious
aspect is the tuning of parameters for AUGUSTUS and EvidenceModeler. Both tools require files that control their behavior. 

Changes to these files will force parts of the pipeline to rerun that rely on stages that make use of these files. EVM and AUGUSTUS are at the very end of the ESGA workflow, so resumeing the pipeline with
modified config files typically does not take too long. You can thus try to improve your annotation iteratively by just playing with these files. More below. 

### Augustus
AUGUSTUS uses extrinsic hints to guide the gene finding. These hints are produced by the ESGA pipeline from provided evidence. How these hints are weighted inside AUGUSTUS however are 
controlled by a config file. 

Please see [here](https://github.com/Gaius-Augustus/Augustus/blob/master/config/extrinsic/extrinsic.cfg) for instructions on how to tweak these parameters for optiomal performance. ESGA uses, by default,
an extrinsic config file that we have set up to work for our typical projects. You can however pass a modified version from the command line instead using the `--aug_config` option together with `-resume`. 

## EVM
ESGA produces several types of inputs (all of them optional) that are then combined by EVM into a consensus gene buid. The weight given to each type of input is controlled by the [weights.txt](../assets/evm/weights.txt) file. 
EVM uses weights between 1 and 10 to determine which evidence to consider first. To pass your own, modified version of the weights file, use the `--evm_weights my_weights.txt` syntax together with the `-resume` option. 
