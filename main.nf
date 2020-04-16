#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include fastaMergeFiles from "./modules/fasta" params(params)
include { repeatmasking_with_lib; repeatmasking_with_species } from "./modules/repeatmasker/main.nf" params(params)
include model_repeats from "./modules/repeatmodeler/main.nf" params(params)
include proteinhint from "./modules/proteins/main.nf" params(params)
include esthint from "./modules/transcripts/main.nf" params(params)
include esthint as trinity_esthint from "./modules/transcripts/main.nf" params(params)
include augustus_prediction from "./modules/augustus/main.nf" params(params)
include merge_hints from "./modules/util" params(params)
include rnaseqhint from "./modules/rnaseq/main.nf" params(params)
include trinity_guided_assembly from "./modules/trinity/main.nf" params(params)
include evm_prediction from "./modules/evm/main.nf" params(params)
include pasa_assembly from "./modules/pasa/main.nf" params(params)
include assembly_preprocessing from "./modules/assembly/main.nf" params(params)

def helpMessage() {
  log.info"""
  =================================================================
   IKMB - de.NBI | ESGA - Evidence-based scalable genome annotation | v${workflow.manifest.version}
  =================================================================
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run ikmb/esga --genome 'Genome.fasta' --proteins 'Proteins.fasta' --reads 'data/*_R{1,2}.fastq'
  Mandatory arguments:
  --genome		Genome reference
      
  At least one of:
  --proteins		Proteins from other species
  --transcripts		ESTs or transcriptome
  --reads		Path to RNA-seq data (must be surrounded with quotes)
  Options:
    -profile            Hardware config to use (optional, will default to 'standard')
    Programs to run:
    --trinity		Run transcriptome assembly with Trinity and produce hints from the transcripts [ true (default) | false ]
    --pasa 		Run the transcriptome-based gene builder PASA (also required when running --training). [ true | false (default) ]. Requires --ESTs and/or --reads with --trinity. 
 	
    Programs parameters:
    --rm_lib		Perform repeatmasking using a library in FASTA format [ default = 'false' ]
    --rm_species	Perform repeatmasking using the built-in DFam library for a given species/taxonomic group (overrides rm_lib) [ default = 'false' ]
    --aug_species	Species model for Augustus [ default = 'human' ]. If "--training true" and you want to do de novo training, give a NEW name to your species
    --aug_config	Location of augustus configuration file [ default = 'bin/augustus_default.cfg' ]
    --max_intron_size	Maximum length of introns to consider for spliced alignments [ default = 20000 ]
    --evm		Whether to run EvicenceModeler at the end to produce a consensus gene build [true | false (default) ]
    --evm_weights	Custom weights file for EvidenceModeler (overrides the internal default)
    --min_contig_size   Discard all contigs from the assembly smaller than this [ default = 5000 ]

    Evidence tuning
    --pri_prot		A positive number between 1 and 5 - the higher, the more important the hint is for gene calling (default: 5)
    --pri_est		A positive number between 1 and 5 - the higher, the more important the hint is for gene calling (default: 3)
    --pri_rnaseq	A positive number between 1 and 5 - the higher, the more important the hint is for gene calling (default: 4)
 	
    How to split programs:
    --nblast		Chunks (# of sequences) to divide genome for blastx jobs [ default = 100 ]
    --nexonerate	Chunks (# of blast hits) to divide Exonerate jobs [ default = 200 ]
    --npart_size	Size in bp to divide RepeatMasker and Augustus jobs [ default = 200000000 ]
    --chunk_size 	Size of sub-regions of the genome on which to run Blastx jobs [ default = 50000 ]

    Other options:
    --singleEnd		Specifies that the RNAseq input is single end reads [ true | false (default) ]
    --rnaseq_stranded	Whether the RNAseq reads were sequenced using a strand-specific method (dUTP) [ true | false (default) ]
    -name		Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}


// Show help message
if (params.help){
	helpMessage()
	exit 0
}

// ****************
// input validation
// ****************

genome = file(params.genome)
if (!genome.exists()) {
	exit 1, "Could not find a genome assembly or file does not exist (---genome)"
}
if (params.proteins) {
	p = file(params.proteins)
	if (!p.exists()) {
		exit 1, "The specified protein file does not exist!"
	}
	proteins = Channel.fromPath(p)
} else {
	proteins = Channel.empty()
}

if (params.transcripts) {
	t = file(params.transcripts)
	if (!t.exists() ) {
		exit 1, "The specified transcript file does not exist!"
	}
	transcripts = Channel.fromPath(t)
} else {
	transcripts = Channel.empty()
}

if (params.rm_species && params.rm_lib) {
	log.warn "Specified both a repeatmasker species and a library - will only use the species!"
}
if (!params.rm_species && !params.rm_lib) {
	log.warn "No repeat data provided, will model repeats de-novo (slow!!!)"
}
if (params.pasa && !params.transcripts && !params.trinity) {
	exit 1, "Cannot run PASA without transcript data (--transcripts and/or --trinity)"
}
// Provide the path to the augustus config folder
// If it's in a container, use the hard-coded path, otherwise the augustus env variable
if (!workflow.containerEngine) {
	Channel.from(file(System.getenv('AUGUSTUS_CONFIG_PATH')))
		.ifEmpty { exit 1; "Looks like the Augustus config path is not set? This shouldn't happen!" }
        	.set { augustus_config_folder }
} else {
// this is a bit dangerous, need to make sure this is updated when we bump to the next release version
	Channel.from(file("/opt/conda/envs/esga-1.0/config"))
        	.set { augustus_config_folder }
}
if (!params.aug_species) {
	exit 1, "Must provide an AUGUSTUS species name to run this pipeline (--aug_species)"
}
if (params.pasa && params.npart_size < 100000000) {
	log.warn "An npart_size smaller than 100MB is risking PASA to fail - increase in case of a crash..."
}

// *******************************
// Print run information
// *******************************
log.info "========================================="
log.info "ESGA Genome Annotation Pipeline v${workflow.manifest.version}"
log.info "d88888b .d8888.  d888b   .d8b.  "
log.info "88'     88'  YP 88' Y8b d8' `8b "
log.info "88ooooo `8bo.   88      88ooo88 "
log.info "88~~~~~   `Y8b. 88  ooo 88~~~88 "
log.info "88.     db   8D 88. ~8~ 88   88 "
log.info "Y88888P `8888Y'  Y888P  YP   YP "
log.info "========================================="
log.info "Genome assembly: 		${params.genome}"
if (params.rm_lib) {
	log.info "Repeatmasker lib:		${params.rm_lib}"
} else if (params.rm_species) {
	log.info "Repeatmasker DFam 2.0 species: ${params.rm_species}"
} else {
	log.info "Repeatmasking:			Compute de-novo"
}
log.info "-----------------------------------------"
log.info "Program settings:"
log.info "AUGUSTUS species:		${params.aug_species}"
log.info "AUGUSTUS config:		${params.aug_config}"
log.info "Maximum intron size:		${params.max_intron_size}"
log.info "Run PASA assembly:		${params.pasa}"
log.info "Run Trinity assembly:		${params.trinity}"
log.info "Run EVM gene building:		${params.evm}"
log.info "EVM weights:			${params.evm_weights}"
log.info "-----------------------------------------"
log.info "Evidences:"
log.info "Proteins:			${params.proteins}"
log.info "Transcripts:			${params.transcripts}"
log.info "RNA-seq:			${params.reads}"
log.info "-----------------------------------------"
log.info "Parallelization settings"
log.info "Size of blastx jobs:			${params.chunk_size} bp"
log.info "# Sequences per Blast job:		${params.nblast}"
log.info "# Sequences per Exonerate job:		${params.nexonerate}"
log.info "Size of genome-level jobs:		${params.npart_size} bp"
log.info "-----------------------------------------"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Command Line:			$workflow.commandLine"
log.info "Run name: 			${params.run_name}"
log.info "========================================="


// ***********************************
// WORKFLOW STARTS HERE
// ***********************************
workflow {

	main:

	// Pre-process the assembly
	// Generate assembly stats and remove small contigs
	assembly_preprocessing(params.genome)
	genome_clean = assembly_preprocessing.out.fasta
	assembly_stats = assembly_preprocessing.out.stats

	// Repeat-mask the assembly
	if (params.rm_species) {
		repeatmasking_with_species(genome_clean,params.rm_species)
		genome_rm = repeatmasking_with_species.out.genome_rm
		repeats = Channel.empty()
	// Use a library provided by the user or compute de-novo
	} else {
		if (params.rm_lib) {
			repeats = Channel.fromPath(params.rm_lib)
		} else {
			model_repeats(params.genome)
			repeats = model_repeats.out.repeats
		}
	        repeatmasking_with_lib(params.genome,repeats)
		genome_rm = repeatmasking_with_lib.out.genome_rm
	}
	
	// Generate hints from proteins (if any)
	if (params.proteins) {
		proteinhint(genome_rm,proteins)
		protein_hints = proteinhint.out.hints
		protein_gff = proteinhint.out.gff
	} else {
		protein_hints = Channel.empty()
		protein_gff = Channel.empty()
	}

	// Generate hints from transcripts (if any)
	if (params.transcripts) {
		esthint(genome_rm,transcripts)
		est_hints = esthint.out.hints
		est_gff = esthint.out.gff
	} else {
		est_hints = Channel.empty()
		est_gff = Channel.empty()
	}

	// Generate hints from RNA-seq (if any)
	if (params.reads) {
		rnaseqhint(genome_clean,params.reads)
		rna_hints = rnaseqhint.out.hints
		rna_bam = rnaseqhint.out.bam
		if (params.trinity) {
			trinity_guided_assembly(rnaseqhint.out.bam)
			trinity_esthint(genome_rm,trinity_guided_assembly.out.assembly)
			trinity_gff = trinity_esthint.out.gff
			trinity_hints = trinity_esthint.out.hints
			trinity_assembly = trinity_guided_assembly.out.assembly
		} else {
			trinity_hints = Channel.empty()
			trinity_gff = Channel.empty()
			trinity_assembly = Channel.empty()
		}	
		
	} else {
		rna_hints = Channel.empty()
		rna_bam = Channel.empty()
		trinity_hints = Channel.empty()
		trinity_gff = Channel.empty()
		trinity_assembly = Channel.empty()
	}

	// Build evidence-based gene models from transcripts
	if (params.pasa) {
		transcript_files = trinity_assembly.concat(transcripts)
		fastaMergeFiles(transcript_files.collect())
		pasa_assembly(genome_rm,fastaMergeFiles.out[0])
		pasa_gff = pasa_assembly.out.gff
	} else {
		pasa_gff = Channel.empty()
	}

	// Merge hints
	hints = protein_hints.concat(est_hints, trinity_hints,rna_hints)
	merge_hints(hints.collect())
	
	// Run AUGUSTUS
	augustus_prediction(genome_rm,merge_hints.out,augustus_config_folder)
	augustus_gff = augustus_prediction.out.gff

	// Combine all inputs into consensus annotation
	if (params.evm) {
		gene_gffs = augustus_gff.concat(pasa_gff).collect()
		transcript_gff = est_gff.concat(trinity_gff).collectFile()
		evm_prediction(genome_rm,protein_gff,transcript_gff,gene_gffs)
		evm_gff = evm_prediction.out.gff
	} else {
		evm_gff = Channel.empty()
	}

	publish:
		genome_rm to: "${params.outdir}/repeatmasking", mode: 'copy'
		assembly_stats to: "${params.outdir}/assembly", mode: 'copy'
		augustus_prediction.out.gff to: "${params.outdir}/annotation/augustus", mode: 'copy'
		est_hints to: "${params.outdir}/evidence/hints", mode: 'copy'
		est_gff to:  "${params.outdir}/evidence/transcripts", mode: 'copy'
		protein_hints to: "${params.outdir}/evidence/hints", mode: 'copy'
		protein_gff to: "${params.outdir}/evidence/proteins", mode: 'copy'
		rna_hints to: "${params.outdir}/evidence/hints", mode: 'copy'
		rna_bam to: "${params.outdir}/evidence/rnaseq", mode: 'copy'
		repeats to: "${params.outdir}/repeatmasking", mode: 'copy'
		evm_gff to: "${params.outdir}/annotation/evm", mode: 'copy'
		pasa_gff to: "${params.outdir}/annotation/pasa", mode: 'copy'
		
}



