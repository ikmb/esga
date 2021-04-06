#!/usr/bin/env nextflow

nextflow.preview.dsl=2
//nextflow.enable.dsl=2

params.version = workflow.manifest.version

// this needs to passed to the imported modules to determine if augustus is run with or without UTR annotation
include fastaMergeFiles from "./modules/fasta" params(params)
include { repeatmasking_with_lib; repeatmasking_with_species } from "./modules/repeatmasker/main.nf" params(params)
include model_repeats from "./modules/repeatmodeler/main.nf" params(params)
include { proteinmodels; proteinhint_spaln } from "./modules/proteins/main.nf" params(params)
include esthint from "./modules/transcripts/main.nf" params(params)
include esthint as trinity_esthint from "./modules/transcripts/main.nf" params(params)
include { augustus_prediction; augustus_prediction_slow; augustus_train_from_spaln; augustus_train_from_pasa  } from "./modules/augustus/main.nf" params(params)
include merge_hints from "./modules/util" params(params)
include rnaseqhint from "./modules/rnaseq/main.nf" params(params)
include trinity_guided_assembly from "./modules/trinity/main.nf" params(params)
include evm_prediction from "./modules/evm/main.nf" params(params)
include { polish_annotation; pasa } from "./modules/pasa/main.nf" params(params)
include assembly_preprocessing from "./modules/assembly/main.nf" params(params)
include rfamsearch from "./modules/infernal/main.nf" params(params)
include map_annotation from "./modules/satsuma/main.nf" params(params)

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
  --proteins_targeted	Exactly one set of proteins from this or a very closely related species
  --transcripts		ESTs or transcriptome
  --reads		Path to RNA-seq data (must be surrounded with quotes)

  Options:
    -profile            Hardware config to use (optional, will default to 'standard')
    Programs to run:
    --pasa 		Run the transcriptome-based gene builder PASA (also required when running --training). [ true | false (default) ]. Requires --ESTs and/or --reads. 
    --evm               Whether to run EvicenceModeler at the end to produce a consensus gene build [true | false (default) ]
    --ncrna		Annotate ncRNAs using RFam
    --trinity		Build transcript models from provided RNAseq data
 	
    Programs parameters:
    --rm_lib		Perform repeatmasking using a library in FASTA format [ default = 'false' ]
    --rm_species	Perform repeatmasking using the built-in DFam library for a given species/taxonomic group (overrides rm_lib) [ default = 'false' ]
    --aug_species	Species model for Augustus [ default = 'human' ]. If "--training true" and you want to do de novo training, give a NEW name to your species
    --aug_config	Location of augustus configuration file [ default = 'bin/augustus_default.cfg' ]
    --max_intron_size	Maximum length of introns to consider for spliced alignments [ default = 20000 ]
    --evm_weights	Custom weights file for EvidenceModeler (overrides the internal default)
    --min_contig_size   Discard all contigs from the assembly smaller than this [ default = 5000 ]
    --min_prot_length 	Length of proteins in AA to consider when building hints [default = 30]

    Evidence tuning
    --pri_prot		A positive number between 1 and 5 - the higher, the more important the hint is for gene calling (default: 4)
    --pri_prot_targeted	A positive number between 1 and 5 - the higher, the more important the hint is for gene calling (default: 5)
    --pri_est		A positive number between 1 and 5 - the higher, the more important the hint is for gene calling (default: 3)
    --pri_rnaseq	A positive number between 1 and 5 - the higher, the more important the hint is for gene calling (default: 4)
    --pri_wiggle	A positive number between 1 and 5 - the higher, the more important the hint is for gene calling (default: 2)
 	
    How to split programs:
    --nproteins		# of sequences to divide protein alignment jobs into [ default = 100 ]
    --npart_size	Size in bp to divide RepeatMasker and Augustus jobs [ default = 200000000, i.e. 200MB ]

    Other options:
    --singleEnd		Specifies that the RNAseq input is single end reads [ true | false (default) ]
    --rnaseq_stranded	Whether the RNAseq reads were sequenced using a strand-specific method (dUTP) [ true | false (default) ]
    --run_name		Name for this pipeline run
    """.stripIndent()
}


// Show help message
if (params.help){
	helpMessage()
	exit 0
}

// ****************
// Collect metadata
// ****************

def summary = [:]

summary['Assembly'] = params.genome
if (params.proteins) {
	summary['Proteins'] = params.proteins
}
if (params.proteins_targeted) {
	summary['ProteinsTargeted'] = params.proteins_targeted
}
if (params.transcripts) {
	summary['Transcripts'] = params.transcripts
}
if (params.reads) {
	summary['Reads'] = params.reads
}
if (params.rm_lib) {
	summary['RepeatLibrary'] = params.rm_lib
}
if (params.rm_species) {
	summary['RepeatSpecies'] = params.rm_species
}
summary['EsgaVersion'] = params.version
summary['EsgaCommitHash'] = workflow.commitId
summary['RunDate'] = workflow.start
summary['RunDir'] = workflow.workDir
summary['Command'] = workflow.commandLine
summary['Container'] = workflow.container

summary['MinContigSize'] = params.min_contig_size
summary['MinProteinLength'] = params.min_prot_length
summary['MaxIntronSize'] = params.max_intron_size

summary['AugustusSpecies'] = params.aug_species
summary['AugustusOptions'] = params.aug_options
summary['AugustusConfig'] = params.aug_config
summary['Priority Proteins'] = params.pri_prot
summary['Priority Transcripts'] = params.pri_est
summary['Priority RNAseq Introns'] = params.pri_rnaseq
summary['Priority RNAseq Exons'] = params.pri_wiggle

summary['PredictUTRs'] = params.utr

summary['PredictNcrna'] = params.ncrna

summary['RunEVM'] = params.evm
summary['EvmWeights'] = params.evm_weights

run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

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

if (params.proteins_targeted) {
	pt = file(params.proteins_targeted)
	if (!pt.exists()) {
		exit 1, "The specified targeted protein file does not exist!"
	}
	proteins_targeted = Channel.from(pt)
} else {
	proteins_targeted = Channel.empty()
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

if (!params.reads && !params.proteins && !params.proteins_targeted && !params.transcripts) {
	exit 1, "Need to specify some kind of annotation evidence for this pipeline to run!"
}

if (params.rm_species && params.rm_lib) {
	log.warn "Specified both a repeatmasker species and a library - will only use the species!"
}
if (!params.rm_species && !params.rm_lib) {
	log.warn "No repeat data provided, will model repeats de-novo (slow!!!)"
}
if (params.aug_training && !params.proteins) {
	// if no proteins are given, we will try to use pasa for training
	params.pasa = true
}
if (params.pasa && !params.transcripts && !params.reads) {
	exit 1, "Cannot run PASA without transcript data (--transcripts and/or --reads)"
}
if (params.pasa && params.reads && !params.trinity) {
	log.info "Will perform de-novo transcriptome assembly from raw reads to inform PASA annotation"
	params.trinity = true
}
if (params.polish && !params.pasa) {
	exit 1, "Cannot polish an annotation without running Pasa (--pasa, --transcripts or --trinity & reads)"
}

if (params.aug_config) {
	aug_config_file = file(params.aug_config)
	if (aug_config_file.exists()) {
		aug_extrinsic_config = Channel.fromPath(params.aug_config)
	} else {
		exit 1, "The specified Augustus extrinsic config file does not seem to exist..."
	}
} else {
	exit 1, "Augustus extrinsic config not defined, cannot proceed..."
}

// Allow input of external hints generated with ProtHint pipelines
if (params.prothint_gff && params.proteins) {
	exit 1, "Only one source of proteins allowed - please choose either --prothint_gff or --proteins"
} else if (params.prothint_gff) {
	prothint_file = file(params.prothint_gff)
	if (!prothint_file.exists()) {
		exit 1, "ProtHint file does not seem to exist..."
	}
	protein_hints= Channel.fromPath(params.prothint_gff)
}

// Path to one or more reference genomes; gtf file is assumed to share the same base name as the genome sequence...
if (params.references) {

	reference_species = Channel.fromPath(params.references)
			.ifEmpty { exit 1, "Could not find the specified reference species..." }
			.map { r -> [ r.getBaseName().toString(), r , file(r.getParent().toString() + "/" + r.getBaseName().toString() + ".gtf") ] }

} else {
	reference_species = Channel.empty()
}

// Provide the path to the augustus config folder
// If it's in a container, use the hard-coded path, otherwise the augustus env variable
if (!workflow.containerEngine) {
	Channel.from(file(System.getenv('AUGUSTUS_CONFIG_PATH')))
		.ifEmpty { exit 1; "Looks like the Augustus config path is not set? This shouldn't happen!" }
        	.set { augustus_config_folder }
} else {
// this is a bit dangerous, need to make sure this is updated when we bump to the next release version
	Channel.from(file("/opt/conda/envs/esga-1.1/config"))
        	.set { augustus_config_folder }
}
if (!params.aug_species) {
	exit 1, "Must provide an AUGUSTUS species name to run this pipeline (--aug_species)"
}

multiqc_report = Channel.empty()

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
	log.info "Repeatmasker DFam 2.0 species:	${params.rm_species}"
} else {
	log.info "Repeatmasking:			Compute de-novo"
}
log.info "-----------------------------------------"
log.info "Program settings:"
log.info "AUGUSTUS species:		${params.aug_species}"
log.info "AUGUSTUS config:		${params.aug_config}"
if (params.aug_training) {
	if (params.proteins) {
		log.info "AUGUSTUS training:		Proteins"
	} else {
		log.info "AUGUSTUS training:		Transcripts"
	}
}
log.info "Predict ncRNAs			${params.ncrna}"
log.info "Run PASA assembly:		${params.pasa}"
log.info "Run Trinity assembly:		${params.trinity}"
log.info "Run EVM gene building:		${params.evm}"
log.info "EVM weights:			${params.evm_weights}"
log.info "Predict UTRs:			${params.utr}"
log.info "-----------------------------------------"
log.info "Evidences:"
log.info "Targeted proteins		${params.proteins_targeted}"
log.info "Other proteins:			${params.proteins}"
log.info "Transcripts:			${params.transcripts}"
log.info "RNA-seq:			${params.reads}"
log.info "References:			${params.references}"
log.info "-----------------------------------------"
log.info "Parallelization settings"
log.info "# Sequences per protein alignment	${params.nproteins}"
log.info "Size of genome-level jobs:		${params.npart_size} bp"
log.info "Max intron length:			${params.max_intron_size}"
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
		repeat_gffs = repeatmasking_with_species.out.genome_rm_gffs
		repeat_hints = repeatmasking_with_species.out.genome_rm_hints
	// Use a library provided by the user or compute de-novo
	} else {
		if (params.rm_lib) {
			repeats = Channel.fromPath(params.rm_lib)
		} else {
			model_repeats(genome_clean)
			repeats = model_repeats.out.repeats
		}
	        repeatmasking_with_lib(genome_clean,repeats)
		genome_rm = repeatmasking_with_lib.out.genome_rm
		repeat_gffs = repeatmasking_with_lib.out.genome_rm_gffs
		repeat_hints = repeatmasking_with_lib.out.genome_rm_hints
	}

	// Annotate ncRNAs using RFam
	if (params.ncrna) {
		rfamsearch(genome_clean)
		ncrna_gff = rfamsearch.out.gff
	} else {
		ncrna_gff = Channel.empty()
	}

	// Generate hints from proteins (if any)
	if (params.proteins) {
		proteinhint_spaln(genome_clean,proteins)
		protein_hints = proteinhint_spaln.out.hints
		protein_evm_align = proteinhint_spaln.out.track
	} else {
		protein_hints = Channel.empty()
		protein_evm_align = Channel.empty()
	}

	// Construct gene models from species specific proteome
	if (params.proteins_targeted) {
		proteinmodels(genome_clean,proteins_targeted)
		protein_targeted_hints = proteinmodels.out.hints
                protein_targeted_gff = proteinmodels.out.gff
		protein_targeted_evm_align = proteinmodels.out.track
	} else {
		protein_targeted_hints = Channel.empty()
		protein_targeted_gff = Channel.empty()
		protein_targeted_evm_align = Channel.empty()
	}

	// Generate hints from transcripts (if any)
	if (params.transcripts) {
		esthint(genome_clean,transcripts)
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
		// Assembly reads into transcripts for PASA
		if (params.trinity || !params.transcripts && params.reads && params.pasa ) {
			trinity_guided_assembly(rnaseqhint.out.bam)
			trinity_esthint(genome_clean,trinity_guided_assembly.out.assembly)
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
		pasa(genome_clean,fastaMergeFiles.out[0])
		pasa_gff = pasa.out.gff
		pasa_fa = pasa.out.fasta
		pasa_db = pasa.out.db
		pasa_transcript_gff = pasa.out.transcript_gff
	} else {
		pasa_gff = Channel.empty()
		pasa_fa = Channel.empty()
		pasa_db = Channel.empty()
		pasa_transcript_gff = Channel.empty()
	}

	if (params.aug_training) {
		// Prefer training from Spaln models since these tend to contain less noise 
		if (params.proteins) {
			augustus_train_from_spaln(genome_rm,protein_gff,augustus_config_folder)
			augustus_conf_folder = augustus_train_from_spaln.out.acf_folder
		} else if (params.transcripts && params.pasa) {
			augustus_train_from_pasa(genome_rm,pasa_gff,augustus_config_folder)
			augustus_conf_folder = augustus_train_from_pasa.out.acf_folder
		} 
	} else {
		augustus_conf_folder = augustus_config_folder
	}

	// map existing gene models from related organisms using Kraken/Satsuma
        if (params.references) {
                map_annotation(reference_species,genome_clean)
                trans_hints = map_annotation.out.hints
		liftovers = map_annotation.out.mapped_gff
        } else {
                trans_hints = Channel.empty()
		liftovers = Channel.empty()
        }

	// Merge hints
	hints = protein_hints.concat(est_hints, trinity_hints,rna_hints,protein_targeted_hints,trans_hints)
	merge_hints(hints.collect())

	// Run AUGUSTUS
	if (!params.fast) {
		augustus_prediction_slow(genome_rm,merge_hints.out,augustus_conf_folder,aug_extrinsic_config)
                augustus_gff = augustus_prediction_slow.out.gff
                augustus_fa = augustus_prediction_slow.out.fasta
		augustus_filtered_gff = augustus_prediction_slow.out.gff_filtered
	} else {	
		augustus_prediction(genome_rm,merge_hints.out,augustus_conf_folder,aug_extrinsic_config)
		augustus_gff = augustus_prediction.out.gff
		augustus_fa = augustus_prediction.out.fasta
		augustus_filtered_gff = augustus_prediction.out.gff_filtered
	}

	// Combine all inputs into consensus annotation
	if (params.evm) {
	
		gene_gffs = augustus_filtered_gff.concat(pasa_gff, protein_targeted_gff, liftovers).collect()
		// Reconcile optional multi-branch transcript evidence into a single channel
		if (params.transcripts && params.reads && params.trinity) {
			transcript_gff = est_gff.concat(trinity_gff)
		} else if (params.transcripts) {
			transcript_gff = est_gff
		} else if (params.reads && params.trinity) {
			transcript_gff = trinity_gff
		} else {
			transcript_gff = Channel.fromPath(params.empty_gff)
		}

		// Reconcile optional multi-branch protein evidence into a single channel
		if (params.proteins && params.proteins_targeted) {
			protein_gff = protein_evm_align.concat(protein_targeted_evm_align).collectFile()
		} else if (params.proteins_targeted) {
			protein_gff = protein_targeted_evm_align
		} else if (params.proteins) {
			protein_gff = protein_evm_align
		} else {
			protein_gff = Channel.empty()
		}

		evm_prediction(genome_rm,protein_gff,transcript_gff,gene_gffs)
		evm_gff = evm_prediction.out.gff
		evm_fa = evm_prediction.out.fasta

		if (params.polish) {
			polish_annotation(genome,evm_gff,pasa_transcript_gff,pasa_db)
			polish_gff = polish_annotation.gff
		} else {
			polish_gff = Channel.empty()
		}
		
		
	} else {
		evm_gff = Channel.empty()
		evm_fa = Channel.empty()
		protein_gff = Channel.empty()
		transcript_gff = Channel.fromPath(params.empty_gff)
		polish_gff = Channel.empty()
	}

	publish:
		genome_rm to: "${params.outdir}/repeatmasking", mode: 'copy'
		repeat_gffs to: "${params.outdir}/repeatmasking", mode: 'copy'
		assembly_stats to: "${params.outdir}/assembly", mode: 'copy'
		augustus_gff to: "${params.outdir}/annotation/augustus", mode: 'copy'
		augustus_fa to: "${params.outdir}/annotation/augustus", mode: 'copy'
		est_hints to: "${params.outdir}/evidence/hints", mode: 'copy'
		est_gff to:  "${params.outdir}/evidence/transcripts", mode: 'copy'
		protein_targeted_gff to: "${params.outdir}/annotation/spaln", mode: 'copy'
		protein_targeted_hints to: "${params.outdir}/evidence/hints", mode: 'copy'
		protein_hints to: "${params.outdir}/evidence/hints", mode: 'copy'
		protein_gff to: "${params.outdir}/evidence/proteins", mode: 'copy'
		rna_hints to: "${params.outdir}/evidence/hints", mode: 'copy'
		rna_bam to: "${params.outdir}/evidence/rnaseq", mode: 'copy'
		repeats to: "${params.outdir}/repeatmasking", mode: 'copy'
		evm_gff to: "${params.outdir}/annotation/evm", mode: 'copy'
		evm_fa to: "${params.outdir}/annotation/evm", mode: 'copy'
		pasa_gff to: "${params.outdir}/annotation/pasa", mode: 'copy'
		pasa_fa to: "${params.outdir}/annotation/pasa", mode: 'copy'
		polish_gff to: "${params.outdir}/annotation/evm", mode: 'copy'
		ncrna_gff to: "${params.outdir}/annotation/ncrna", mode: 'copy'
		protein_targeted_evm_align to: "${params.outdir}/evidence_modeler", mode: 'copy'
		protein_evm_align to: "${params.outdir}/evidence_modeler", mode: 'copy'
		transcript_gff to: "${params.outdir}/evidence_modeler", mode: 'copy'
}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="

  def email_fields = [:]
  email_fields['version'] = workflow.manifest.version
  email_fields['session'] = workflow.sessionId
  email_fields['runName'] = run_name
  email_fields['success'] = workflow.success
  email_fields['dateStarted'] = workflow.start
  email_fields['dateComplete'] = workflow.complete
  email_fields['duration'] = workflow.duration
  email_fields['exitStatus'] = workflow.exitStatus
  email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
  email_fields['errorReport'] = (workflow.errorReport ?: 'None')
  email_fields['commandLine'] = workflow.commandLine
  email_fields['projectDir'] = workflow.projectDir
  email_fields['script_file'] = workflow.scriptFile
  email_fields['launchDir'] = workflow.launchDir
  email_fields['user'] = workflow.userName
  email_fields['Pipeline script hash ID'] = workflow.scriptId
  email_fields['genome'] = params.genome
  email_fields['manifest'] = workflow.manifest
  email_fields['summary'] = summary

  email_info = ""
  for (s in email_fields) {
	email_info += "\n${s.key}: ${s.value}"
  }

  def output_d = new File( "${params.outdir}/pipeline_info/" )
  if( !output_d.exists() ) {
      output_d.mkdirs()
  }

  def output_tf = new File( output_d, "pipeline_report.txt" )
  output_tf.withWriter { w -> w << email_info }	

 // make txt template
  def engine = new groovy.text.GStringTemplateEngine()

  def tf = new File("$baseDir/assets/email_template.txt")
  def txt_template = engine.createTemplate(tf).make(email_fields)
  def email_txt = txt_template.toString()

  // make email template
  def hf = new File("$baseDir/assets/email_template.html")
  def html_template = engine.createTemplate(hf).make(email_fields)
  def email_html = html_template.toString()
  
  def subject = "ESGA annotation finished ($run_name)."

  if (params.email) {

  	def mqc_report = null
  	try {
        	if (workflow.success && !params.skip_multiqc) {
            		mqc_report = multiqc_report.getVal()
            		if (mqc_report.getClass() == ArrayList){
                		log.warn "[ESGA] Found multiple reports from process 'multiqc', will use only one"
                		mqc_report = mqc_report[0]
                	}
        	}
    	} catch (all) {
        	log.warn "[ESGA] Could not attach MultiQC report to summary email"
  	}

	def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
	def sf = new File("$baseDir/assets/sendmail_template.txt")	
    	def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    	def sendmail_html = sendmail_template.toString()

	try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
        }

  }

}

