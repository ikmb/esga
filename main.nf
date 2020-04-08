#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include repeatmasking from "./modules/repeatmasker/repeatmasking" params(params)
include proteinhint from "./modules/proteins/proteinevidence" params(params)
include esthint from "./modules/transcripts/estevidence" params(params)
include augustus_prediction from "./modules/augustus/augustus" params(params)
include merge_hints from "./modules/util" params(params)
//include rnaseqhint from "./modules/rnaseq" params(params)
//include evm from "./modules/evidencemodeler" params(params)
//include pasa from "./modules/pasa" params(params)

// ****************
// input validation
// ****************

if (params.proteins) {
	proteins = file(params.proteins)
	if (!proteins.exists()) {
		exit 1, "The specified protein file does not exist!"
	}
}
if (params.ESTs) {
	ests = file(params.ESTs)
	if (!ests.exists() ) {
		exit 1, "The specified EST file does not exist!"
	}
}

// Provide the path to the augustus config folder
// If it's in a container, use the hard-coded path, otherwise the augustus env variable
if (!workflow.containerEngine) {
	Channel.from(file(System.getenv('AUGUSTUS_CONFIG_PATH')))
		.ifEmpty { exit 1; "Looks like the Augustus config path is not set? This shouldn't happen!" }
        	.set { augustus_config_folder }
} else {
// this is a bit dangerous, need to make sure this is updated when we bump to the next release version
	Channel.from(file("/opt/conda/envs/genome-annotation-1.0/config"))
        	.set { augustus_config_folder }
}

workflow {

        repeatmasking(params.genome,params.rm_lib)
	
	if (params.proteins) {
		proteinhint(repeatmasking.out.genome_rm,proteins)
		protein_hints = proteinhint.out.hints
	} else {
		protein_hints = Channel.empty()
	}

	if (params.ESTs) {
		esthint(repeatmasking.out.genome_rm,ests)
		est_hints = esthint.out.hints
	} else {
		est_hints = Channel.empty()
	}

	if (params.reads) {
		rnaseqhint(repeatmasking.out.genome_rm,params.reads)
	}

	merge_hints(protein_hints,est_hints)
	
	augustus_prediction(repeatmasking.out.genome_rm,merge_hints.out,augustus_config_folder)
}



