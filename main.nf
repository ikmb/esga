#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include repeatmasking from "./modules/repeatmasker/main.nf" params(params)
include proteinhint from "./modules/proteins/main.nf" params(params)
include esthint from "./modules/transcripts/main.nf" params(params)
include esthint as trinity_esthint from "./modules/transcripts/main.nf" params(params)
include augustus_prediction from "./modules/augustus/main.nf" params(params)
include merge_hints from "./modules/util" params(params)
include rnaseqhint from "./modules/rnaseq/main.nf" params(params)
include trinity_guided_assembly from "./modules/trinity/main.nf" params(params)
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

	if (params.rm_lib) {
		repeats = Channel.fromPath(params.rm_lib)
	} else {
		model_repeats(params.genome)
		repeats = model_repeats.out.repeats
	}

        repeatmasking(params.genome,repeats)
	
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
		rnaseqhint(params.genome,params.reads)
		rna_hints = rnaseqhint.out.hints
		if (params.trinity) {
			trinity_guided_assembly(rnaseqhint.out.bam)
			trinity_esthint(repeatmasking.out.genome_rm,trinity_guided_assembly.out.assembly)
			trinity_hints = trinity_esthint.out.hints
		} else {
			trinity_hints = Channel.empty()
		}	
		
	} else {
		rna_hints = Channel.empty()
		trinity_hints = Channel.empty()
	}

	hints = protein_hints.merge(est_hints).merge(trinity_hints).merge(rna_hints)
	merge_hints(hints)
	
	augustus_prediction(repeatmasking.out.genome_rm,merge_hints.out,augustus_config_folder)

}



