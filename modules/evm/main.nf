include GffToFasta from "./../util" params(params)

workflow evm_prediction {

	take:
		genome
		protein_gff
		transcript_gff
		gene_models
	main:
		evmMergeGenes(gene_models)
		evmPartition(genome,evmMergeGenes.out[0],transcript_gff,protein_gff)
		runEvm(evmPartition.out[0].splitText(by: params.nevm, file: true))
		evmMerge(runEvm.out.collect(),evmPartition.out[1],genome)
		evmToGff(evmMerge.out[0].collect())
		GffToFasta(evmToGff.out[0],genome)

	emit:
		gff = evmToGff.out
		fasta = GffToFasta.out[0]
}

process evmMergeGenes {

	label 'short_running'

	input:
	path gene_models

	output:
	path merged_models

	script:
	merged_models = "genes.merged.gff"

	"""
		cat $gene_models | grep -v "^#" >> $merged_models
	"""

}

process evmPartition {

	label 'long_running'

	// DO NOT USE SCRATCH WITH THIS!

	input:
	path genome_rm
	path gene_models
	path est_gff
	path protein_gff

	output:
	path evm_commands
	path partitions

	script:

	partitions = "partitions_list.out"
	evm_commands = "commands.list"

	protein_options = ""
	transcript_options = ""
	if (protein_gff) {
		protein_options = "--protein_alignments $protein_gff "
	}
	if (est_gff) {
		transcript_options = "--transcript_alignments $est_gff "
	}

	"""
		\$EVM_HOME/EvmUtils/partition_EVM_inputs.pl --genome $genome_rm \
			--gene_predictions $gene_models \
			--segmentSize 2000000 --overlapSize 200000 --partition_listing $partitions \
			$protein_options $transcript_options
				
		\$EVM_HOME/EvmUtils/write_EVM_commands.pl --genome $genome_rm \
			--weights ${params.evm_weights} \
			--gene_predictions $gene_models \
			--output_file_name evm.out \
			--partitions $partitions > $evm_commands
	"""	
			
}

// The outputs don't do nothing here, EVM combines the chunks based on the original partition file
process runEvm {

	//scratch true

	label 'long_running'


	input:
	path evm_chunk	

	output:
	path log_file

	script:
	log_file = evm_chunk.getBaseName() + ".log"
	"""
		\$EVM_HOME/EvmUtils/execute_EVM_commands.pl $evm_chunk | tee $log_file
	"""
		
}
	
// Merge all the separate partition outputs into a final gff file
process evmMerge {
	
	label 'long_running'

	input:
	path evm_logs
	path partitions
	path genome_rm
	
	output:
	path partitions
	path done

	script:
	evm_final = "evm.out"
	done = "done.txt"
	"""
		\$EVM_HOME/EvmUtils/recombine_EVM_partial_outputs.pl --partitions $partitions --output_file_name evm.out
		\$EVM_HOME/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions $partitions --output evm.out --genome $genome_rm
		touch $done
	"""

}

// We merge the partial gffs from the partitions with a perl script, 
// since the output folders are not transferred between processes
process evmToGff {

	label 'medium_running'

	//publishDir "${params.outdir}/annotation/evm", mode: 'copy'

	input:
	path partitions

	output:
	path evm_gff

	script:
	evm_gff = "annotations.evm.gff"

	"""
		merge_evm_gff.pl --partitions $partitions --gff $evm_gff
	"""

}
