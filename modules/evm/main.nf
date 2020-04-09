workflow evm {

	take:
		genome
		gene_models
		protein_gff
		transcript_gff
	main:

	emit:
}

process predEvmPartition {

	label 'long_running'

	// CANNOT USE SCRATCH FOR THIS!
	//scratch true

	//publishDir "${OUTDIR}/annotation/evm/jobs", mode: 'copy'

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
	transcripts = "transcripts.merged.gff"
        gene_models = "gene_models.gff"

	protein_options = ""
	transcript_options = ""
	if (protein_gff) {
		protein_options = "--protein_alignments $protein_gff "
	}
	if ( est_gff ) {
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

	// run n commands per chunk
	evm_command_chunks = inputToEvm.splitText(by: params.nevm, file: true)

	// The outputs doesn't do nothing here, EVM combines the chunks based on the original partition file
	process predEvm {

		scratch true

		label 'long_running'

                //publishDir "${OUTDIR}/annotation/evm/chunks", mode: 'copy'

		input:
		file(evm_chunk) from evm_command_chunks

		output:
		file(log_file) into EvmOut

		script:
		log_file = evm_chunk.getBaseName() + ".log"
		"""
			\$EVM_HOME/EvmUtils/execute_EVM_commands.pl $evm_chunk | tee $log_file
		"""
		
	}
	
	// Merge all the separate partition outputs into a final gff file
	process predEvmMerge {
	
		label 'long_running'

		publishDir "${OUTDIR}/annotation/evm", mode: 'copy'

		input:

		file(evm_logs) from EvmOut.collect()
		file(partitions) from EvmPartition
		set file(genome_rm),file(genome_index) from genome_to_evm_merge

		output:
		file(partitions) into EvmResult
		file(done)

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
	process predEvmToGff {

		label 'medium_running'

		publishDir "${OUTDIR}/annotation/evm", mode: 'copy'

		input:
		file(partitions) from EvmResult

		output:
		file(evm_gff) into EvmGFF

		script:
		evm_gff = "annotations.evm.gff"

		"""
			merge_evm_gff.pl --partitions $partitions --gff $evm_gff
		"""

	}
