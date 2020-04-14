include prepHintsToBed from "./../util" params(params)
include fastaSplitChunks from "./../fasta" params(params)

// Predict gene models on a genome sequence using hints
// We need the augustus_config_dir so we can create custom-trained models and pass them to Augustus in the context of a container
workflow augustus_prediction {

	take:
		genome
		hints
		augustus_config_dir

	main:
		fastaSplitChunks(genome,params.nchunks)
		prepHintsToBed(hints)
		prepAugustusConfig(augustus_config_dir)
		runAugustusBatch(fastaSplitChunks.out.flatMap(),prepHintsToBed.out,prepAugustusConfig.out.collect().map{ it[0].toString() } )
		mergeAugustusGff(runAugustusBatch.out.collect())

	emit:
		gff = mergeAugustusGff.out
		config = prepAugustusConfig.out

}

process runAugustusBatch {

	input:
	path genome_chunk
	path hints
	path regions
	env AUGUSTUS_CONFIG_PATH

	output:
	path augustus_result

	script:
	chunk_name = genome_chunk.getName().tokenize("_")[-1]
	augustus_result = "augustus.${chunk_name}.out.gff"
	genome_fai = genome_chunk.getName() + ".fai"
	command_file = "commands." + chunk_name + ".txt"

	"""
		samtools faidx $genome_chunk
		fastaexplode -f $genome_chunk -d . 
		augustus_from_regions.pl --genome_fai $genome_fai --model $params.aug_species --utr ${params.utr} --isof false --aug_conf ${params.aug_config} --hints $hints --bed $regions > $command_file
		parallel -j ${task.cpus} < $command_file
		cat *augustus.gff > $augustus_result
		rm *augustus.gff
	"""
}

// Merge all the chunk GFF files into one file
process mergeAugustusGff {

	label 'short_running'

	input:
	path augustus_gffs

	output:
	path augustus_merged_gff

	script:
	augustus_merged_gff = "augustus.merged.out.gff"
	
	"""	
		cat $augustus_gffs >> merged.gff
		create_gff_ids.pl --gff merged.gff > $augustus_merged_gff
	"""
}

// containerized AUGUSTUS_CONFIG_PATH does not work, need to move into work/
process prepAugustusConfig {

	input:
	path augustus_config_dir

	output:
	path copied_dir

	script:
	copied_dir = "augustus_config"

	"""
		mkdir -p $copied_dir
		cp -R $augustus_config_dir/* $copied_dir/
	"""

}
