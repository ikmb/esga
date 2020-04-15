include fastaSplitChunks from "./../fasta"  params(params)
include estMinimap from "./../transcripts/main.nf" params(params)

workflow pasa_assembly {

	take:
		genome
		transcripts

	main:
		fastaSplitChunks(genome, params.nchunks)
		runSeqClean(transcripts)
		estMinimap(runSeqClean.out[0],genome)
		runMinimapSplit(fastaSplitChunks.out.flatMap(),runSeqClean.out.collect(),estMinimap.out.collect())
		runPasaFromCustom(runMinimapSplit.out[0],runMinimapSplit.out[1],runMinimapSplit.out[2])
		PasaToModels(runPasaFromCustom.out[0].collect(),runPasaFromCustom.out[1].collect())

	emit:
		gff = PasaToModels.out[1]
		alignments = PasaToModels.out[2]

}

// Currently does not work in singularity/conda so we just copy the input until we can fix this
process runSeqClean {

	label 'short_running'

	input:
	path transcripts

	output:
	path transcripts_clean

	script:
	transcripts_clean = transcripts.getName() + ".clean"

	"""
		cp $transcripts $transcripts_clean
	"""
}
		
// We parallelize PASA by filtering the minimap alignments per genome chunk
// meaning we select only those alignments present in this part of the assembly
process runMinimapSplit {
	
	label 'short_running'

	input:
	path genome_chunk
	path transcripts
	path minimap_gff

	output:
	path genome_chunk
	path transcripts_minimap
	path minimap_chunk
			
	script:
	minimap_chunk = genome_chunk.getBaseName() + ".minimap.gff"
	transcripts_minimap = genome_chunk.getBaseName() + ".transcripts.fasta"
	genome_chunk_index = genome_chunk + ".fai"
	// filter the gff file to only contain entries for our scaffolds of interest
	// then make a list of all transcript ids and extract them from the full transcript fasta
	"""
		samtools faidx $genome_chunk
		minimap_filter_gff_by_genome_index.pl --index $genome_chunk_index --gff $minimap_gff --outfile  $minimap_chunk
		minimap_gff_to_accs.pl --gff $minimap_chunk | sort -u > list.txt
		faSomeRecords $transcripts list.txt $transcripts_minimap
	"""
}

// Run the PASA pipeline
process runPasaFromCustom {
			
	scratch true
			
	input:
	path genome_rm
	path transcripts
	path custom_gff
	
	output:
	path pasa_assemblies_fasta
	path pasa_assemblies_gff

	script:
	trunk = genome_rm.getBaseName()
	pasa_assemblies_fasta = "pasa_DB_${trunk}.sqlite.assemblies.fasta"
	pasa_assemblies_gff = "pasa_DB_${trunk}.sqlite.pasa_assemblies.gff3"
			
	// optional MySQL support
	// create a config file with credentials, add config file to pasa execute
	// and somehow check that the DB doesn't already exists
	mysql_create_options = ""
	mysql_config_option = ""
	mysql_db_name = ""
	if (params.pasa_mysql_user) {
		mysql_options = "make_pasa_mysql_config.pl --infile \$PASAHOME/pasa_conf/conf.txt --outfile pasa_mysql_conf.txt --user ${params.pasa_mysql_user} --pass ${params.pasa_mysql_pass} --host ${params.pasa_mysql_host} --port ${params.pasa_mysql_port}"	
		mysql_config_option = "-C pasa_mysql_conf.txt"
		mysql_db_name = "--mysql $run_name"
	}

	// The pasa sqlite file must have a fully qualified path, the script is a workaround as this seems difficult to do inside 
	// the script statement

	"""
		make_pasa_config.pl --infile $pasa_config --trunk $trunk --outfile pasa_DB.config $mysql_db_name
		$mysql_create_options
		\$PASAHOME/Launch_PASA_pipeline.pl \
			-c pasa_DB.config -C -R \
			-t $transcripts \
			-I $params.max_intron_size \
			-g $genome_rm \
			--IMPORT_CUSTOM_ALIGNMENTS_GFF3 $custom_gff \
			--CPU ${task.cpus} \
			$mysql_config_option
	"""	
}

// Extract gene models from PASA database
// All chunks are merged using a perl script into the base name pasa_db_merged
// this is not...ideal. 
process PasaToModels {

	scratch true

	label 'long_running'

	input:
	path pasa_assemblies_fasta
	path pasa_assemblies_gff

	output:
	path pasa_transdecoder_fasta
	path pasa_transdecoder_gff
	path merged_gff 

	script:
	base_name = "pasa_db_merged"
	merged_fasta = base_name + ".assemblies.fasta"
	merged_gff = base_name + ".pasa_assemblies.gff"
	pasa_transdecoder_fasta = merged_fasta + ".transdecoder.pep"
	pasa_transdecoder_gff = merged_fasta + ".transdecoder.genome.gff3"

	script:

	"""
		pasa_merge_chunks.pl --base $base_name
		\$PASAHOME/scripts/pasa_asmbls_to_training_set.dbi \
		--pasa_transcripts_fasta $merged_fasta \
		--pasa_transcripts_gff3 $merged_gff \
                     	
	"""
}
