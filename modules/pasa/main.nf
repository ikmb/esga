include fastaSplitSize from "./../fasta"  params(params)
include estMinimap from "./../transcripts/main.nf" params(params)
include GffToFasta from "./../util" params(params)

workflow pasa {

	take:
		genome
		transcripts

	main:
		runSeqClean(transcripts)
		estMinimap(runSeqClean.out[0],genome)
		runPasa(genome,transcripts,estMinimap.out[0])
		PasaToModels(runPasa.out[0],runPasa.out[1])
		GffToFasta(PasaToModels.out[1],genome)	

	emit:
		gff = PasaToModels.out[1]
		alignments = runPasa.out[1]
		fasta = GffToFasta.out[0]

}

// Currently does not work in singularity/conda so we just copy the input until we can fix this
process runSeqClean {

	label 'pasa'

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

process runPasa {

	label 'pasa'

	input:
	path genome
	path transcripts
	path minimap_gff

	output:
	path pasa_assemblies_fasta
        path pasa_assemblies_gff

	script:

        trunk = genome.getName()

	pasa_assemblies_fasta = "pasa_DB_${trunk}.sqlite.assemblies.fasta"
        pasa_assemblies_gff = "pasa_DB_${trunk}.sqlite.pasa_assemblies.gff3"

	mysql_create_options = ""
        mysql_config_option = ""
        mysql_db_name = ""
        if (params.pasa_mysql_user) {
                mysql_options = "make_pasa_mysql_config.pl --infile \$PASAHOME/pasa_conf/conf.txt --outfile pasa_mysql_conf.txt --user ${params.pasa_mysql_user} --pass ${params.pasa_mysql_pass} --host ${params.pasa_mysql_host} --port ${params.pasa_mysql_port}"
                mysql_config_option = "-C pasa_mysql_conf.txt"
                mysql_db_name = "--mysql $run_name"
        }

	"""
		make_pasa_config.pl --infile ${params.pasa_config} --trunk $trunk --outfile pasa_DB.config $mysql_db_name

		\$PASAHOME/Launch_PASA_pipeline.pl \
                        -c pasa_DB.config -C -R \
                        -t $transcripts \
                        -I $params.max_intron_size \
                        -g $genome \
                        --IMPORT_CUSTOM_ALIGNMENTS_GFF3 $minimap_gff \
                        --CPU ${task.cpus} \
	"""

}

process PasaToModels {


        label 'pasa'

        //publishDir "${params.outdir}/annotation/pasa", mode: 'copy'

        input:
        path pasa_assemblies_fasta
        path pasa_assemblies_gff

        output:
        path pasa_transdecoder_fasta
        path pasa_transdecoder_gff

        script:
        pasa_transdecoder_fasta = pasa_assemblies_fasta.getName() + ".transdecoder.pep"
        pasa_transdecoder_gff = pasa_assemblies_fasta.getName() + ".transdecoder.genome.gff3"

        script:

        """
                \$PASAHOME/scripts/pasa_asmbls_to_training_set.dbi \
                --pasa_transcripts_fasta $pasa_assemblies_fasta \
                --pasa_transcripts_gff3 $pasa_assemblies_gff \

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
