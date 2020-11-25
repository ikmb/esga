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
		db = runPasa.out[2]

}

workflow polish_annotation {

	take:
		genome
		genes
		transcripts
		db

	main:
		runPasaPolish(genome,genes,transcripts)

	emit:
		gff = runPasaPolish.out
	
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
	path db_name

	script:

        trunk = genome.getBaseName()

	pasa_assemblies_fasta = "pasa_DB_${trunk}.sqlite.assemblies.fasta"
        pasa_assemblies_gff = "pasa_DB_${trunk}.sqlite.pasa_assemblies.gff3"
	db_name = "pasa_DB_" + trunk + ".sqlite"

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

        publishDir "${params.outdir}/logs/pasa", mode: 'copy'

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
		
process runPasaPolish {

	publishDir "${params.outdir}/annotation/pasa", mode: 'copy'

	input:
	path genome
	path genes_gff
	path minimap_gff
	path pasa_db

	output:
        path "*gene_structures_post_PASA_updates*.gff3"

	script:
	
	trunk = genome.getBaseName() + "_pasa"

	"""

		make_pasa_config.pl --infile ${params.pasa_config} --trunk $trunk --outfile pasa_DB.assemble.config
	
		\$PASAHOME/scripts/Load_Current_Gene_Annotations.dbi \
			-P $genes_gff \
			-g $genome \
			-c pasa_DB.assemble.config

		\$PASAHOME/Launch_PASA_pipeline.pl \
		        -c pasa_DB.config -A \
		        -g $genome \
                        --IMPORT_CUSTOM_ALIGNMENTS_GFF3 $minimap_gff \
                        --CPU ${task.cpus} \
                        -I $params.max_intron_size \
		
	"""

}
