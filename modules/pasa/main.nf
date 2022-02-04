// ****************
// PASA transcriptome assembly
// ****************

include { fastaCleanProteins; fastaSplitSize ; fastaCleanNames } from "./../fasta"  params(params)
include { estMinimap; estMinimapToGff } from "./../transcripts/main.nf" params(params)
include { GffToFasta } from "./../util" addParams(folder: "${params.outdir}/annotation/pasa")

// Run the PASA pipeline
workflow pasa {

	take:
		genome
		transcripts

	main:
		fastaCleanProteins(transcripts)	
		fastaCleanNames(fastaCleanProteins.out)
		seqclean(fastaCleanNames.out)
		estMinimap(seqclean.out[0],genome)
		estMinimapToGff(estMinimap.out)
		pasa_assembly(genome,seqclean.out[0],seqclean.out[1],fastaCleanNames.out)
		PasaToModels(pasa_assembly.out[0],pasa_assembly.out[1])
		GffToFasta(PasaToModels.out[1],genome)	

	emit:
		gff = PasaToModels.out[1]
		alignments = pasa_assembly.out[1]
		fasta = GffToFasta.out[0]
		db = pasa_assembly.out[2]
		transcript_gff = estMinimap.out[0]

}

// Add UTRs do a gene build
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
process seqclean {

	label 'pasa'

	input:
	path transcripts

	output:
	path transcripts_clean
	path transcripts_cln
	script:
	transcripts_clean = transcripts.getName() + ".clean"
	transcripts_cln = transcripts.getName() + ".cln"

	"""
		export USER=${workflow.userName}
		seqclean $transcripts -c ${task.cpus}
	"""
}

// Run the PASA pipeline from pre-aligned sequences (estMinimap)
// Using the built-in alignment is way too slow!
process pasa_assembly {

	label 'pasa'

	scratch true 

        publishDir "${params.outdir}/logs/pasa", mode: 'copy'

	input:
	path genome
	path transcripts_clean
	path transcripts_cln
	path transcripts_untrimmed

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
        mysql_db_name = params.pasa_mysql_db_name
        if (params.pasa_mysql_user) {
                mysql_options = "make_pasa_mysql_config.pl --infile \$PASAHOME/pasa_conf/conf.txt --outfile pasa_mysql_conf.txt --user ${params.pasa_mysql_user} --pass ${params.pasa_mysql_pass} --host ${params.pasa_mysql_host} --port ${params.pasa_mysql_port}"
                mysql_config_option = "-C pasa_mysql_conf.txt"
                mysql_db_name = "--mysql pasa_db"
        }

	"""
		make_pasa_config.pl --infile ${params.pasa_config} --trunk $trunk --outfile pasa_DB.config $mysql_db_name

		\$PASAHOME/Launch_PASA_pipeline.pl \
			--ALIGNERS ${params.pasa_aligner} \
                        -c pasa_DB.config -C -R \
                        -t $transcripts_clean \
                        -I $params.max_intron_size \
			--transcribed_is_aligned_orient \
                        -g $genome \
                        --CPU ${task.cpus} \
	"""

}

// Turn the pasa results into full gff3 file
process PasaToModels {

	label 'pasa'

	//scratch true 
        publishDir "${params.outdir}/annotation/pasa", mode: 'copy'

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
		
// Use Pasa to polish a gene build with transcript data
// requires a gene build in gff3 format, aligned transcripts and the genome
process runPasaPolish {

	label 'pasa'

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
