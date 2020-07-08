include { prepHintsToBed ; GffToFasta } from "./../util" params(params)
include fastaSplitSize from "./../fasta" params(params)

// Predict gene models on a genome sequence using hints
// We need the augustus_config_dir so we can create custom-trained models and pass them to Augustus in the context of a container
workflow augustus_prediction {

	take:
		genome
		hints
		augustus_config_dir

	main:
		fastaSplitSize(genome,params.npart_size)
                prepHintsToBed(hints)
		prepAugustusConfig(augustus_config_dir)
		runAugustusBatch(fastaSplitSize.out.flatMap(),prepHintsToBed.out,prepAugustusConfig.out.collect().map{ it[0].toString() } )
		mergeAugustusGff(runAugustusBatch.out.collect())
		GffToFasta(mergeAugustusGff.out[0],genome)

	emit:
		gff = mergeAugustusGff.out
		config = prepAugustusConfig.out
		fasta = GffToFasta.out[0]

}

workflow augustus_prediction_slow {

	take:
                genome
                hints
                augustus_config_dir

        main:
                fastaSplitSize(genome,params.npart_size)
                prepAugustusConfig(augustus_config_dir)
                runAugustus(fastaSplitSize.out.flatMap(),hints,prepAugustusConfig.out.collect().map{ it[0].toString() } )
                mergeAugustusGff(runAugustus.out.collect())
                GffToFasta(mergeAugustusGff.out[0],genome)

        emit:
                gff = mergeAugustusGff.out
                config = prepAugustusConfig.out
                fasta = GffToFasta.out[0]

}

workflow augustus_training_from_pasa {

	take:
		genome
		pasa_gff
		model
		augustus_config

	main:

		PasaGffToTraining(pasa_gff)
		trainAugustus(genome,PasaGffToTraining.out[0],augustus_config,augustus_config)

	emit:
		acf_folder = trainAugustus.out[0]
		stats = trainAugustus.out[1]

}

process PasaGffToTraining {

	label 'short_running'

	input:
	path pasa_gff

	output:
	path training_gff

	script:
	training_gff = "transdec.complete.gff3"

	"""
		pasa_select_training_models.pl --nmodels $params.aug_training_models --infile $pasa_transdecoder_gff >> $training_gff
        """
}

// Run one of two training routines for model training
process trainAugustus {

	label 'extra_long_running'

	scratch true

	//publishDir "${params.outdir}/augustus/training/", mode: 'copy'

	input:
	path genome
	path complete_models
	env AUGUSTUS_CONFIG_PATH
	path acf_folder

        output:
	path acf_folder
	path training_stats

	script:
        complete_gb = "complete_peptides.raw.gb"
        train_gb = "complete_peptides.raw.gb.train"
        test_gb = "complete_peptides.raw.gb.test"
        training_stats = "training_accuracy.out"

        // If the model already exists, do not run new_species.pl
        //model_path = "${acf_training_path}/species/${params.aug_species}"
        model_file = file("${acf_folder}/species/${params.aug_species}")
	options = ""
	if (!model_file.exists()) {
		options = "new_species.pl --species=${params.aug_species}"
	}

	"""
		echo ${acf_folder.toString()} >> training.txt
		gff2gbSmallDNA.pl $complete_models $genome 1000 $complete_gb
		split_training.pl --infile $complete_gb --percent 90
		$options
		etraining --species=$params.aug_species --stopCodonExcludedFromCDS=false $train_gb
		optimize_augustus.pl --species=$params.aug_species $train_gb --cpus=${task.cpus} --UTR=off 
		augustus --stopCodonExcludedFromCDS=false --species=$params.aug_species $test_gb | tee $training_stats
	"""
}

process runAugustus {

	input:
        path genome_chunk
        path hints
        env AUGUSTUS_CONFIG_PATH

        output:
        path augustus_result

        script:
        chunk_name = genome_chunk.getName().tokenize("_")[-1]
        augustus_result = "augustus.${chunk_name}.out.gff"

        """
		augustus --species=${params.aug_species} --alternatives-from-sampling=false --alternatives-from-evidence=false --hintsfile=$hints --gff3=on --UTR=${params.utr} --extrinsicCfgFile=${params.aug_config} --uniqueGeneId=true $genome_chunk > $augustus_result
 
        """

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

	//publishDir "${params.outdir}/annotation/augustus", mode: 'copy'

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
