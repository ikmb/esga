include { prepHintsToBed ; GffToFasta } from "./../util" params(params)
include fastaSplitSize from "./../fasta" params(params)

// Predict gene models on a genome sequence using hints
// We need the augustus_config_dir so we can create custom-trained models and pass them to Augustus in the context of a container
// Defunct, consider using augustus_parallel
workflow augustus_prediction {

	take:
		genome
		hints
		augustus_config_dir
		aug_extrinsic_config

	main:
		fastaSplitSize(genome,params.npart_size)
                prepHintsToBed(hints)
		runAugustusBatch(fastaSplitSize.out.flatMap(),prepHintsToBed.out,augustus_config_dir.out.collect().map{ it[0].toString() },aug_extrinsic_config.collect() )
		mergeAugustusGff(runAugustusBatch.out.collect())
		GffToFasta(mergeAugustusGff.out[0],genome)
		AugustusFilterModels(mergeAugustusGff.out[0],genome)
	emit:
		gff = mergeAugustusGff.out
		fasta = GffToFasta.out[0]
		gff_filtered = AugustusFilterModels.out[0]
}

// Run an AUGUSTUS annotation using the most exhaustive settings
// This is the slowest option!
workflow augustus_prediction_slow {

	take:
                genome
                hints
                augustus_config_dir
		aug_extrinsic_config

        main:
                fastaSplitSize(genome,params.npart_size)
                runAugustus(fastaSplitSize.out.flatMap(),hints,augustus_config_dir.collect().map{ it[0].toString() },aug_extrinsic_config.collect() )
                mergeAugustusGff(runAugustus.out.collect())
                GffToFasta(mergeAugustusGff.out[0],genome)
                AugustusFilterModels(mergeAugustusGff.out[0],genome)		

        emit:
                gff = mergeAugustusGff.out
                fasta = GffToFasta.out[0]
		gff_filtered = AugustusFilterModels.out[0]

}

// Run AUGUSTUS in parallel, splitting contigs into chunks and merging them later
// This is the fastest option!
workflow augustus_parallel {

	take:
		genome
		hints
		augustus_config_dir
		aug_extrinsic_config

	main:
		fastaSplitSize(genome,params.npart_size)
                runAugustusChunks(fastaSplitSize.out.flatMap(),hints.collect(),augustus_config_dir.collect().map{ it[0].toString() },aug_extrinsic_config.collect() )
                joinAugustusChunks(runAugustusChunks.out)
		mergeAugustusGff(joinAugustusChunks.out.collect())
                GffToFasta(mergeAugustusGff.out[0],genome)
                AugustusFilterModels(mergeAugustusGff.out[0],genome)
	emit:
		gff = mergeAugustusGff.out
                fasta = GffToFasta.out[0]
                gff_filtered = mergeAugustusGff.out

}

// Stage the AUGUSTUS config so we can modify it later
workflow augustus_prep_config {

	take:
		augustus_config_dir

	main:
                prepAugustusConfig(augustus_config_dir)

	emit:
		config = prepAugustusConfig.out

}

// Search for putative genic regions with AUGUSTUS
// Does not require evidences, just an existing profile
workflow augustus_prescan {

	take:
		genome
		augustus_config_dir

	main:
                fastaSplitSize(genome,params.npart_size)
                runAugustusScan(fastaSplitSize.out.flatMap(),augustus_config_dir.collect().map{ it[0].toString() },aug_extrinsic_config.collect() )
                mergeAugustusGff(runAugustusScan.out.collect())
                GffToFasta(mergeAugustusGff.out[0],genome)

        emit:
                gff = mergeAugustusGff.out
                fasta = GffToFasta.out[0]
		
}

// Train a novel AUGUSTUS profile from SPALN alignments
workflow augustus_train_from_spaln {

	take:
		genome
		spaln_gff
		augustus_config_dir

	main:
		SpalnGffToTraining(spaln_gff)
		trainAugustus(genome,SpalnGffToTraining.out,augustus_config_dir.collect().map{ it[0].toString() },augustus_config_dir)

	emit:
		acf_folder = trainAugustus.out[0]
		stats = trainAugustus.out[1]

}

// Train novel AUGUSTUS profile from PASA gene models
workflow augustus_train_from_pasa {

	take:
		genome
		pasa_gff
		augustus_config_dir

	main:
		PasaGffToTraining(pasa_gff)
		trainAugustus(genome,PasaGffToTraining.out[0],augustus_config_dir.collect().map{ it[0].toString() },augustus_config_dir)

	emit:
		acf_folder = trainAugustus.out[0]
		stats = trainAugustus.out[1]

}

// Filter a set of augustus predictions to remove spurious models
process AugustusFilterModels {

        publishDir "${params.outdir}/logs/augustus", mode: 'copy'

	input:
	path gff
	path genome

	output:
	path gff_good
	path gff_bad
	path proteins

	script:
	gff_good = gff.getName() + ".good.gff"
	gff_bad = gff.getName() + ".bad.gff"
	proteins = "augustus.protein_supported_models.proteins.fa"
	def support = 1.0
	// FIX THIS
	// Filtering makes no sense if we do not have evidences
	def m = "PE"
	if (params.evidence) {
		if (!params.proteins && !params.proteins_targeted) {
			m = "E"
		}
	}
	if (params.evidence) {
		"""
			augustus_filter_models.pl --infile $gff --support $support --mode $m
			gffread -g $genome -y augustus.protein_supported_models.proteins.fa $gff_good
		"""
	} else {
		"""
			gffread -g $genome -y $proteins $gff
			cp $gff $gff_good
			touch $gff_bad
		"""
	}
}

// Convert a SPALN alignment to full GFF format
// Can then be used for training AUGUSTUS
process SpalnGffToTraining {

	label 'short_running'

        publishDir "${params.outdir}/logs/augustus", mode: 'copy'

	input:
	path spaln_gff

	output:
	path models

	script:
	models = spaln_gff.getBaseName() + ".training.gff"

	"""
		spaln_add_exons.pl --infile $spaln_gff > $models
	"""

}

// Extract full length models from PASA gff
// Can then be used for training AUGUSTUS
process PasaGffToTraining {

	label 'short_running'

	input:
	path pasa_gff

	output:
	path training_gff

	script:
	training_gff = "transdec.complete.gff3"

	"""
		pasa_select_training_models.pl --nmodels $params.aug_training_models --infile $pasa_gff >> $training_gff
        """
}

// Run one of two training routines for model training
process trainAugustus {

	label 'augustus'

	//scratch true

	publishDir "${params.outdir}/augustus/training/", mode: 'copy'

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
	// Need to use bash for this as the config folder does not yet exist at this point of the process so can't use Groovy...
	aug_folder = "augustus_config/species/${params.aug_species}"
	options = "new_species.pl --species=${params.aug_species}"

	// Re-training and optimization is quite slow and yields minimal gains. Disable if a fast annotation is requested. 
	retrain_options = ""
	if (!params.fast) {
		retrain_options = "optimize_augustus.pl --species=${params.aug_species} ${train_gb} --cpus=${task.cpus} --UTR=off && etraining --species=${params.aug_species} --stopCodonExcludedFromCDS=true ${train_gb}"
	}
		
	"""
		echo $aug_folder > test.txt
		gff2gbSmallDNA.pl $complete_models $genome 1000 $complete_gb
		randomSplit.pl $complete_gb 250
		if [ ! -d $aug_folder ]; then
			$options
		fi
		etraining --species=$params.aug_species --stopCodonExcludedFromCDS=true $train_gb
		$retrain_options
		augustus --stopCodonExcludedFromCDS=true --species=$params.aug_species $test_gb | tee $training_stats
	"""
}

// runs a naked augustus prediction on the genome 
// can be used to find regions of potential interest for e.g. targeted alignments
process runAugustusScan {

	label 'augustus'

        publishDir "${params.outdir}/logs/augustus", mode: 'copy'

	scratch true

        input:
        path genome_chunk
        env AUGUSTUS_CONFIG_PATH
        path aug_extrinsic_config

        output:
        path augustus_result

        script:
        chunk_name = genome_chunk.getName().tokenize("_")[-1]
        augustus_result = "augustus.${chunk_name}.prescan.out.gff"
        config_file = file(params.aug_config)

        """
                augustus --species=${params.aug_species} ${params.aug_options} --softmasking=on --gff3=on --uniqueGeneId=true $genome_chunk > $augustus_result

        """

}

// Runs AUGUSTUS on a genome or chunk thereof
// Exposes several parameters and accepts hints etc
process runAugustus {

	label 'augustus'

	publishDir "${params.outdir}/logs/augustus", mode: 'copy'

	scratch true

	input:
        path genome_chunk
        path hints
        env AUGUSTUS_CONFIG_PATH
	path aug_extrinsic_config

        output:
        path augustus_result

        script:
        chunk_name = genome_chunk.getName().tokenize("_")[-1]
        augustus_result = "augustus.${chunk_name}.out.gff"
	config_file = file(params.aug_config)

	utr = (params.utr) ? "on" : "off"	

        """
		augustus --species=${params.aug_species} ${params.aug_options} --softmasking=on --hintsfile=$hints --gff3=on --UTR=${utr} --extrinsicCfgFile=${config_file} --uniqueGeneId=true $genome_chunk > $augustus_result
 
        """

}

// Runs AUGUSTUS on a genome or chunk thereof
// Speeds up  the search by limiting AUGUSTUS to regions that hold evidences +/- flanks
process runAugustusBatch {

	label 'augustus'

	input:
	path genome_chunk
	path hints
	path regions
	env AUGUSTUS_CONFIG_PATH
        path aug_extrinsic_config

	output:
	path augustus_result

	script:
	chunk_name = genome_chunk.getName().tokenize("_")[-1]
	augustus_result = "augustus.${chunk_name}.out.gff"
	genome_fai = genome_chunk.getName() + ".fai"
	command_file = "commands." + chunk_name + ".txt"
        utr = (params.utr) ? "on" : "off"
	"""
		samtools faidx $genome_chunk
		fastaexplode -f $genome_chunk -d . 		
		augustus_from_regions.pl --genome_fai $genome_fai --model $params.aug_species --utr ${utr} --options '${params.aug_options}' --aug_conf ${params.aug_config} --hints $hints --bed $regions > $command_file
		parallel -j ${task.cpus} < $command_file
		cat *augustus.gff > $augustus_result
		rm *augustus.gff
	"""
}

// Implements are more parallelized version of the Braker parallelization scheme
// Assemblies are pre-split and each chunk is then broken up and parallized
process runAugustusChunks {

	label 'augustus'

        publishDir "${params.outdir}/logs/augustus", mode: 'copy'

	input:
	path genome_chunk
	path hints
	env AUGUSTUS_CONFIG_PATH
	path augustus_extrinsic_config

	output:
	path augustus_result

	script:
	chunk_name = genome_chunk.getName().tokenize("_")[-1]
	augustus_result = "augustus.${chunk_name}.out.gff"
	command_file = "commands." + chunk_name + ".txt"
	utr = (params.utr) ? "on" : "off"
	
	def support = 1.0
	def m = "PE"
	def options = ""
        if (params.evidence) {
                if (!params.proteins && !params.proteins_targeted) {
                        m = "E"
                }
	}

	"""
		samtools faidx $genome_chunk
		fastaexplode -f $genome_chunk -d .
		augustus_from_chunks.pl --chunk_length $params.aug_chunk_length --genome_fai ${genome_chunk}.fai --model $params.aug_species --utr ${utr} --options '${params.aug_options}' --aug_conf ${params.aug_config} --hints $hints > $command_file
		parallel -j ${task.cpus} < $command_file
		for i in \$(ls *.out | sort -n); do echo \$i >> files.txt ; done;
		joingenes -f files.txt -o ${augustus_result}
		rm files.txt
	"""
	
}

process joinAugustusChunks {

	label 'augustus'

        publishDir "${params.outdir}/logs/augustus/chunks", mode: 'copy'

	input:
	path chunks

	output:
	path augustus_models

	script:
	chunk_name = chunks.getBaseName() 
	augustus_models = chunk_name + ".merged.gff"

	"""
		fix_joingenes_gtf.pl < $chunks > $augustus_models
	"""
}

// Merge all the chunk GFF files into one file and make new IDs
process mergeAugustusGff {

	label 'short_running'

	publishDir "${params.outdir}/logs/augustus", mode: 'copy'

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
// folder must be editable!
process prepAugustusConfig {

	label 'augustus'

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
