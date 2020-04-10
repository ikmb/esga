include { fastaSplitChunks; fastaMergeFiles } from "./../fasta" params(params)

workflow repeatmasking_with_lib {

	take:
		genome
		rm_lib
	
	main:
	        fastaSplitChunks(genome,params.nchunks)
		repeatLib(rm_lib)
		repeatMaskLib(fastaSplitChunks.out.flatMap(),repeatLib.out.collect().map{it[0].toString()},rm_lib.collect())
		fastaMergeFiles(repeatMaskLib.out[0].collect())

	emit:
		genome_rm = fastaMergeFiles.out[0]
		genome_rm_gffs = repeatMaskLib.out[1].collectFile()

}

workflow repeatmasking_with_species {

	take:
		genome
		species

	main:
		fastaSplitChunks(genome,params.nchunks)
                repeatLibSpecies(species)
                repeatMaskSpecies(fastaSplitChunks.out.flatMap(),repeatLibSpecies.out.collect().map{it[0].toString()},species)
                fastaMergeFiles(repeatMaskSpecies.out[0].collect())

	emit:
		genome_rm = fastaMergeFiles.out[0]
                genome_rm_gffs = repeatMaskSpecies.out[1].collectFile()

}
// RepeatMasker library needs ot be writable. Need to do this so we can work with locked containers
// Solution: we trigger initial library formatting and then pass the resulting folder as REPEATMASKER_LIB_DIR
process repeatLib {

	label 'short_running'

	input:
	path repeat_lib

	output:
	path "Library"

	script:
	
	"""
		cp ${baseDir}/assets/repeatmasker/my_genome.fa .
		cp ${baseDir}/assets/repeatmasker/repeats.fa .
		mkdir -p Library
                cp ${baseDir}/assets/repeatmasker/DfamConsensus.embl Library/
                gunzip -c ${baseDir}/assets/repeatmasker/taxonomy.dat.gz > Library/taxonomy.dat
		export REPEATMASKER_LIB_DIR=\$PWD/Library
		RepeatMasker -lib $repeat_lib my_genome.fa > out
	"""	
}

process repeatLibSpecies {

	input:
	val species

	output:
	path "Library"

	script:

	"""
		cp ${baseDir}/assets/repeatmasker/my_genome.fa .
                cp ${baseDir}/assets/repeatmasker/repeats.fa .

		mkdir -p Library
		cp ${baseDir}/assets/repeatmasker/DfamConsensus.embl Library/
		gunzip -c ${baseDir}/assets/repeatmasker/taxonomy.dat.gz > Library/taxonomy.dat
                export REPEATMASKER_LIB_DIR=\$PWD/Library

		RepeatMasker -species $species my_genome.fa > out

	"""
	
}

// generate a soft-masked sequence for each assembly chunk
// if nothing was masked, return the original genome sequence instead and an empty gff file. 
process repeatMaskLib {

	scratch true

	input: 
	path genome
	env REPEATMASKER_LIB_DIR
	path repeats

	output:
	path genome_rm
	path rm_gff

	script:

	def options = ""
	options = "-lib $repeats"
	base_name = genome.getName()
	genome_rm = base_name + ".masked"
	rm_gff = base_name + ".out.gff"
	rm_tbl = base_name + ".tbl"
	rm_out = base_name + ".out"
	
	"""
		echo \$REPEATMASKER_LIB_DIR > lib_dir.txt
		RepeatMasker $options -gff -xsmall -q -pa ${task.cpus} $genome
		test -f ${genome_rm} || cp $genome $genome_rm && touch $rm_gff
	"""
}


process repeatMaskSpecies {

        scratch true

        input:
        path genome
        env REPEATMASKER_LIB_DIR
        val species

        output:
        path genome_rm
        path rm_gff

        script:

        def options = ""
        options = "-species $species"
        base_name = genome.getName()
        genome_rm = base_name + ".masked"
        rm_gff = base_name + ".out.gff"
        rm_tbl = base_name + ".tbl"
        rm_out = base_name + ".out"

        """
                echo \$REPEATMASKER_LIB_DIR > lib_dir.txt
                RepeatMasker $options -gff -xsmall -q -pa ${task.cpus} $genome
                test -f ${genome_rm} || cp $genome $genome_rm && touch $rm_gff
        """
}

