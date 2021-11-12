// Repeatmasking a genome assembly

include { fastaSplitSize; fastaMergeChunks } from "./../fasta" params(params)

// Repeatmasking using a library in Fasta format
workflow repeatmasking_with_lib {

	take:
		genome
		rm_lib

	main:
		assembly_name = Channel.value("genome.rm.fa")

		fastaSplitSize(genome,params.npart_size)
		repeatLib(rm_lib)
		repeatMaskLib(fastaSplitSize.out.flatMap(),repeatLib.out.collect().map{it[0].toString()},rm_lib.collect())
		fastaMergeChunks(repeatMaskLib.out[0].collect(),assembly_name)
		repeats_to_hints(repeatMaskLib.out[1].collect())
		repeats_to_gmod(repeatMaskLib.out[1].collect())
	emit:
		genome_rm = fastaMergeChunks.out[0]
		genome_rm_gffs = repeatMaskLib.out[1].collectFile()
		genome_rm_hints = repeats_to_hints.out

}

//Repeatmasking providing a defined taxonomic group (use with care!)
workflow repeatmasking_with_species {

	take:
		genome
		species

	main:
                assembly_name = Channel.value("genome.rm.fa")

		fastaSplitSize(genome,params.npart_size)
		trigger_library(species)
                repeatMaskSpecies(fastaSplitSize.out.flatMap(),trigger_library.out.collect())
                fastaMergeChunks(repeatMaskSpecies.out[0].collect(),assembly_name)
		repeats_to_hints(repeatMaskSpecies.out[1].collect())
		repeats_to_gmod(repeatMaskLib.out[1].collect())
	emit:
		genome_rm = fastaMergeChunks.out[0]
                genome_rm_gffs = repeatMaskSpecies.out[1].collectFile()
		genome_rm_hints = repeats_to_hints.out

}

// RepeatMasker library needs ot be writable. Need to do this so we can work with locked containers
// Solution: we trigger initial library formatting and then pass the resulting folder as REPEATMASKER_LIB_DIR
process repeatLib {

	label 'repeatmasker'

	publishDir "${params.outdir}/logs/repeatmasker", mode: 'copy'

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

process trigger_library {

	label 'repeatmasker'

	input:
	val species

	output:
	val species

	script:

	"""
		cp ${baseDir}/assets/repeatmasker/my_genome.fa .
		RepeatMasker -species $species my_genome.fa > out
	"""
}

// trigger one-time repeat database formatting before running parallel masking steps
process repeatLibSpecies {

	label 'repeatmasker'

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

	label 'repeatmasker'

	publishDir "${params.outdir}/repeatmasker", mode: 'copy'

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
		RepeatMasker $options -gff -xsmall -q -nolow -pa ${task.cpus} $genome
		test -f ${genome_rm} || cp $genome $genome_rm && touch $rm_gff
	"""
}


process repeatMaskSpecies {

        scratch true

	label 'repeatmasker'

	publishDir "${params.outdir}/repeatmasker", mode: 'copy'	

        input:
        path genome
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
                RepeatMasker $options -gff -xsmall -nolow -q -pa ${task.cpus} $genome
                test -f ${genome_rm} || cp $genome $genome_rm && touch $rm_gff
        """
}

// Convert repeat annotations to AUGUSTUS hints
process repeats_to_hints {

	publishDir "${params.outdir}/logs/repeatmasker", mode: 'copy'

	scratch true

	input:
	path gffs

	output:
	path hints

	script:
	hints = "hints.repeats.gff"

	"""
		cat $gffs > merged.gff
		gff_repeats2hints.pl --infile merged.gff > $hints
	"""

}

process repeats_to_gmod {

	publishDir "${params.outdir}/gmod", mode: 'copy'

	input:
	path gffs

	output:
	path gff

	script:
	
	"""
		cat $gffs > $gff
	"""

}
