include fastaSplitChunks from "./../fasta" params(params)

workflow repeatmasking {

	take:
		genome
		rm_lib
	
	main:
	        fastaSplitChunks(genome,params.nchunks)
		repeatLib(rm_lib)
		repeatMask(fastaSplitChunks.out.flatMap(),repeatLib.out.map{it.toString()},rm_lib)
		repeatMerge(repeatMask.out[0].collect())

	emit:
		genome_rm = repeatMerge.out[0]
		genome_rm_gffs = repeatMask.out[1].collectFile()

}

// RepeatMasker library needs ot be writable. Need to do this so we can work with locked containers
// Solution: we trigger initial library formatting and then pass the resulting folder as REPEATMASKER_LIB_DIR
process repeatLib {

	label 'short_running'

	publishDir "${params.outdir}/repeatmasker/", mode: 'copy'

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

// generate a soft-masked sequence for each assembly chunk
// if nothing was masked, return the original genome sequence instead and an empty gff file. 
process repeatMask {

	//scratch true

	publishDir "${params.outdir}/repeatmasker/chunks"

	input: 
	path(genome)
	env(REPEATMASKER_LIB_DIR)
	path(repeats)

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

// Merge the repeat-masked assembly chunks
process repeatMerge {

	label 'short_running'

        publishDir "${params.outdir}/repeatmasker", mode: 'copy'

	input:
	path genome_chunks

	output:
	path masked_genome

	script:
	
	masked_genome = "assembly.rm.fa"
	masked_genome_index = masked_genome + ".fai"

	"""
		cat $genome_chunks >> merged.fa
		fastasort -f merged.fa > $masked_genome
		samtools faidx $masked_genome
		rm merged.fa
	"""	
}
