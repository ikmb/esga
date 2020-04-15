include { AssemblyStats; AssemblyFilterSize  } from "./../fasta" params(params)

workflow assembly_preprocessing {

	take:
		genome

	main:

		AssemblyStats(genome)
		AssemblyFilterSize(genome,params.min_contig_size)

	emit:
		stats = AssemblyStats.out[0]
		fasta = AssemblyFilterSize.out[0]
}
