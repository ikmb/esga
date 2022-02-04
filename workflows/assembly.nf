// ****************
// Module for assembly processing and validation
// ****************

include { fastaCleanNames } from "./../modules/fasta" addParams(results: "${params.outdir}/logs/fasta")
include { AssemblyStats ; AssemblyFilterSize } from "./../modules/assembly/main.nf" params(params)

// Clean up the assembly and generate some basic stats
workflow ASSEMBLY_PREPROCESS {

	take:
		genome

	main:

		AssemblyStats(genome)
		AssemblyFilterSize(genome,params.min_contig_size)
		fastaCleanNames(AssemblyFilterSize.out)

	emit:
		stats = AssemblyStats.out[0]
		fasta = fastaCleanNames.out[0]
}

