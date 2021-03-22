
include { fastaCleanNames } from "./../fasta" addParams(results: "${params.outdir}/logs/fasta")

// Clean up the assembly and generate some basic stats
workflow assembly_preprocessing {

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

// Get basic stats for the assembly
process AssemblyStats {

	//publishDir "${params.outdir}/assembly", mode: 'copy'

	label 'gaas'

        input:
        path fasta

        output:
        path stats_dir

        script:
	stats_dir = "stats"

        """
                gaas_fasta_statistics.pl -f $fasta -o $stats_dir
        """

}

// Remove contigs below a certain length
process AssemblyFilterSize {

        //publishDir "${params.outdir}/assembly", mode: 'copy'

	label 'gaas'

        input:
        path fasta
        val min_size

        output:
        path fasta_filtered

        script:
        fasta_filtered = fasta.getBaseName() + ".filtered.fa"

        """
                gaas_fasta_filter_by_size.pl -f $fasta -s $min_size -o $fasta_filtered
        """

}

