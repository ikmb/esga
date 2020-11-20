include { gtf2hints } from "./../gtf" params(params)
include { fastaSplitSize } from "./../fasta" params(params)

// note that ref below is [ name, genome, gtf ]
workflow map_annotation {

	take:
		ref
		query_genome

	main:
		fastaSplitSize(query_genome,params.npart_size)
		fastaCleanNames(ref)
		align_genomes(fastaCleanNames.out.combine(fastaSplitSize.out.flatMap()) )
		
		grouped_chains = align_genomes.out.groupTuple()

		grouped_input = fastaCleanNames.out.join(grouped_chains)

		map_gtf(grouped_input,query_genome)
		gtf2hints(map_gtf.out.collectFile())

	emit:
		mapped_gtf = map_gtf.out
		hints = gtf2hints.out

}

process fastaCleanNames {

        input:
        tuple val(name),path(fasta),path(gtf)

        output:
        tuple val(name),path(fasta_clean),path(gtf)

        script:
        fasta_clean = fasta.getBaseName() + ".clean.fa"

        """
                sed 's/ .*//' $fasta > $fasta_clean
        """

}

// run Satsuma2, redirect output to null because Satsuma2 is very chatty
process align_genomes {

	label 'satsuma'
	tag "${ref_genome}"

	input:
	tuple val(ref_name),path(ref_genome),path(gtf),path(query_genome)

	output:
	tuple val(ref_name),path(satsuma_chain_chunk)

	script:
	query_chunk = query_genome.getBaseName()
	satsuma_chain_chunk = "satsuma_summary.chained.out_" + query_chunk

	"""
		SatsumaSynteny2 -q $query_genome -t $ref_genome -threads ${task.cpus} -o align 2>&1 >/dev/null
		cp align/satsuma_summary.chained.out $satsuma_chain_chunk
	"""
}

process map_gtf {


	label 'satsuma'

	input:
	tuple val(ref_name), path(ref_genome),path(ref_gtf),path(satsuma_chunks)
	path query_genome

	output:
	path mapped_gtf

	script:

	mapped_gtf = ref_gtf.getBaseName() + "." + query_genome.getBaseName() + ".mapped.gtf"

	"""
		cat $satsuma_chunks > satsuma_summary.chained.out
		kraken_build_config.pl --ref_fa $ref_genome --query_fa $query_genome --chain satsuma_summary.chained.out > kraken.config
		RunKraken -c kraken.config -T QUERY -S REF -s $ref_gtf -o $mapped_gtf
	"""

}


