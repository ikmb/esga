include { gtf2hints } from "./../gtf" params(params)
include { fastaSplitSize } from "./../fasta" params(params)

// note that ref below is [ genome, gtf ]
workflow map_annotation {

	take:
		ref
		query_genome

	main:
		fastaSplitSize(query_genome,params.npart_size)
		align_genomes(ref.collect(),fastaSplitSize.out.flatMap())
		map_gtf(align_genomes.out.collectFile(),query_genome,ref)
		gtf2hints(map_gtf.out.collectFile())

	emit:
		mapped_gtf = map_gtf.out
		hints = gtf2hints.out

}

// run Satsuma2, redirect output to null because Satsuma2 is very chatty
process align_genomes {

	label 'satsuma'
	tag "${ref_genome}"

	input:
	tuple path(ref_genome),path(gtf)
	path query_genome

	output:
	path satsuma_chain

	script:
	satsuma_chain = "satsuma_summary.chained.out"

	"""
		SatsumaSynteny2 -q $query_genome -t $ref_genome -threads ${task.cpus} -o align 2>&1 >/dev/null
		mv align/satsuma_summary.chained.out satsuma_summary.chained.out
	"""
}

process map_gtf {

	label 'satsuma'

	input:
	path satsuma_chain
	path query_genome
	tuple path(ref_genome),path(ref_gtf)

	output:
	path mapped_gtf

	script:

	mapped_gtf = ref_gtf.getBaseName() + "." + ref_genome.getBaseName() + ".mapped.gtf"

	"""
		kraken_build_config.pl --ref_fa $ref_genome --query_fa $query_genome --chain $satsuma_chain > kraken.config
		runKraken -c kraken.config -T QUERY -S REF -s $ref_gtf
	"""

}


