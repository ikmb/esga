include { gtf2hints; kraken2gff } from "./../gtf" params(params)
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
		gtf2hints(map_gtf.out[0].collectFile())
		kraken2gff(map_gtf.out[0])

	emit:
		mapped_gtf = map_gtf.out[0].collect()
		mapped_gff = kraken2gff.out.collect()
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
		rm -Rf align
	"""
}

// run Kraken to map the annotation from the reference to our genome of interest
process map_gtf {

	publishDir "${params.outdir}/logs/satsuma", mode: 'copy'

	label 'satsuma'

	input:
	tuple val(ref_name), path(ref_genome),path(ref_gtf),path(satsuma_chunks)
	path query_genome
	
	output:
	path mapped_gtf
	path chain_file

	script:
	chain_file = "satsuma_summary.chained." + ref_name + ".out"
	mapped_gtf = ref_gtf.getBaseName() + "." + query_genome.getBaseName() + ".mapped.gtf"

	"""
		cat $satsuma_chunks > $chain_file
		kraken_build_config.pl --ref_fa $ref_genome --query_fa $query_genome --chain $chain_file > kraken.config
		RunKraken -c kraken.config -T QUERY -S REF -s $ref_gtf -o $mapped_gtf -f gene,transcript,mRNA,CDS,exon
	"""

}


