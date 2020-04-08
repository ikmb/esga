// ******************
// Generic functions for FASTA processing
// ******************

// split a fasta file into a number of chunks
// will delete all empty chunks since no checking of the absolute minimum number of chunks is performed.
process fastaSplitChunks {

        publishDir "${params.outdir}/fasta/chunks"

        input:
        path fasta
	val nchunks

        output:
        path "*_chunk_*"

        script:
        ref_name = fasta.getBaseName() + "_chunk_%3.3d"

        """
                fastasplitn -in $fasta -n $nchunks -t $ref_name
		find \$PWD -size -10c -print -delete
        """
}

process assemblySplit {

	publishDir "${params.outdir}/databases/genome/", mode: 'copy'

	label 'short_running'

	input:
	path genome
	val chunk_size

	output:
	path genome_chunks
	path genome_agp

	script:

	genome_chunks = genome + ".chunk"
	genome_agp = genome + ".agp"

	"""
		chromosome_chunk.pl -fasta_file $genome -size $chunk_size
	"""
}
