// ******************
// Generic functions for FASTA processing
// ******************

// split a fasta file into a number of chunks
// will delete all empty chunks since no checking of the absolute minimum number of chunks is performed.
process fastaSplitChunks {

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

process fastaSplitSize {

	input:
	path fasta
	val fsize

	output:
	path "*.part-*.*"

	script:

	"""
		fasta-splitter.pl -part-sequence-size $fsize $fasta
	"""

}

process assemblySplit {

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

process fastaMergeChunks {

        label 'short_running'

	//publishDir "${params.outdir}/fasta", mode: 'copy'

        input:
        path chunks
	val assembly_name

        output:
        path merged_chunks

        script:

	merged_chunks = assembly_name

        """
                cat $chunks >> merged.fa
                fastasort -f merged.fa > $merged_chunks
                rm merged.fa
        """
}

process fastaMergeFiles {

	label 'short_running'

        input:
        path chunks

        output:
        path merged_fasta

        script:

        merged_fasta = "sequences.merged.fa"

        """
                cat $chunks >> merged.fa
                fastasort -f merged.fa > $merged_fasta
                rm merged.fa
        """

}

process fastaRemoveShort {

	input:
	path fasta
	val min_len

	output:
	path fasta_filtered

	script:

	fasta_filtered = fasta.getBaseName() + "." + min_len + "." + fasta.getExtension()

	"""
		gaas_fasta_filter_by_size.pl -f $fasta -s $min_len -o $fasta_filtered
	"""

}

process fastaCleanProteins {

	input:
	path fasta

	output:
	path fasta_clean

	script:

	fasta_clean = fasta.getBaseName() + ".clean.fa"

	"""

		gaas_fasta_cleaner.pl -f $fasta > tmp
		fastaclean -f tmp -p | sed 's/:filter(clean)//' > $fasta_clean
		rm tmp
	"""


}
