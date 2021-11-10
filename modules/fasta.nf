// ******************
// Generic functions for FASTA processing
// ******************

process fastaSplitSize {

	label 'medium_running'

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

	publishDir "${params.results}", mode: 'copy'
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

	publishDir "${params.outdir}/logs/fasta", mode: 'copy'

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
                cat $chunks >> $merged_fasta
                #fastasort -f merged.fa > $merged_fasta
                #rm merged.fa
        """

}

process fastaRemoveShort {

	label 'gaas'
	
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

process fastaGaasClean {

	label 'gaas'

	input:
	path fasta

	output:
	path fasta_clean

	script:

	fasta_clean = fasta.getBaseName() + ".gaas_cleaned.fa"

	"""
		sed 's/[.]\$//' $fasta > cleaned.fa
                gaas_fasta_cleaner.pl -f cleaned.fa -o $fasta_clean
		rm cleaned.fa
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
		fastaclean -f $fasta -p | sed 's/:filter(clean)//' | sed 's/ pep .*//' > $fasta_clean
	"""
}

process fastaCleanNames {

	input:
	path(fasta)
	
	output:
	path(fasta_clean)

	script:
	fasta_clean = fasta.getBaseName() + ".clean.fa"

	"""
		sed 's/ .*//' $fasta > $fasta_clean		
	"""

}
