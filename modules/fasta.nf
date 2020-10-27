// ******************
// Generic functions for FASTA processing
// ******************


process fastaToBlastnDBMasked {

        input:
        path genome_fa

        output:
        path "${dbName}*.n*"

        script:
        dbName = genome_fa.getBaseName()
        mask = dbName + ".asnb"
        """

                ${params.makeblastdb} -in $genome_fa -parse_seqids -dbtype nucl -out $dbName

                convert2blastmask -in $genome_fa -parse_seqids -masking_algorithm repeat \
                        -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin \
                        -out $mask

                ${params.makeblastdb} -in $genome_fa -parse_seqids -dbtype nucl -mask_data $mask -out $dbName
        """
}

process fastaToBlastnDB {

        input:
        path genome_fa

        output:
        path "${dbName}*.n*"

        script:
        dbName = genome_fa.getBaseName()

        """
                ${params.makeblastdb} -in $genome_fa -parse_seqids -dbtype nucl -out $dbName
        """
}

// Get Accession numbers from FASTA file

process fastaToList {

	label 'short_running'

	input:
	path fasta

	output:
	path accessions

	script:
	accessions = fasta.getBaseName() + ".accs.txt"

	"""
		grep "^>.*" $fasta | sed 's/^>//'  > $accessions
	"""

}

// create a cdbtools compatible  index
// we need this to do very focused exonerate searches later
process fastaToCdbindex {

        label 'short_running'

        input:
        path fasta

        output:
        path fasta_index

        script:
        fasta_index = fasta.getName()+ ".cidx"

        """
                cdbfasta $fasta
        """
}

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

		gaas_fasta_cleaner.pl -f $fasta -o tmp
		fastaclean -f tmp -p | sed 's/:filter(clean)//' | sed 's/ pep .*//' > $fasta_clean
		rm tmp
	"""

}
