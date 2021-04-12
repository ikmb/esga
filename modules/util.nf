// this is currently super simplistic and should be made a bit smarter down the line
process merge_hints {

	input:
	path all_hints

	output:
	path hints

	script:
	hints = "hints.merged.gff"

	"""
		cat $all_hints > merged.txt
		cat merged.txt | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -k 3,3 | sort -s -k 1,1 | join_mult_hints.pl > $hints
		rm merged.txt
	"""

}

// Convert a set of hints into hint-covered regions in BED format
process prepHintsToBed {

        label 'short_running'

        input:
        path hints

        output:
	path hints
        path bed

        script:

        bed = "regions.bed"

        """
                grep -v "#" $hints | grep -v "false" | sort -k1,1 -k4,4n -k5,5n -t\$'\t' > hints.sorted
                gff_by_strand.pl --infile hints.sorted
                gff2clusters.pl --infile hints.plus.gff --max_intron $params.max_intron_size > $bed
                gff2clusters.pl --infile hints.minus.gff --max_intron $params.max_intron_size >> $bed
        """
}

// Use StringTies gffread tool to extract protein sequences from genomes and gtf files
process GffToFasta {

	label 'short_running'

	input:
	path gff
	path genome

	output:
	path fasta

	script:
	fasta = gff.getBaseName() + ".proteins.fasta"

	"""
		gffread $gff -g $genome -y $fasta
	"""
	
}

