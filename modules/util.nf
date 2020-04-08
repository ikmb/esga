process merge_hints {


	input:
	path protein_hints
	path est_hints

	output:
	path hints

	script:
	hints = "hints.merged.gff"

	"""

		cat $protein_hints $est_hints >> $hints
	"""

}

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

