process gtf2hints {

	input:
	path gtf

	output:
	path hints
	
	script:
	hints = gtf.getBaseName() + ".hints.gff"

	"""
		gtf2hints.pl --gtf $gtf --pri $params.pri_trans --source T > $hints
	"""
}

// The sed command here is to fix an issue with malformed kraken output
process kraken2gff {

	publishDir "${params.outdir}/logs/kraken", mode: 'copy'

	label 'gaas'

	input:
	path gtf

	output:
	path gff

	script:
	gff = gtf.getBaseName() + ".gff3"

	"""
		sed -i.bak 's/;\"/\"/g' $gtf
		sed -i.bak2 's/\t\t/\tensembl\t/' $gtf
		kraken2gff.pl --infile $gtf > $gff

		rm *.bak*
	"""

}

