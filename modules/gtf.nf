process gtf2hints {

	input:
	path gtf

	output:
	path hints

	script:
	hints = gtf.getBaseName() + ".hints.gff"

	"""
		gtf2hints.pl --gtf $gtf --pri 4 --source T > $hints
	"""
}
