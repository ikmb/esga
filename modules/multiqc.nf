// MultiQC reporting

process multiqc {


	input:
	path reports

	output:
	path 'multiqc_report.html'

	script:

	"""
		multiqc .
	"""
}
