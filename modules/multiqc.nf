// MultiQC reporting

workflow multiqc {

	take:
		software_versions

	main:
		report(software_versions)

	emit:

		html = report.out

}

process report {

	label 'multiqc'

	publishDir "${params.outdir}/MultiQC", mode: 'copy'

	input:
	path reports

	output:
	path 'multiqc_report.html'

	script:

	"""
		multiqc .
	"""
}
