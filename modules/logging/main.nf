// Module to log relevant information about the pipeline

workflow get_software_versions {

	main:
		base_container_versions()
		combine_versions(base_container_versions.out.collect()) 

	emit:
		yaml = combine_versions.out
		
}

process base_container_versions {

	publishDir "${params.outdir}/logs/versions", mode: 'copy'

	output:
	path "*.txt"

	script:

	"""
		minimap2 --version &> v_minimap2.txt || true
		RepeatMasker -v &> v_repeatmasker.txt || true
		augustus | head -n1 &> v_augustus.txt || true
                spaln &> v_spaln.txt || true

	"""
}

process trinity_version {

	output:
	path "v_trinity.txt"

	script:

	"""
		Trinity --version | grep "Trinity version" > v_trinity.txt
	"""
	
}

process fastp_version {

	output:
	path "v_fastp.txt"

	script:

	"""
		fastp -v 2> v_fastp.txt
	"""

}

process combine_versions {

	input:
	path versions

	output:
	path yaml

	script:
	yaml = "software_versions_mqc.yaml"

	"""
		parse_versions.pl > $yaml
	"""
}
