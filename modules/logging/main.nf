// Module to log relevant information about the pipeline

workflow get_software_versions {

	main:
		base_container_versions
		trinity_version
		fastp_version		
		combine_versions(base_container_versions.out.concat(trinity_version,fastp_version) )

	emit:

		versions = combine_versions.out
		

}

process base_container_versions {

	output:
	path "v*.txt"

	script:

	"""
		exonerate -version > v_exonerate.txt
		minimap2 --version > v_minimap2.txt
		RepeatMasker -v > v_repeatmasker.txt
		spaln 2> v_spaln.txt
		augustus 2> v_augustus.txt
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

process satsuma_version {

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
