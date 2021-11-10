// Module to log relevant information about the pipeline

workflow get_software_versions {

	take:
		genome
	main:
		base_container_versions()
		trinity_version()
		augustus_version()
		fastp_version()
		minimap2_version()
		versions = base_container_versions.out.collect().concat( trinity_version.out.collect() , augustus_version.out.collect() , fastp_version.out.collect() , minimap2_version.out.collect() )
		combine_versions( versions.collect() )

	emit:
		yaml = combine_versions.out

}

process base_container_versions {

	publishDir "${params.outdir}/logs/versions", mode: 'copy'

	output:
	path "*.txt"

	script:

	"""
		RepeatMasker -v &> v_repeatmasker.txt || true
                spaln &> v_spaln.txt || true
		samtools --version | head -n 1 &> v_samtools.txt || true
		gffread --version &> v_gffread.txt || true
		exonerate --version | head -n1 > v_exonerate.txt || true 
		gth -version | head -n1 > v_gth.txt || true 
	"""
}

process augustus_version {

        publishDir "${params.outdir}/logs/versions", mode: 'copy'

	label 'augustus'

	output:
	path "v_augustus.txt"

	script:

	"""
		augustus | head -n1 &> v_augustus.txt || true 
	"""

}

process trinity_version {

        publishDir "${params.outdir}/logs/versions", mode: 'copy'

	label 'trinity'

	output:
	path "v_trinity.txt"

	script:

	"""
		Trinity --version | grep "Trinity version" &> v_trinity.txt || true
	"""
	
}

process fastp_version {

        publishDir "${params.outdir}/logs/versions", mode: 'copy'

	label 'fastp'

	output:
	path "v_fastp.txt"

	script:

	"""
		fastp -v 2> v_fastp.txt || true
	"""

}

process minimap2_version {

        publishDir "${params.outdir}/logs/versions", mode: 'copy'

	label 'minimap'

	output:
	path "v_minimap2.txt"

	script:
	
	"""
	 	minimap2 --version &> v_minimap2.txt || true
	"""
}


process combine_versions {

	publishDir "${params.outdir}/MultiQC", mode: 'copy'

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
