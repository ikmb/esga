// De-novo transcriptome assembly using Trinity

workflow trinity_guided_assembly {
	
	take:
		bam

	main:
		runTrinityGuided(bam)				

	emit:
		assembly = runTrinityGuided.out
}


process runTrinityGuided {
	
	label 'trinity'

	scratch true 

	//publishDir "${params.outdir}/transcripts/trinity", mode: 'copy'
	
	input:
	path bam	

	output:
	path "transcriptome_trinity/Trinity-GG.fasta" 
	
	script:

	trinity_option = ( params.rnaseq_stranded == true ) ? "--SS_lib_type RF" : ""

	"""
		Trinity --genome_guided_bam $bam \
		--genome_guided_max_intron ${params.max_intron_size} \
		--CPU ${task.cpus} \
		--max_memory ${task.memory.toGiga()-1}G \
		--output transcriptome_trinity \
		$trinity_option
	"""
}
