// **************************
// Production of annotation hints from EST/transcript FASTA sequences
// **************************

workflow esthint {

	take:
		genome_rm
		est

	main:
		estMinimap(est,genome_rm)
		estMinimapToHints(estMinimap.out)
	emit:
		gff = estMinimap.out
		hints = estMinimapToHints.out
}

process estMinimap {

	scratch true

	input:
	path est
	path genome_rm	

	output:
	path minimap_gff	

	script:
	minimap_gff = est.getBaseName() + ".minimap.gff"
	minimap_bam = est.getBaseName() + ".minimap.bam"

	"""
		samtools faidx $genome_rm
		minimap2 -t ${task.cpus} -ax splice:hq -c -G ${params.max_intron_size}  $genome_rm $est | samtools sort -O BAM -o $minimap_bam
		minimap2_bam2gff.pl $minimap_bam > $minimap_gff
	"""	

}

// Combine exonerate hits and generate hints
process estMinimapToHints {

	label 'short_running'

	input:
	path minimap_gff
	
	output:
	path minimap_hints
	
	script:
	minimap_hints = minimap_gff.getBaseName() + ".hints.gff"
			
	"""
		minimap2hints.pl --source est2genome --pri ${params.pri_est} --infile $minimap_gff --outfile $minimap_hints
	"""
}

