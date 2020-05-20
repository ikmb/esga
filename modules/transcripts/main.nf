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

	publishDir "${params.outdir}/transcripts", mode: 'copy'

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

process blastTranscripts {


	input:
	path est_chunk
	path blastdb_files
		
	output:
	path blast_report

	script:
	db_name = blastdb_files[0].baseName
	chunk_name = est_chunk.getName().tokenize('.')[-2]
	blast_report = "${est_chunk.baseName}.${db_name}.est.blast"

	"""
		blastn -db $db_name -evalue $params.blast_evalue -query $est_chunk -outfmt "${params.blast_options}" -num_threads ${task.cpus} > $blast_report
	"""
}

process blastToTargets {

	input:
	path blast_reports

	output:
	path targets

	script:
	targets = "EST.blast.targets.txt"

	"""
		cat $blast_reports >> merged.out
		blast2exonerate_targets.pl --infile merged.out --max_intron_size $params.max_intron_size > $targets
	"""
}

// from a protein database and genome sequence to run exonerate
process transcriptExonerateBatch {

        scratch true

	publishDir "${params.outdir}/logs/transcripts/", mode: 'copy'

        input:
        path hits_chunk
        path protein_db
        path protein_db_index
        path genome

        output:
        path exonerate_chunk

        script:
        genome_faidx = genome.getName() + ".fai"
        query_tag = protein_db.baseName
        chunk_name = hits_chunk.getName().tokenize('.')[-2]
        commands = "commands." + chunk_name + ".txt"
        exonerate_chunk = "${hits_chunk.baseName}.${query_tag}.exonerate.out"

        // get the protein fasta sequences, produce the exonerate command and genomic target interval fasta, run the whole thing,
        // merge it all down to one file and translate back to genomic coordinates
        // remove all the untracked intermediate files

        """
                samtools faidx $genome
                extractMatchTargetsFromIndex.pl --matches $hits_chunk --db $protein_db_index
                exonerate_from_blast_hits.pl --matches $hits_chunk --assembly_index $genome --max_intron_size $params.max_intron_size --query_index $protein_db_index --analysis protein2genome --outfile $commands
                parallel -j ${task.cpus} < $commands
                cat *.exonerate.align | grep -v '#' | grep 'exonerate:protein2genome:local' >> merged.${chunk_name}.exonerate.out 2>/dev/null
                exonerate_offset2genomic.pl --infile merged.${chunk_name}.exonerate.out --outfile $exonerate_chunk
                test -f $exonerate_chunk || cat "#" > $exonerate_chunk
		rm *.align
		rm *._target_.fa*
		rm *._query_.fa*
        """
}

process transcriptExonerate2Hints {

	input:
	path exonerate_results
	
	output:
	path exonerate_hints
	
	script:
	exonerate_hints = "ESTs.exonerate.hints.gff"
			
	"""
		cat $exonerate_results >> merged.out
		exonerate2gff.pl --infile merged.out --source est --outfile $exonerate_hints
	"""
}
