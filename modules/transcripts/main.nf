// **************************
// Production of annotation hints from EST/transcript FASTA sequences
// **************************
include { fastaToBlastnDBMasked; fastaToCdbindex } from "./../fasta" params(params)

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

workflow esthint_slow {

	take:
		genome_rm
		est

	main:
		fastaToBlastnDBMasked(genome_rm)
		fastaToCdbindex(est)
		estBlastNMasked( est.splitFasta(by: params.nblast, file: true), fastaToBlastnDBMasked.out.collect() )
		blastnToTargets(estBlastNMasked.out.collect())
		estExonerate(blastnToTargets.out.splitText(by: params.nexonerate, file:true),genome_rm.collect(),fastaToCdbindex.out.collect(),est.collect())
		estExonerateToHints( estExonerate.out.collect() )
	emit:
		gff = estExonerate.out
		hints = estExonerateToHints.out

}

process estBlastN {

	input:
	path est_chunk
	path blast_db_files

	output:
	path blast_result

	script:
        db_name = blast_db_files[0].getBaseName()
        chunk_name = est_chunk.getName().tokenize('.')[-2]
        blast_result = "${est_chunk.baseName}.blast"

        """
                blastn -num_threads ${task.cpus} -evalue ${params.blast_evalue} -outfmt \"${params.blast_options}\" -db $db_name -query $est_chunk > $blast_result
        """

}

process estBlastNMasked {

	input:
        path est_chunk
        path blast_db_files

        output:
        path blast_result

        script:
        db_name = blast_db_files[0].getBaseName()
        chunk_name = est_chunk.getName().tokenize('.')[-2]
        blast_result = "${est_chunk.baseName}.blast"

        """
                blastn -num_threads ${task.cpus} -db_soft_mask 40 -evalue ${params.blast_evalue} -outfmt \"${params.blast_options}\" -db $db_name -query $est_chunk > $blast_result
        """
}

process blastnToTargets {

	input:
	path reports

	output:
	path targets

	script:

	targets = reports[0].getBaseName() + ".targets.txt"

	"""
		cat $reports > reports.txt
		tblastn2exonerate_targets.pl --infile reports.txt --max_intron_size ${params.max_intron_size} > $targets
	"""

}

process estExonerate {

	input:
	path hits_chunk
	path genome
	path est_db_index
	path est

	output:
	path exonerate_report

	script:
	commands = "commands.txt"
	exonerate_report = hits_chunk.getBaseName() + ".est.exonerate.out"
	chunk_name = hits_chunk.getBaseName()
	
	"""
		samtools faidx $genome
                extractMatchTargetsFromIndex.pl --matches $hits_chunk --db $est_db_index
                exonerate_from_blast_hits.pl --matches $hits_chunk --assembly_index $genome --max_intron_size $params.max_intron_size  --analysis est2genome --outfile $commands
                parallel -j ${task.cpus} < $commands
                cat *.exonerate.align | grep -v '#' | grep 'exonerate:est2genome' >> merged.${chunk_name}.exonerate.out 2>grep.log
                exonerate_offset2genomic.pl --infile merged.${chunk_name}.exonerate.out --outfile $exonerate_report
                [ -s $exonerate_report ] || cat "#" > $exonerate_report

                rm *.align
                rm *._target_.fa*
                rm *._query_.fa*

	"""
}

// merge the exonerate hits and create the hints
process estExonerateToHints {

        label 'medium_running'

        publishDir "${params.outdir}/logs/exonerate", mode: 'copy'

        input:
        path chunks

        output:
        path exonerate_gff

        script:
        exonerate_gff = "transcripts.exonerate.hints.gff"

        """
                cat $chunks > all_chunks.out
                exonerate2gff.pl --infile all_chunks.out --pri ${params.pri_est} --source est --outfile $exonerate_gff
        """
}

process estMinimap {

	//scratch true

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
		minimap2hints.pl --src $params.t_est --source est2genome --pri ${params.pri_est} --infile $minimap_gff --outfile $minimap_hints
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
