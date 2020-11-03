// **************************
// Production of annotation hints from protein FASTA sequences
// **************************

include { fastaToBlastnDBMasked; fastaToCdbindex; fastaCleanProteins; fastaRemoveShort; assemblySplit } from "./../fasta" params(params)

workflow proteinhint_slow {

	take:
                genome_rm
                protein_fa

        main:
                fastaToBlastnDBMasked(genome_rm)
		fastaCleanProteins(protein_fa)
		fastaRemoveShort(fastaCleanProteins.out,params.min_prot_length)
                fastaToCdbindex(fastaRemoveShort.out)
		tblastnMasked(fastaRemoveShort.out.splitFasta(by: params.nblast, file: true) , fastaToBlastnDBMasked.out)
		tblastnToTargets(tblastnMasked.out[0].collect())
		TargetsFindMissing(tblastnToTargets.out,fastaRemoveShort.out)
		protExonerateFromList(TargetsFindMissing.out.splitText(by: params.nexonerate_exhaustive, file: true),fastaRemoveShort.out.collect(),fastaToCdbindex.out.collect(),genome_rm.collect())	
                protExonerateBatch(tblastnToTargets.out.splitText(by: params.nexonerate, file: true),fastaRemoveShort.out.collect(),fastaToCdbindex.out.collect(),genome_rm.collect())
                protExonerateToHints(protExonerateBatch.out.concat(protExonerateFromList.out).collect())

        emit:
                hints = protExonerateToHints.out
                gff = protExonerateBatch.out[0].collectFile()


}

workflow proteinhint_spaln {

	take:
		genome_rm
		protein_fa

	main:
                fastaCleanProteins(protein_fa)
                fastaRemoveShort(fastaCleanProteins.out,params.min_prot_length)
		spalnMakeIndex(genome_rm)
		spalnAlign(fastaRemoveShort.out.splitFasta(by: params.nblast, file: true),spalnMakeIndex.out)
		spalnMerge(spalnAlign.out.collect(),spalnMakeIndex.out,60)
		spalnToHints(spalnMerge.out[0], params.pri_prot)
	
	emit:
		hints = spalnToHints.out
		gff = spalnMerge.out[0]
		track = spalnMerge.out[1]

}

workflow proteinmodels {

	take:
                genome_rm
                protein_fa

        main:
                fastaCleanProteins(protein_fa)
                fastaRemoveShort(fastaCleanProteins.out,params.min_prot_length)
                spalnMakeIndex(genome_rm)
                spalnAlign(fastaRemoveShort.out.splitFasta(by: params.nblast, file: true),spalnMakeIndex.out)
                spalnMerge(spalnAlign.out.collect(),spalnMakeIndex.out,90)
                spalnToHints(spalnMerge.out[0], params.pri_prot_target)
		spaln2evm(spalnMerge.out[0])

        emit:
                hints = spalnToHints.out
                gff = spalnMerge.out[0]
		track = spalnMerge.out[1]

}

process tblastn {

        publishDir "${params.outdir}/logs/tblastn" , mode: 'copy'

	label 'short_running'
	
	input:
	path protein_chunk
	path blastdb_files

	output:
	path protein_blast_report

	script:
	db_name = blastdb_files[0].baseName
	chunk_name = protein_chunk.getName().tokenize('.')[-2]
	protein_blast_report = "${protein_chunk.baseName}.blast"

	"""
		${params.tblastn} -num_threads ${task.cpus} -evalue ${params.blast_evalue} -max_intron_length ${params.max_intron_size} -outfmt \"${params.blast_options}\" -db $db_name -query $protein_chunk > $protein_blast_report
	"""

}

process tblastnMasked {

        publishDir "${params.outdir}/logs/tblastn" , mode: 'copy'

	label 'short_running'
	
	input:
	path protein_chunk
	path blastdb_files

	output:
	path protein_blast_report

	script:
	db_name = blastdb_files[0].baseName
	chunk_name = protein_chunk.getName().tokenize('.')[-2]
	protein_blast_report = "${protein_chunk.baseName}.blast"

	"""
		${params.tblastn} -num_threads ${task.cpus} -evalue ${params.blast_evalue} -db_soft_mask 40 -max_intron_length ${params.max_intron_size} -outfmt \"${params.blast_options}\" -db $db_name -query $protein_chunk > $protein_blast_report
	"""

}

process tblastnToTargets {

        publishDir "${params.outdir}/logs/tblastn" , mode: 'copy'

	input:
	path blast_reports

	output:
	path targets
	
	script:
	query_tag = "ProteinDB"
	targets = "${query_tag}.targets"
	
	"""
		cat $blast_reports | sort -k1 -k2 -k3 >> merged.txt
		tblastn2exonerate_targets.pl --infile merged.txt --max_intron_size $params.max_intron_size --length_percent $params.blast_length_percent --min_id $params.blast_pident > $targets
	"""
}

// Compare a target file with a protein fasta file and find the proteins missing from the target file
process TargetsFindMissing {

	label 'short_running'

	input:
	path targets
	path protein_fa

	output:
	path missing

	script:
	missing = file(protein_fa).getBaseName() + ".missing.txt"

	"""
		blast_find_missing_targets.pl --targets $targets --proteins $protein_fa > $missing
	"""

}

process spalnMakeIndex {

	label 'spaln'

	input:
	path genome

	output:
	path("genome_spaln*")

	script:

	"""
		cp $genome genome_spaln.fa
		spaln -W -KP -E -t${task.cpus} genome_spaln.fa
	"""
	
}

process spalnAlign {

	label 'spaln'

	publishDir "${params.outdir}/logs/spaln", mode: 'copy'

	input:
	path proteins
	path spaln_files

	output:
	path ("${chunk_name}.*")

	script:
	chunk_name = proteins.getBaseName() 
	spaln_gff = chunk_name + ".gff"
	spaln_grd = chunk_name + ".grd"

	"""
		spaln -o $chunk_name -Q${params.spaln_q} -O12 -t${task.cpus} -dgenome_spaln $proteins
		
		if [ -s $spaln_grd ] || [ ! -f $spaln_grd ] 
		then
			spaln -o $chunk_name -Q6 -O12 -t${task.cpus} -dgenome_spaln $proteins
		fi
			
	"""

}

process spalnMerge {

	label 'spaln'

	publishDir "${params.outdir}/logs/spaln" , mode: 'copy'

	input:
	path spaln_reports
	path spaln_files
	val similarity

	output:
	path spaln_final
	path spaln_track 

	script:
	spaln_final = spaln_reports[0].getBaseName() + ".merged.${similarity}.final.gff"
	spaln_track = spaln_reports[0].getBaseName() + ".merged.${similarity}.final.webapollo.gff"
	
	"""
		sortgrcd -C${similarity} -F2 -I${similarity} -J180 -O0 -n0 *.grd > $spaln_final
		sortgrcd -C${similarity} -F2 -J180 -O2 -n0 *.grd > $spaln_track
	"""

}

process spaln2evm {

	input:
	path spaln_models

	output:
	path spaln_evm

	script:
	spaln_evm = spaln_models.getBaseName() + ".evm.gff"

	"""
		sed 's/ALN/spaln2genome/' $spaln_models > tmp
		gff_add_exons.pl --infile tmp > $spaln_evm
		rm tmp
	"""
}

process spalnToHints {

	publishDir "${params.outdir}/logs/spaln", mode: 'copy'

	input:
	path gff
	val priority

	output:
	path hints

	script:
	hints = "spaln.proteins.${priority}.hints.gff"
	"""
		align2hints.pl --in=$gff --maxintronlen=${params.max_intron_size} --prg=spaln --priority=${priority} --out=$hints
	"""
}


// Run exonerate on full genomes for select proteins unable to be located using heuristics
process protExonerateFromList {

        publishDir "${params.outdir}/logs/exonerate_exhaustive", mode: 'copy'

	input:
	path accessions
	path protein_db
	path protein_db_index
	path genome

	output:
	path exonerate_aligns

	script:
	chunk_name = accessions.getBaseName()
	exonerate_aligns = chunk_name +  ".exonerate.out"

	"""
		for i in \$(cat $accessions); do cdbyank -a \$i $protein_db_index > \$i.fasta ; done;
		exonerate_from_list.pl --accessions $accessions --db $protein_db_index --genome $genome --max_intron_size ${params.max_intron_size} > commands.txt
		echo "Starting Exonerate run..."
		parallel -j ${task.cpus} < commands.txt 2>/dev/null
		echo "Finished exonerate run"
		echo '# exonerate alignments for ${accessions}' > $exonerate_aligns
                cat *.exonerate.align | grep -v '#' | grep 'exonerate:protein2genome:local' >> $exonerate_aligns  2>/dev/null  || true
	"""
}

// Run Exonerate on the blast regions
// Takes a list of blast matches and will extract protein sequences and potential target regions
// from a protein database and genome sequence to run exonerate 
process protExonerateBatch {

	scratch true

        publishDir "${params.outdir}/logs/exonerate", mode: 'copy'

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
		exonerate_from_blast_hits.pl --matches $hits_chunk --assembly_index $genome --max_intron_size $params.max_intron_size --analysis protein2genome --outfile $commands
		parallel -j ${task.cpus} < $commands
		cat *.exonerate.align | grep -v '#' | grep 'exonerate:protein2genome:local' >> merged.${chunk_name}.exonerate.out 2>/dev/null
		exonerate_offset2genomic.pl --infile merged.${chunk_name}.exonerate.out --outfile $exonerate_chunk
		[ -s $exonerate_chunk ] || cat "#" > $exonerate_chunk

		rm *.align
		rm *._target_.fa*
		rm *._query_.fa*
	"""
}

// merge the exonerate hits and create the hints
process protExonerateToHints {

	label 'medium_running'

        publishDir "${params.outdir}/logs/exonerate", mode: 'copy'

	input:
	path chunks

	output:
	path exonerate_gff

	script:
	exonerate_gff = "proteins.exonerate.hints.gff"
	"""
		cat $chunks > all_chunks.out
		exonerate2gff.pl --infile all_chunks.out --pri ${params.pri_prot} --source protein --outfile $exonerate_gff
	"""
}
