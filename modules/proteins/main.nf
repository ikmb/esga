// **************************
// Production of annotation hints from protein FASTA sequences
// **************************

include assemblySplit from "./../fasta" params(params)

workflow proteinhint {

	take:
		genome_rm
		protein_fa

	main:
		assemblySplit(genome_rm,params.chunk_size)
		fastaToDiamondDB(protein_fa)
		fastaToCdbindex(protein_fa)		
		runDiamondx(assemblySplit.out[0].splitFasta(by: params.nblast, file: true),fastaToDiamondDB.out.collect())
		diamondxToTargets(runDiamondx.out.collect(),assemblySplit.out[1])
		protExonerateBatch(diamondxToTargets.out.splitText(by: params.nexonerate, file: true),protein_fa,fastaToCdbindex.out,genome_rm)
		protExonerateToHints(protExonerateBatch.out.collect())

	emit:
		hints = protExonerateToHints.out
		gff = protExonerateBatch.out[0].collectFile()
}

process fastaToDiamondDB {

	label 'medium_running'

        input:
	path protein_fa   

        output:
       	path "${dbName}.dmnd"

        script:
       	dbName = protein_fa.getBaseName()
        """
		diamond makedb --in $protein_fa --db $dbName
       	"""
}

// create a cdbtools compatible  index
// we need this to do very focused exonerate searches later
process fastaToCdbindex {

	label 'short_running'

	input:
	path fasta

	output:
	path protein_index

	script:
	protein_index = fasta.getName()+ ".cidx"

	"""
		cdbfasta $fasta 
	"""
}

// Blast each genome chunk against the protein database
// This is used to define targets for exhaustive exonerate alignments
process runDiamondx {

	//scratch true

	input:
	path genome_chunk
	path db_files

	output:
	path protein_blast_report

	script:
	db_name = db_files[0].baseName
	chunk_name = genome_chunk.getName().tokenize('.')[-2]
	protein_blast_report = "${genome_chunk.baseName}.blast"
	"""
		diamond blastx --sensitive --threads ${task.cpus} --evalue ${params.blast_evalue} --outfmt ${params.blast_options} --db $db_name --query $genome_chunk --out $protein_blast_report
	"""
}

// Parse Protein Blast output for exonerate processing
process diamondxToTargets {

	label 'short_running'

	input:
	file blast_reports
	path genome_agp

	output:
	path targets
	
	script:
	query_tag = "proteinDB"
	targets = "${query_tag}.targets"
	
	"""
		cat $blast_reports > merged.txt
		blast_chunk_to_toplevel.pl --blast merged.txt --agp $genome_agp > merged.translated.txt
		blast2exonerate_targets.pl --infile merged.translated.txt --max_intron_size $params.max_intron_size > $targets
		rm merged.*.txt
	"""
}

// Run Exonerate on the blast regions
// Takes a list of blast matches and will extract protein sequences and potential target regions
// from a protein database and genome sequence to run exonerate 
process protExonerateBatch {

	scratch true

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
	"""
}

// merge the exonerate hits and create the hints
process protExonerateToHints {

	label 'medium_running'

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
