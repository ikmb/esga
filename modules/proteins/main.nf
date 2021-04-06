// **************************
// Production of annotation hints from protein FASTA sequences
// **************************

include { fastaCleanProteins; fastaRemoveShort; assemblySplit } from "./../fasta" params(params)

// Used for related protein sequences
workflow proteinhint_spaln {

	take:
		genome_rm
		protein_fa

	main:
                fastaCleanProteins(protein_fa)
                fastaRemoveShort(fastaCleanProteins.out,params.min_prot_length)
		spalnMakeIndex(genome_rm)
		spalnAlign(fastaRemoveShort.out.splitFasta(by: params.nproteins, file: true),spalnMakeIndex.out,1)
		spalnMerge(spalnAlign.out.collect(),spalnMakeIndex.out,60)
		spalnToHints(spalnMerge.out, params.pri_prot)
		spaln2evm(spalnMerge.out)
	
	emit:
		hints = spalnToHints.out
		gff = spalnMerge.out
		track = spaln2evm.out
}

// Used for targeted proteins to make a gene build from
workflow proteinmodels {

	take:
                genome_rm
                protein_fa

        main:
                fastaCleanProteins(protein_fa)
                fastaRemoveShort(fastaCleanProteins.out,params.min_prot_length)
                spalnMakeIndex(genome_rm)
                spalnAlign(fastaRemoveShort.out.splitFasta(by: params.nproteins, file: true),spalnMakeIndex.out,1)
                spalnMerge(spalnAlign.out.collect(),spalnMakeIndex.out,90)
                spalnToHints(spalnMerge.out, params.pri_prot_target)
		spaln2evm(spalnMerge.out)

        emit:
                hints = spalnToHints.out
                gff = spalnMerge.out
		track = spaln2evm.out

}

// Create a genome index for spaln
process spalnMakeIndex {

	publishDir "${params.outdir}/logs/spaln", mode: 'copy'

	input:
	path genome

	output:
	path("genome_spaln*")

	script:

	"""
		cp $genome genome_spaln.fa
		spaln -W -KP -t${task.cpus} genome_spaln.fa
	"""
	
}

// the comp_para parameter defines whether an alignment is intra (0) or inter-species (1). For targeted proteins, we use 0
// output can be empty if no alignments are found - hence it is optional
process spalnAlign {

	scratch true

	publishDir "${params.outdir}/logs/spaln", mode: 'copy'

	input:
	path proteins
	path spaln_files
	val comp_para

	output:
	path ("${chunk_name}.*") optional true

	script:
	chunk_name = proteins.getBaseName() 
	spaln_gff = chunk_name + ".gff"
	spaln_grd = chunk_name + ".grd"

	"""
		spaln -o $chunk_name -Q${params.spaln_q} -T${params.spaln_taxon} ${params.spaln_options} -O12 -t${task.cpus} -Dgenome_spaln $proteins

	"""

}

// Merge spaln output across chunks using the companion tool sortgrcd
process spalnMerge {

	publishDir "${params.outdir}/logs/spaln" , mode: 'copy'

	input:
	path spaln_reports
	path spaln_files
	val similarity

	output:
	path spaln_final

	script:
	spaln_final = spaln_reports[0].getBaseName() + ".merged.${similarity}.final.gff"
	
	"""
		sortgrcd  -I${similarity} -O0 -n0 *.grd > merged.gff
		spaln_add_exons.pl --infile merged.gff > $spaln_final
		rm merged.gff
	"""

}

// Convert spaln models into EVM compatible format
process spaln2evm {

        publishDir "${params.outdir}/logs/spaln", mode: 'copy'

	input:
	path spaln_models

	output:
	path spaln_evm

	script:
	spaln_evm = spaln_models.getBaseName() + ".evm.gff"

	"""
		spaln2evm.pl --infile $spaln_models > $spaln_evm
	"""
}

// Convert spaln models into AUGUSTUS compatible hint format
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

// make a blast index 
process blast_index {

	publishDir "${params.outdir}/blast/", mode: 'copy'

	input:
	path genome_fa
	
	output:
	path "${dbName}*.n*"

	script:
	dbName = genome_fa.getBaseName()

	"""
		makeblastdb -in $fasta -dbtype nucl -parse_seqids -out $dbName
	"""
}

process blast_proteins {

	input:
	path protein_chunk
	path blast_files

	output:
	path blast_results

	script:
	blast_results = protein_chunk.getBaseName() + ".blast"
	db_name = blastdb_files[0].baseName

	"""
		tblastn -num_threads ${task.cpus} \
			-evalue ${params.blast_evalue} \
			-max_intron_length ${params.max_intron_size} \
			-db $db_name -query $protein_chunk > $blast_results
	"""

}

process blast2targets {

	input:
	path blast_results

	output:
	path targets

	script:
	targets = "blast_targets.bed"

	"""
		tblastn2exonerate_targets.pl --infile merged.txt > $targets
	"""
}	

process blast2exonerate {

	input:
	path genome
	path targets
	
	output:
	path exonerate_out

	script:
	exonerate_out = blast_result.getBaseName() + ".exonerate.out"

	"""
		exonerate_from_blast_hits.pl --matches $targets \
			--genome_index $genome \
			--max_intron_size $params.max_intron_size \
			--protein_index $protein_db_index \
			--outfile $commands
	"""
}
