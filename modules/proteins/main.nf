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
		spaln2gmod(spalnMerge.out)
	
	emit:
		hints = spalnToHints.out
		gff = spalnMerge.out
		evm = spaln2evm.out
		track = spaln2gmod.out
}

// Used for targeted proteins to make a gene build from
workflow proteinmodels_spaln {

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
		spaln2gmod(spalnMerge.out)

        emit:
                hints = spalnToHints.out
                gff = spalnMerge.out
		evm = spaln2evm.out
		track = spaln2gmod.out

}

workflow proteinhint_gth {

	take:
		genome
		protein_fa

	main:
		blast_index(genome)
		fastaCleanProteins(protein_fa)
                fastaRemoveShort(fastaCleanProteins.out,params.min_prot_length)
		blast_proteins(fastaRemoveShort.out.splitFasta(by: params.nproteins, file: true),blast_index.out.collect())
		blast2targets(blast_proteins.out.collect())
		targets2gth(genome.collect(),protein_fa.collect(),blast2targets.out.splitText(by: params.nproteins, file: true))	
		gthToHints(targets2gth.out.collect(),params.pri_prot)
		spaln2evm(targets2gth.out.collect().collectFile())
	emit:
		hints = gthToHints.out
		gff = targets2gth.out.collectFile()
		track = spaln2evm.out

}

workflow proteinmodels_gth {

	take:
                genome
                protein_fa

        main:
                blast_index(genome)
                fastaCleanProteins(protein_fa)
                fastaRemoveShort(fastaCleanProteins.out,params.min_prot_length)
                blast_proteins(fastaRemoveShort.out.splitFasta(by: params.nproteins, file: true),blast_index.out.collect())
                blast2targets(blast_proteins.out.collect())
                targets2gth(genome.collect(),protein_fa.collect(),blast2targets.out.splitText(by: params.nproteins, file: true))
                gthToHints(targets2gth.out.collect(),params.pri_prot)
		spaln2evm(targets2gth.out.collect().collectFile())
        emit:
                hints = gthToHints.out
                gff = targets2gth.out.collectFile()
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

	//publishDir "${params.outdir}/logs/spaln", mode: 'copy'

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

process spaln2gmod {

	publishDir "${params.outdir}/tracks", mode: 'copy'

	input:
	path spaln_models

	output:
	path spaln_track

	script:
	spaln_track = spaln_models.getBaseName() + ".gmod.gff"

	"""
		spaln2gmod.pl --infile $spaln_models > $spaln_track
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

process gthToHints {

	publishDir "${params.outdir}/logs/gth", mode: 'copy'
	
	input:
	path gffs
	val priority

	output:
	path hints

	script:
	hints = file(gffs[0]).getBaseName() + ".proteins.${priority}.hints.gff"

	"""
		cat $gffs > gff
		align2hints.pl --in=gff --maxintronlen=${params.max_intron_size} --prg=gth --priority=${priority} --out=$hints
		rm gff
	"""

}

// make a blast index 
process blast_index {

	label 'blast'

	publishDir "${params.outdir}/blast/", mode: 'copy'

	input:
	path genome_fa
	
	output:
	path "${dbName}*.n*"

	script:
	dbName = genome_fa.getBaseName()

	"""
		makeblastdb -in $genome_fa -dbtype nucl -parse_seqids -out $dbName
	"""
}

process blast_proteins {

	scratch true

	label 'blast'

	publishDir "${params.outdir}/logs/tblastn", mode: 'copy'

	input:
	path protein_chunk
	path blastdb_files

	output:
	path blast_results

	script:
	blast_results = protein_chunk.getBaseName() + ".blast"
	db_name = blastdb_files[0].baseName

	"""
		tblastn -num_threads ${task.cpus} \
			-evalue ${params.blast_evalue} \
			-max_intron_length ${params.max_intron_size} \
			-outfmt 6 \
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
		cat $blast_results > blast.out
		cut -f1,2 blast.out | sort -u | sort -k2 > $targets
		rm blast.out
	"""
}	

process targets2gth {

	scratch true

	input:
	path genome
	path protein_fa
	path target_chunk

	output:
	path alignments

	script:
	alignments = target_chunk + ".gth"

	def options = ""	
	if (params.gth_options) {
		options = params.gth_options
	}

	"""
		cut -f1 $target_chunk | sort -u > proteins.txt
		cut -f2 $target_chunk | sort -u > dna.txt
		mkdir -p protein_db && cd protein_db && fasta_extract_from_list.pl --list ../proteins.txt --fasta ../${protein_fa} > jobs.sh && bash jobs.sh && cd ..
		mkdir -p genome_db && cd genome_db && fasta_extract_from_list.pl --list ../dna.txt --fasta ../${genome} > jobs.sh && bash jobs.sh && cd ..
		gth_from_targets.pl --targets $target_chunk --options $options > all_commands.txt
		grep create all_commands.txt > indexing.sh
		bash indexing.sh
		rm *.gth.out
		grep -v create all_commands.txt > commands.txt
		parallel -j ${task.cpus} < commands.txt
		cat *.gth.out > $alignments
		rm -Rf genome_db protein_db *.gth.out
	"""

}

