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
process spalnAlign {

	scratch true

	publishDir "${params.outdir}/logs/spaln", mode: 'copy'

	input:
	path proteins
	path spaln_files
	val comp_para

	output:
	path ("${chunk_name}.*")

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
