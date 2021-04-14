process make_blast_nt_index {

	input:
	path fasta

	output:
	path ("${fasta}*")

	script:

	"""
		makeblastdb -in $fasta -dbtype nucl 
	"""

}

process align_transcripts {


	input:
	path fasta
	path blast_index

	output:
	path blast_results

	script:
	index_name = blast_index[0]
	base_name = blast_index[0].getBaseName()
	blast_result = fasta.getBaseName() + "." + base_name + ".blast"

	"""
		blastn -query $fasta -db $index_name -num_threads ${task.cpus} -evalue 0.0001 -
	"""

}
