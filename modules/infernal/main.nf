include fastaSplitSize from "./../fasta" params(params)

workflow rfamsearch {

	take:
		genome
		
	main:
                fastaSplitSize(genome,params.npart_size)		
		download_rfam(genome)
		infernal_press(download_rfam.out[0])
		infernal_search(infernal_press.out.collect(),fastaSplitSize.out.flatMap())
		infernal2gff(infernal_search.out.collectFile(),download_rfam.out[1].collect())
	emit:
		gff = infernal2gff.out
	
}

process download_rfam {

	executor 'local'

	input:
	path genome

	output:
	path "Rfam.cm"
	path "family.txt"

	script:

	"""
		wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.2/Rfam.cm.gz
		gunzip Rfam.cm.gz
		wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.2/database_files/family.txt.gz
		gunzip family.txt.gz
	"""

}

process infernal_press {

	label 'infernal'

	input:
	path rfam

	output:
	path "Rfam.cm*"

	script:

	"""
		cmpress Rfam.cm
	"""
}

process infernal_search {

	publishDir "${params.outdir}/logs/rfam", mode: 'copy'

	label 'infernal'

	input:
	path rfam_db
	path genome_chunk

	output:
	path rfam_tbl

	script:

	rfam_tbl = genome_chunk.getBaseName() + ".rfam.tbl"
	rfam_txt = genome_chunk.getBaseName() + ".rfam.out"

	"""
		cmsearch --rfam --cpu ${task.cpus} --cut_tc --tblout $rfam_tbl -o $rfam_txt Rfam.cm $genome_chunk
	"""
}

process infernal2gff {

	input:
	path rfam_report
	path families

	output:
	path gff

	script:
	gff = "rfam.gff"

	"""
		rfam2gff.pl --infile $rfam_report --family $families > $gff
	"""
}
