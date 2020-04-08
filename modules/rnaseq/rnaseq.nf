workflow rnaseqhint {

	take:
		genome
		reads

	main:
		rseqMakeDB(genome)
		rseqTrim(Channel.fromFilePairs(reads).ifEmpty { exit 1, "Did not find any matching read files" } )

}


// trim reads
process rseqTrim {

	publishDir "${params.outdir}/evidence/rnaseq/fastp", mode: 'copy'

	scratch true 

	input:
	val name
	path reads
		
	output:
	path "*_trimmed.fastq.gz"
	path json
	path html

	script:
	prefix = reads[0].toString().split("_R1")[0]
	json = prefix + ".fastp.json"
	html = prefix + ".fastp.html"

	if (params.singleEnd) {
		left = file(reads[0]).getBaseName() + "_trimmed.fastq.gz"
		"""
              		fastp -i ${reads[0]} --out1 ${left} -w ${task.cpus} -j $json -h $html
                """
	} else {
		left = file(reads[0]).getBaseName() + "_trimmed.fastq.gz"
		right = file(reads[1]).getBaseName() + "_trimmed.fastq.gz"
		"""
			fastp --detect_adapter_for_pe --in1 ${reads[0]} --in2 ${reads[1]} --out1 $left --out2 $right -w ${task.cpus} -j $json -h $html
		"""
	}
}

// Generate an alignment index from the genome sequence
process rseqMakeDB {

	label 'long_running'

	publishDir "${params.outdir}/databases/HisatDB", mode: 'copy'

	input:
	path genome
	
	output:
	path "${dbName}.*.ht2"
	
	script:
	dbName = genome.baseName
	dbName_1 = dbName + ".1.ht2"
		
	prefix = dbName
	"""
	hisat2-build $genome $dbName -p ${task.cpus}
	"""
}

/*
 * STEP RNAseq.3 - Hisat2
 */

process rseqMap {

	publishDir "${OUTDIR}/evidence/rnaseq/Hisat2/libraries", mode: 'copy'
	
	scratch true

	input:
	path reads
	path hs2_indices
	
	output:
	path "*accepted_hits.bam"
	
	script:
	indexBase = hs2_indices[0].toString() - ~/.\d.ht2/
	ReadsBase = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/

	prefix = ReadsBase

	if (params.singleEnd) {
		"""
		hisat2 -x $indexBase -U $reads -p ${task.cpus} | samtools view -bS - | samtools sort -m 2G -@ 4 - > ${prefix}_accepted_hits.bam
		"""
	} else {
		"""
		hisat2 -x $indexBase -1 ${reads[0]} -2 ${reads[1]} -p ${task.cpus} | samtools view -bS - | samtools sort -m 2G -@4 - > ${prefix}_accepted_hits.bam
		"""
	}
}

// Combine all BAM files for hint generation
process rseqMergeBams {

	publishDir "${params.outdir}/evidence/rnaseq/Hisat2", mode: 'copy'

	scratch true 

	input:
	path hisat_bams

	output:
	path bam	

	script:
	bam = "hisat2.merged.bam"
	avail_ram_per_core = (task.memory/task.cpus).toGiga()-1
	
	"""
		samtools merge - $hisat_bams | samtools sort -@ ${task.cpus} -m${avail_ram_per_core}G - > $bam
	"""
}

/*
 * STEP RNAseq.4 - Hisat2 into Hints
 */	
process rseqHints {

	publishDir "${params.outdir}/evidence/rnaseq/hints", mode: 'copy'

	input:
	path bam
	
	output:
	path hisat_hints

	script:
	hisat_hints = "RNAseq.hisat.hints.gff"

	"""
		bam2hints --intronsonly 0 -p ${params.pri_rnaseq} -s 'E' --in=$bam --out=$hisat_hints
	"""
}
