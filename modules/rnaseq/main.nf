// Workflow to perform RNA-seq mapping against a reference
// requires: genome sequence in FASTA format (unmasked is fine) and one or more RNA-seq read sets

include {bam_index} from "../util.nf" params(params)

workflow rnaseqhint_hisat {

	take:
		genome
		reads

	main:
		HisatMakeDB(genome)
		runFastp(reads )
		HisatMap(runFastp.out[0],HisatMakeDB.out.collect())
		makeBigWig(HisatMap.out[0],HisatMap.out[1])
		mergeBams(HisatMap.out[0].collect())
		BamToExonHint(mergeBams.out)
		BamToIntronHint(mergeBams.out)
		filterRseqHints(BamToIntronHint.out.concat(BamToExonHint.out).collectFile())		
			
	emit:
		bam = mergeBams.out[0]
		hints = filterRseqHints.out[0]		
}

workflow rnaseqhint_star {

	take:
                genome
                reads

        main:
                STARmakeDB(genome)
                runFastp(reads)
                STARalign(runFastp.out[0],STARmakeDB.out.collect())
		STARalignTwo(runFastp.out[0],STARmakeDB.out.collect(),STARalign.out.collect())
		bam_index(STARalignTwo.out[0])
		makeBigWig(bam_index.out[0],bam_index.out[1])
                mergeBams(STARalignTwo.out[0].collect())
                BamToExonHint(mergeBams.out)
                BamToIntronHint(mergeBams.out)
                filterRseqHints(BamToIntronHint.out.concat(BamToExonHint.out).collectFile())

        emit:
                bam = mergeBams.out[0]
                hints = filterRseqHints.out[0]


}


// trim RNAseq reads
process runFastp {

	label 'fastp'

	//publishDir "${params.outdir}/evidence/rnaseq/fastp", mode: 'copy'

	scratch true 

	input:
	tuple name, path(reads)
		
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

process STARmakeDB {

	label 'star'

	publishDir "${params.outdir}/databases/STAR", mode: 'copy'

	input:
	path genome

	output:
	path star_dir

	script:
	star_dir = "STAR_INDEX"
	"""
		mkdir -p $star_dir
		STAR --runThreadN ${task.cpus} \
			--runMode genomeGenerate \
			--genomeDir $star_dir \
			--genomeFastaFiles $genome \
			--limitGenomeGenerateRAM ${task.memory.toBytes()}
	"""
}

process STARalign {

	scratch true

	label 'star'

	input:
	path reads
	path star_index

	output:
	path junctions

	script:
	options = "--outFilterType BySJout --outFilterMultimapNmax 5 --outSAMstrandField intronMotif"
	ReadsBase = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
	junctions = ReadsBase + ".SJ.out.tab"
	bam = ReadsBase + ".aligned.bam"
	"""
		STAR --runThreadN ${task.cpus} \
			--genomeDir $star_index \
			--readFilesCommand zcat \
			--limitBAMsortRAM ${task.memory.toBytes()/4} \
			--readFilesIn $reads \
			--alignIntronMin 20 \
			--alignIntronMax $params.max_intron_size \
			--outSAMtype BAM SortedByCoordinate \
			$options

			cp SJ.out.tab $junctions
	"""
}

process STARalignTwo {

        scratch true

        label 'star'

        publishDir "${params.outdir}/logs/rnaseq", mode: 'copy'

        input:
        path reads
        path star_index
	path all_juncs

        output:
        path bam
        path junctions

        script:
        options = "--outFilterType BySJout --outFilterMultimapNmax 5 --outSAMstrandField intronMotif"
        ReadsBase = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
        junctions = ReadsBase + ".SJ.out.tab"
        bam = ReadsBase + ".aligned.bam"
        """
                STAR --runThreadN ${task.cpus} \
                        --genomeDir $star_index \
                        --readFilesCommand zcat \
                        --readFilesIn $reads \
                        --alignIntronMin 20 \
			--limitBAMsortRAM ${task.memory.toBytes()/4} \
			--sjdbFileChrStartEnd $all_juncs \
                        --alignIntronMax $params.max_intron_size \
                        --outSAMtype BAM SortedByCoordinate \
                        $options
                        cp Aligned.sortedByCoord.out.bam $bam
                        cp SJ.out.tab $junctions
        """
}

// Generate an alignment index from the genome sequence
// We use HiSat so we don't need to know the read length during index construction, unlike STAR
process HisatMakeDB {

	label 'hisat'

	//publishDir "${params.outdir}/databases/HisatDB", mode: 'copy'

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

process HisatMap {

	label 'hisat'

	publishDir "${params.outdir}/logs/rnaseq", mode: 'copy'
	
	scratch true

	input:
	path reads
	path hs2_indices
	
	output:
	path "*accepted_hits.bam"
	path "*.bai"

	script:
	indexBase = hs2_indices[0].toString() - ~/.\d.ht2/
	ReadsBase = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/

	prefix = ReadsBase

	if (params.singleEnd) {
		"""
		hisat2 -x $indexBase -U $reads -p ${task.cpus} | samtools view -bS - | samtools sort -m 2G -@ 4 - > ${prefix}_accepted_hits.bam
		samtools index  ${prefix}_accepted_hits.bam
		"""
	} else {
		"""
		hisat2 -x $indexBase -1 ${reads[0]} -2 ${reads[1]} -p ${task.cpus} | samtools view -bS - | samtools sort -m 2G -@4 - > ${prefix}_accepted_hits.bam
		samtools index  ${prefix}_accepted_hits.bam
		"""
	}
}

process makeBigWig {

	scratch true

	label 'deeptools'

	publishDir "${params.outdir}/tracks", mode: 'copy'

	input:
	path bam
	path bai

	output:
	path "*.bw"

	script:

	bigwig = bam.getBaseName() + ".bw"
	bigwig_fw = bam.getBaseName() + ".fw.bw"
	bigwig_rev = bam.getBaseName() + ".rev.bw"

	if (params.rnaseq_stranded) {
		"""
			bamCoverage -of bigwig -bs 10 --filterRNAstrand forward --ignoreDuplicates  -p ${task.cpus} -b $bam -o $bigwig_fw
        	        bamCoverage -of bigwig -bs 10 --filterRNAstrand reverse --ignoreDuplicates  -p ${task.cpus} -b $bam -o $bigwig_rev
		"""

	} else {
		"""
			bamCoverage -of bigwig -bs 10 --ignoreDuplicates -p ${task.cpus} -b $bam -o $bigwig
		"""
	}

}

// Combine all BAM files for hint generation
process mergeBams {

	//publishDir "${params.outdir}/evidence/rnaseq/", mode: 'copy'

	scratch true 

	input:
	path bams

	output:
	path bam	

	script:
	bam = "rnaseq.merged.bam"
	avail_ram_per_core = (task.memory/task.cpus).toGiga()-1
	
	"""
		samtools merge - $bams | samtools sort -@ ${task.cpus} -m${avail_ram_per_core}G - > $bam
	"""
}

// Only make ep hints if utr annotation is active, else we get tons of small genes in UTR regions. 
process BamToExonHint {

	scratch true

        publishDir "${params.outdir}/evidence/rnaseq/", mode: 'copy'

	input:
	path bam

	output:
	path hints

	script:
	hints = "rnaseq.merged.exons.gff"

	if (params.utr) {
	 
		"""
		samtools depth $bam | perl -ne 'BEGIN{ print "track type= print wiggle_0 name=merged_reads description=merged_reads\n"}; (\$c, \$start, \$depth) = split;if (\$c ne \$lastC) { print "variableStep chrom=\$c span=10\n"; };\$lastC=\$c;next unless \$. % 10 ==0;print "\$start\t\$depth\n" unless \$depth<3' | wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --radius=4.5 --pri=${params.pri_wiggle} --strand='.' > $hints
		"""
	} else {
		"""
		echo '##gff-version 3' >> $hints
		"""
	}
}

/*
 * STEP RNAseq.4 - Hisat2 into Hints
 */	
process BamToIntronHint {

	//publishDir "${params.outdir}/evidence/rnaseq/hints", mode: 'copy'
	label 'augustus'

	input:
	path bam
	
	output:
	path hisat_hints

	script:
	hisat_hints = "rnaseq.merged.introns.gff"

	"""
		bam2hints --intronsonly 0 -p ${params.pri_rnaseq} -s 'E' --in=$bam --out=$hisat_hints
	"""
}

// Remove  hints that only have little multi support (spurious!)
process filterRseqHints {

	input:
	path hints

	output:
	path filtered_hints

	script:
	filtered_hints = hints.getBaseName() + ".filtered.gff"

	"""
		gff_filter_by_multi.pl --infile $hints > $filtered_hints
	"""

}
