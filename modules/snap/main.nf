include { fastaSplitSize } from "./../fasta" params(params)

workflow snap_train_from_spaln {

	take:
	genome
	models

	main:
	spaln_to_snap(models)
	train_snap(genome,spaln_to_snap.out)

	emit:
	hmm = train_snap.out

}

workflow snap_train_from_pasa {

        take:
        genome
        models

        main:
        pasa_to_snap(models)
        train_snap(genome,pasa_to_snap.out)

        emit:
        hmm = train_snap.out

}


workflow snap {

	take:
	genome
	hmm

	main:
	fastaSplitSize(genome,params.npart_size)
	run_snap(fastaSplitSize.out,hmm.collect())

	emit:
	annotation = run_snap.out.collectFile()

}

process pasa_to_snap {

	input:
	path gff

	output:
	path zff

	script:
	zff = gff.getBaseName() + ".zff"

	"""
		touch $zff"
	"""

}

process spaln_to_snap {

	input:
	path gff

	output:
	path zff

	script:
	zff = gff.getBaseName() + ".zff"

	"""
		gff2snap.pl < $gff > $zff
	"""

}


process train_snap {

	label 'snap'

	publishDir "${params.outdir}/ab-initio/snap/", mode: 'copy'

	input:
	path genome
	path models

	output:
	path hmm

	script:
	base_name = genome.getBaseName() + ".snap"
	hmm = base_name + ".hmm"

	"""
		fathom -categorize 1000 $models $genome
	        fathom -export 1000 uni.ann uni.dna
        	forge export.ann export.dna
	        hmm-assembler.pl $base_name . > $hmm

	"""
}

process run_snap {

	label 'snap'

	publishDir "${params.outdir}/ab-initio/snap/", mode: 'copy'

	input:
	path genome
	path hmm

	output:
	path annotation

	script:
	annotation = genome.getBaseName() + "snap.gff3"

	"""
		snap -lcmask $hmm $genome > models.zff
		zff2gff3.pl < models.zff > $annotation
	"""
}
