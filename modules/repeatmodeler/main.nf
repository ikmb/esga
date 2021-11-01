// model repeats denovo

workflow model_repeats {

	take:
	genome

	main:
	repeatModel(genome)

	emit:
	repeats = repeatModel.out[0]

}

process repeatModel {

	label 'repeatmodeler'

        publishDir "${params.outdir}/repeatmodeler/", mode: 'copy'

        scratch true

        input:
        path genome_fa 

        output:
	path repeats

        script:

        repeats = "consensi.fa"
        """
                BuildDatabase -name genome_source -engine ncbi $genome_fa
                RepeatModeler -engine ncbi -pa ${task.cpus} -database genome_source
        	cp RM_*/consensi.fa .
	"""
}

