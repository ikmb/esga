process repeatModel {

        publishDir "${params.outdir}/repeatmodeler/", mode: 'copy'

        scratch true

        input:
        file genome_fa 

        output:
	file repeats

        script:

        repeats = "consensi.fa"
        """
                BuildDatabase -name genome_source -engine ncbi $genome_fa
                RepeatModeler -engine ncbi -pa ${task.cpus} -database genome_source
        	cp RM_*/consensi.fa .
	"""
}

