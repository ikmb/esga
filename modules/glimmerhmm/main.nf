workflow glimmerhmm_train {


	take:
	genome
	models

	main:
	models2exons(models)
	train_glimmer(genome,models2exons.out)
	
	emit:
	glimmer_model = train_glimmer.out

}


workflow glimmerhmm {

	take:
	genome
	glimmer_model

	main:
	glimmerhmm(genome,glimmermodel)

	emit:
	annotation = glimmerhmm.out

}
