![](images/deNBI_logo.jpg) ![](images/ikmb_bfx_logo.png)

# ESGA - Genome Annotation (v2)

[![Nextflow](https://img.shields.io/badge/nextflow-20.01.0-brightgreen)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](http://singularity.lbl.gov)

This pipeline peforms annotation of novel genomes using a combination of evidence alignments, evidence-based gene building and ab-initio gene building.

### Pipeline main steps

The minimum requirements are a genome file and at least one type of evidence.

From this, the pipeline can run the following processing steps:

* align proteins against a genome and generate annotation hints
* align transcripts against a genome and generate annotation hints
* align RNA-seq reads against a genome and generate annotation hints
* assemble transcripts from aligned RNA-seq reads and generate annotation hints (Trinity pipeline)
* Produce evidence-based gene models from aligned transcript sequences (PASA pipeline)
* Produce ab-initio, hint-supported gene models (AUGUSTUS pipeline)
* Produce consensus annotation from all of the above (EvidenceModeler pipeline)

Optional:

* Train a novel ab-initio prediction profile for AUGUSTUS (using PASA transcripts)

### Test data

A simple test data set can be downloaded [here](https://drive.google.com/open?id=1VFqLnRJiuj5Vhj2KCOdY58jwxZKkkMVU)

### Documentation

Documentation about the pipeline can be found in the `docs/` directory or under the links below:

1. [What happens in this pipeline?](docs/pipeline.md)
2. [Recommendations](docs/recommendations.md)
3. [Installation and configuration](docs/installation.md)
4. [Running the pipeline](docs/usage.md)
5. [Output](docs/output.md)
6. [Troubleshooting](docs/troubleshooting.md)
7. [What's next](docs/whatsnext.md)

### Credits

This pipeline was written by Dr. Montserrat Torres ([MontseTor](https://github.com/MontseTor)) and Dr. Marc Höppner ([marchoeppner](https://github.com/marchoeppner)) at [IKMB](http://www.ikmb.uni-kiel.de).
The authors gratefully acknowledge inspiration, fruitful discussions and a few useful code snippets from the [nf-core](https://www.nf-co.re).


