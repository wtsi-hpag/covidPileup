# covidPileup v1.0
Pipeline for building SNP and GAP pileup from multiple COVID-19 whole genome sequences.

### Download and Compile:

    $ git clone  https://github.com/wtsi-hpag/covidPileup.git 
    $ cd covidPileup 
    $ ./install.sh
		
If everything compiled successfully you must see the final comment: 
		"Congrats: installation successful!"		


#### External packages
The genome aligner BWA (http://bio-bwa.sourceforge.net) is downloaded and compiled by covidPileup.

### Run the pipelines

#### Run covidPileup:
           $ /full/path/to/covidPileup/src/covidPileup -nodes <nodes> -SNP <plot> -GAP <plot> \
		 reference.fasta COVID-19_genomes.fasta sample-pileup
           

	       Parameters:
             nodes:        number of CPUs requested  [ default = 30 ]
             length:       reference genome length [ default = 40000 ]
             SNP:          output image file on SNP pileup [ default = not plot ]
             GAP:          output image file on GAP pileup [ default = not plot ]
	    
#### Run example with reference.fasta COVID-19_genomes.fasta 

           $ /full/path/to/covidPileup/src/covidPileup -nodes 30 -SNP plot -GAP plot \
		 reference.fasta COVID-19_genomes.fasta sample-pileup
	    


