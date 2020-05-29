# covidPileup v1.1
Pipeline for building SNP and GAP pileup from multiple COVID-19 whole genome sequences.

### Download and Compile:

    $ git clone  https://github.com/wtsi-hpag/covidPileup.git 
    $ cd covidPileup 
    $ bash install.sh
		
If everything compiled successfully you must see the final comment: 
		"Congrats: installation successful!"		


#### External packages
The genome aligner BWA (http://bio-bwa.sourceforge.net) and SMALT (http://www.sanger.ac.uk/science/tools/smalt-0) are downloaded and compiled by covidPileup.

### Run the pipelines

#### Run covidPileup:
           $ /full/path/to/covidPileup/src/covidPileup -nodes <nodes> -SNP <plot> -GAP <plot> \
		 -country UK reference.fasta COVID-19_genomes.fasta sample-pileup
           

	       Parameters:
             nodes:        number of CPUs requested  [ default = 30 ]
             length:       reference genome length [ default = 40000 ]
             SNP:          output 4 image files on SNP pileup [ default = not plot ]
                           1. sample-pileup.pileupSNP.png; 2. sample-pileup.frequeSNP.png
                           3. sample-pileup.uniqueSNP.png; 4. sample-pileup.regionSNP.png

             GAP:          output 2 image files on GAP pileup [ default = not plot ]
                           5. sample-pileup.pileupGAP.png; 6. sample-pileup.frequeGAP.png
	         country       specific SNPs in the country [ default = UK ]
                           "EU" and Europen countries like "France", "Italy" et al should work;
                           "UK" and "England", "Scotland", "Wales" et al should work;

               Other output files:
             sample-pileup.name - country names and sample number for earch country
             sample-pileup.snps - all the SNPs detected by the pipeline
             sample-pileup.spec - country specific SNPs 
                  
#### Run example with reference.fasta COVID-19_genomes.fasta 

           $ /full/path/to/covidPileup/src/covidPileup -nodes 30 -SNP plot -GAP plot -country UK\
		 reference.fasta COVID-19_genomes.fasta sample-pileup
	   
or using SMALT 
           $ /full/path/to/covidPileup/src/covidPileup -nodes 30 -SNP plot -GAP plot -country UK\
		 -align smalt reference.fasta COVID-19_genomes.fasta sample-pileup
	    
#### Data input fasta files 
     The names in the fasta file have to follow the format used by Gisaid

     Australia/VIC110/2020
     USA/WA-UW-4032/2020
     England/20124095202/2020
 


