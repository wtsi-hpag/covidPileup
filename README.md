# covidPileup v1.1
Pipeline for COVID-19 variation analysis from whole genome sequences.

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
		 -country <country_name> -cover <n_cover> -length <reference_length> \
                 reference.fasta COVID-19_genomes.fasta sample-pileup

	       Parameters:
             nodes:        Number of CPUs requested  [ default = 30 ]
             length:       Reference genome length [ default = 40000 ]
             cover:        Threshold coverage number to report specific SNPs  [ default = 5 o] \
	         country:      Specific SNPs in the country [ default = UK ]
                           "EU" and Europen countries like "France", "Italy" et al should work;
                           "UK" and "England", "Scotland", "Wales" et al should work;
             SNP:          Output 4 image files on SNP pileup [ default not plot ]
                           1. sample-pileup.pileupSNP.png; 2. sample-pileup.frequeSNP.png
                           3. sample-pileup.uniqueSNP.png; 4. sample-pileup.regionSNP.png

             GAP:          Output 2 image files on GAP pileup [ default not plot ]
                           5. sample-pileup.pileupGAP.png; 6. sample-pileup.frequeGAP.png

               Other output files:
             sample-pileup.name - country names and sample number for earch country
             sample-pileup.snps - all the SNPs detected by the pipeline
             sample-pileup.spec - country specific (unique) SNPs 
                  
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
 
#### Start a new run without doing the alignment 
       The aligner SMALT is much slower than BWA, but it offers better alignments for sequences \
       (1) with a lot of "N"s; (2) lower similarity with the reference such as bats or pangolins \
       Personally, I would recommend SMALT, which takes 2-3 hours with 60 CPUs for 50K genomes.

       If you had made a run previously with a temporary directory:       \ 

          /lustre/team117/zn1/project/covid/tmp_rununik_71885/            \

       You may use the already existing alignment file to save time:      \
  
           $ /full/path/to/covidPileup/src/covidPileup -nodes <nodes> -SNP <plot> -GAP <plot> \
		 -country <country_name> -cover <n_cover> -length <reference_length> \
                 -data /lustre/team117/zn1/project/covid/tmp_rununik_71885/align.dat \
                 reference.fasta COVID-19_genomes.fasta sample-pileup

       Note: "-data ?" expects a full path.

