/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2020  Genome Research Ltd.                                *
 *                                                                          *
 *  Author: Zemin Ning (zn1@sanger.ac.uk)                                   *
 *                                                                          *
 *  This file is part of covidPileup pipeline.                              *
 *                                                                          *
 *  covidPileup is a free software: you can redistribute it and/or modify it*
 *  under the terms of the GNU General Public License as published by the   *
 *  Free Software Foundation, either version 3 of the License, or (at your  *
 *  option) any later version.                                              *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful, but     *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        *
 *  General Public License for more details.                                *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License along *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                          *
 ****************************************************************************
 ****************************************************************************/
/****************************************************************************/

 
#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include "fasta.h"

#define GT '>'
#define GT4 (((((GT<<8)+GT)<<8)+GT)<<8)+GT

#define ENDS_EXTRA 0
#define PADCHAR '-'
#define MAX_N_BRG 50000 
#define MAX_N_ROW 40000 
#define Max_N_NameBase 400 
#define Max_N_NameBase2 400 
#define Max_N_Pair 100
static char bindir[2000];
static char tmpdir[2000];
static char **S_Name;
static int *insert_siz,*insert_dev,*core_list,*ctg_list,*ctg_head,*read2contig;
static int *readIndex;

/* SSAS default parameters   */
static int n_group=0;
static int seq_len = 40000;
static int file_tag = 0;
static int run_align = 1;
static int plot_SNP = 0;
static int plot_GAP = 0;
static int min_len = 3000;
static int sam_flag = 0;
static int ances_flag = 0;
static int nation_flag = 0;
static int n_cover = 5;

typedef struct
{
       int foffset;
       int fsindex;
} SIO;

void
RunSystemCommand(char *cmd)
{
    int ret;
    if ((ret = system(cmd)) != 0)
    {
        fprintf(stderr, "Error running command: %s\n", cmd);
        exit(EXIT_FAILURE);
    }
}


int main(int argc, char **argv)
{
    int i,nSeq,args;
    char *st,*ed;
    char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
    void ArraySort_String2(int n,char **Pair_Name,int *brr);
    fasta *seq;
    FILE *fp,*namef,*namef2;
    int size_file;
    int m_score,n_nodes,n_debug,num_sigma;
    void Memory_Allocate(int arr);
    char tempa[2000],tempc[2000],syscmd[2000],workdir[2000];
    char file_tarseq[2000],file_pileup[2000],file_snpfrq[2000],file_datas[2000],file_ctname[2000],file_ctspec[2000];
    char file_genome[2000],file_read2[2000],samname[500],bamname[500],toolname[500],snpname[500],gapname[500],datname[500];
    char countryname[100];
    char file_pileupSNP[2000],file_pileupGAP[2000],file_frequeSNP[2000],file_frequeGAP[2000],file_uniqueSNP[2000],file_regionSNP[2000];
    int systemRet = system (syscmd);
    int systemChd = chdir(tmpdir);
    pid_t pid;

    seq=NULL;
    
    if(argc < 2)
    {
         printf("Program: covidPileup -  Pileup pipeline for COVID-19 SNPs\n");
         printf("Version: 1.0\n");
         printf("\n");
         
         printf("Usage: %s -nodes 30 -length 40000 -cover 5 -country UK <Input_Reference> <Input_Covid-Genomes> <Output_Pileup-file>\n",argv[0]);
         printf("       nodes    (30)    - Number of CPUs requested\n");
         printf("       length   (40000) - Genome sequence length\n");
         printf("       cover    (5)     - Threshold coverage number to report specific SNPs \n");
         printf("       country  (UK)    - Specific SNPs in the country \n");
         printf("       SNP      (plot)  - Output image file on SNP pileup information\n");
         printf("       GAP      (plot)  - Output image file on GAP pileup information\n");
         exit(1);
    }

    m_score = 50;
    n_nodes = 20;
    n_debug = 1;

    strcpy(toolname,"bwa");
    strcpy(countryname,"UK");
    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-nodes"))
       {
         sscanf(argv[++i],"%d",&n_nodes);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-align"))
       {
         memset(toolname,'\0',500);
         sscanf(argv[++i],"%s",toolname);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-ancestry"))
       {
         sscanf(argv[++i],"%d",&ances_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-data"))
       {
         run_align = 0;
         sscanf(argv[++i],"%s",datname);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-SNP"))
       {
         sscanf(argv[++i],"%s",snpname);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-GAP"))
       {
         sscanf(argv[++i],"%s",gapname);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-country"))
       {
         sscanf(argv[++i],"%s",countryname);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-length"))
       {
         sscanf(argv[++i],"%d",&seq_len);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-cover"))
       {
         sscanf(argv[++i],"%d",&n_cover);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-help"))
       {
         printf("Program: covidPileup -  Pileup pipeline for COVID-19 SNPs\n");
         printf("Version: 1.0\n");
         printf("\n");
         
         printf("Usage: %s -nodes 30 -length 40000 -cover 5 -country UK <Input_Reference> <Input_Covid-Genomes> <Output_Pileup-file>\n",argv[0]);
         printf("       nodes    (30)    - Number of CPUs requested\n");
         printf("       length   (40000) - Genome sequence length\n");
         printf("       cover    (5)     - Threshold coverage number to report specific SNPs \n");
         printf("       country  (UK)    - Specific SNPs in the country \n");
         printf("       SNP      (plot)  - Output image file on SNP pileup information\n");
         printf("       GAP      (plot)  - Output image file on GAP pileup information\n");
         exit(1);
       }
       else if(!strcmp(argv[i],"-debug"))
       {
         sscanf(argv[++i],"%d",&n_debug);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-file"))
       {
         sscanf(argv[++i],"%d",&file_tag); 
         args=args+2;
       }
    }

    plot_SNP = 0;
    plot_GAP = 0;
    if(strncmp(snpname,"plot",4) == 0)
      plot_SNP = 1;
    if(strncmp(gapname,"plot",4) == 0)
      plot_GAP = 1;

    nation_flag = 1;
    if((strcmp(countryname,"UK")==0)||(strcmp(countryname,"USA")==0)||(strcmp(countryname,"EU")==0))
      nation_flag = 0;
 
    pid = getpid();
    memset(tempa,'\0',2000);
    if (!getcwd(tempa, sizeof(tempa)))
    {
      exit(1);
    } 
    memset(tmpdir,'\0',2000);
    memset(workdir,'\0',2000);
    sprintf(tmpdir,"%s/",tempa);
    sprintf(workdir,"%s/tmp_rununik_%d/",tempa,pid);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"mkdir %s",workdir);

    RunSystemCommand(syscmd);

    if(chdir(workdir) == -1)
    {
      printf("System command error: chdir\n");
      exit(EXIT_FAILURE);
    }
     
    st = argv[0];
    ed = strrchr(argv[0],'/');
    memset(tempc,'\0',2000);
    strncpy(tempc,argv[0],ed-st);
    memset(bindir,'\0',2000);
    sprintf(bindir,"%s/covid-bin",tempc);

    memset(file_tarseq,'\0',2000);
    memset(file_genome,'\0',2000);
    memset(file_pileup,'\0',2000);
    memset(file_snpfrq,'\0',2000);
    memset(file_ctspec,'\0',2000);
    memset(file_ctname,'\0',2000);
    memset(file_datas,'\0',2000);

    sprintf(file_tarseq,"%s/%s",tempa,argv[args]);
    sprintf(file_genome,"%s/%s",tempa,argv[args+1]);
    sprintf(file_pileup,"%s/%s",tempa,argv[args+2]);
    sprintf(file_snpfrq,"%s/%s.snps",tempa,argv[args+2]);
    sprintf(file_ctname,"%s/%s.name",tempa,argv[args+2]);
    sprintf(file_ctspec,"%s/%s.spec",tempa,argv[args+2]);
    sprintf(file_pileupSNP,"%s/%s.pileupSNP.png",tempa,argv[args+2]);
    sprintf(file_pileupGAP,"%s/%s.pileupGAP.png",tempa,argv[args+2]);
    sprintf(file_frequeSNP,"%s/%s.frequeSNP.png",tempa,argv[args+2]);
    sprintf(file_frequeGAP,"%s/%s.frequeGAP.png",tempa,argv[args+2]);
    sprintf(file_uniqueSNP,"%s/%s.uniqueSNP.png",tempa,argv[args+2]);
    sprintf(file_regionSNP,"%s/%s.regionSNP.png",tempa,argv[args+2]);

    if((namef = fopen(file_tarseq,"r")) == NULL)
    {
      printf("File not in the working directory!\n");
      if((namef = fopen(argv[args],"r")) == NULL)
      {
        printf("File %s not found and please copy it to your working directory!\n",argv[args]);
        exit(1);
      }
      else
      {
        memset(file_tarseq,'\0',2000);
        strcpy(file_tarseq,argv[args]);
        printf("Input target assembly file1: %s\n",file_tarseq);
      }
    }
    else
    {
      printf("Input target assembly file2: %s\n",file_tarseq);
    } 

    if((namef = fopen(file_genome,"r")) == NULL)
    {
      printf("File not in the working directory!\n");
      if((namef = fopen(argv[args+1],"r")) == NULL)
      {
        printf("File %s not found and please copy it to your working directory!\n",argv[args+1]);
        exit(1);
      }
      else
      {
        memset(file_genome,'\0',2000);
        strcpy(file_genome,argv[args+1]);
        printf("Input read1 file: %s\n",file_genome);
      }
    }
    else
    {
      printf("Input read1 file: %s\n",file_genome);
    }

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/covid_fastq -name tarseq -len 20000 %s tarseq.fastq tarseq.tag > try.out",bindir,file_tarseq);
    RunSystemCommand(syscmd);
   
    if(run_align)
    {
      if((strcmp(toolname,"bwa") == 0))
      { 
        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/bwa index tarseq.fastq > try.out",bindir);
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/bwa mem -t %d tarseq.fastq %s | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat",bindir,n_nodes,file_genome,"($2<100){print $1,$2,$3,$4,$5,$6,$10}");
        RunSystemCommand(syscmd);
      }
      else if((strcmp(toolname,"smalt") == 0))
      {
        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/smalt index -k 13 -k 13 hash_genome tarseq.fastq > try.out",bindir);
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/smalt map -n %d -m 100 -f samsoft -O hash_genome %s | egrep tarseq_ | awk '%s' | egrep -v '^@' > align.dat",bindir,n_nodes,file_genome,"($2<100){print $1,$2,$3,$4,$5,$6,$10}");
        RunSystemCommand(syscmd);
      }
      else
      {
        printf("Give an alignment tool! \n");
        exit(1);
      }
    }
    else
    {
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"cp %s align.dat",datname);
      RunSystemCommand(syscmd);
    } 

    memset(syscmd,'\0',2000);
    printf("%s/covidSNP -length %d tarseq.fastq align.dat pileup.dat > try.out",bindir,seq_len);
    sprintf(syscmd,"%s/covidSNP -length %d -ancestry %d tarseq.fastq align.dat pileup.out > try.out",bindir,seq_len,ances_flag);
    RunSystemCommand(syscmd);

/*  Process specific SNP patterns   */
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"cat align.dat | awk '{print $1}' | sort > name.dat");
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/covid_names name.dat country.dat > try.out",bindir);
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"egrep SNP pileup.out | sort -k 4,4n -k 5,5 > snp-sort.dat");
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/covid_nation -country %s -cover %d country.dat snp-sort.dat country.spec country.snps > try.out",bindir,countryname,n_cover);
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"egrep SNP pileup.out | sort -n -k 4 | awk '{print $2,$4}' > pileupSNP.dat");
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"egrep GAP pileup.out | sort -n -k 4 | awk '{print $2,$4}' > pileupGAP.dat");
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"egrep NUMsnp pileup.out | sort -n -k 5  > frequeSNP-all.dat");
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"egrep NUMsnp pileup.out | sort -n -k 5 | awk '{print $2,$5}' > frequeSNP.dat");
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"egrep NUMgap pileup.out | sort -n -k 5 | awk '{print $2,$5}' > frequeGAP.dat");
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/covid_frequ pileupSNP.dat sample-pileupSNP.freq  > try.out",bindir);
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/covid_frequ -gap 0 pileupGAP.dat sample-pileupGAP.freq  > try.out",bindir);
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/covid_frequ frequeSNP.dat sample-frequeSNP.freq  > try.out",bindir);
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/covid_frequ -gap 0 frequeGAP.dat sample-frequeGAP.freq  > try.out",bindir);
    RunSystemCommand(syscmd);
    
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"cp snp-sort.dat %s",file_snpfrq);
    RunSystemCommand(syscmd);

/*
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"cp frequeSNP-all.dat %s",file_snpfrq);
    RunSystemCommand(syscmd);
                                    */
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"cp country.dat %s",file_ctname);
    RunSystemCommand(syscmd);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"cp country.spec %s",file_ctspec);
    RunSystemCommand(syscmd);

    if(plot_SNP == 1)
    {
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"cp %s/public-pileupSNP.freq public-pileupSNP.freq",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/covid_comms -plot 1 frequeSNP.dat plot-pileupSNP.sh > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"bash plot-pileupSNP.sh");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"mv data.png %s ",file_pileupSNP);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"cp %s/public-frequeSNP.freq public-frequeSNP.freq",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/covid_comms -plot 2 frequeSNP.dat plot-frequeSNP.sh > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"bash plot-frequeSNP.sh");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"mv data.png %s ",file_frequeSNP);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/covid_nation -country USA -cover %d country.dat snp-sort.dat snp-usa-uniq.dat snp-usa-all.dat > try.out",bindir,n_cover);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/covid_nation -country UK -cover %d country.dat snp-sort.dat snp-uk-uniq.dat snp-uk-all.dat > try.out",bindir,n_cover);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/covid_nation -country EU -cover %d country.dat snp-sort.dat snp-eu-uniq.dat snp-eu-all.dat > try.out",bindir,n_cover);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"egrep Uniq_ snp-eu-uniq.dat | awk '{print $3,$5}' > uniqsnp-eu.dat");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"egrep All_ snp-eu-all.dat | awk '{print $3,$5}' > country-eu.dat");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"egrep Uniq_ snp-uk-uniq.dat | awk '{print $3,$5}' > uniqsnp-uk.dat");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"egrep All_ snp-uk-all.dat | awk '{print $3,$5}' > country-uk.dat");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"egrep Uniq_ snp-usa-uniq.dat | awk '{print $3,$5}' > uniqsnp-usa.dat");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"egrep All_ snp-usa-all.dat | awk '{print $3,$5}' > country-usa.dat");
      RunSystemCommand(syscmd);

      if(nation_flag == 0)
      {
        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/covid_comms -plot 5 frequeSNP.dat plot-uniqueSNP.sh > try.out",bindir);
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"bash plot-uniqueSNP.sh");
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"mv data.png %s ",file_uniqueSNP);
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/covid_comms -plot 6 frequeSNP.dat plot-regionSNP.sh > try.out",bindir);
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"bash plot-regionSNP.sh");
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"mv data.png %s ",file_regionSNP);
        RunSystemCommand(syscmd);
      }
      else
      {
        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/covid_nation -country %s -cover %d country.dat snp-sort.dat snp-sm-uniq.dat snp-sm-all.dat > try.out",bindir,countryname,n_cover);
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"egrep Uniq_ snp-sm-uniq.dat | awk '{print $3,$5}' > uniqsnp-sm.dat");
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"egrep All_ snp-sm-all.dat | awk '{print $3,$5}' > country-sm.dat");
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/covid_comms -plot 7 -country %s frequeSNP.dat plot-uniqueSNP.sh > try.out",bindir,countryname);
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"bash plot-uniqueSNP.sh");
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"mv data.png %s ",file_uniqueSNP);
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"%s/covid_comms -plot 8 -country %s frequeSNP.dat plot-regionSNP.sh > try.out",bindir,countryname);
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"bash plot-regionSNP.sh");
        RunSystemCommand(syscmd);

        memset(syscmd,'\0',2000);
        sprintf(syscmd,"mv data.png %s ",file_regionSNP);
        RunSystemCommand(syscmd);
      }
    }
    
    if(plot_GAP == 1)
    {
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"cp %s/public-pileupGAP.freq public-pileupGAP.freq",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/covid_comms -plot 3 frequeSNP.dat plot-pileupGAP.sh > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"bash plot-pileupGAP.sh");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"mv data.png %s ",file_pileupGAP);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"cp %s/public-frequeGAP.freq public-frequeGAP.freq",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/covid_comms -plot 4 frequeSNP.dat plot-frequeGAP.sh > try.out",bindir);
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"bash plot-frequeGAP.sh");
      RunSystemCommand(syscmd);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"mv data.png %s ",file_frequeGAP);
      RunSystemCommand(syscmd);
    }
    
    if(n_debug == 0)
    {
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf * > /dev/null");
    
      RunSystemCommand(syscmd);
      if(chdir(tmpdir) == -1)
      {
        printf("System command error: chdir\n");
        exit(EXIT_FAILURE);
      }

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf %s > /dev/null",workdir);
    
      RunSystemCommand(syscmd);
    }
    return EXIT_SUCCESS;
    return EXIT_FAILURE;

}
/* end of the main */



/*   subroutine to sort out read pairs    */
/* =============================== */
void Read_Index(int nSeq, char *namefile)
/* =============================== */
{
     int i,j,nseq;
     int i_reads,n_reads,c_reads,insertsize=0;
     FILE *namef;
     char *ptr;
     char line[500];
     char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);


     if((namef = fopen(namefile,"r")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }
     nseq = 0;
     i_reads = 0;
     n_reads = 0;
     c_reads = 0;
     while(!feof(namef))
     {
       if(fgets(line,500,namef) == NULL)
            printf("Data input file problem!\n");
       if(feof(namef)) break;
       nseq++;
     }
     fclose(namef);

     nseq = 2*nseq;
     if((ctg_head = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_head\n");
       exit(1);
     }
     if((ctg_list = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_list\n");
       exit(1);
     }
     if((core_list = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - core_list\n");
       exit(1);
     }
     if((read2contig = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - read2contig\n");
       exit(1);
     }

     nseq = nseq*3;
     S_Name=cmatrix(0,nseq,0,Max_N_NameBase);
     if((insert_siz = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - insert\n");
       exit(1);
     }
     if((insert_dev = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - insert\n");
       exit(1);
     }
 
     if((namef = fopen(namefile,"r")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }

     j = 0;
     insertsize = 0;
     while(!feof(namef))
     {
       int nPair=0,len=0;
       char line2[500],line3[500],base[500];

       if(fgets(line,500,namef) == NULL)
            printf("Data input file problem!\n");
       if(feof(namef)) break;
       strcpy(line2,line);
       strcpy(line3,line);
       insertsize = 0;
       if((strncmp(line,"readnames",9))==0)
       {
         i=0;
         c_reads = 0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==5)
            {
                memset(base,'\0',500);
//                len=strlen(ptr);
//                strncpy(base,ptr,len-1);
                strcat(base,ptr);
                c_reads = atoi(base);
            }
         }
//       printf("creads: %d %d\n",c_reads,n_reads);
         if(n_group>0)
           ctg_list[n_group-1]=n_reads;
         n_group++;
         n_reads = 0;
       }
       else
       {      
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
         i=0;
         for(ptr=strtok(line3," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(nPair>1)
            {
              if(i==(nPair-2))
              {
                memset(base,'\0',500);
                strcat(base,ptr);
                insertsize = atoi(base);
              }
            }
         }
         i=0;
         j=0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==0)
            {

              len=strlen(ptr);
              if(nPair==1)
                strncpy(S_Name[i_reads],ptr,len-1);
              else
                strncpy(S_Name[i_reads],ptr,len);
//       printf("reads: %d %d %d %d %s\n",j,i_reads,insertsize,c_reads,S_Name[i_reads]);
              i_reads++;
              j++;
            }
//            else if(insertsize<50000)
            else if((insertsize<50000)&&(c_reads>4))
            {
              if(i==1)
              {
                len=strlen(ptr);
                strncpy(S_Name[i_reads],ptr,len);
                i_reads++;
                j++;
              }
              else if((i==2)&&(i<(nPair-2)))
              {
                len=strlen(ptr);
                strncpy(S_Name[i_reads],ptr,len);
                i_reads++;
                j++;
              }
/*              else if(i==(nPair-2))
              {
                memset(base,'\0',500);
                strcat(base,ptr);
                for(k=0;k<j;k++)
                   insert_siz[i_reads-k-1] = atoi(base);
              }
              else if(i==(nPair-1))
              {
                memset(base,'\0',500);
                strcat(base,ptr);
                for(k=0;k<j;k++)
                   insert_dev[i_reads-k-1] = atoi(base);
              }   */
            }
         }
         n_reads = n_reads+j;
       }
       c_reads++;
     }
     fclose(namef);
     printf("contig: %d %d\n",n_reads,i_reads);
     ctg_list[n_group-1]=n_reads;
     ctg_head[0]=0;
     for(i=1;i<n_group;i++)
        ctg_head[i]=ctg_head[i-1]+ctg_list[i-1];
      
}

#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  

/* =============================== */
void ArraySort_Long(int n, B64_long *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* =============================== */
void ArraySort_Int(int n, int *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* =============================== */
void ArraySort_Mix(int n, B64_long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/*   function to sort an array into a decreasing order:  a>b>c>....    */  
/* =============================== */
void ArraySort2_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]>=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]<arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]<arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]<arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]>a);
             do j--; while (arr[j]<a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Mix3(int n, B64_long *arr, int *brr, int *crr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,c,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             c=crr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
                crr[i+1]=crr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
             crr[i+1]=c;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);
          SWAP(crr[k],crr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
            SWAP(crr[m],crr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
            SWAP(crr[m+1],crr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
            SWAP(crr[m],crr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          c=crr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
             SWAP(crr[i],crr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          crr[m+1]=crr[j];
          crr[j]=c;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/*   to swap the string arrays           */
/* ============================================= */
void s_swap(char **Pair_Name, int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String(int n, char **Pair_Name, int *brr)
/* ============================================= */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int temp,MIN=7;
     char p[Max_N_NameBase];

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             strcpy(p,Pair_Name[j]);
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(strcmp(Pair_Name[i],p)<=0) break;
                strcpy(Pair_Name[i+1],Pair_Name[i]);
                brr[i+1]=brr[i];
             }
             strcpy(Pair_Name[i+1],p);
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          s_swap(Pair_Name,k,m+1);
          SWAP(brr[k],brr[m+1]);

          if(strcmp(Pair_Name[m],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m,ir);
            SWAP(brr[m],brr[ir]);
          }

          if(strcmp(Pair_Name[m+1],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m+1,ir);
            SWAP(brr[m+1],brr[ir]);
          }

          if(strcmp(Pair_Name[m],Pair_Name[m+1])>0)
          {
            s_swap(Pair_Name,m,m+1);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          strcpy(p,Pair_Name[m+1]);
          b=brr[m+1];
          for(;;)
          {
             do i++; while (strcmp(Pair_Name[i],p)<0);
             do j--; while (strcmp(Pair_Name[j],p)>0);
             if(j<i) break;
             s_swap(Pair_Name,i,j);
             SWAP(brr[i],brr[j]);
          }
          strcpy(Pair_Name[m+1],Pair_Name[j]);
          strcpy(Pair_Name[j],p);
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/*   to swap the string arrays           */
/* ============================================= */
void s_swap2(char **Pair_Name, int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase2];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String2(int n, char **Pair_Name, int *brr)
/* ============================================= */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int temp,MIN=7;
     char p[Max_N_NameBase2];

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             strcpy(p,Pair_Name[j]);
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(strcmp(Pair_Name[i],p)<=0) break;
                strcpy(Pair_Name[i+1],Pair_Name[i]);
                brr[i+1]=brr[i];
             }
             strcpy(Pair_Name[i+1],p);
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          s_swap2(Pair_Name,k,m+1);
          SWAP(brr[k],brr[m+1]);

          if(strcmp(Pair_Name[m],Pair_Name[ir])>0)
          {
            s_swap2(Pair_Name,m,ir);
            SWAP(brr[m],brr[ir]);
          }

          if(strcmp(Pair_Name[m+1],Pair_Name[ir])>0)
          {
            s_swap2(Pair_Name,m+1,ir);
            SWAP(brr[m+1],brr[ir]);
          }

          if(strcmp(Pair_Name[m],Pair_Name[m+1])>0)
          {
            s_swap2(Pair_Name,m,m+1);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          strcpy(p,Pair_Name[m+1]);
          b=brr[m+1];
          for(;;)
          {
             do i++; while (strcmp(Pair_Name[i],p)<0);
             do j--; while (strcmp(Pair_Name[j],p)>0);
             if(j<i) break;
             s_swap2(Pair_Name,i,j);
             SWAP(brr[i],brr[j]);
          }
          strcpy(Pair_Name[m+1],Pair_Name[j]);
          strcpy(Pair_Name[j],p);
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int  **m;

        /* allocate pointers to rows        */
        if((m=(int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(int *)calloc(nrow*ncol,sizeof(int)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}

/* creat char matrix with subscript ange cm[nrl...nrh][ncl...nch]  */
char    **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char **cm;

        /* allocate pointers to rows        */
        if((cm=(char **)calloc(nrow,sizeof(char*)))==NULL)
        {
           printf("error cmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        cm+=0;
        cm-=nrl;

        /* allocate rows and set pointers to them        */
        if((cm[nrl]=(char *)calloc(nrow*ncol,sizeof(char)))==NULL)
        {
           printf("error cmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        cm[nrl]+=0;
        cm[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           cm[i]=cm[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return cm;
}

