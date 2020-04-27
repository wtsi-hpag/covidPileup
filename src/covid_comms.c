/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2020  Genome Research Ltd.                                *
 *                                                                          *
 *  Author: Zemin Ning (zn1@sanger.ac.uk)                                   *
 *                                                                          *
 *  This file is part of covidPileup pipeline.                              *
 *                                                                          *
 *  Scaff10x is a free software: you can redistribute it and/or modify it   *
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

#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>

static int IMOD=0;
static int plot_flag = 1;

int main(int argc, char **argv)
{
    int i=0,j=0,k,args=0,num_steps,nSeq,nRead,num_samples;
    int *s_len,BAR = 0,nstep = 0,stopflag;
    char line[500],tempc1[100];
    char KKK1[100],KKK2[100],KKK3[100],KKK4[100],KKK5[100],KKK6[100],KKK0[100];
    FILE *namef;
    long num_base,base;
    double rate;

    if(argc < 2)
    {
      printf("Usage: %s -gap 0 <input frequence file> <output frequency file>\n",argv[0]);
      exit(1);
    }

    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-mod"))
       {
         sscanf(argv[++i],"%d",&IMOD); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-plot"))
       {
         sscanf(argv[++i],"%d",&plot_flag);
         args=args+2;
       }
    }

    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file: %s \n",argv[args]);
      exit(1);
    }

    nRead = 0;
    while(!feof(namef))
    {
      if(fgets(line,500,namef) == NULL)
        printf("Data input file problem! %s\n",argv[args]);
      if(feof(namef)) break;
      nRead++;
    }
    fclose(namef); 

    nSeq = nRead;
    if((namef = fopen(argv[args+1],"w")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

    if(plot_flag == 1)
    {
      strcpy(KKK0,"set logscale y"); 
      strcpy(KKK1,"set terminal svg");     
      strcpy(KKK2,"set style line 1 lt 1 lw 2 pt 2 linecolor rgb \\\"red\\\"");     
      strcpy(KKK3,"set style line 2 lt 1 lw 2 pt 2 linecolor rgb \\\"blue\\\"");     
      strcpy(KKK4,"set xlabel \\\"COVID-19 genome coordinates\\\"");     
      strcpy(KKK5,"set ylabel \\\"SNP Pileup Frequency / 100 \\\"");     
      strcpy(KKK6,"plot [ 0 to 30000 ] [ 1 to 100 ]");
 
      fprintf(namef,"#!/bin/tcsh\n");
      fprintf(namef,"\n"); 
      fprintf(namef,"function plotcmd\n");
      fprintf(namef,"\{\n"); 
      fprintf(namef,"printf \"%s\\n\"\n",KKK0);
      fprintf(namef,"printf \"%s\\n\"\n",KKK1);
      fprintf(namef,"printf \"%s\\n\"\n",KKK2);
      fprintf(namef,"printf \"%s\\n\"\n",KKK3);
      fprintf(namef,"printf \"%s\\n\"\n",KKK4);
      fprintf(namef,"printf \"%s\\n\"\n",KKK5);
      fprintf(namef,"printf \"%s \\\"%s\\\" title \\\"%s \\\" %s\\\"%s\\\" title \\\"%s %d %s \\\" %s\" \n",KKK6,"public-pileupSNP.freq","COVID-19 Public: 368 genomes","with lines ls 1, ","sample-pileupSNP.freq","COVID-19 Sample: ",nSeq," genomes "," with lines ls 2");
      fprintf(namef,"}\n");
      fprintf(namef,"plotcmd | gnuplot > data.svg\n");
      fprintf(namef,"inkscape -z --export-text-to-path --export-pdf data.pdf data.svg\n");
      fprintf(namef,"gs -r600 -dNOPAUSE -dBATCH -sDEVICE=png256 -sOutputFile=%s.png data.pdf\n","data");
      fprintf(namef,"\n"); 
      fprintf(namef,"\n");
    }
    else if(plot_flag == 2)
    {
      strcpy(KKK0,"set logscale y"); 
      strcpy(KKK1,"set terminal svg");     
      strcpy(KKK2,"set style line 1 lt 1 lw 2 pt 2 linecolor rgb \\\"red\\\"");     
      strcpy(KKK3,"set style line 2 lt 1 lw 2 pt 2 linecolor rgb \\\"blue\\\"");     
      strcpy(KKK4,"set xlabel \\\"Number of SNPs\\\"");     
      strcpy(KKK5,"set ylabel \\\"SNP Frequency / 100\\\"");     
      strcpy(KKK6,"plot [ 0 to 20 ] [ 1 to 100 ]");
 
      fprintf(namef,"#!/bin/tcsh\n");
      fprintf(namef,"\n"); 
      fprintf(namef,"function plotcmd\n");
      fprintf(namef,"\{\n"); 
      fprintf(namef,"printf \"%s\\n\"\n",KKK0);
      fprintf(namef,"printf \"%s\\n\"\n",KKK1);
      fprintf(namef,"printf \"%s\\n\"\n",KKK2);
      fprintf(namef,"printf \"%s\\n\"\n",KKK3);
      fprintf(namef,"printf \"%s\\n\"\n",KKK4);
      fprintf(namef,"printf \"%s\\n\"\n",KKK5);
      fprintf(namef,"printf \"%s \\\"%s\\\" title \\\"%s \\\" %s\\\"%s\\\" title \\\"%s %d %s \\\" %s\" \n",KKK6,"public-frequeSNP.freq","COVID-19 Public: 368 genomes","with lines ls 1, ","sample-frequeSNP.freq","COVID-19 Sample: ",nSeq," genomes "," with lines ls 2");
      fprintf(namef,"}\n");
      fprintf(namef,"plotcmd | gnuplot > data.svg\n");
      fprintf(namef,"inkscape -z --export-text-to-path --export-pdf data.pdf data.svg\n");
      fprintf(namef,"gs -r600 -dNOPAUSE -dBATCH -sDEVICE=png256 -sOutputFile=%s.png data.pdf\n","data");
      fprintf(namef,"\n"); 
      fprintf(namef,"\n");
    } 
    else if(plot_flag == 3)
    {
      strcpy(KKK0,"set logscale y"); 
      strcpy(KKK1,"set terminal svg");     
      strcpy(KKK2,"set style line 1 lt 1 lw 2 pt 2 linecolor rgb \\\"red\\\"");     
      strcpy(KKK3,"set style line 2 lt 1 lw 2 pt 2 linecolor rgb \\\"blue\\\"");     
      strcpy(KKK4,"set xlabel \\\"COVID-19 genome coordinates\\\"");     
      strcpy(KKK5,"set ylabel \\\"GAP Pileup Frequency / 100 \\\"");     
      strcpy(KKK6,"plot [ 0 to 30000 ] [ 1 to 100 ]");
 
      fprintf(namef,"#!/bin/tcsh\n");
      fprintf(namef,"\n"); 
      fprintf(namef,"function plotcmd\n");
      fprintf(namef,"\{\n"); 
      fprintf(namef,"printf \"%s\\n\"\n",KKK0);
      fprintf(namef,"printf \"%s\\n\"\n",KKK1);
      fprintf(namef,"printf \"%s\\n\"\n",KKK2);
      fprintf(namef,"printf \"%s\\n\"\n",KKK3);
      fprintf(namef,"printf \"%s\\n\"\n",KKK4);
      fprintf(namef,"printf \"%s\\n\"\n",KKK5);
      fprintf(namef,"printf \"%s \\\"%s\\\" title \\\"%s \\\" %s\\\"%s\\\" title \\\"%s %d %s \\\" %s\" \n",KKK6,"public-pileupGAP.freq","COVID-19 Public: 368 genomes","with lines ls 1, ","sample-pileupGAP.freq","COVID-19 Sample: ",nSeq," genomes "," with lines ls 2");
      fprintf(namef,"}\n");
      fprintf(namef,"plotcmd | gnuplot > data.svg\n");
      fprintf(namef,"inkscape -z --export-text-to-path --export-pdf data.pdf data.svg\n");
      fprintf(namef,"gs -r600 -dNOPAUSE -dBATCH -sDEVICE=png256 -sOutputFile=%s.png data.pdf\n","data");
      fprintf(namef,"\n"); 
      fprintf(namef,"\n");
    } 
    else if(plot_flag == 4)
    {
      strcpy(KKK0,"set logscale y"); 
      strcpy(KKK1,"set terminal svg");     
      strcpy(KKK2,"set style line 1 lt 1 lw 2 pt 2 linecolor rgb \\\"red\\\"");     
      strcpy(KKK3,"set style line 2 lt 1 lw 2 pt 2 linecolor rgb \\\"blue\\\"");     
      strcpy(KKK4,"set xlabel \\\"Number of GAPs\\\"");     
      strcpy(KKK5,"set ylabel \\\"GAP Frequency / 100\\\"");     
      strcpy(KKK6,"plot [ 0 to 20 ] [ 1 to 100 ]");
 
      fprintf(namef,"#!/bin/tcsh\n");
      fprintf(namef,"\n"); 
      fprintf(namef,"function plotcmd\n");
      fprintf(namef,"\{\n"); 
      fprintf(namef,"printf \"%s\\n\"\n",KKK0);
      fprintf(namef,"printf \"%s\\n\"\n",KKK1);
      fprintf(namef,"printf \"%s\\n\"\n",KKK2);
      fprintf(namef,"printf \"%s\\n\"\n",KKK3);
      fprintf(namef,"printf \"%s\\n\"\n",KKK4);
      fprintf(namef,"printf \"%s\\n\"\n",KKK5);
      fprintf(namef,"printf \"%s \\\"%s\\\" title \\\"%s \\\" %s\\\"%s\\\" title \\\"%s %d %s \\\" %s\" \n",KKK6,"public-frequeGAP.freq","COVID-19 Public: 368 genomes","with lines ls 1, ","sample-frequeGAP.freq","COVID-19 Sample: ",nSeq," genomes "," with lines ls 2");
      fprintf(namef,"}\n");
      fprintf(namef,"plotcmd | gnuplot > data.svg\n");
      fprintf(namef,"inkscape -z --export-text-to-path --export-pdf data.pdf data.svg\n");
      fprintf(namef,"gs -r600 -dNOPAUSE -dBATCH -sDEVICE=png256 -sOutputFile=%s.png data.pdf\n","data");
      fprintf(namef,"\n"); 
      fprintf(namef,"\n");
    }
    else
    {
      printf("Wrong plot index!\n");
      exit(1);
    } 
    fclose(namef);
    return EXIT_SUCCESS;  
}


