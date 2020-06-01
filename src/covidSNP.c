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
#include "fasta.h"

#define GT '>'
#define GT4 (((((GT<<8)+GT)<<8)+GT)<<8)+GT

#define ENDS_EXTRA 0
#define PADCHAR '-'
#define MAX_N_BRG 50000 
#define MAX_N_ROW 40000 
#define nfm 800000
#define nfm_sub 500000
#define Max_N_NameBase 60
#define Max_N_Pair 100
static long *h_dna;

/* SSAS default parameters   */
static int IMOD=0;
static int len_covid=40000;
static int set_qual=15;
static int set_len=10;
static int ances_flag =0;

int main(int argc, char **argv)
{
    FILE *fp,*namef,*namef2;
    long dataSize,totalBases;
    int i,j,k,m,n,nSeq,nRead,args,num_SNPs=0,num_GAPs=0;
    int proce_flag;
    int nseq,num_align,*ctg_index,*ctg_hitst,*cig_base,*cig_code,*snp_offset;
    fasta *seq,*seqp;
    char *ptr,**ctgname,line[100000],cigarline[10000],*seqline,*read_base,*refe_base;
    char score1[60],score2[60];
    char **cmatrix(long nrl,long nrh,long ncl,long nch);
    fasta *segg;
    long Size_pdata,Size_q_pdata;
    unsigned int *pidata;
    int num_seqque;
    char *pdata,*st;

    seq=NULL;
    fflush(stdout);
    if(system("ps aux | grep covidSNP; date") == -1)
    {
//        printf("System command error:\n);
    }

    ances_flag = 0;
    proce_flag = 1;
    if(argc < 2)
    {
      printf("Usage: %s <-length 40000> <-ancestry 0> <input_reference_fasta> <alignment_file> <output_SNP_file>\n",argv[0]);
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
       else if(!strcmp(argv[i],"-ancestry"))
       {
         sscanf(argv[++i],"%d",&ances_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-length"))
       {
         sscanf(argv[++i],"%d",&len_covid);
         args=args+2;
       }
    }

    if((fp=fopen(argv[args],"rb"))==NULL) error("Cannot open file\n");
    fseek(fp, 0, SEEK_END);
    Size_q_pdata = ftell(fp) + 1;
    fclose(fp);
    if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL)
      error("calloc pdata\n");
    num_seqque = extractFastq(argv[args],pdata,Size_q_pdata);
    if((segg=(fasta*)calloc((num_seqque+1),sizeof(fasta)))==NULL)
      error("calloc segg\n");
    if((seq=decodeFastq(argv[args],&num_seqque,&totalBases,pdata,Size_q_pdata,segg))==NULL)
      error("no query data found.\n");
    nseq=0;
    nSeq = num_seqque;
    printf("Number of shotgun reads  %d \n",nSeq);

    if(totalBases == 0)
      return EXIT_SUCCESS;

/*  input read alignment info line   */
    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: reads group file: %s \n",argv[args+1]);
      exit(1);
    }

    nRead = 0;
    while(!feof(namef))
    {
      if(fgets(line,100000,namef) == NULL)
        printf("Data input file problem! %s\n",argv[args+1]);
      if(feof(namef)) break;
      nRead++;
    }
    fclose(namef); 

    printf("Number of genomes  %d %d\n",nSeq,nRead);

    ctgname=cmatrix(0,nRead,0,500);
    if((cig_base= (int *)calloc(10000,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - ctg_index\n");
      exit(1);
    }
    if((cig_code= (int *)calloc(10000,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - ctg_index\n");
      exit(1);
    }
    if((snp_offset= (int *)calloc(len_covid,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - ctg_index\n");
      exit(1);
    }
    if((seqline= (char *)calloc(len_covid,sizeof(char))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - seqline\n");
      exit(1);
    }
    if((read_base= (char *)calloc(len_covid,sizeof(char))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - read_base\n");
      exit(1);
    }
    if((refe_base= (char *)calloc(len_covid,sizeof(char))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - read_base\n");
      exit(1);
    }
    if((ctg_index= (int *)calloc(nRead,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - ctg_index\n");
      exit(1);
    }
    if((ctg_hitst= (int *)calloc(nRead,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - ctg_index\n");
      exit(1);
    }

    dataSize=1000;
/*  process contigs one by one   */
    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    if((namef2 = fopen(argv[args+2],"w")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the SNP output file         */
    num_align=0;
    while(!feof(namef))
    {
      int nPair=0,seq_len,hitst;
      char line2[100000],line3[100000],base[500],cbase[20];

      if(fgets(line,100000,namef) == NULL)
        printf("Data input file problem2!\n");
      if(feof(namef)) break;
      strcpy(line2,line);
      strcpy(line3,line);
      proce_flag = 1;
      if(ances_flag == 0)
      {
        if(((strncmp(line,"bat",3))==0)||((strncmp(line,"pan",3))==0))
          proce_flag = 0;
        else
          proce_flag = 1; 
      }
      else
        proce_flag = 1; 
      if(proce_flag)
      { 
        for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
        {
        }
        i=0;
        for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
        {
           if(i==0)
           {
             memset(base,'\0',500);
             strcat(base,ptr);
             strcpy(ctgname[num_align],base);
           }
           else if(i==3)
           {
             memset(base,'\0',500);
             strcat(base,ptr);
             ctg_hitst[num_align] = atoi(base);
             hitst = atoi(base)-1;
           }
           else if(i==5)
           {
             int c_len,r_len,stopflag; 
             int num_hits,offset,r_offset,s_offset; 
             memset(cigarline,'\0',10000);
             strcat(cigarline,ptr);
           }
           else if(i==6)
           {
             int c_len,r_len,stopflag; 
             int num_hits,offset,r_offset,s_offset; 
             int gbase = 0;

             memset(seqline,'\0',len_covid);
             memset(snp_offset,'\0',4*len_covid);
             strcat(seqline,ptr);
             seq_len = strlen(seqline);
             st = cigarline;
             c_len = strlen(cigarline);
	     memset(cig_code,0,4*c_len);
	     memset(cig_base,0,4*c_len);
             c_len = strlen(cigarline)-1;
             for(k=0;k<c_len;k++)
             {
		if(cigarline[k] == 'M')
		  cig_code[k] = 1;
		else if(cigarline[k] == 'I')
		  cig_code[k] = 2;
		else if(cigarline[k] == 'D')
		  cig_code[k] = 3;
		else if(cigarline[k] == 'S')
		  cig_code[k] = 4;
		else if(cigarline[k] == 'H')
		  cig_code[k] = 5;
		else
		  cig_code[k] = 0;
             }
             num_hits = 0;
             offset = 0;
             r_offset = hitst;
             s_offset = 0;
	     for(k=0;k<c_len;k++)
             {
		stopflag = 0;
		j = k+1;
		while((j < c_len)&&(stopflag == 0))
                {
                   if(cig_code[j]==0)
                   {
                     j++;
                   }
                   else
                     stopflag=1;
                }
                memset(cbase,'\0',20);
                if(num_hits == 0)
                {
                  for(m=k;m<j;m++)
		    cbase[m-k] = cigarline[m];
                }
                else
                {
                  for(m=(k+1);m<j;m++)
		    cbase[m-k-1] = cigarline[m];
                } 
                r_len = atoi(cbase);
//                printf("yyy: %d %d %d %s %s %c || %d %d %d\n",k,c_len,r_len,cbase,cigarline,cigarline[j],offset,r_offset,s_offset);
                if(cigarline[j] == 'S')
                {
                  s_offset = r_len;
                }
                else if(cigarline[j] == 'M')
                {
//                printf("yyy2: %d %d %d %s %c %s || %d %d %d\n",j,c_len,r_len,cbase,cigarline[j],seqline,offset,r_offset,s_offset);
                  for(n=0;n<r_len;n++)
                  {
                     refe_base[offset+n] = seq->data[r_offset+n]; 
                     read_base[offset+n] = seqline[s_offset+n];
                     snp_offset[offset+n] = r_offset+n+1; 
//                  printf("Base0: %d %d %d %c %c %c\n",num_align,n,offset,refe_base[n],read_base[n],seqline[n]);
                  }
                  offset = offset + r_len;
                  r_offset = r_offset + r_len;
                  s_offset = s_offset + r_len;
                }
                else if(cigarline[j] == 'I')
                {
                  for(n=0;n<r_len;n++)
                  {
                     refe_base[offset+n] = '-'; 
                     read_base[offset+n] = seqline[s_offset+n]; 
                     snp_offset[offset+n] = offset+1; 
             if(num_align == 167)
                  printf("Base: %d %d %d %c %c\n",n,r_len,offset+n,refe_base[n+offset],read_base[n+offset]);
                  }
                  s_offset = s_offset+r_len;
                  offset = offset + r_len;
                }
                else if(cigarline[j] == 'D')
                {
                  for(n=0;n<r_len;n++)
                  {
                     refe_base[offset+n] = seq->data[r_offset+n]; 
                     read_base[offset+n] = '-'; 
                     snp_offset[offset+n] = r_offset+n+1; 
//             if(num_align == 66)
//                  printf("Base: %d %d %c %c\n",n,offset+n,refe_base[n+offset],read_base[n+offset]);
                  }
                  offset = offset + r_len;
                  r_offset = r_offset+r_len;
                  
                }
                num_hits++;
                k = j-1;
             }
         
             num_SNPs = 0; 
             num_GAPs = 0; 
             if(seq_len > 20000)
               printf("Seq %d %s %s %s\n",num_align,seq->name,ctgname[num_align],cigarline);
             for(n=0;n<offset;n++)
             {
                gbase = 0;
                if((read_base[n] == 'A')||(read_base[n] == 'C')||(read_base[n] == 'G')||(read_base[n] == 'T'))
                  gbase = 1;
                if(read_base[n] == 'N')
                  gbase = 2;
                if(read_base[n] == '-')
                  gbase = 0;
                if(refe_base[n] == '-')
                  gbase = 0;
//	          printf("    PPP %d %d %d %s %c %c %d\n",num_align,n,snp_offset[n],ctgname[num_align],refe_base[n],read_base[n],ctg_hitst[num_align]);
                if((refe_base[n]!=read_base[n])&&(gbase == 1)&&(seq_len > 20000)&&(snp_offset[n] > ctg_hitst[num_align]))
                {
	          printf("    SNP %d %d %d %s %c %c %d\n",num_align,n,snp_offset[n],ctgname[num_align],refe_base[n],read_base[n],ctg_hitst[num_align]);
	          fprintf(namef2,"SNP %d %d %d %s %c %c\n",num_align,n,snp_offset[n],ctgname[num_align],refe_base[n],read_base[n]);
                  num_SNPs++; 
                }
                if(gbase == 2)
                {
	          printf("GAP %d %d %d %s %c %c %d\n",num_align,n,snp_offset[n],ctgname[num_align],refe_base[n],read_base[n],ctg_hitst[num_align]);
	          fprintf(namef2,"GAP %d %d %d %s %c %c\n",num_align,n,snp_offset[n],ctgname[num_align],refe_base[n],read_base[n]);
                  num_GAPs++; 
                }
             }
             if(seq_len > 20000)
	     {
               printf("NUMsnp %d %s %s %d\n",num_align,seq->name,ctgname[num_align],num_SNPs);
               fprintf(namef2,"NUMsnp %d %s %s %d\n",num_align,seq->name,ctgname[num_align],num_SNPs);
               printf("NUMgap %d %s %s %d\n",num_align,seq->name,ctgname[num_align],num_GAPs);
               fprintf(namef2,"NUMgap %d %s %s %d\n",num_align,seq->name,ctgname[num_align],num_GAPs);
               num_align++;
	     }	     
           }
        }
      }
    }
    fclose(namef);
    fclose(namef2);


    if(seq){
        free(seq->name);
        free(seq);
        seq = NULL;
    }    
    printf("Job finished for %d contigs!\n",nSeq);
    return EXIT_SUCCESS;

}
/* end of the main */

#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  

/* =============================== */
void ArraySort_Long(int n, long *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

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
void ArraySort_Mix(int n, long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

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
void ArraySort_Mix3(int n, long *arr, int *brr, int *crr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,c,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

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
void s_swap(char Pair_Name[][Max_N_NameBase], int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String(int n, char Pair_Name[][Max_N_NameBase], int *brr)
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


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **imatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
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
char    **cmatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
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

