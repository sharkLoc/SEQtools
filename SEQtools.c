#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

#include "func.h"

void usage(int lowq, int meanQ, double qRate, double nR)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: SEQtools [options]\n\n");
	fprintf(stderr, "\t-i\t:<char>   fq1 file\n");
	fprintf(stderr, "\t-I\t:<char>   fq2 file\n");
	fprintf(stderr, "\t-l\t:<int>    low quality threshold (default: [%d])\n",lowq);
	fprintf(stderr, "\t-m\t:<int>    filter reads with low average quality less than [%d]\n",meanQ);
	fprintf(stderr, "\t-q\t:<float>  low quality rate (default: [%.2f])\n",qRate);
	fprintf(stderr, "\t-n\t:<float>  N rate threshold (default: [%.2f])\n",nR);
	fprintf(stderr, "\t-o\t:<char>   clean fq1 file name\n");
	fprintf(stderr, "\t-O\t:<char>   clean fq2 file name\n");
	fprintf(stderr, "\t-h\t:<char>   show this help\n");
	fprintf(stderr, "\n");
	return ;
}

int main(int argc, char *argv[])
{
	gzFile fp, fd, fo1,fo2;
	int result, lowq=10, meanQ=30; double qRate = 0.5, nR = 0.05;

	while((result = getopt(argc, argv, "i:I:l:m:q:n:o:O:h")) != -1)
	{
		switch(result)
		{
			case 'i': fp = gzopen(optarg,"r"); break;
			case 'I': fd = gzopen(optarg,"r"); break;
			case 'l': lowq = atoi(optarg); break;
			case 'm': meanQ = atoi(optarg); break;
			case 'q': qRate = atof(optarg); break;
			case 'n': nR = atof(optarg); break;
			case 'o': fo1 = gzopen(optarg,"w"); break;
			case 'O': fo2 = gzopen(optarg,"w"); break;
			case 'h': usage(lowq,meanQ,qRate,nR); break;
			case '?': usage(lowq,meanQ,qRate,nR); break;
		}
	}
	if(argc==1) usage(lowq,meanQ,qRate,nR);
	
	if(!fp || !fd || !fo1 || !fo2)  exit(1);
	stat stat_1={0.,0.,0,0,0,0,0,0,0,0,0,0}; stat stat_2={0.,0.,0,0,0,0,0,0,0,0,0,0}; int flag=1;

	while(flag)
	{
		rinfo read1={{NULL,NULL,NULL,NULL},0,0,0,0.,0.};
		rinfo read2={{NULL,NULL,NULL,NULL},0,0,0,0.,0.};
		for(int i=0; i<4; i++)
		{
			 char *r1=readline(fp);
			 char *r2=readline(fd);
			 if(r1 && r2)
			 {
			 	size_t len1 = strlen(r1);
				size_t len2 = strlen(r2);
				read1.line[i]=(char *)malloc(len1+1);
				read2.line[i]=(char *)malloc(len2+1);
				strncpy(read1.line[i],r1,len1+1);
				strncpy(read2.line[i],r2,len2+1);
				if(i==1){read1.len=len1; read2.len=len2;}
			 }
			 else { flag=0; }
			 free(r1); free(r2);
		}
		if(flag)
		{
			for(int i=0; i<4; i++)
			{
				if(i==1)
				{
					int k=0; stat_1.readCount++;
					while(read1.line[i][k] != '\0')
					{
						if(read1.line[i][k] == 'A') stat_1.aCount++;
						if(read1.line[i][k] == 'T') stat_1.tCount++;
						if(read1.line[i][k] == 'G') stat_1.gCount++;
						if(read1.line[i][k] == 'C') stat_1.cCount++;
						if(read1.line[i][k] == 'N') {stat_1.nCount++; read1.nCount++;}
						k++; stat_1.baseCount++;
					}
					read1.nRate = nRate_calc(read1.len, read1.nCount);

					int s=0; stat_2.readCount++;
					while(read2.line[i][s] != '\0')
					{
						if(read2.line[i][s] == 'A') stat_2.aCount++;
						if(read2.line[i][s] == 'T') stat_2.tCount++;
						if(read2.line[i][s] == 'G') stat_2.gCount++;
						if(read2.line[i][s] == 'C') stat_2.cCount++;
						if(read2.line[i][s] == 'N') {stat_2.nCount++; read2.nCount++;}
						s++; stat_2.baseCount++;
					}
					read2.nRate = nRate_calc(read2.len, read2.nCount);
				}
				if(i==3)
				{
					read1.mean_qval = meanQvalue(read1.len,read1.line[i]);
					read1.lowQval = lowQvrate(lowq,read1.line[i]);
					stat_1.q20 += quality20(read1.line[i]);
					stat_1.q30 += quality30(read1.line[i]);
					read2.mean_qval = meanQvalue(read2.len,read2.line[i]);
					read2.lowQval = lowQvrate(lowq,read2.line[i]);
					stat_2.q20 +=quality20(read2.line[i]);
					stat_2.q30 +=quality30(read2.line[i]);
				}
				//fprintf(stdout,"%s\t%zu\t%zu\t%zu\t%f\t%f\n",read1.line[i],read1.len,read1.nCount,read1.mean_qval,read1.lowQval,read1.nRate);
				//fprintf(stdout,"%s\t%zu\t%zu\t%zu\t%f\t%f\n",read2.line[i],read2.len,read2.nCount,read2.mean_qval,read2.lowQval,read2.nRate);
			}
			
			if(read1.lowQval<qRate && read2.lowQval<qRate && read1.nRate<nR && read2.nRate<nR)
			{
				for(int i=0; i<4; i++)
				{
					gzprintf(fo1,"%s\n",read1.line[i]);
					gzprintf(fo2,"%s\n",read2.line[i]);
				}
			}
			else
			{
				stat_1.filter++;
				stat_2.filter++;
			}
			freeRead(read1); freeRead(read2);
		}
	}
	gzclose(fp); gzclose(fd);
	gzclose(fo1); gzclose(fo2);
	
	stat_1.averageLen = (float)stat_1.baseCount / stat_1.readCount;
	stat_2.averageLen = (float)stat_2.baseCount / stat_2.readCount;
	stat_1.gc = (float)(stat_1.gCount + stat_1.cCount) / stat_1.baseCount * 100;
	stat_2.gc = (float)(stat_2.gCount + stat_2.cCount) / stat_2.baseCount * 100;

	fprintf(stdout,"Iterm\t%s\t%s\n",argv[1],argv[2]);
	fprintf(stdout,"read average length:\t%d\t%d\n",(int)stat_1.averageLen,(int)stat_2.averageLen);
	fprintf(stdout,"read GC content(%%):\t%.2f\t%.2f\n",stat_1.gc,stat_2.gc);
	fprintf(stdout,"total read Count:\t%lu\t%lu\n",stat_1.readCount,stat_2.readCount);
	fprintf(stdout,"total base Count:\t%lu\t%lu\n",stat_1.baseCount,stat_2.baseCount);
	fprintf(stdout,"\n");
	fprintf(stdout,"base A Count:\t%lu(%.2f%%)\t%lu(%.2f%%)\n",stat_1.aCount,(float)stat_1.aCount/stat_1.baseCount*100.,stat_2.aCount,(float)stat_2.aCount/stat_2.baseCount*100.);
	fprintf(stdout,"base C Count:\t%lu(%.2f%%)\t%lu(%.2f%%)\n",stat_1.cCount,(float)stat_1.cCount/stat_1.baseCount*100.,stat_2.cCount,(float)stat_2.cCount/stat_2.baseCount*100.);
	fprintf(stdout,"base G Count:\t%lu(%.2f%%)\t%lu(%.2f%%)\n",stat_1.gCount,(float)stat_1.gCount/stat_1.baseCount*100.,stat_2.gCount,(float)stat_2.gCount/stat_2.baseCount*100.);
	fprintf(stdout,"base T Count:\t%lu(%.2f%%)\t%lu(%.2f%%)\n",stat_1.tCount,(float)stat_1.tCount/stat_1.baseCount*100.,stat_2.tCount,(float)stat_2.tCount/stat_2.baseCount*100.);
	fprintf(stdout,"base N Count:\t%lu(%.2f%%)\t%lu(%.2f%%)\n",stat_1.nCount,(float)stat_1.nCount/stat_1.baseCount*100.,stat_2.nCount,(float)stat_2.nCount/stat_2.baseCount*100.);
	fprintf(stdout,"\n");
	fprintf(stdout,"filtered read Count:\t%lu\t%lu\n",stat_1.filter,stat_2.filter);
	fprintf(stdout,"Number of base calls with quality value of 20 or higher (Q20+) (%%)\t%lu(%.2f%%)\t%lu(%.2f%%)\n",\
		stat_1.q20,(float)stat_1.q20/stat_1.baseCount*100.,stat_2.q20,(float)stat_2.q20/stat_2.baseCount*100.);
	fprintf(stdout,"Number of base calls with quality value of 30 or higher (Q30+) (%%)\t%lu(%.2f%%)\t%lu(%.2f%%)\n",\
		stat_1.q30,(float)stat_1.q30/stat_1.baseCount*100.,stat_2.q30,(float)stat_2.q30/stat_2.baseCount*100.);

	exit(0);
}
