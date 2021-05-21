#include <zlib.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "func.h"

char *readline(gzFile file)
{
	size_t baselen = 256;
	char *line = (char *)calloc(baselen,sizeof(char));
	if(!line) exit(1);

	int ch;
	size_t index = 0;
	while((ch=gzgetc(file)) != -1 && ch != 10)
	{
		line[index] = ch;
		index++;
		if(index == baselen)
		{
			baselen += 128;
			line=(char *)realloc(line,baselen);
			if(!line)
			{
				free(line);
				exit(1);
			}
		}
	}
	line[index] = '\0';  // tail '\0' 
	if(ch == -1)  return NULL;// end of file
	return line;
}

uint64_t quality20(char *line)
{
	int i=0; uint64_t qCount=0;
	while(line[i]!='\0')
	{
		int q=line[i]-33;
		if(q>=20) qCount++;
		i++;
	}
	return qCount;
}

uint64_t quality30(char *line)
{
	int i=0; uint64_t qCount=0;
	while(line[i]!='\0')
	{
		int q=line[i]-33;
		if(q>=30) qCount++;
		i++;
	}
	return qCount;
}

size_t meanQvalue(size_t len, char *read)
{
	size_t sum=0;
	for(int i=0; read[i]!='\0'; i++)
	{
		sum += read[i]-33;
	}
	return (size_t)sum / len;
}

float lowQvrate(int q,char *read)
{
	size_t len=strlen(read);
	size_t lowCount=0;
	for(int i=0; read[i]!='\0'; i++)
	{
		if(read[i]-33<=q) lowCount++;	
	}
	return lowCount / len ;
}

float nRate_calc(size_t len, size_t count)
{
	return (float)count / len ;
}

void freeRead(rinfo read)
{
	for(int i=0; i<4; i++)
	{
		free(read.line[i]);
		read.line[i]=NULL;
	}
	return ;
}

