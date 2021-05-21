#ifndef _H_SEQ
#define _H_SEQ

#include <zlib.h>
#include <stdint.h>

typedef struct Stat {
	float averageLen,gc;
	uint64_t readCount,baseCount,nCount,aCount,tCount,gCount,cCount,filter,q20,q30;
} stat;

typedef struct readInfo {
	char *line[4];
	size_t len, nCount,mean_qval;
	double lowQval, nRate;
} rinfo;

char *readline(gzFile file);

uint64_t quality20(char *line);

uint64_t quality30(char *line);

size_t meanQvalue(size_t len, char *read);

float lowQvrate(int q,char *read);

float nRate_calc(size_t len, size_t count);

void freeRead(rinfo read);

#endif
