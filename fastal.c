/* Fastal - get sequence length distribution from fasta file
** Copyright (C) 2014 Louis Ranjard (louis.ranjard@gmail.com)
**
** This program is free software; you can redistribute it and/or
** modify it under the terms of the GNU General Public License
** as published by the Free Software Foundation; either version 2
** of the License, or (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
******************************************************************************/

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include "kseq/kseq.h"
#include <getopt.h>

// get distribution length of fasta file

// compilation:
// gcc -O3 fastal.c -o fastal -lz
// ./fastal sequences.fasta
// gcc -g -O2 fastal.c -o fastal -lz
// valgrind --tool=memcheck --leak-check=yes --track-origins=yes ./fastal sequences.fasta


/* PROTOTYPES */
int cmpfunc(const void * a, const void * b);


// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
  extern char *optarg;
	extern int optind;
	gzFile fp;
	kseq_t *seq;
	int l,m,n=0,t=0,id05,id5,id95;
	double mean; 
	int median, p05, p95;
  int pos, *pos1, npos1, nchunk1, memchunk=1000000;
	static char usage[] = "Usage: %s <in.fasta>\n";
  
	if (argc<1) {
		fprintf(stderr, usage, argv[0]);
		exit(1);
	}
  
	fp = gzopen(argv[1], "r"); // STEP 1: open the file handler
  if (fp == NULL) {fprintf(stderr, "\n***Error: cannot open file %s\n", argv[1]);exit(1);}
  
	seq = kseq_init(fp); // STEP 2: initialize seq
	while ((l = kseq_read(seq)) >= 0) { // STEP 3: get total number of sequences
	  n++;
	}
	int lengths[n];
	
	kseq_destroy(seq); // STEP 4: get length of each sequence, needs to first reinitialise data stream
	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp); 
  for (m=0; m<n; m++) {
    l = kseq_read(seq);
	  lengths[m] = strlen(seq->seq.s);
	  //printf("%d\n",lengths[m]);
	  t = t+lengths[m];
	}
	
	qsort(lengths, n, sizeof(int), cmpfunc);
	id5 = (int) (n-1)*0.5;
	median = lengths[id5];
	id05 = (int) (n-1)*0.05;
	p05 = lengths[id05];
	id95 = (int) (n-1)*0.95;
	p95 = lengths[id95];
	
	printf("Number of sequences: %d\n", n);
	printf("Total length: %d\n", t);
	mean = (double)t/n;
	printf("Average length: %.3f\n", mean);
	printf("min \t\t 5%% \t\t med. \t\t 95%% \t\t max\n");
	printf("%d \t\t %d \t\t %d \t\t %d \t\t %d \n", lengths[0], p05, median, p95, lengths[n-1]);

  //for (m=0; m<n; m++) {
	//  printf("%d\n",lengths[m]);
	//}
	
	kseq_destroy(seq); // STEP 5: destroy seq
	gzclose(fp); // STEP 6: close the file handler
	
	return 0;
}

int cmpfunc(const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

