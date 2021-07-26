/*
	ITERATED GREEDY COLORING OF A GRAPH
*/

/*
        Title: Iterated Greedy Sepup and Call Source File.
        file: igmain.c
        does: igmain routine
        Source: Joseph Culberson's Coloring Page
                http://web.cs.ualberta.ca/~joe/Coloring/index.html
        Author: Joseph Culberson
        email: joe@cs.ualberta.ca
	Copyright (c) 1997 Joseph Culberson. All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions
	are met:
	1. Redistributions of source code must retain the above copyright
   	notice, this list of conditions and the following disclaimer.
	2. Redistributions in binary form must reproduce the above copyright
   	notice, this list of conditions and the following disclaimer in the
   	documentation and/or other materials provided with the distribution.
	3. All advertising materials mentioning features or use of this 
	software must display the following acknowledgement:
     	This product includes software developed by J. Culberson at the
     	University of Alberta, Edmonton.
	4. Neither the name of the University nor the names of its contributors
   	may be used to endorse or promote products derived from this software
   	without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE CONTRIBUTORS ``AS IS'' AND ANY
	EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
	THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
	PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE
	CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
	SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
	NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
	HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
	OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
	EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

	THIS SOFTWARE IS SUPPLIED WITHOUT ANY SUPPORT SERVICES. 

*/

#include "mysys.h"
#include "itrgrdy.h"
#include "greedy.h"
#include "graph.h"
#include "colorrtns.h"


/* Global Time info */
long seconds,microsecs,i;
struct rusage tmp;


void colorsearch(char *name)
/* 
	This version uses the greedy, this routine 
	asks for the various parameters.
*/

{

	/* Iterated Greedy requires a lot of parameters */
	popmembertype m;

	int maxiter; /* number of iterations of greedy without improvement
			before terminating */
	int out_freq; /* how often to print progress line */

	int i,j,k;

	/* Reordering control information (probablistic variation) */
	int controlselect[CONTROLSIZE], /* frequency info */
	    controlset[CONTROLSIZE],  /* which heuristic */
	    controltot; /* total for computing random selection */


	int switchback,retrylimit;

	colortype targetclr; /* terminate if this color achieved */

	greedytype greedyWt[MAXGRDY];
	int totalgreedywt;
	int kempe_freq;

	char info[256];
	long lseconds, lmicrosecs;
	int pid;

	/* process id */
	pid = getpid();
	printf("Process pid = %d\n",pid);

	/* get the target color */
	printf("Enter target color ");
	scanf("%hu",&targetclr);
	printf("%d\n",targetclr);


	printf("GREEDY TYPE SELECTION\n");
	printf("\t%1d\tSimple Greedy\n",SIMPLEG);
	printf("\t%1d\tLargest First Greedy\n",LARGESTG);
	printf("\t%1d\tSmallest First Greedy\n",SMALLESTG);
	printf("\t%1d\tRandom Sequence Greedy\n",RANDSEQG);
	printf("\t%1d\tReverse Order Greedy\n",REVERSEG);
	printf("\t%1d\tStir Color Greedy\n",STIRGRDY);

	printf("For each greedy type, enter its weight for selection:\n");
	totalgreedywt = 0;
	for(i=SIMPLEG;i<MAXGRDY;i++) {
		printf("Enter weight for %d ",i);
		scanf("%d",&greedyWt[i]);
		printf("%d\n",greedyWt[i]);
		totalgreedywt += greedyWt[i];
	}
	if (totalgreedywt <= 0) {
		printf("Must have a positive total\n");
		exit(1);
	}

	printf("Enter the frequency for Kempe reductions: ");
	scanf("%d",&kempe_freq);
	printf("%d\n",kempe_freq);

	printf("Enter number of iterations before quitting (after last ");
	printf("improvement) ");
	scanf("%d",&maxiter);
	printf("%d\n",maxiter);

	printf("REORDER CONTROL INFORMATION\n");

	printf("Heuristic Controls\n");
	printf("\t1\treverse order (else in order)\n");
	printf("\t2\tshuffle\n");
	printf("\t4\tlargest first\n");
	printf("\t8\tsmallest first\n");
	printf("\t16\tincreasing total degree\n");
	printf("\t32\tdecreasing total degree\n");
	printf("Sum any subset, sorts [size[degree[shuffle,reverse,order]]]\n");

	printf("\nFrequency, control - terminate by 0 0\n");

	controltot = 0;
	k= 0;
	while (k < CONTROLSIZE) {
		scanf("%d %d",&i, &j);
		printf("%d %d\n",i,j);
		if (i == 0) {
			if (k < 1 ) {
				printf("ERROR: must have at least one control\n");
				exit(1);
			}
			controlset[k] = j;
			controlselect[k] = controltot+5;
			k = CONTROLSIZE+5;
		}
		else {
			controlset[k] = j;
			controltot += i;
			controlselect[k] = controltot;
			k++;
			if (k >= CONTROLSIZE) {
				printf("ERROR: too many controls\n");
				exit(1);
			}
		}
	}

	/* use default values instead */
	/* printf("Switchback, retry limit\n");
	   scanf("%d %d",&switchback, &retrylimit);
	   printf("%d %d\n",switchback,retrylimit);

	   printf("Output Frequency ");
	   scanf("%d",&out_freq);
	   printf("%d\n",out_freq); */

	switchback = 0; retrylimit = 0;
        out_freq = 100;

	/* compute degree sequences */
	computedeg();


	getacoloring(&m, name, &j); /* open the color file for this graph
			and get the desired coloring */

	verifycolor(&m);

	itrgrdy(&m, maxiter, out_freq, controlset, controlselect, controltot,
			switchback, retrylimit, targetclr, greedyWt,
			totalgreedywt, kempe_freq);

	getcolorinfo(&m);
	
	printinfo(&m);

	verifycolor(&m);

        i = getrusage(RUSAGE_SELF,&tmp);
        lseconds = tmp.ru_utime.tv_sec-seconds;
        lmicrosecs = tmp.ru_utime.tv_usec-microsecs;

	sprintf(info, "ITERATED GREEDY[%d] cpu = %5.2f pid = %d",
                j,lseconds+(lmicrosecs/1000000.0) , pid);

	fileres(name, &m,info);

	/* printcoloring(&m); */
	

}

int main(int argc, char *argv[])
{
	int seed;
        struct rlimit rlp;


	if (argc < 2) {
		printf("\tUsage: itrgrdy file [cpulimit]\n");
		exit(1);
	}

	about("ITERATED GREEDY");

	if (argc > 2 ) {
                getrlimit(RLIMIT_CPU,&rlp);
                rlp.rlim_cur =atoi(argv[2]);
#ifndef linux
                printf("Cpulimit set to %d\n", rlp.rlim_cur);
#else
                printf("Cpulimit set to %ld\n", rlp.rlim_cur);
#endif
                if ( setrlimit(RLIMIT_CPU, &rlp)) {
                        printf("setrlimit failed\n");
#ifndef linux
                        printf("errno = %d\n",errno);
#endif
                        exit(1);
                }
                signal(SIGXCPU,cpulimig);
        }

        signal(SIGINT,intfcig);


	/* Set up graph */
	i = getrusage(RUSAGE_SELF,&tmp);
	seconds = tmp.ru_utime.tv_sec;
	microsecs = tmp.ru_utime.tv_usec;
	getgraph(argv[1]);
	i = getrusage(RUSAGE_SELF,&tmp);
	seconds = tmp.ru_utime.tv_sec-seconds;
	microsecs = tmp.ru_utime.tv_usec-microsecs;
	printf("GRAPH SETUP cpu = %5.2f\n",
		seconds+(microsecs/1000000.0));

	/* perform color search */

	/* Random initialization for search */
	printf("Enter seed for search randomization: ");
	scanf("%d",&seed);
	printf(" %d\n",seed);
	srandom(seed);

	i = getrusage(RUSAGE_SELF,&tmp);
	seconds = tmp.ru_utime.tv_sec;
	microsecs = tmp.ru_utime.tv_usec;
	colorsearch(argv[1]);
	i = getrusage(RUSAGE_SELF,&tmp);
	seconds = tmp.ru_utime.tv_sec-seconds;
	microsecs = tmp.ru_utime.tv_usec-microsecs;
	printf("Coloring time cpu = %5.2f\n",
		seconds+(microsecs/1000000.0));

	return(0);
}

