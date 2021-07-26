/*
	TABU COLORING OF A GRAPH
*/

/*
        Title: TABU setup and Call Source File.
        file: tbumain.c
        does: tbumain routine
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
#include "tabu.h"
#include "greedy.h"
#include "itrgrdy.h"
#include "graph.h"
#include "colorrtns.h"
#include "partition.h"

/* Global Time info */
long seconds,microsecs,i;
struct rusage tmp;


void colorsearch(char *name)

{

	popmembertype m;

	int tabunbr,tabuminnbr,tabusize;

	int maxiter, tabusteps;

        char info[256];
        long lseconds, lmicrosecs;
        int j,pid;

        /* process id */
        pid = getpid();
        printf("Process pid = %d\n",pid);


	printf("TABU CONTROL INFORMATION\n");

	printf("Maximum Iterations for tabu: ");
	scanf("%d",&maxiter);

	printf("Enter number of neighbors ");
	scanf("%d",&tabunbr);
	printf("%d\n",tabunbr);
	if (tabunbr > MAXTABULIST) {
		printf("ERROR: Cannot exceed %d\n",MAXTABULIST);
		exit(1);
	}

	printf("Enter minimum number of neighbors ");
	scanf("%d",&tabuminnbr);
	printf("%d\n",tabuminnbr);

	printf("Enter tabu list size ");
	scanf("%d",&tabusize);
	printf("%d\n",tabusize);
	if (tabusize > MAXTABU) {
		printf("Tabu size must be no more than %d\n",MAXTABU);
		exit(1);
	}
	if (tabusize < 1) {
		printf("Tabu size must be > 0\n");
		exit(1);
	}

	printf("Enter the target color you  want to try for: ");
	scanf("%hu",&targetK);
	printf("%hu\n",targetK);

	printf("If target not achieved, how many increases allowed: ");
	scanf("%d",&tabusteps);
	printf("%d\n",tabusteps);


	/* compute degree sequences */
	computedeg();


	getacoloring(&m, name, &j); /* open the color file for this graph
			and get the desired coloring */

	verifycolor(&m);

	/* NEEDED target color  as parameter; allow trivial start */

	getcolorinfo(&m);

	memtopart(&m);

	tabucol(maxiter, tabunbr, tabuminnbr, tabusize, tabusteps);

	parttomem(&m);

	greedy(&m, 0, order, MAXVERTEX, SIMPLEG,1);

	getcolorinfo(&m);
	
	printinfo(&m);

	verifycolor(&m);

        i = getrusage(RUSAGE_SELF,&tmp);
        lseconds = tmp.ru_utime.tv_sec-seconds;
        lmicrosecs = tmp.ru_utime.tv_usec-microsecs;

        sprintf(info, "TABU[%d] cpu = %5.2f pid = %d",
                j,lseconds+(lmicrosecs/1000000.0) , pid);

        fileres(name, &m,info);


	/* printcoloring(&m); */
	

}

int main(int argc, char *argv[])
{
	int seed;

	if (argc!=2) {
		printf("\tUsage: tabu file\n");
		exit(1);
	}

	about("TABU");


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

