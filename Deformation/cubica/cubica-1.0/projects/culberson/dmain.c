/*
	DSATUR ALGORITHM APPLIED TO GRAPH
*/

/*
        Title: DSATUR Source File.
        file: dsatur.c
        does: dsatur routine
        Source: Joseph Culberson's Coloring Page
                http://web.cs.ualberta.ca/~joe/Coloring/index.html
        Author: (Of this source code) Joseph Culberson
	Based On: ``New Methods to Color the Vertices of a Graph''
		Daniel Br\'{e}laz, CACM 22, 251--256, 1979.
		(Some modifications may apply).
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
#include "dsatur.h"
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
	popmembertype m;
	int	vordering;
	int i;

        char info[256];
        long lseconds, lmicrosecs;
        int pid;

        /* process id */
        pid = getpid();
        printf("Process pid = %d\n",pid);

  /*
	printf("\n DSATUR COLORING\n");
	printf("Initial Vertex Ordering:\n\t1 -- inorder\n\t2 -- random\n");
	printf("\t3 -- decreasing degree\n\t4 -- increasing degree\n");

	printf("Using: ");
	scanf("%d",&vordering);
	printf("%d\n",vordering);
  */

  vordering = 2;

	for(i=0;i<order;i++) m.vc[i].vertex = i;
	switch(vordering) {
		case 1 : break;
		case 2 : permute(&m,0,order); 
			break;
		case 3 : 
			computedeg();
			qsort((char *) m.vc,(int) order,
                          sizeof(struct vrtxandclr), (compfunc)decdeg);
			break; 
		case 4: 
			computedeg();
			qsort((char *) m.vc,(int) order,
                          sizeof(struct vrtxandclr),(compfunc)incdeg);
			break;
		default : {
			printf("DSATUR: Illegal Vertex Ordering\n");
			exit(1);
		}
	}

	dsatur(&m);

	getcolorinfo(&m);
	
	printinfo(&m);

	verifycolor(&m);

  i = getrusage(RUSAGE_SELF,&tmp);
  lseconds = tmp.ru_utime.tv_sec-seconds;
  lmicrosecs = tmp.ru_utime.tv_usec-microsecs;

  sprintf(info, "DSATUR cpu = %5.2f pid = %d",
          lseconds+(lmicrosecs/1000000.0) , pid);

  fileres(name, &m,info);
	//printcoloring(&m);
}

int main(int argc, char *argv[])
{
	int seed;

	if (argc!=3) {
		printf("\tUsage: dsatur file seed\n");
		exit(1);
	}

	about("DSATUR");

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
  seed = atoi(argv[2]);
  printf("Randomization seed =  ");
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

