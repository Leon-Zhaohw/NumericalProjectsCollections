/*
	BACKTRACK DSATUR ALGORITHM APPLIED TO GRAPH
*/
/*
        Title: Backtrack Dsatur Source File.
        file: bkdmain.c
        does: Start up for the bktdsat routine
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
#include "bktdsat.h"
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
	int i;

	int branch, min, max;
	colortype targetcolor;



        char info[256];
        long lseconds, lmicrosecs;
        int pid;

        /* process id */
        pid = getpid();
        printf("Process pid = %d\n",pid);



	printf("\n BACKTRACK DSATUR COLORING\n");

	for(i=0;i<order;i++) m.vc[i].vertex = i;
	computedeg();
	qsort((char *) m.vc,(int) order,
		sizeof(struct vrtxandclr), (compfunc)decdeg);

	printf("Target Color - terminates if achieved: ");
	scanf("%hu",&targetcolor);
	printf("%d\n",targetcolor);

	printf("Maximum branching factor: ");
	scanf("%d",&branch);
	printf("%d\n",branch);

	printf("Block off which area\n");
	printf("  ('x y' means no branching in depth range x to y): ");
	scanf("%d %d",&min,&max);
	printf("%d %d\n",min,max);

	bktdsat(&m,branch,targetcolor,min,max);

	getcolorinfo(&m);
	
	printinfo(&m);

	verifycolor(&m);

        i = getrusage(RUSAGE_SELF,&tmp);
        lseconds = tmp.ru_utime.tv_sec-seconds;
        lmicrosecs = tmp.ru_utime.tv_usec-microsecs;

        sprintf(info, "BACKTRACK DSATUR cpu = %5.2f pid = %d",
                lseconds+(lmicrosecs/1000000.0) , pid);

        fileres(name, &m,info);


	/* printcoloring(&m); */

}

int main(int argc, char *argv[])
{
	int seed;
        struct rlimit rlp;

	if (argc < 2) {
		printf("\tUsage: bktdsat [cpulimit]file\n");
		exit(1);
	}
	
	about("BACKTRACK DSATUR");

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
                signal(SIGXCPU,cpulimbk);
        }

        signal(SIGINT,intfcbk);


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

