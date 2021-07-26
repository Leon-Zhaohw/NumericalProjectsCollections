/*
	GREEDY COLORING OF A GRAPH
*/

/*
        Title: Greedy Set up and Call Source File.
        file: gmain.c
        does: Takes user parameters and calls greedy.
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
#include "greedy.h"
#include "graph.h"
#include "colorrtns.h"

/* Global Time info */
long seconds,microsecs,i;
struct rusage tmp;

void lbfs(popmembertype *m, int *secondary)
/*
	generate a lexicographic breadthfirst search order
	using secondary to break ties (max first)
*/
{
	int i,oldi,j,t, first;
	vertextype vec[MAXVERTEX+1]; 
	int next[MAXVERTEX+1], start[MAXVERTEX+1];
	int current,zfree;
	vertextype temp;

	for(i=0;i<order;i++) vec[i] = i;

	/* start[i] is the first vertex of the i+1st class */
	/* first indicates the first lexico classs still active */
	/* next is a linked list of classes in lexico order */
	first = 0;
	next[first]=order;
	start[first] = order; /* initially all vertices are in the null 
				label class */
	zfree = 1;
	current = 0;
	while(current < order) {
		/* select vertex from front class with maximum secondary 
		   note current will always be pointing at the initial vertex
		   of the first availble class
		*/

		for(i=current+1;i<start[first]; i++)
			if (secondary[vec[i]] > secondary[vec[current]] ) {
				temp = vec[i];
				vec[i] = vec[current];
				vec[current] = temp;
			}

		i= current+1;
		t = first;

		/* split the classes */
		while (i < order) { 
			j = start[t]-1;
			oldi = i;
			while (i<=j) {
				if (edge(vec[i],vec[current]) ) i++;
				else {
					temp = vec[i];
					vec[i] = vec[j];
					vec[j] = temp;
					j--;
				}
			}
			/* i has passed j; if both i and j moved,
			   then i is now first of new sub-class 
			   insert it.
			*/
			if ( (oldi != i)  && (j != start[t]-1) ) {
			   next[zfree] = next[t];
			   start[zfree] = start[t];
			   next[t] = zfree;
			   start[t] = i;
			   zfree++;
			   t = next[t]; /* so we will skip below */
			}

			i = start[t];
			t = next[t];
		}

		/*advance current and determine whether first should advance*/
		current++;
		if (current >= start[first]) first = next[first];
	}

	/* copy vec to the popmember */
	for(i=0;i<order;i++) m->vc[i].vertex = vec[i];

#ifdef DEBUG
	for(i=0;i<30;i++) {
		printf("%d: ",i);
		for(j=0;j<i;j++) {
			if (edge(m->vc[i].vertex,m->vc[j].vertex))
				printf("%d ", j);
		}
		printf("\n");
	}
#endif

}


void colorsearch(char *name)
/* 
	This version uses the greedy, this routine 
	asks for the various parameters.
*/

{
	popmembertype m;
	greedytype initGreedy;
	int	vordering,kempewh;
	char	usekempe[5];
	int i;

	int *secondary;

        char info[256];
        long lseconds, lmicrosecs;
        int pid;

	/* process id */
	pid = getpid();
	printf("Process pid = %d\n",pid);

	printf("GREEDY TYPE SELECTION\n");
	printf("\t%1d\tSimple Greedy\n",SIMPLEG);
	printf("\t%1d\tLargest First Greedy\n",LARGESTG);
	printf("\t%1d\tSmallest First Greedy\n",SMALLESTG);
	printf("\t%1d\tRandom Sequence Greedy\n",RANDSEQG);
	printf("\t%1d\tReverse Order Greedy\n",REVERSEG);
	printf("\t%1d\tStir Color Greedy\n",STIRGRDY);
	
	printf("Which for this program ");
	scanf("%d",&initGreedy);
	printf("%d\n",initGreedy);
	if ((initGreedy >= MAXGRDY) || (initGreedy < 1)){
		printf("GREEDY: Illegal Greedy type\n");
		exit(1);
	}

	printf("Initial Vertex Ordering:\n\t1 -- inorder\n\t2 -- random\n");
	printf("\t3 -- decreasing degree\n\t4 -- increasing degree\n");
	printf("\t5 -- LBFS random\n\t6 -- LBFS decreasing degree\n");
	printf("\t7 -- LBFS increasing degree\n");

	printf("Using: ");
	scanf("%d",&vordering);
	printf("%d\n",vordering);

	switch(vordering) {
		case 1 :
			for(i=0;i<order;i++) m.vc[i].vertex = i;
			break;
		case 2 : 
			for(i=0;i<order;i++) m.vc[i].vertex = i;
			permute(&m,0,order); 
			break;
		case 3 : 
			for(i=0;i<order;i++) m.vc[i].vertex = i;
			computedeg();
			qsort((char *) m.vc,(int) order,
                          sizeof(struct vrtxandclr), (compfunc)decdeg);
			break; 
		case 4: 
			for(i=0;i<order;i++) m.vc[i].vertex = i;
			computedeg();
			qsort((char *) m.vc,(int) order,
                          sizeof(struct vrtxandclr),(compfunc)incdeg);
			break;
		case 5: 
			secondary = (int *)malloc(order*sizeof(int));
			for(i=0;i<order;i++) secondary[i] = random();
			lbfs(&m,secondary);
			free(secondary);
			break;
		case 6:
			secondary = (int *) malloc(order*sizeof(int));
			computedeg();
			for(i=0;i<order;i++) secondary[i] = degseq[i];
			lbfs(&m,secondary);
			free(secondary);
			break;
		case 7: 
			secondary = (int *) malloc(order*sizeof(int));
			computedeg();
			for(i=0;i<order;i++) secondary[i] = order - degseq[i];
			lbfs(&m,secondary);
			free(secondary);
			break;
		default : {
			printf("GREEDY: Illegal Vertex Ordering\n");
			exit(1);
		}
	}

	printf("Use kempe reductions y/n ");
	scanf("%s",usekempe);
	printf("%s\n",usekempe);
	if  ( usekempe[0] == 'y') {
		kempewh = 0;
	} else if (usekempe[0] == 'n' ) {
		kempewh = 1;
	} else {
		printf("GREEDY: illegal input for Kempe selection\n");
		exit(1);
	}

	greedy(&m, 0, order, MAXVERTEX, initGreedy, kempewh);

	getcolorinfo(&m);
	
	printinfo(&m);

	verifycolor(&m);

        i = getrusage(RUSAGE_SELF,&tmp);
        lseconds = tmp.ru_utime.tv_sec-seconds;
        lmicrosecs = tmp.ru_utime.tv_usec-microsecs;

        sprintf(info, "GREEDY cpu = %5.2f pid = %d",
                lseconds+(lmicrosecs/1000000.0) , pid);

        fileres(name, &m,info);



	/* printcoloring(&m); */
	

}

int main(int argc, char *argv[])
{
	int seed;

	if (argc!=2) {
		printf("\tUsage: greedy file\n");
		exit(1);
	}

	about("GREEDY");

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

