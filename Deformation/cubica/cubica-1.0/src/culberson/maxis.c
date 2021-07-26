
/*
        Title: MAXIS Source File.
        file: maxis.c
        does: maxis routine
        Source: Joseph Culberson's Coloring Page
                http://web.cs.ualberta.ca/~joe/Coloring/index.html
        Author: Joseph Culberson
	Based on:
	        Bollobas and Thomason "Random Graphs of Small Order" in
	        Random Graphs 83, Annals of Disc. Math. 28, pp47-97.
	        See in particular section 6. "Colouring large random graphs"
	        pp 86-90. Modifications may have been made.

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

/*
	Bollobas and Thomason "Random Graphs of Small Order" in
	Random Graphs 83, Annals of Disc. Math. 28, pp47-97.
	See in particular section 6. "Colouring large random graphs"
	pp 86-90
*/


#include "mysys.h"
#include <math.h>
#include <stdlib.h>
#include "colorrtns.h"
#include "graph.h"

#define MAXIS 500
#define MAXCUT 100


typedef struct isinfo {
	vertextype	possible, /*vertex not incident to any in IS */
			degree;   /*degree of vertex*/
} istype;

int udgcmp( istype *a, istype *b)
/*
	Comparison routine for sorting by degree downwards
*/
{
	if (a->degree < b->degree) return(1);
	else if (a->degree == b->degree) return(0);
	else return(-1);
}

int dgcmp( istype *a, istype *b)
/*
	Comparison routine for sorting by degree
*/
{
	if (a->degree > b->degree) return(1);
	else if (a->degree == b->degree) return(0);
	else return(-1);
}

int middeg;

int midcmp( istype *a, istype *b)
/*
	comparison routine for sorting mid values to front 
*/
{
	if (abs(middeg - (int) (a->degree)) > abs(middeg - (int) (b->degree)))
		return(1);
	else if (a->degree == b->degree) return(0);
	else	return(-1);
}

typedef istype twoDarray[MAXIS][MAXVERTEX];
typedef twoDarray *twoDarrayp;

void indset(
/*
	brute force indset finding 
*/
int	freesize,	/* number of remaining vertices */
int	*retsize,	/* size of best independent set found */
vertextype bestset[],	/* returned independent set */
vertextype okaylist[], 	/* list of vertices we start with */
int	cutlimit[],
int 	cutvertex[],
int	limit,
int	usort,
int	msort)
{
/*
	A complete rewrite of program for selecting maximum IS
	Data Structure:
*/
	vertextype curset[MAXIS]; /* current independent set */
 	twoDarrayp indsetinf;

	istype *nptr, *pptr;
	adjacencytype  *x,*y; /* speed up pointers*/

	int nextv[MAXIS]; /* pointer to next vertex to try in possible */
	int numposs[MAXIS]; /* number in possible */

	int degseq[MAXVERTEX];

	int nextis; /* stack pointer to curset */
	int bestis, bestdeg; /* size counts for best set */
	int lcldeg; /* for computing degree of is */

	int i,j,k; /* for loop controls etc.*/
	int prev; /* temporary for speed */

	int firsttime; /* depth control */
	int cutoff[MAXCUT];
	int usortlimit,msortlimit;

	int degtot;

	/* coefficients for quadratic cutoff rates */

	indsetinf = (twoDarrayp) malloc (sizeof(twoDarray));
	if (indsetinf == NULL) printf("ERROR: indset: not enough memory.\n");

	/* initialize 0th to the initial set of vertices */
	numposs[0] = freesize;
	pptr = (*indsetinf)[0];
	for(i=0;i<freesize;i++) {
		pptr[i].possible = okaylist[i];
		pptr[i].degree = 0;
		initnbr(x,pptr[i].possible);
		for(k=0;k<i;k++)
			if (isnbr(x,(pptr[k].possible))) {
				pptr[i].degree++;
				pptr[k].degree++;
			}
	}	
	nextv[0] = 0;

	usortlimit = (freesize * usort) /100;
	msortlimit = (freesize * msort) /100;
#ifdef DBG
	printf("freesize = %d, usortlimit = %d, msortlimit = %d\n",
		freesize,usortlimit, msortlimit);
#endif
	k=0;
	while (cutvertex[k] > numposs[0]) k++;
	if (cutlimit[k] < numposs[0]) 
		cutoff[0] = cutlimit[k];
	else	cutoff[0] = numposs[0];

	/* set degree sequence */
	degtot = 0;
	for(i=0;i<freesize;i++) {
		degseq[pptr[i].possible] = pptr[i].degree;
		degtot += pptr[i].degree;
	}
	
	if (freesize >=1) middeg = degtot / freesize;
	else middeg = degtot;

	/* sort middles to front */
	if (numposs[0] >= usortlimit) {
#ifdef DBG
	printf("First Sorting decreasing degree\n");
#endif
		qsort((char *) pptr, (int) numposs[0], sizeof(istype), (compfunc)udgcmp);
	} else if (numposs[0] >= msortlimit) {
#ifdef DBG
	printf("First Sorting medium degree\n");
#endif
		qsort((char *) pptr, (int) numposs[0], sizeof(istype), (compfunc)midcmp);
	}
	else {
#ifdef DBG
	printf("First Sorting increasing degree\n");
#endif
		qsort((char *) pptr, (int) numposs[0], sizeof(istype), (compfunc)dgcmp);
	}


#ifdef DSEQ
	printf("FIRST SORT\n");
        printf("numposs[0] = %d\n",numposs[0]);
        for(i=0;i<numposs[0];i++)
                printf("(%d,%d) ",(*indsetinf)[0][i].possible,
                        (*indsetinf)[0][i].degree);
        printf("\n");
        fflush(stdout);
#endif

	bestis = 0;
	bestdeg = 0;
	nextis = 1;
	firsttime = 1; 
	while ((nextis > limit) || firsttime) { 
		if (nextis >= limit) firsttime = 0;
		/* select next vertex */
		prev = nextis -1;

#ifdef DBG
                printf("deglim of prev = %d\n",cutoff[prev]);
#endif

		if  ( nextv[prev] >= cutoff[prev] ) {
			/* BACKTRACK */
			nextis--;
		}
		else if (bestis > (prev+(numposs[prev]-nextv[prev]))) {
			/* BOUNDED BACKTRACK  - there are too few vertices
					left to build a better set
					this is most useful on k-colorable Graph
			*/
			nextis--;
		}
		else { 
			/* use some speed up variables */
			nptr = (*indsetinf)[nextis];
			pptr = (*indsetinf)[prev];

			/* select the next vertex */
			curset[nextis] = pptr[nextv[prev]].possible;
			initnbr(x,(curset[nextis]));

			/* reset previous next */
			nextv[prev]++;

			/* create the possible list */
			nextv[nextis] = 0;
			j=0;
			/* note: we consider only the remaining vertices
			   of previous possible list */
			for(i=nextv[prev];i<numposs[prev];i++) {
				if ( !(isnbr(x,(pptr[i].possible))) ) {
					nptr[j].possible = pptr[i].possible;
					nptr[j].degree = 0;
					initnbr(y, (nptr[j].possible));
					for(k=0;k<j;k++) 
					   if (isnbr(y,(nptr[k].possible)) ) {
						nptr[k].degree++;
						nptr[j].degree++;
					}
					j++;
				}
			}
			numposs[nextis]=j;
			degtot = 0; /* mindeg = order; */
			for(i=0;i<j;i++) {
				degtot += nptr[i].degree;
			}
			if (j > 0) middeg = degtot / j;
			else middeg = degtot;

			if (numposs[nextis] >= usortlimit) {
#ifdef DBG
			printf("Sorting decreasing degree\n");
#endif
				qsort((char *) nptr, (int) numposs[nextis], 
				sizeof(istype), (compfunc)udgcmp);
			} else if (numposs[nextis] >= msortlimit) {
#ifdef DBG
			printf("Sorting median degree\n");
#endif
				qsort((char *) nptr, (int) numposs[nextis], 
				sizeof(istype), (compfunc)midcmp);
			} else {
#ifdef DBG
			printf("Sorting increasing degree\n");
#endif
				qsort((char *) nptr, (int) numposs[nextis], 
				sizeof(istype), (compfunc)dgcmp); 
			}

#ifdef DSEQ
	printf("INTERNAL SORT nextis=%d\n",nextis);
        printf("numposs[nextis] = %d\n",numposs[nextis]);
        for(i=0;i<numposs[nextis];i++)
                printf("(%d,%d) ",(*indsetinf)[nextis][i].possible,
                        (*indsetinf)[nextis][i].degree);
        printf("\n");
        fflush(stdout);
#endif
			
			k=0;
			while(cutvertex[k] > numposs[nextis]) k++;
			if (cutlimit[k] < numposs[nextis])
				cutoff[nextis] = cutlimit[k];
			else cutoff[nextis] = numposs[nextis];

#ifdef DBG
        printf("numposs[%d] = %d\n",nextis,numposs[nextis]);
        if (DBG & 2) {
        for(i=0;i<numposs[nextis];i++)
                printf("(%d,%d) ",(*indsetinf)[nextis][i].possible,
                        (*indsetinf)[nextis][i].degree);
        printf("\n");
        }
        fflush(stdout);
#endif


			/* keep track of the best so far */
			if (bestis < nextis ) {
				/* copy the set */
				bestdeg = 0;
				for(i=1;i<=nextis;i++) {
					bestset[i] = curset[i];
					bestdeg += degseq[curset[i]];
				}
				bestis = nextis;
			}
                        else if (bestis == nextis) {
                                /* compute degree */
                                lcldeg =0;
                                for(i=1;i<=nextis;i++)
                                        lcldeg += degseq[curset[i]];
                                if (bestdeg < lcldeg){
                                 for(i=1;i<=nextis;i++)
                                        bestset[i] = curset[i];
                                 bestdeg = lcldeg;
                                }
                        }


			/* next iteration */
			nextis++;
		}
	}
	*retsize = bestis;
	free(indsetinf);
}

static	int firstrmaxis = 1;
static	int limit,usort,msort;
static	int cutlimit[MAXCUT], cutvertex[MAXCUT]; /* branch & bound values */


void rmaxis( popmembertype *m, int respect)
/* respect a previous coloring 1==yes */
{
	int i,j,k,idx,flag,numremaining;
	colortype clr,oldclr;

	adjacencytype *x;

	long tsecs, tmicros;

	vertextype okaylist[MAXVERTEX];
	vertextype bc[MAXIS];
	vertextype v;
	int bcindex,baseidx;
	colortype used[MAXVERTEX];

	/* initialize the used set */
	for(i=0;i<order; i++) {
		used[i] = 0;
	}

	if (firstrmaxis) {
		firstrmaxis = 0;
		printf("Vertex Num, Cutlim in decreasing order to 0\n");
		i= -1;
		do {
			i++;
			scanf("%d %d",&cutvertex[i],&cutlimit[i]);
			printf("%d %d\n",cutvertex[i],cutlimit[i]);
		} while (cutvertex[i] > 0);

		printf("Backtrack limit = 0 means no limit to backtrack\n");
		printf("Backtrack limit = k means do not backtrack over first k\n");
		printf("Upsort limit(U) and Midsort limit(M) with |G| =N\n\
\tif NumRem > (U*N)/100 then sort by decreasing degree\n\
\telse if NumRem > (M*N)/100 then sort by medium degree\n\
\telse sort by increasing degree\n");

		printf("Enter Backtrack Limit, UpSort Limit, Midsort Limit");
		scanf("%d %d %d",&limit,&usort,&msort);
		printf("%d %d %d\n",limit,usort,msort);

		fflush(stdout);
	}

	clr = 0;
	numremaining = order;
	while (numremaining > 0) {
		i = getrusage(RUSAGE_SELF,&tmp);
		tsecs = tmp.ru_utime.tv_sec;
		tmicros = tmp.ru_utime.tv_usec;

		clr++;

		if (respect != 1) {
			/* put all unused in the list */
			j=0;
			for(i=0;i<order;i++) {
				if (0 == used[i]) {
					okaylist[j] = i;
					j++;
				}
			}
		} else {
			/* find the subset non-adjacent to next color */

			/* scan past previously colored sets */
			i=0;
			while (0 != used[m->vc[i].vertex]) i++;

			/* recolor remaining vertices of this color */
			oldclr = m->vc[i].color;
#ifdef DEBUG
	printf("i= %d oldclr= %d\n",i,oldclr);
#endif

			k = 0;
			while ( (i < order) && (oldclr == m->vc[i].color) ) {
				/* if unused then can color it */
				if (0 == used[(v= m->vc[i].vertex)]) {
					used[v] = clr;
					numremaining--;
					bc[k] = v;
					k++;
				}
				i++;
			}

			baseidx = k;
			/* compute the vertices not adjacent to this
			   independent set and to the right of it in 
			   the ordering */

			j = 0;
			for(idx=i; idx < order; idx++) {
				v = m->vc[idx].vertex;
				if (0 == used[v] ) {
					k=0; 
					flag = 1;
					initnbr(x,v);
					while ( (k<baseidx) && ( flag)) {
						if ( isnbr(x,(bc[k])) )
							flag = 0;
						k++;
					}
					if (flag) {
						okaylist[j] = v;
						j++;
					}
				}
			}
			printf("Base set for clr %d baseidx = %d\n",
					clr,baseidx);
			fflush(stdout);
#ifdef DEBUG
			for(i=0;i<baseidx;i++)
				printf("%d ",bc[i]);
			printf("\n");
			fflush(stdout);
#endif
		
		}

		bcindex = 0;

		if (j > 0) 
			indset((int) j,&bcindex,bc,okaylist,cutlimit,
				cutvertex,limit,usort,msort);

		i = getrusage(RUSAGE_SELF,&tmp);
		tsecs = tmp.ru_utime.tv_sec - tsecs;;
		tmicros = tmp.ru_utime.tv_usec - tmicros;

		/* preserve indset by coloring it */
		for(i=1;i<=bcindex;i++) used[bc[i]] = clr;
		numremaining -= bcindex;
		printf("Independent Set for clr %d bcindex = %d ",clr,bcindex);
		printf("CPU = %5.2f\n", tsecs+(tmicros/1000000.0));
		fflush(stdout);
#ifdef DEBUG
		for(i=1;i<=bcindex;i++)
			printf("%d ",bc[i]);
		printf("\n");
		fflush(stdout);
#endif
	}
	for(i=0;i<order;i++) {
		m->vc[i].vertex = i;
		m->vc[i].color = used[i];
	}
	m->clrdata.numcolors = clr;
}
