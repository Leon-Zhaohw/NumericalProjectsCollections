
/*
        Title: Maximal Clique  Source File.
        file: maxclique.c
        does: maxclique routine
        Source: Joseph Culberson's Coloring Page
                http://web.cs.ualberta.ca/~joe/Coloring/index.html
        Author: Joseph Culberson, Denis Papp
	Based On: Bollobas and Thomason "Random Graphs of Small Order" in
        Random Graphs 83, Annals of Disc. Math. 28, pp47-97.
        See in particular section 6. "Colouring large random graphs"
        pp 86-90

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
	Based on the maxis routines: Modification by D. Papp.
	For use by bktdsat.

	Bollobas and Thomason "Random Graphs of Small Order" in
	Random Graphs 83, Annals of Disc. Math. 28, pp47-97.
	See in particular section 6. "Colouring large random graphs"
	pp 86-90
*/

/*	Finds a large clique then permutes the vertices with the clique
	first and then sorts the rest sorted by degree
 */

#include "mysys.h"
#include <math.h>
#include <stdlib.h>
#include "colorrtns.h"
#include "maxclique.h"

#define MAXIS 500
#define MAXCUT 100

#define isnotnbr(x,i)   !(((x)[(i) >> SHIFT]) & (1 << ((i) & MASK)))

typedef struct isinfo {
	vertextype	possible, /*vertex not incident to any in IS */
			degree;   /*degree of vertex*/
} istype;


typedef istype twoDarray[MAXIS][MAXVERTEX];
typedef twoDarray *twoDarrayp;

/* globals */
int degree[MAXVERTEX];

int middeg;

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

int midcmp( istype *a, istype *b)
/*
        comparison routine for sorting mid values to front
*/
{
        if (abs(middeg - (int) (a->degree)) > abs(middeg - (int) (b->degree)))
                return(1);
        else if (a->degree == b->degree) return(0);
        else    return(-1);
}

/*
	brute force clique finding
	using indset() function from maxis - so variable names are
	not necessarily appropriate
*/
void clique(
int	freesize,	/* number of remaining vertices */
int	*retsize,	/* size of best independent set found */
vertextype bestset[],	/* returned independent set */
vertextype okaylist[], 	/* list of vertices we start with */
int	cutlimit[],
int	cutvertex[],
int	limit,
int	usort,
int	msort)
{
/*
	A complete rewrite of program for selecting maximum clique
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
			if (isnotnbr(x,(pptr[k].possible))) {
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
		qsort((char *) pptr, (int) numposs[0], sizeof(istype), 
							(compfunc)udgcmp);
	} else if (numposs[0] >= msortlimit) {
#ifdef DBG
	printf("First Sorting medium degree\n");
#endif
		qsort((char *) pptr, (int) numposs[0], sizeof(istype), 
							(compfunc)midcmp);
	}
	else {
#ifdef DBG
	printf("First Sorting increasing degree\n");
#endif
		qsort((char *) pptr, (int) numposs[0], sizeof(istype), 
							(compfunc)dgcmp);
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
				if ( !(isnotnbr(x,(pptr[i].possible))) ) {
					nptr[j].possible = pptr[i].possible;
					nptr[j].degree = 0;
					initnbr(y, (nptr[j].possible));
					for(k=0;k<j;k++) 
					   if (isnotnbr(y,(nptr[k].possible))) {
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

static	int limit,usort,msort;
static	int cutlimit[MAXCUT], cutvertex[MAXCUT]; /* branch & bound values */

int maxclique( popmembertype *m)
{
	int i,j;
	vertextype okaylist[MAXVERTEX];
	vertextype bc[MAXIS];
	int bcindex;

	istype info[MAXVERTEX];

	int findclique;
	int updown;

	printf("Find maxclique? (0-no,1-yes) ");
	scanf("%d",&findclique);
	printf("%d\n",findclique);

	printf("Sort vertices by decreasing (0) or increasing (1) degree first? ");
	scanf("%d",&updown);
	printf("%d\n",updown);

	/* ### set defaults for the input - not sure what to use */
    if (findclique) {
	printf("Vertex Num, Cutlim in decreasing order to 0\n");
	i = -1;
	do {
		i++;
		scanf("%d %d",&cutvertex[i],&cutlimit[i]);
		printf("%d %d\n",cutvertex[i],cutlimit[i]);
	} while (cutvertex[i] > 0);

	printf("Backtrack limit = 0 means no limit to backtrack\n");
	printf("Backtrack limit = k means do not backtrack over first k\n");
	printf("Upsort limit(U) and Midsort limit(M) with |G| =N\n\
\tif NumRem > (U*N)/100 then sort by decreasing degree\n\
\telse if NumRem > (U*M)/100 then sort by medium degree\n\
\telse sort by increasing degree\n");

	printf("Enter Backtrack Limit, UpSort Limit, Midsort Limit");
	scanf("%d %d %d",&limit,&usort,&msort);
	printf("%d %d %d\n",limit,usort,msort);

	fflush(stdout);

	for(i=0;i<order;i++) okaylist[i] = i;

	bcindex = 0;

	clique(order,&bcindex,bc,okaylist,cutlimit,
			cutvertex,limit,usort,msort);

	printf("Max Clique: ");
	for(i=1;i<=bcindex;i++) printf("%d ",bc[i]);
	printf("\n");

	printf("Clique found size bcindex = %d ",bcindex);
	fflush(stdout);
     }	


	/* calculate degrees */
	for (i=0;i<order;i++) {
		info[i].degree=0;
		info[i].possible=i;
	}

	for (i=1;i<order;i++)
	for (j=0;j<i;j++)
	if(edge(i,j)) {
		info[i].degree++;
		info[j].degree++;
	}

     if (findclique) {

	/* place clique at beginning of permutation by giving max/min degree*/
	if (updown==0) {
		for (i=1;i<=bcindex;i++) 
			info[bc[i]].degree = order+info[bc[i]].degree;
	} else { 
		for (i=1;i<=bcindex;i++) 
			info[bc[i]].degree = -order+info[bc[i]].degree;
	}
     } else bcindex=1;


	/* sort vertices by degree */
	if (updown==0)  /* decreasing */
		qsort( info,order,sizeof(istype), (compfunc)udgcmp);
	else            /* increasing */
		qsort( info,order,sizeof(istype), (compfunc)dgcmp);
		
	for (i=0;i<order;i++) m->vc[i].vertex = info[i].possible;

	return(bcindex);
}
