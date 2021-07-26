
/*
        Title: name Source File.
        file: bktdsat.c
        does: bktdsat routine
        Source: Joseph Culberson's Coloring Page
                http://web.cs.ualberta.ca/~joe/Coloring/lclindex.html
        Author: Joseph Culberson, Denis Papp
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
//#include <values.h>
#include "colorrtns.h"
#include "graph.h"
#include "maxclique.h"
#include "bktdsat.h"

/*choose minsat */

/*#define DEBUG */ 
/*#define IMPACT */

#define MAXCLR 640 
#define LONG_BITS	BITS(long)
#define MOD_MASK	(LONG_BITS-1)

#ifdef ANZAC
#define  	SHIFT_VALUE	6
#else
#define 	SHIFT_VALUE	5	
#endif

/* prototypes */
void BlockColor(vertextype v, colortype c, colortype maxclr, int branch,
		popmembertype *m);
void FindPair(colortype maxclr,vertextype *v, colortype *c,int *impval);
void move(vertextype v, int newsatur);
void fix(void);
void ApplyColor(vertextype v,colortype c, colortype maxclr, 
		int branch, popmembertype *m);
void Color(colortype maxclr, int branch, popmembertype *m);
#ifndef IMPACT
int impact(vertextype v, colortype c);
#endif


/* global variables */

vertextype nextv[MAXCLR+MAXVERTEX];
vertextype prev[MAXCLR+MAXVERTEX];
vertextype lclindex[MAXCLR+MAXVERTEX];

/* how many colors are conflicting */
vertextype satur[MAXVERTEX];
/* pointer to position in lists */
vertextype current[MAXCLR];

/* total of each adjacent color to vertex */
short clrset[MAXVERTEX][MAXCLR];

#ifdef IMPACT
/* color impacts */
short impact[MAXVERTEX][MAXCLR];
#endif

colortype bestcolor,maxsat,minsat;
vertextype numcolored;
popmembertype bestmember;
int fixed;

colortype target;
int maxbranch;
int minlimit,maxlimit;

/* choose min or max sat? */
int minmax;

#ifdef IMPACT
/* code to print out impact array */
#define USECLR 10
void pimpact(void)
{
	int i,j;

	printf("## ");
	for (i=0;i<USECLR;i++)
		printf("%4d ",i);
	printf("\n");

	for (i=0;i<order;i++) {
	printf("%2d ",i);
	for (j=0;j<USECLR;j++) 
		printf("%4d ",impact[i][j]);
	printf("\n");
	}
}
#endif

/* interrupt control */
int stopflag = 0;
void cpulimbk()
{
	printf("CPU TIME EXCEEDED -- let me clean up\n");
	stopflag = 1;
}
void intfcbk()
{
	printf("INTERRUPT DETECTED -- cleaning up\n");
	stopflag = 1;
}




/* m is an ineffecient way to organize data  NOTE */
void bktdsat(
popmembertype *m,
int branch,
colortype targetclr,
int min,int max)
{
	/* vertices of equal saturation are kept in a circular doubly linked
	   list. There are MAXCLR such lists. Vertex i is represented 0 <= i <
	   MAXVERTEX by the ith position. nextv and prev indicate the nextv and
	   previous vertices in the list. 
	   lclindex indicates the lclindex number of the vertex
	   in the permutation handed to brelaz. Each list is kept in its 
	   lclindex order. The positions MAXVERTEX <= i < MAXCLR represent the
	   reference positions in the circular lists; lclindex for these
	   positions is ENDLIST to make it easy to detect end of list, and
	   by keeping ENDLIST large makes it easy to detect insertion
	   conditions at end of list.
	*/
	vertextype v,w;
	vertextype i,j;
	int cliquesize;

	stopflag = 0;

	printf("At each iteration choose vertex of minimum (0) or maximum (1) saturation? ");
	scanf("%d",&minmax);
	printf("%d\n",minmax);

	/* find a large clique, like MAXIS */
	cliquesize = maxclique(m);

	if (cliquesize > targetclr) {
		printf("WARNING: a clique > target color was found\n");
		target = cliquesize;
	} else {
		target = targetclr;
	}
	/** if target is hard, bestcolor=target+1 else =n+1 */
	bestcolor = order+1;

	/* initialize */

	maxsat=0;
	minsat=0;
	numcolored=0;
	bestmember=*m;
	fixed = 0;
	maxbranch = order;
	minlimit = min;
	maxlimit = max;

#ifdef IMPACT
	/* initialize impact values to degree */
	for (i=0;i<order;i++) 
	for (j=0;j<i;j++)
		if (edge(i,j)) {
			impact[i][0]++;
			impact[j][0]++;
		}
	for (i=0;i<order;i++)
	for (j=1;j<MAXCLR;j++)
		impact[i][j] = impact[i][0];
		
	printf("initially\n");
	#ifdef DEBUG
	pimpact();
	#endif
#endif

	/* initially no color saturation */
	for(i=0;i<order;i++) {
		for(j=0;j<MAXCLR;j++) clrset[i][j]=0;
		satur[i]=0;
	}

	/* all vertices are on zero conflict list */
	for(i=MAXVERTEX;i<MAXCLR+MAXVERTEX;i++) {
		nextv[i] = prev[i] = i;
		lclindex[i] = ENDLIST;
	}

	/* the 0th conflict list is anchored in MAXVERTEX position of array */
 	w = MAXVERTEX;
	for(i=0;i<order;i++) {
		lclindex[v = m->vc[i].vertex] = i;
		/* vertices are coming in from smallest to largest,
		   so insert at end of list  thus never changing w.  
		 */
		nextv[v] = w;
	 	prev[v] = prev[w];
	 	nextv[prev[w]] = v;
  		prev[w] = v;
	}

	Color(0,branch,m);

	*m = bestmember;
	m->clrdata.numcolors = bestcolor;
	printf("Best Coloring found: %d\n",bestcolor);
}

#ifndef IMPACT
int impact(vertextype v, colortype c)
{
	adjacencytype *x;
	vertextype w;
	int impval=0;

	initnbr(x,v);
	for (w=0; w<order; w++) {
		if (isnbr(x,w) && lclindex[w] != ENDLIST && clrset[w][c]==0)
			impval++;
	}
	return (impval);
}
#endif

void fix(void)
{
	int j;

	/* initialize pointers to reference positions */
	for(j=0;j<MAXCLR;j++ ) current[j] = j+MAXVERTEX;

	/* scan for maximum saturation list */
	while (nextv[current[maxsat]] == current[maxsat] && maxsat>0) maxsat--;
	/* scan for min saturation list */
	while (nextv[current[minsat]] == current[minsat] && minsat<maxsat) minsat++;

	fixed = 1;
}

void Color(colortype maxclr, int branch, popmembertype *m)
{
	/* maxcnf is maxsat */
	vertextype v;
	colortype c;
	int impval;

#ifdef DEBUG
printf("maxsat =%d\n",maxsat);
#endif

	if (stopflag) return;

	if (numcolored>=order) {
		if (maxclr<bestcolor) {
			printf("\nCOLORED %d\n",maxclr);
			bestcolor=maxclr;
			bestmember=*m;
		} else 
			printf("End of branch, no better coloring\n");
	} else if (bestcolor<=target) {
		printf("Found target color\n");
	} else if (maxclr>=bestcolor) {
		printf("Worse or equal coloring, returning\n");
	} else {
		fix();
		
		if (maxsat==maxclr) {
			/* some vertex is completely saturated */
			v=nextv[current[maxsat]];	
			if (maxclr+1 < bestcolor) {
			/*printf("##COLOR %d at DEPTH %d\n",maxclr+1,numcolored);*/
				ApplyColor(v,maxclr+1,maxclr,branch,m);
			} 
			#ifdef DEBUG
			else 
				printf("%d:%d: Worse or equal coloring\n",
					     bestcolor, numcolored);
			#endif
		} else {
			FindPair(maxclr,&v,&c,&impval);
			ApplyColor(v,c,maxclr,branch,m);
			if (maxclr>=bestcolor) return;
			/* if impact==0 then this color was not at fault */
			if (branch>0 && impval>0  
			    && (numcolored<minlimit || numcolored>maxlimit)){
				if (numcolored <= maxbranch) {
					maxbranch = numcolored;
					printf("Branched at %d (%d)\n",
						maxbranch,bestcolor);
					fflush(stdout);	
				}
				BlockColor(v,c,maxclr,branch-1,m);
			}
		}
	}
}

/* #### use impact size array to save list of nbrs? */
void ApplyColor(vertextype v,colortype c, colortype maxclr, 
		int branch, popmembertype *m)
{
	vertextype oldlclindex,w;
	int oldmaxsat,oldminsat;
	int j;
	adjacencytype *x;

	oldmaxsat = maxsat;
	oldminsat = minsat;

	/* pull v off its list */
	nextv[prev[v]] = nextv[v];
	prev[nextv[v]] = prev[v];

	if (c>maxclr) maxclr = c;

	numcolored++;
	m->vc[lclindex[v]].color = c;
	oldlclindex = lclindex[v];
 	lclindex[v] = ENDLIST;	/* no longer on any list */

	/*#### use impact size array to save list of nbrs? */

	/* update saturation and impact lists */
	initnbr(x,v);	
	for (j=0;j<order;j++) {
		w = m->vc[j].vertex;	
		if (isnbr(x,w) && lclindex[w] != ENDLIST) {
			#ifdef IMPACT
			/* do impact */
			for (i=1;i<MAXCLR;i++)
			if (clrset[v][i]==0)
				impact[w][i]--;		
			#endif

		 	/* mark color in colorset and check if */
			/* color not previously adjacent to w */
			if (0==(clrset[w][c]++)) {
				#ifdef IMPACT
				/* do impact */
				initnbr(x2,w);
				for (i=0;i<order;i++) {
				  w2 = m->vc[i].vertex;
				  if (isnbr(x2,w2) && lclindex[w2] != ENDLIST) 
				    impact[w2][c]--;
				}
				#endif

				/* move vertex to nextv list */
				move(w,satur[w]+1);
				satur[w]++;
				if (maxsat < satur[w])
					maxsat = satur[w];
#ifdef DEBUG
printf("newmaxsat=%d\n",maxsat);
#endif
			}
		}
	}
	fixed=0;

	#ifdef IMPACT	
	#ifdef DEBUG
	pimpact();
	#endif
	#endif

	Color(maxclr,branch,m);

	if (!fixed) { 
		fix(); 
#ifdef DEBUG
printf("ERROR: current[] was not fixed\n");
#endif
	}

	/* restore saturation and impact lists */
	for (j=0;j<order;j++) {
		w = m->vc[j].vertex;
		if (isnbr(x,w) && lclindex[w] != ENDLIST) {
			#ifdef IMPACT
			/* do impact */
			for (i=0;i<MAXCLR;i++)
			if (clrset[v][i]==0)
				impact[w][i]++;
			#endif

			/* unmark color in colorset and check if */
			/* color now not adjacent to w */
			if (0==(--clrset[w][c])) {
				#ifdef IMPACT
				/* do impact */
				initnbr(x2,w);
				for (i=0;i<order;i++) {
				  w2 = m->vc[i].vertex;
				  if (isnbr(x2,w2) && lclindex[w2] != ENDLIST) 
				    impact[w2][c]++;
				}
				#endif

				/* assume satur[w]>0 */
				/* return vertex to prev list */
				move(w,satur[w]-1);
				satur[w]--;
				if (minsat > satur[w]) 
					minsat = satur[w];
			}
		}
	}
	fixed=0;

	#ifdef IMPACT
	#ifdef DEBUG
	pimpact();
	#endif
	#endif

	/* put v back on its list */
	if (prev[nextv[v]]!= prev[v])
		printf("ERROR: prev nextv %d != prev %d\n",v,v);
	if (nextv[prev[v]]!= nextv[v]) 
		printf("ERROR: nextv prev %d != nextv %d\n",v,v);

	prev[nextv[v]] = v;
	nextv[prev[v]] = v;

	numcolored--;
	lclindex[v] = oldlclindex;
	m->vc[lclindex[v]].color = 0;

	maxsat = oldmaxsat;
	minsat = oldminsat;
}


void move(vertextype v, int newsatur)
{
	vertextype z,zp;

	nextv[prev[v]] = nextv[v];
	prev[nextv[v]] = prev[v];
	if (current[satur[v]] == v) 
		current[satur[v]] = prev[v];

	/* insert v into new list */
	/* note use of maximal ENDLIST lclindex */
	z = current[newsatur];
	while (lclindex[zp=nextv[z]] < lclindex[v]) z=zp;

	nextv[v] = zp;
	prev[v] = z;
	nextv[z] = v;
	prev[zp]= v;
	current[newsatur] = v;
	
#ifdef DEBUG
	if (current[newsatur]==nextv[current[newsatur]] && current[newsatur]<MAXVERTEX) {
		printf("ERROR: move(): current==nextv[current]\n");
	}
#endif

}

void FindPair(colortype maxclr,vertextype *v, colortype *c, int *impval)
{
	int w,i,t;

	*impval=order; *c=1; *v=0;

	if (minmax==0) 
		w = nextv[current[minsat]];
	else
		w = nextv[current[maxsat]];

	for (i=1;i<=maxclr;i++) {
		if (clrset[w][i]==0) {
			#ifdef IMPACT
			  t = impact[w][i];
			#else
			  t = impact(w,i);
			#endif
			if (t < *impval) { 
				*impval = t; 
				*c = i;
				*v = w;
			}
		}
		if (*impval==0) break;
	}

/* CHECK ALL VERTICES FOR MIN IMPACT*/
	/*for (w = 0; w<order; w++) { 
		if (lclindex[w]!=ENDLIST) {
			for (i=1;i<=maxclr;i++) { 
				if (clrset[w][i]==0) {
					#ifdef IMPACT
  					  t = impact[w][i];
					#else
					  t = impact(w,i);
					#endif
					if (t < *impval) {
						*impval = t;
						*c = i;
						*v = w;
#ifdef  DEBUG
printf("choosing v=%d c=%d, imp=%d\n",*v,*c,*impval);
#endif
					}
				}
				if (*impval==0) break;
			}
		}
		if (*impval==0) break;
	}
*/	
}

void BlockColor(vertextype v, colortype c, colortype maxclr, int branch,
		popmembertype *m)
{
	clrset[v][c] = order;

	fix();

	move(v,satur[v]+1);
	satur[v]++;
	if (maxsat < satur[v])
		maxsat = satur[v];
	fixed=0;

	Color(maxclr,branch,m);

	if (!fixed) { 
		fix(); 
#ifdef DEBUG
fprintf(stderr,"ERROR: BlockColor: current[] was not fixed\n");
#endif
	} 

	move(v,satur[v]-1);
	satur[v]--;
	if (minsat > satur[v])
		minsat = satur[v];
	fixed=0;

	clrset[v][c] = 0;
}

