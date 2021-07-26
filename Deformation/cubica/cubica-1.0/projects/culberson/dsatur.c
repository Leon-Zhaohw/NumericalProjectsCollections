
/*
        Title: DSATUR Source File.
        file: dsatur.c
        does: dsatur routine
        Source: Joseph Culberson's Coloring Page
                http://web.cs.ualberta.ca/~joe/Coloring/index.html
        Author: ((Of this source code) Joseph Culberson
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
//#include <values.h>
#include "graph.h"
#include "colorrtns.h"
#include "dsatur.h"

#define MAXCLR 640 
#ifndef linux
#define LONG_BITS	BITS(long)
#else 
#define LONG_BITS       32
#endif
#define MOD_MASK	(LONG_BITS-1)

/* ANZAC is a DEC alpha and requires special consideration */
#ifdef ANZAC
#define  	SHIFT_VALUE	6
#else
#define		SHIFT_VALUE	5	
#endif

/* m is an ineffecient way to organize data  NOTE */
void dsatur( popmembertype *m)
/*
	Ref: New Methods to Color the Vertices of a Graph,
	Daniel Brelaz, CACM 22, 251--256, 1979.
	Like Brelaz's DSATUR, except preference on vertices of equal
	restriction will be decided on the permutation weights.
	These weights will act as the permutation order does in
	greedy: but only secondarily to the constraint that they
	must be of maximal conflict. These comments are subject to
	several interpretations.
*/
{
	/* we really need the graph as an adjacency list now ? */

	vertextype i,j,maxcnf,v,w;
	colortype c,clr;
	adjacencytype *x;
	vertextype z,zp;

	long clrset[MAXVERTEX][MAXCLR/LONG_BITS]; /* use bits to set colors used */

	/* vertices of equal saturation are kept in a circular doubly linked
	   list. There are MAXCLR such lists. Vertex i is represented 0 <= i <
	   MAXVERTEX by the ith position. next and prev indicate the next and
	   previous vertices in the list. 
	   index indicates the index number of the vertex
	   in the permutation handed to brelaz. Each list is kept in its 
	   index order. The positions MAXVERTEX <= i < MAXCLR represent the
	   reference positions in the circular lists; index for these
	   positions is ENDLIST to make it easy to detect end of list, and
	   by keeping ENDLIST large makes it easy to detect insertion
	   conditions at end of list.
	*/

	vertextype next[MAXCLR+MAXVERTEX];
	vertextype prev[MAXCLR+MAXVERTEX];
	vertextype index[MAXCLR+MAXVERTEX];

	/* how many colors are conflicting */
	vertextype satur[MAXVERTEX];
	/* pointer to position in lists */
	vertextype current[MAXCLR];

	clr = 0;

	/* initially no color saturation */
	for(i=0;i<order;i++) {
		for(j=0;j<(MAXCLR/LONG_BITS);j++) clrset[i][j]=0;
		satur[i] =0;
	}
	
	/* all vertices are on zero conflict list */
	for(i=MAXVERTEX;i<MAXCLR+MAXVERTEX;i++) {
		next[i] = prev[i] = i;
		index[i] = ENDLIST;
	}
	/* the 0th conflict list is anchored in MAXVERTEX position of array */
	w = MAXVERTEX;
	for(i=0;i<order;i++) {
		index[v = m->vc[i].vertex] = i;
		/* vertices are coming in from smallest to largest,
		so insert at end of list  thus never changing w.
		*/
		next[v] = w;
		prev[v] = prev[w];
		next[prev[w]] = v;
		prev[w] = v;
	}
		
	/* now do the actual coloring */
	maxcnf = 0;
	for(i=0;i<order;i++) { /* color each vertex */
		/* initialize pointers to reference positions */
		for(j=0;j<MAXCLR;j++ ) current[j] = j+MAXVERTEX;

		/* scan for maximum saturation list */
		while (next[current[maxcnf]] == current[maxcnf]) maxcnf--;
#ifdef DEBUG
printf("maxcnf =%d\n",maxcnf);
#endif


		/* v is vertex to color, remove from list */
		v = next[current[maxcnf]];
		prev[next[v]] = current[maxcnf];
		next[current[maxcnf]]= next[v];

		/* assign minimum color not adjacent to v */
		/* recall floating point trick FUTURE */
		for (c=0;c<MAXCLR;c++)
			/* divide by LONG_BITS and mod LONG_BITS */
			if (0 == (clrset[v][c >> SHIFT_VALUE] 
				 & (1L << ( c & MOD_MASK )))) {
				/* note for consistency we make color 1..nclrs*/
				m->vc[index[v]].color = 1+(clr = c);
				break;
			}
		index[v] = ENDLIST; /* no longer on any list */
#ifdef DEBUG
printf("v=%d clr =%d\n",v,clr+1);
#endif

		/* scan list of adjacencies in permutation order */
		initnbr(x,v);
		for(j=0;j<order;j++) {
			w = m->vc[j].vertex;
			if (isnbr(x,w) )
			  if (index[w] != ENDLIST) {
				/* color not previously adjacent to w */
				if (0== (clrset[w][clr >> SHIFT_VALUE] 
					& (1L << ( clr & MOD_MASK )))) {
					/* mark color in colorset */
					clrset[w][clr >> SHIFT_VALUE] |= 
							1L << (clr & MOD_MASK );

					/* remove vertex from current list */
					next[prev[w]] = next[w];
					prev[next[w]] = prev[w];
					if (current[satur[w]] == w)
						current[satur[w]] = prev[w];

					/* increase saturation and check max */
					satur[w]++;
					if (maxcnf < satur[w])
						maxcnf = satur[w];
#ifdef DEBUG
printf("newmaxcnf=%d\n",maxcnf);
#endif

					/* insert w into new list */
					/* Note use of maximal ENDLIST index*/
					z = current[satur[w]];
					while (index[zp=next[z]] < index[w])
						z=zp;
					next[w] = zp;
					prev[w] = z;
					next[z] = w;
					prev[zp] =w;
					current[satur[w]] = w;
				}
			}
		}
	}
}
