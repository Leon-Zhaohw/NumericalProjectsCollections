
/*
        Title: Partition Source File.
        file: partition.c
        does: partition routines used in TABU coloring procedures. 
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
#include "partition.h"
#include "brute.h"
#include "itrgrdy.h"

/* Global variables */
/* PARTITION */
parttype  partlist[MAXVERTEX];
colortype targetK; /* number of partition elements */

listtype  vertexlist[MAXVERTEX];


/* conflict info needed globally */
long numInConflict; /* number of vertices in conflict */
vertextype conflictList[MAXVERTEX]; /* actual vertices in conflict */
long cedges; /* number of edges between conflicting vertices */

void countvertexconf(vertextype v, colortype p, int *ecnt)
/* 
	ecnt = the number of conflicts between v and the members of p.
*/
{
	vertextype w;
	adjacencytype *x;

	w = partlist[p].first;
	*ecnt = 0;
	initnbr(x,v);
	while (w != ENDLIST) { /* v not in p */
		if ( isnbr(x,w) ) {
			(*ecnt)++;
		}
		w = vertexlist[w].next;
	}
}

void deleteconf(vertextype v)
{
	vertextype w;
	/* sanity check */
	if (vertexlist[v].where == ENDLIST) {
		printf("Error: deleting element not in conflict : deleteconf\n");
		exit(1);
	}

	/* move last element to fill position */
	numInConflict--;
	w = conflictList[numInConflict];
	vertexlist[w].where = vertexlist[v].where;
	conflictList[vertexlist[w].where] = w;
	vertexlist[v].where = ENDLIST;
}

void addconf(vertextype v)
/*
	ADD V to conflict list
*/
{
	/* sanity check */
	if (vertexlist[v].where != ENDLIST) {
		printf("Error: vertex already in the list : addconf\n");
		return;
	}
	vertexlist[v].where = numInConflict;
	conflictList[numInConflict] = v;
	numInConflict++;
}

void insert(vertextype v, colortype p) /* at head of list p */
/* vertextype v; which vertex */
/* colortype p;  which list */
{
	vertextype w;
	/* sanity check */
	if (vertexlist[v].part != ENDLIST) {
		printf("Error: double insertion : insert\n");
		exit(1);
	}

	/* put v in the partition list p */
	vertexlist[v].part = p;
	vertexlist[v].next = partlist[p].first;
	partlist[p].first = v;

	/* update the conflict list */
	vertexlist[v].nbconf = 0;
	w = partlist[p].first; 
	while (w != ENDLIST) {
		if (edge(v,w)) {
			vertexlist[v].nbconf++;
			vertexlist[w].nbconf++;
			if (vertexlist[w].nbconf == 1)
				addconf(w);
		}
		w = vertexlist[w].next;
	}
	if (vertexlist[v].nbconf > 0) addconf(v);
	cedges += vertexlist[v].nbconf;
}

void vdelete(vertextype v)
/*
	DELETE v from whatever partition it is in
	and update conflict sets
*/
{
	colortype p;
	vertextype w;

	/* sanity check */
	if (vertexlist[v].part == ENDLIST ) {
		printf("Error: deleting vertex not in any partition : delete\n");
		exit(1);
	}

	/* remove v from partlist */
	p = vertexlist[v].part;
	w = partlist[p].first;
	if (w == v) {
		partlist[p].first = vertexlist[v].next;
	}
	else {
		while (v != vertexlist[w].next) 
			w = vertexlist[w].next;
		vertexlist[w].next = vertexlist[v].next;
	}

	/* update conflicts */
	if (vertexlist[v].where != ENDLIST) { /* if no conflicts do nothing */
		cedges -= vertexlist[v].nbconf;
		w = partlist[p].first;
		while (w != ENDLIST) {
			if (edge(v,w)) {
				vertexlist[w].nbconf--;
				if (vertexlist[w].nbconf == 0) 
					deleteconf(w);
			}
			w = vertexlist[w].next;
		}
		vertexlist[v].nbconf = 0;
		deleteconf(v);
	}
	vertexlist[v].part = ENDLIST;	/* not on any list */
}

void movevertex(vertextype v, colortype p)
/*
	MOVE vertex v to partition p : must work even if v in p
*/
{
	/* if v in p do nothing */
	if (vertexlist[v].part != p) {
		vdelete(v);
		insert(v,p);
	}
}

void structinit()
{
	vertextype i;
	colortype j;

	/* initialize partition list to null */
	for(j=0;j<targetK;j++) {
		partlist[j].first = ENDLIST;
	}

	/* initialize conflict list and values */
	cedges = 0;
	numInConflict = 0; 
	for(i=0;i<order;i++) {
		vertexlist[i].nbconf = 0;
		vertexlist[i].where = ENDLIST;
		vertexlist[i].next = ENDLIST;
		vertexlist[i].part = ENDLIST;
	}
}

/* take the current partition and put in order into a popmember*/
void parttomem( popmembertype *m)
{
	int i,cnt;
	vertextype v;
	colortype p;
	char Testvec[MAXVERTEX];

	i = -1;
	v = partlist[targetK-1].first;
	while( v != ENDLIST) {
		i++;
		m->vc[i].vertex = v;
		v = vertexlist[v].next;
	}
	for(p=0;p<targetK-1;p++) {
		v = partlist[p].first;
		while (v != ENDLIST) {
			i++;
			m->vc[i].vertex = v;
			v = vertexlist[v].next;
		}
	}

	/* Test the vertex set for errors */
	for(i=0;i<order;i++) Testvec[i] = 0;
	for(i=0;i<order;i++) Testvec[m->vc[i].vertex] = 1;
	cnt = 0;
	for(i=0;i<order;i++) 
		if (Testvec[i] == 0) cnt++;
	if (cnt > 0 ){
		printf("Error %d missing vertices in parttomem\n",cnt);
		exit(1);
	}
}

/* take colored vertex set and seed a partition for tabu*/
void memtopart( popmembertype *m)
{
	colortype p,clr,minp;
	vertextype i,start;
	int x,mincnf;

	structinit();

	/* help tabu by putting largest first*/
	setvec((int) m->clrdata.numcolors);
	revblckcount(m);
	qsort((char *) m->vc,(int) order,
		sizeof(struct vrtxandclr),(compfunc)ccompare);

	clr = m->vc[0].color;
	p = 0;
	start = order;
	for(i=0;((i<order) && (p < targetK));i++)  {
		if (clr != m->vc[i].color) {
			clr = m->vc[i].color;
			p++;
		}
		if ( p < targetK) insert(m->vc[i].vertex,p);
		else start = i;
	}

	/* try to minimize impact of remaining vertices */
	printf("Vertices to be inserted = %d\n",order-start);
	for(i=start;i<order;i++) {
		minp =0;
		mincnf = MAXVERTEX;
		for(p=0;p<targetK;p++) {
			countvertexconf(m->vc[i].vertex,p,&x);
			if (x < mincnf){
				mincnf = x;
				minp = p;
			}
		}
		insert(m->vc[i].vertex,minp);
	}
}
