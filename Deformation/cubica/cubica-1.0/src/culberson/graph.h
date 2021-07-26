/*
	GRAPH FORMAT DEFINITIONS FILE
*/

/*
        Title: Graph Definitions File.
        file: graph.h
        does: 
	        Graph definition and input
	        Macros for access to graph
	        cheat Info from Graph input

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

/*
	Graph definition and input
	Macros for access to graph
	Cheat Info from Graph input
*/
#ifndef GRAPHDEFS
#define GRAPHDEFS

/*
	!!!WARNING!!!!
	MAXVERTEX must be divisible by 8
*/
#define MAXVERTEX  2000

#define SHIFT 3
#define MASK 7
#define ROWSIZE ((MAXVERTEX >> SHIFT) +1 )
#define GRAPHSIZE  (MAXVERTEX * ROWSIZE)


/* Definitions useful for checking and setting edges */

/*                     ***NOTE***
    set and clear are asymmetric - use setedge(i,j) setedge(j,i) etc.
*/
#define setedge(i,j)  graph[((i)*ROWSIZE) + ((j) >> SHIFT)] |= (1 << ((j) & MASK))
#define clearedge(i,j) graph[((i)*ROWSIZE) + ((j) >> SHIFT)] &= ~(1 << ((j) & MASK))

#define edge(i,j)    (graph[((i)*ROWSIZE)+((j) >> SHIFT)] & (1 << ((j) & MASK)) )

/* for loops involving potential neighbors */
#define initnbr(x,i)  (x) = graph + ((i)*ROWSIZE)
#define isnbr(x,i)  (((x)[(i) >> SHIFT]) & (1 << ((i) & MASK)))

typedef int vertextype;
typedef unsigned char adjacencytype;

extern adjacencytype graph[GRAPHSIZE];
extern vertextype order;

/* CHEAT INFORMATION FROM GRAPH */
extern int partset[MAXVERTEX];
extern int partitionflag;
extern int partitionnumber;

extern void printgraph();
extern void getgraph(char a[]);

#endif
