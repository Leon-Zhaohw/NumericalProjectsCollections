/* GREEDY.H Definitions */

/*
        Title: Greedy Definitions File.
        file: greedy.h
        does: prototype to greedy routine
		Types of Greedy routines, Kempe interface.
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

#ifndef GREEDYDEF
#define GREEDYDEF

#include "colorrtns.h"
#include "graph.h"

/* The greedy types allowed for the greedy algorithm */

#define SIMPLEG 1
#define LARGESTG 2
#define SMALLESTG 3
#define RANDSEQG 4
#define REVERSEG 5
#define STIRGRDY 6

/* if any added above, change the following to act as a bound */
#define MAXGRDY 7

typedef int greedytype;

extern void greedy(
/* the generalized greedy algorithm */
	popmembertype *member, /* defines initial permutation of vertices,
				and assignes the colors */
	vertextype first, 
	vertextype last, /* color vertices in positions i, first <= i < last */
	colortype oldcolor, /* previous coloring: MAXVERTEX if first coloring */
	greedytype which, /* which greedy type to use */
	int KEMPE); /* if 0, then do Kempe reduction, else not */

#endif
