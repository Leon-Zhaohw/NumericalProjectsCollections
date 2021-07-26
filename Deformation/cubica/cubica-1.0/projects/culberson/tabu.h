/*
	 TABU DEFINITIONS
*/

/*
        Title: TABU Definitions File.
        file: tabu.h
        does: prototype to tabucol routine
        Source: Joseph Culberson's Coloring Page
                http://web.cs.ualberta.ca/~joe/Coloring/index.html
        Author: Code by Joseph Culberson, based on the description
		by Hertz and de Werra,  ``Using Tabu Search Techniques
		for Graph Coloring'' in Computing 39, 345--351, 1987,
		but with several modifications.
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

#ifndef TABUDEF
#define TABUDEF

#include "colorrtns.h"
#include "graph.h"
/* the maximum number of tabu elements */
#define MAXTABU 100

/* For the tabu neighbors information list */
#define MAXTABULIST 1000

typedef struct tbinf {
	vertextype vertex; /* which vertex is being moved */
	colortype to; /* to which partition  */
	int cedges; /* conflict edges if used */
	
} tabuinfotype;

extern void tabucol(
                /* Note the presence of global targetK in partition.c */
int     nbmax,  /* maximum iterations before improvement */
int     rep,    /* number of representative neighbors */
int     minrep, /* minimum nbrs to generate before quitting on new min fnc*/
int     tabusize,       /* size of tabu list */
int     steps); /* number of increments to target color to make if failure */

#endif
