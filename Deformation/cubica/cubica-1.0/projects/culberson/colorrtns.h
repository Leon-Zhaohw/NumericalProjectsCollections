/* COLOR ROUTINES AND STRUCTURES DEFINITIONS */

/*
        Title: Color routines Definitions File.
        file: colorrtns.h
        does:  Color Storage Structures.
		Prototypes to routines for manipulating the color storage
		structures, output of colors, graph information (degree seq.) 
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

#ifndef COLORRTNDEF
#define COLORRTNDEF

#include "graph.h"

/* COLOR STORAGE STRUCTURES */


typedef unsigned short colortype;

struct clrinfo {
        colortype       numcolors; /* number of colors used */
        int     total;  /* sum of colors (weight)  used */
};
typedef struct clrinfo  clrinfotype;

struct vrtxandclr {
        vertextype vertex;
        colortype color;
};
typedef struct vrtxandclr vrtxandclrtype;

struct popmember {
        clrinfotype clrdata;
        vrtxandclrtype vc[MAXVERTEX];
};
typedef struct popmember popmembertype;

/* COLOR MANIPULATION ROUTINES */
extern void printinfo( popmembertype *member);
/* print number of colors and color sum */

extern void printcoloring( popmembertype *member);
/* print the colors of vertices in order */

extern void printpurity( popmembertype *member);
/* the purity computation for closeness to hidden color */

extern void getcolorinfo( popmembertype *member);
/* computes the maximum and sum of colors in member */

extern void permute( popmembertype *member, vertextype first, vertextype last);
/* randomly permute the vertex ordering of member in the
   range first to last -1 */

extern void trivial_color( popmembertype *m);
/* apply the colors 1 to order to the vertices of m,
   whatever order they are in */

extern void verifycolor(popmembertype *m);
/* verify the coloring of m as to correctness and print
   and appropriate message */

/* print a permutation member - for debugging mostly*/
extern void printperm(popmembertype *m);

extern void getacoloring(popmembertype *m, char *name, int *which);
/* open the indicated file and obtain the coloring asked for by the user */

/* PROPERTIES OF GRAPH */

extern int degseq[];

extern int decdeg( vrtxandclrtype *a, vrtxandclrtype *b);
/* comparison routine for decreasing sort by degree */

extern int incdeg(vrtxandclrtype *a, vrtxandclrtype *b);
/* comparison routine for increasing sort by degree */

extern void computedeg();
/* Compute degree sequence.  */

/* FINAL OUTPUT TO RESULTS FILE */
extern void fileres(char *name, popmembertype *m, char *prgm);
/* name will be appended with .res, coloring data will be appended to 
		file name.res, prgm is name of program  */

extern void about(char *pgrmname);
#endif
