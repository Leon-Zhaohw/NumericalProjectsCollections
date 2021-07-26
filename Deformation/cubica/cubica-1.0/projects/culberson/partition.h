/*
	PARTITION DEFINITIONS
*/

/*
        Title: Partition Definitions File.
        file: partition.h
        does: prototypes to routines used in TABU search
		to partition and re-partiion vertex sets.
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


#ifndef PARTITIONDEF
#define PARTITIONDEF

#include "colorrtns.h"
#define ENDLIST 65535

typedef struct vlistelem {
	colortype part; /* which partition is vertex in */
	vertextype nbconf; /* how many conflicts between this vertex and others
				in its partition */
	vertextype where; /* position in conflictList if it has > 0 conflicts*/
	vertextype next; /* next vertex in partition list */
} listtype;

typedef struct plistelem {
	vertextype first;
} parttype;

extern colortype targetK; /* number of partition elements */
extern long cedges;
extern long numInConflict;
extern vertextype conflictList[MAXVERTEX];
extern listtype vertexlist[MAXVERTEX];
extern parttype  partlist[MAXVERTEX];

extern void initpart();
extern void countvertexconf(vertextype v, colortype p, int *ecnt);
extern void movevertex(vertextype v, colortype p);
extern void insert(vertextype v, colortype p);
extern void structinit();
extern void parttomem( popmembertype *m);
extern void memtopart( popmembertype *m);

#endif
