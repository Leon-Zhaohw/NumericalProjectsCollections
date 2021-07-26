/*
	ITERATED GREEDY DEFINITION FILE
*/
/*
        Title: Iterated Greedy Definitions File.
        file: itrgrdy.h
        does: prototype to itrgrdy routine
		various sorting routines
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

#ifndef ITERGREEDYDEF
#define ITERGREEDYDEF

#include "colorrtns.h"
#include "greedy.h"

/* control array size */
#define CONTROLSIZE 50

/* FUNCTIONS AND VARIABLES */
/*
	These routines use permvec to sort a vector through qsort.
	permvec is a sequence of values, one for each color.
	When a new permutation is desired, permvec is given
	a set of values that will cause the vertices to be sorted
	according to the desired criteria.
	
	The following functions initialize and use permvec.

	These may be useful in other programs, e.g. iterated maxis
*/

/* Since CNTR is largest, block size when used takes precedence,
   degree is second and then previous order (or reverse) is final
*/
#define CNTR 1000000
#define DEGCNT 200


extern int permvec[MAXVERTEX+1];

/* comparison function to be passed as arg to qsort - 
	sorting is by color group using permvec */
extern int ccompare( struct vrtxandclr *a, struct vrtxandclr *b);

/* initialize permvec to sort in ascending order by color */
extern void setvec(int size);

/* initialize permvec to sort in descending order by color */
extern void revvec(int size);

/* Sort first by increasing set size */
/* Adds to current value CNTR for each vertex of color c;
   so if initialized by revvec or setvec then secondary sorting by that.
*/
extern void blckcount( popmembertype *mem);

/* Sort first by decreasing set size */
/* like blockcount, except CNTR is subtracted */
extern void revblckcount( popmembertype *mem);

/* sorting by increasing total degree */
/* total degree is the total degree if the color set */
extern void updegcnt( popmembertype *mem);

/* sorting by decreasing total degree */
extern void dwndegcnt( popmembertype *mem);

/* randomly shuffle the portion of the permvec from 
   positions start to num. ASSUMES that permvec has
   been initialized with at least one of the above
   (usually setvec) */
extern void shufflevec( int start,int num, int permvec[]);

/* THE ACTUAL ITERATED GREEDY ROUTINE */
extern void itrgrdy(
popmembertype *member,  /* contains vertex permutation and coloring */
int     maxiter,        /* maximum iterations without improvement */
int     out_freq,       /* how often to output data */
int     controlset[],   /* which control */
int     controlselect[], /* density vector for controls */
int     controltot,     /* total for selection limit */
int     switchback, 
int     retrylimit, /* retry controls from best so far */
colortype targetclr,    /* if this coloring achieved stop */
greedytype greedyWt[], /* randomized selection of greedy type */
int totalgreedywt,
int kempe_freq); /* kempe reduction applied once every kempe_freq iters */

/* interrupt cutoffs */
void cpulimig();
void intfcig();

#endif

