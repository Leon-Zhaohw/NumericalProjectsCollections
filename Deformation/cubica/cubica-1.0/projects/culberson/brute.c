
/*
        Title: Brute Source File.
        file: brute.c
        does: brute subroutine - does brute force extensions of
		TABU search lines.
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
#include "brute.h"
#include "colorrtns.h"
#include "graph.h"
#include "partition.h"

#define MAXSIZE 20

int cretrn[3]; /* for counting recursion stats */

int search( int depth, int conflim)
{
	colortype p,q;
	vertextype i;
	vertextype list[MAXSIZE], cnf;
	
	/* make copy of vertices  as conflictList will change*/
	for(i=0;i<numInConflict;i++){
		list[i] = conflictList[i];
	}
		
	cnf = numInConflict;

	/* try each v in each possible group */
	for(i=0;i<cnf;i++) {
		q = vertexlist[list[i]].part;
		for (p=0;p<targetK;p++) {
		   if (p != q) {
			cretrn[0]++;
			movevertex(list[i],p);
			/* terminate successfully when no more conflicts */
			if (numInConflict <=0 ) {
				return(1);
			}

			/* NOTE: We could reduce conflim to constrain search*/
			if ((depth >1) && (numInConflict < conflim))  {
				cretrn[1]++;
				/* terminate when recursion succeeds */
				if (search(depth-1,conflim-1) == 1) {
					return(1);
				}
			}
		   }
		}
		/* restore vertex to starting position */
		movevertex(list[i],q);
	}
	cretrn[2]++;
	return(0); /* failure */
}
			
int brute(int depth,int conflim)
/*
	Do a brute force search to depth to find a solution
*/
{
	int i,temp;
	for(i=0;i<3;i++) cretrn[i] = 0;

	if (conflim >= MAXSIZE) {
		printf(" ERROR: conflim being reset to MAXSIZE = %d\n"
			,MAXSIZE);
		conflim = MAXSIZE-1;
		return(0);
	}

	temp = numInConflict;

	if (1==search(depth,conflim)) {
		printf("Brute succeeds on cnf = %d : STATS :",temp);
		printf("%d %d %d\n",cretrn[0],cretrn[1],cretrn[2]);
		return(1);
	} else  {
/*
		printf("Brute fails : STATS :");
		printf("%d %d %d\n",cretrn[0],cretrn[1],cretrn[2]);
*/
		return(0);
	}
}

