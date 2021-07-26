
/*
        Title: TABU Source File.
        file: tabu.c
        does: tabucol routine
        Source: Joseph Culberson's Coloring Page
                http://web.cs.ualberta.ca/~joe/Coloring/index.html
        Author(source Code): Joseph Culberson
	Based On:  A. Hertz and D. de Werra ``Using Tabu Search
	Techniques for Graph Coloring'' in Computing(39) 345-351,
	1987. Several modifications have been made to the original algorithm.
	No claim as to equivalent behavior with respect to the original is
	made.

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
#include "colorrtns.h"
#include "graph.h"
#include "tabu.h"
#include "partition.h"
#include "brute.h"

/* the hash table size  must be power of two */
#define MAXHASH 131072

/*
	This tabu version uses a hash code and if the element in hash
	is not zero, then it assumes it is tabu. This may occasionally
	give false tabus, but is should be much faster. Moves are added
	and deleted from the has table through a circular list of keys
*/

#define HASH(V,P)  ((MAXHASH-1) & (((V) << 7)  | (P)))
char hash[MAXHASH];
int key[MAXTABU];

void selectrand( tabuinfotype *info)
/* 
	select vertex and partition, subject to tabu list etc.
*/
{
	vertextype v;
	colortype p;
	int i;
	char flag;

	flag = 1;
	while (flag) {
		v = conflictList[ random() % numInConflict];
		p = (colortype) (random() % targetK);
		if (vertexlist[v].part != p) {
			flag =  hash[HASH(v,p)];
		}
	}
        /* fill in info */
        info->vertex = v;
        info->to = p;

	/* compute changes in conflicts if this is used */
	info->cedges = cedges - vertexlist[v].nbconf;
        countvertexconf(v,info->to,&i); /* number added */
	info->cedges += i;
}

/*
	More or less like the algorithm in the paper
	by Herz and de Werra
*/
void tabucol(
		/* Note the presence of global targetK in partition.c */
int	nbmax,	/* maximum iterations before improvement */
int	rep,	/* number of representative neighbors */
int	minrep,	/* minimum nbrs to generate before quitting on new min fnc*/
int	tabusize,	/* size of tabu list */
int	steps) /* number of increments to target color to make if failure */
{
	tabuinfotype infolist[MAXTABULIST];

	colortype j;
	int i,tabunext,nbiter,min,lim;
	int maxTarget;
	int bestiter,bestfnc;
	int brutefail;

	tabunext = 0;
	nbiter = 0;
	for(i=0;i<tabusize;i++) key[i] = 0;
	for(i=0;i<MAXHASH;i++) hash[i] = 0;

	maxTarget = targetK + steps;

	bestfnc = cedges; bestiter = 0;
	brutefail = 0;
	while ((cedges >0) && (targetK <= maxTarget)) {
		printf("Trying for target %d\n",targetK);
		printf("bestiter %d bestfnc %d\n",bestiter,bestfnc);
	        while ((cedges>0) && ( (nbiter - bestiter) < nbmax) ) {

	                /* generate neighbors */
			lim = rep;
			for(i=0;i<rep;i++) {
				selectrand(&infolist[i]);
				/* immediate selection */
				if (infolist[i].cedges < cedges) {
					if (i >= minrep) {
						lim = i+1;
						break;
					}
				}
			}
	
	                /* record best */
			min = 0;
			for(i=1;i<lim;i++) 
				if (infolist[i].cedges < infolist[min].cedges)
					 min = i;
	
	
	                /* update tabu list */
			/* - remove old conflict */
			tabunext = (1+tabunext) % tabusize;
			hash[key[tabunext]] = 0;
			hash[(key[tabunext] = HASH((infolist[min].vertex),
				(vertexlist[infolist[min].vertex].part)))] = 1;
	
			/* make the move */
			movevertex(infolist[min].vertex, infolist[min].to);
	
			if (numInConflict == 2)  {
				if (1== brute(2,3)) cedges = 0;
				else brutefail++;
			}
			else if ((numInConflict == 3) && (cedges == 2)) {
				if (1 == brute(2,3)) cedges = 0;
				else brutefail++;
			}
			
	
	                nbiter++;
			if (cedges < bestfnc) {
				bestfnc = cedges;
				bestiter = nbiter;
				printf("bestiter = %d bestfnc = %d\n",
					bestiter,bestfnc);
				fflush(stdout);
			}
			if (0==nbiter%1000) {
			  printf("nbiter %d cedges %ld nbconf %ld\n", 
				nbiter,cedges, numInConflict);
			  fflush(stdout);
			}
		}
		if (cedges > 0) {
			partlist[targetK].first = ENDLIST;
			targetK++;
			printf("Brute failed %d times\n",brutefail);
			brutefail = 0;
			nbiter = 0;
		}
        }

	if (cedges == 0) printf("Success nbiter = %d brutefail = %d\n",
			nbiter,brutefail);
	else printf("Failure nbiter = %d brutefail = %d\n",nbiter,brutefail);

	/* Create a special partition for all of the conflicting vertices */
	/* this will be used by reorder */
	j = targetK;
	targetK++;
	partlist[j].first = ENDLIST;
	for(i=0;i<numInConflict;i++) 
		movevertex(conflictList[i],j);
}
