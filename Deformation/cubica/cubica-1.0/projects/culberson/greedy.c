
/*
        Title: Greedy Routines Source File.
        file: greedy.c
        does: greedy routine
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
#include "greedy.h"
#include "graph.h"
#include "colorrtns.h"

/* 
	GREEDY ALGORITHM
	assign a coloring to a subrange of member
	using the permutation order in the member.
*/
/*
	MODIFICATIONS May 1994
	Change Greedy to select the order of testing based on
	various criteria.
*/

static colortype maxindex;
static colortype clrorder[MAXVERTEX+1];
static int clrcnt[MAXVERTEX+1];
static colortype clrindex;

void stirclr(
popmembertype *member,
vertextype first, vertextype last)
		 /* color vertices in positions i, first <= i < last */
{
	/* PROTOTYPE _ NOT TUNED */
	/* for each vertex select the color sets that minimize
	   the number of vertices that are in conflict.
	   Insert the new vertex there, and recolor the 
	   ejected ones using greedy.
	*/

	int v, cfl[MAXVERTEX][MAXVERTEX];
	int cnfcnt[MAXVERTEX];
	colortype clr,j,minc;
	int x,i,k,h,mincnf;
	int flag;

	clr = 1;
	for(i=first;i<last;i++) {
		v = member->vc[i].vertex;
		for(j=1;j<=clr;j++) cnfcnt[j] =0;
		for(k=first;k<i;k++) 
		     if (edge(v,member->vc[k].vertex) ) {
			x = (cnfcnt[member->vc[k].color]++);
			cfl[member->vc[k].color][x] =k;
		}
		mincnf = order;
		minc = clr;
		for(j=1;j<=clr;j++) {
			if (cnfcnt[j] < mincnf) {
				mincnf = cnfcnt[j];
				minc = j;
			}
		}
		/* use the minimum for this vertex */
		member->vc[i].color = minc;
		if (mincnf != 0) {
			/* we must recolor the conflicts */
			flag = 0;
			for(h=1;h<=mincnf;h++){
				v = cfl[minc][h];
				for(j=1;j<=clr+1;j++) cnfcnt[j] =0;
				for(k=first;k<=i;k++) 
                     			if (edge(v,member->vc[k].vertex) ) {
                        			cnfcnt[member->vc[k].color]++;
                			}       
				k = clr+1;
				for(j=1;j<=k;j++)
					if (cnfcnt[j] == 0) {
						if (j<=clr) 
						  member->vc[i].color = j;
						else flag = 1;
					}
				}

			if (flag == 1) {
				 member->vc[i].color = k;
				 clr = k;
			}
		}
	}
}

void initializeClrSelect(
colortype oldcolor, /* we should never use more than the previous coloring
			-- Use MAXVERTEX if not previously colored */
greedytype	which)
{
	int i;
	maxindex = 1; /* we always have to use at least one */
	for(i=0;i<=oldcolor;i++) clrorder[i] =i;
	switch (which) {
		case SIMPLEG:	
				break;
		case LARGESTG: 	
				for(i=1;i<=oldcolor;i++) clrcnt[i] = 0;
				clrcnt[0] = MAXVERTEX;
				break;
		case SMALLESTG:	
				for(i=1;i<=oldcolor;i++) clrcnt[i] = 0;
				clrcnt[0] = MAXVERTEX;
				break;
		case RANDSEQG:	
				break;
		case REVERSEG: 
				break;
		default: printf("initializeClrSelect:Illegal greedy type %d\n",
				which );
			exit(1);
			break;
	}
}

colortype selectFirstClr(
greedytype which)
{
	switch (which) {
                case SIMPLEG:   
				clrindex = 1; 
				break;
                case LARGESTG:	
				clrindex = 1; 
				break;
		case SMALLESTG:	
				clrindex = maxindex; 
				break;
                case RANDSEQG:  
				clrindex = 1; 
				break;
                case REVERSEG:  
				clrindex = maxindex; 
				break;
		default: printf("selectFirstClr:Illegal greedy type %d\n",
				which );
                        exit(1);
                        break;
        }
	return(clrorder[clrindex]); 
}

colortype nextClr(
greedytype which)
{
        switch (which) {
                case SIMPLEG:   
				clrindex++;
				break;
                case LARGESTG:  
				clrindex++;
				break;
		case SMALLESTG:	
				 if (clrindex>1) clrindex--;
				else clrindex = maxindex+1;
				break;
                case RANDSEQG:  
				clrindex++;
				break;
                case REVERSEG: 
				 if (clrindex>1) clrindex--;
				else clrindex = maxindex+1;
				break;
		default: printf("nextClr:Illegal greedy type %d\n",
				which );
                        exit(1);
                        break;
        }
	return(clrorder[clrindex]);
}

void updateClrSelect(
greedytype which)
{
	int i, temp;
        switch (which) {
                case SIMPLEG: 
				if (clrindex > maxindex) maxindex = clrindex;
				break;
                case LARGESTG:
				if (clrindex <= maxindex) {
				  clrcnt[clrorder[clrindex]]++;
				  i = clrindex;
				  /* move larger forward */
				  while (clrcnt[clrorder[i-1]] <
					        clrcnt[clrorder[i]] ) {
					temp = clrorder[i];
					clrorder[i] = clrorder[i-1];
					i--;
					clrorder[i] = temp;
				  }
				} else { /* start new color */
				  maxindex = clrindex;
				  clrcnt[clrorder[clrindex]] = 1;
				}
				break;
				
/* NOTE: As implemented, SMALLESTG and LARGESTG are not symmetric with respect
	to equal sized sets. LARGESTG favors previous equals, SMALLESTG
	favors new equals */
                case SMALLESTG:
                                if (clrindex <= maxindex) {
                                  clrcnt[clrorder[clrindex]]++;
                                  i = clrindex;
				  /* move larger forward */
                                  while (clrcnt[clrorder[i-1]] <
                                                clrcnt[clrorder[i]] ) {
                                        temp = clrorder[i];
                                        clrorder[i] = clrorder[i-1];
                                        i--;
                                        clrorder[i] = temp;
                                  }     
                                } else {
                                  maxindex = clrindex;
                                  clrcnt[clrorder[clrindex]] = 1;
                                }
                                break;

                case RANDSEQG:
				if (clrindex > maxindex) {
				/* NOTE: This ensures that the permutation
				   is equally probable, but there are biases:
				   if maxindex does not change, then the
				   permutation is the same for the next vertex
				   colored, otherwise it has a minor change */
					maxindex++;
					i = (random() % maxindex)+1;
					temp = clrorder[i];
					clrorder[i] = clrorder[maxindex];
					clrorder[maxindex] = temp;
				}
                                break;
                case REVERSEG:
				if (clrindex > maxindex) maxindex++;
                                break;
		default: printf("updateClrSelect:Illegal greedy type %d\n",
				which );
                        exit(1);
                        break;
        }
}


/* pointers for coloring sets */
vertextype start[MAXVERTEX + 1], next[MAXVERTEX+1];
colortype bestclr,fixedclr;

#ifdef DFSGREEDY
#define BRUTE 15
popmembertype best;
void dfsclr(
popmembertype *member;
vertextype first,last;
vertextype vert,crnt;
colortype maxclr)
{
	adjacencytype *x;
	vertextype v;
	colortype c;
	vertextype i,k;
	int flag;

	if (maxclr > bestclr) return;

	/* we assume throughout that fixedclr <= maxclr */

	v = member->vc[crnt].vertex;
	initnbr(x,v);
	for(c=1;c<=maxclr+1;c++) {
		k = start[c];
		flag = 1;
		/* check colors from fixed part */
		if (c <= fixedclr) while(k <= order) {
			if ( isnbr(x,k) ) 
				{flag = 0;break;}
			else 
				k = next[k];
		}
		/* we can use this color wrt fixed */
		if (flag) for(k=vert;k<crnt;k++)
			if ((member->vc[k].color == c) &&
				(isnbr(x,member->vc[k].vertex))) {
				flag = 0; break;
		}
		/* okay to use this color */
		if ( flag) {
			member->vc[crnt].color = c;
			if (crnt+1 < last) {
				if (c <=maxclr) 
				 dfsclr(member,first,last,vert,crnt+1, maxclr);
				else
				 dfsclr(member,first,last,vert,crnt+1, c);
			} else {
				if (maxclr < bestclr) {
					for(i=first;i<last;i++)
						best.vc[i] = member->vc[i];
					bestclr = maxclr;
				}
			}
		}
	}
}
#endif

extern void kempe(popmembertype *member);

void greedy(
popmembertype *member,
vertextype first,
vertextype last, /* color vertices in positions i, first <= i < last */
colortype oldcolor, /* previous coloring: MAXVERTEX if first coloring */
greedytype which,
int KEMPE) /* if 0 then kempe chain reduction invoked */
{
	colortype clr,maxclr;
	vertextype i,k,v;
	/* To increase efficiency, we use a data structure to store
	   a list of vertices of each color.  If p = 0.5,
	   then with prob. p we get a conflict after first test,
	  etc. On average we look at 2 vertices per color before rejecting.
	*/
	adjacencytype *x;


	if (which == STIRGRDY) {
		stirclr(member,first,last);
		return;
	}

	/* All color sets are empty initially */
	for(i=1;i<=order;i++) start[i] =order+1;

	initializeClrSelect(oldcolor,which);

	maxclr = 1;
#ifdef DFSCOLOR
	for(i=first;i<last-BRUTE;i++) 
#else
	for(i=first;i<last;i++) 
#endif
	{
		clr = selectFirstClr(which); 
		k = start[clr];
		v = member->vc[i].vertex;

		initnbr(x,v);
		while (k <= order) {
			/* no conflicts with clr but vertex at k untested */
			if (isnbr(x, k )) {
				/* conflict with clr */
				clr = nextClr(which); /* next clr */
				k = start[clr]; 
			}
			else k = next[k];  /* next vertex position to test */
		}
		/* k> order means clr had no conflicts so use it */
		member->vc[i].color = clr;
		if (clr >maxclr) maxclr = clr;
		
		/* add vertex to appropriate color list */
		next[v] = start[clr];
		start[clr] = v;
		updateClrSelect(which);
	}

#ifdef DFSCOLOR
	bestclr = maxclr+BRUTE;
#else
	bestclr=maxclr;
#endif
	if (bestclr > oldcolor) bestclr = oldcolor;
	fixedclr = maxclr;
#ifdef DFSGREEDY
	dfsclr(member,first,last,i,i,maxclr);
	for(i=first;i<last;i++)
		member->vc[i] = best.vc[i];
#endif

	/* At this point we have the lists of colors, so easy to
	   implement a Kempe chain reduction.
	*/

	if (KEMPE == 0) {
		kempe(member);
	}

}

int  bfscomponent(
vertextype v,vertextype clr1,vertextype clr2,vertextype flagset[],
vertextype visit,
int compnum)
/*
	find the compoent on v (must be in clr2) and return
	1 if more of clr2 than clr1, else return 0
	flagset will be set to the color of the vertex
	for each vertex in the component.
*/
{
	adjacencytype *x;

	vertextype cnt1,cnt2;
	vertextype queue1[MAXVERTEX+1];
	vertextype queue2[MAXVERTEX+1];
	int  first1,last1;
	int  first2,last2;
	vertextype k;

	/* We can assume that the color classes are much smaller
		than nbrhoods, so traverse the color classes */

	first2 = 0; last2 = 0; queue2[first2] = v;
	flagset[v] = compnum;

	cnt1 = 0;
	cnt2 = 0;

	while(first2 <= last2) {
		/* assume queue 1 is empty */
		last1 = -1;
		first1 = 0;
		
		while(first2 <= last2) {
			v = queue2[first2]; 
			first2++;
			cnt2++;
	
			initnbr(x,v);
			k = start[clr1];
			while (k <= order) {
				if (isnbr(x,k) ) {
					if (flagset[k] == visit) {
						flagset[k] = compnum;
						last1++;
						queue1[last1] = k;
					}
				}
				k = next[k];
			}
		}

		/* now empty queue2 and search through queue1 */
		first2 = 0;
		last2 = -1;
	
		while (first1 <= last1) {
			v = queue1[first1]; 
			first1++;
			cnt1++;
	
			initnbr(x,v);
			k = start[clr2];
			while (k <= order) {
				if (isnbr(x,k) ) {
					if (flagset[k] == visit) {
						flagset[k] = compnum;
						last2++;
						queue2[last2] = k;
					}
				}
				k = next[k];
			}
		}
	}
	if (cnt2 > cnt1) return(1); else return(0);
}

#ifdef DEBUGKEMPE
void printsets(
vertextype clr1, vertextype clr2, vertextype flagset[],
int components[], int compnum)
{
	vertextype k;

	printf("Components compnum=%d:\n",compnum);
	for(k=0;k<compnum;k++) {
		printf("\tk = %d components[k]=%1d\n",k,components[k]);
	}
	printf("Colorset 1 = %d\n",clr1);
	printf("vertex\t flagset\n");
	k = start[clr1];
	while (k <=order) {
		printf("%5d\t%5d\n",k,flagset[k]);
		k = next[k];
	}
	printf("Colorset 2 = %d\n",clr2);
	printf("vertex\t flagset\n");
	k = start[clr2];
	while (k <=order) {
		printf("%5d\t%5d\n",k,flagset[k]);
		k = next[k];
	}
}
#endif

	 

void Swap(
vertextype clr1, vertextype clr2, vertextype flagset[],
int components[])
{
	vertextype i,k;
	vertextype vset1[MAXVERTEX+1],l1;
	vertextype vset2[MAXVERTEX+1],l2;
	int c1swap,c2swap,c1same,c2same;

	/* swap the colors that match the flag set */
	
	l1 = 0; l2 = 0;

	/* okay its dumb, but I just delete them all off then put them
	all back on the appropriate lists */

	c1swap = 0;
	c1same = 0;
	k=start[clr1];
	while (k <= order) {
		if (components[flagset[k]] ) {
			c1swap++;
			vset2[l2] = k;
			l2++;
		} else {
			c1same++;
			vset1[l1] =k;
			l1++;
		}
		start[clr1] = next[k];
		k = start[clr1];
	} 
	
	c2swap =0;
	c2same=0;
	k=start[clr2];
	while (k <= order) {
		if (components[flagset[k]] ) {
			c2swap++;
			vset1[l1] = k;
			l1++;
		} else {
			c2same++;
			vset2[l2] =k;
			l2++;
		}
		start[clr2] = next[k];
		k = start[clr2];
	} 

#ifdef SWAPOUT
	if ( ((c1swap > 1)  && (c1same >1)) &&
		((c2swap > 1)  && (c2same >1)) ) {
		printf("Real swap clr1 c1swap = %d c1same = %d\n",
			c1swap,c1same);
		printf("Real swap clr2 c2swap = %d c2same = %d\n",
			c2swap,c2same);
	}
#endif

	for(i=0;i<l1;i++) {
		next[vset1[i]] = start[clr1];
		start[clr1] = vset1[i];
	}
	for(i=0;i<l2;i++) {
		next[vset2[i]] = start[clr2];
		start[clr2] = vset2[i];
	}
}
	
			

int CheckComponents( vertextype clr1,vertextype clr2)
/* check for components between color classes 1 and 2;
	if any clr1 is smaller than clr2, then swap the class 
	and return a flag value of one, else f no swaps return zero
*/
{
	vertextype flagset[MAXVERTEX+1];
	int components[MAXVERTEX+1];
	vertextype k;
	vertextype visit;
	int compnum,swaps;
	int cnt1,cnt2;
	
	/* initialize flagset as not visited */
	visit = order+1;

	/* sometimes color 1 will not all be visited in bfs */
	components[visit] = 0;

	cnt1 = 0;
	k = start[clr1];
	while (k <= order) {
		cnt1++; /* for debugging */
		flagset[k] = visit;
		k = next[k];
	}

	cnt2 = 0;
	k = start[clr2];
	while (k <= order) {
		cnt2++; /* for debugging */
		flagset[k] = visit;
		k = next[k];
	}
	

	compnum = 0;
	swaps = 0;
	k= start[clr2];
	while (k <= order) {
		if (flagset[k] == visit) {
			/* not visited yet */
			if ( (components[compnum] = 
			       bfscomponent(k,clr1,clr2,flagset,visit,compnum)))
					swaps = 1;
			compnum++;
		} 
		k = next[k];
	}
#ifdef DEBUGKEMPE
	if (compnum > 1) {
		printf("CheckComp: compnum=%d swaps=%d\n",compnum,swaps);
		printf("\tclr1 = %d clr2 = %d\n",clr1,clr2);
	}
#endif
	if (swaps) {
#ifdef DEBUGKEMPE
		if (compnum >1) {
		printf("call swap c1,c2 %d, %d\n",clr1,clr2);
		printf("Before:\n");
		printsets(clr1,clr2,flagset,components,compnum);
		printf("before swap cnt1= %d cnt2 = %d\n",cnt1,cnt2);
		}
#endif
		Swap(clr1,clr2,flagset,components);
#ifdef DEBUGKEMPE
		if (compnum >1) {
		printf("After:\n");
		printsets(clr1,clr2,flagset,components,compnum);
		fflush(stdout);
		}
#endif
		return(1);
	} else  return(0);
}

void kempe( popmembertype *member)
{

	vertextype c,i,j;

	if (bestclr>30) j= bestclr-30;
	else 
		j=1;
	for(c=j;c<=bestclr;c++){
		i = 1;
		while (i<c ){
			if (CheckComponents(i,c)) {
			}
			i++;
		}
	}
	/*
	printf("Component checking finished\n");
	fflush(stdout);
	*/
	/* now collect colorsets into member */
	j=0;
	for(c=1;c<=bestclr;c++) {
		i = start[c];
		while (i<=order) {
			member->vc[j].vertex = i;
			member->vc[j].color = c;
			j++;
			i = next[i];
		}
	}
}
