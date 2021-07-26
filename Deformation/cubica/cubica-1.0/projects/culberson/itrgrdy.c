/*
	ITERATED GREEDY
*/

/*
        Title: Iterated Greedy  Source File.
        file: itrgrdy.c
        does: itrgrdy routine
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
#include "colorrtns.h"
#include "graph.h"
#include "itrgrdy.h"


/*
	If two colors get the same total, then color number may increase

*/
	
/*
	ROLLOVER WARNING
	if there are more than about 2000 vertices of the same color
	then sorting by size may cause int roll-over on 32 bit
	machines
*/


/* 
	Reorder a colored permutation so that the vertices 
	of a color are grouped. 
*/

/* keep the best seen so far from reorder */
popmembertype bestmem,overallbest;

int permvec[MAXVERTEX+1];


/* interrupt processes */
int stopflag = 0; /* reset to terminate */

void cpulimig()
{
	printf("CPU TIME LIMIT EXCEEDED -- let me clean up\n");
	fflush(stdout);
	stopflag = 1;
}

void intfcig()
{
	printf("INTERRUPT DETECTED -- let me clean up\n");
	fflush(stdout);
	stopflag = 1;
}

/* FUNCTIONS */

/*
	See itrgrdy.h for a description of these functions
*/

/* comparison function - sort by color group using permvec */
int ccompare( struct vrtxandclr *a, struct vrtxandclr *b)
{
	if (permvec[a->color] > permvec[b->color]) return(1);
	else if (permvec[a->color] == permvec[b->color]) return(0);
	else return(-1);
}

void setvec(int size)
{
	int c;
	for(c=0;c<=size;c++)
		permvec[c] = c;
}

void revvec(int size)
{
	int c;
	
	for(c=0;c<=size;c++)
		permvec[c] = 1+size -c;
}

void blckcount( popmembertype *mem)
{
	vertextype i;
	
	/* if the following is set then sorting is first by group size */
	for(i=0;i<order;i++)
		permvec[mem->vc[i].color] += CNTR;

}

void revblckcount( popmembertype *mem)
{
	vertextype i;
	
	/* if the following is set then sorting is first by largest group */
	for(i=0;i<order;i++)
		permvec[mem->vc[i].color] -= CNTR;

}

void updegcnt( popmembertype *mem)
{
	vertextype i;
	
	/* if the following is set then sorting is by increasing total deg*/
	for(i=0;i<order;i++)
		permvec[mem->vc[i].color] += degseq[mem->vc[i].vertex]*DEGCNT;

}

void dwndegcnt( popmembertype *mem)
{
	vertextype i;
	
	/* if the following is set then sorting is by decreasing total deg*/
	for(i=0;i<order;i++)
		permvec[mem->vc[i].color] -= degseq[mem->vc[i].vertex]*DEGCNT;

}

void shufflevec( int start,int num, int permvec[])
{
	int i,j,k;
	for(i=start;i<num;i++) {
		j = i + (random() % (1+num-i));
		k = permvec[i];
		permvec[i] = permvec[j];
		permvec[j] = k;
	}
}


void itrgrdy(
popmembertype *member,	/* contains vertex permutation and coloring */
int	maxiter,	/* maximum iterations without improvement */
int	out_freq,	/* how often to output data */
int	controlset[],	/* which control */
int	controlselect[], /* density vector for controls */
int	controltot,	/* total for selection limit */
int	switchback, 
int	retrylimit, /* retry controls from best so far */
colortype targetclr,	/* if this coloring achieved stop */
greedytype greedyWt[],
int totalgreedywt,
int kempe_freq)
{
	int mintotal;
	colortype minclr;
	int nbiter,bestiter;
	int tsecs, tmicros;
	int cntrlidx,control;
	int retry;
	int i;
	int findWt;
	greedytype greedyindex;
	int kempe_val;

#ifdef DEBUG
	if(DEBUG & 1)
		printinfo(member);
	if (DEBUG & 2) {
		printperm(member);
		printcoloring(member);
	}
#endif

	stopflag = 0;
	getcolorinfo(member);
	minclr = member->clrdata.numcolors;
	mintotal = member->clrdata.total;

	for(i=0;i<order;i++) bestmem.vc[i] = member->vc[i];
	bestmem.clrdata = member->clrdata;

	printf("\n---------\n      ST: ");
	printinfo(&bestmem);

	printpurity(member);


	nbiter=0; bestiter =0;
	retry = 0;

	kempe_val = 0;

	while (((nbiter - bestiter) < maxiter) && (minclr > targetclr)
			&& (!stopflag) ) {
		nbiter++;
		cntrlidx = random() % controltot;
		i = 0;
		while (controlselect[i] < cntrlidx) i++;
		control = controlset[i];
		
		if (control & 1)  revvec((int) member->clrdata.numcolors);
		else setvec((int) member->clrdata.numcolors);

		if (control & 2) shufflevec((int) 1,
			(int) member->clrdata.numcolors, permvec);

		if (control & 4) revblckcount(member);

		if (control & 8) blckcount(member);

		if (control & 16) updegcnt(member);

		if (control & 32) dwndegcnt(member);


		qsort((char *) member->vc,(int) order,
			sizeof(struct vrtxandclr),(compfunc)ccompare);

#ifdef DEBUG
		if (DEBUG & 2) {
			printperm(member);
		}
#endif
		/* select a greedy type by proportional selection */
		findWt = random() % totalgreedywt;
		greedyindex =SIMPLEG;
		while (0 <= (findWt -= greedyWt[greedyindex])) greedyindex++;


		kempe_val = (kempe_val +1 ) % kempe_freq;
		greedy(member, 0, order, MAXVERTEX, greedyindex, kempe_val);
		/* rmaxis(member,1); maximum set respecting coloring*/


		getcolorinfo(member);

		if (minclr > member->clrdata.numcolors) {
			mintotal = member->clrdata.total;
			minclr = member->clrdata.numcolors;
			bestiter = nbiter;
			retry = 0;
			for(i=0;i<order;i++) bestmem.vc[i] = member->vc[i];
			bestmem.clrdata = member->clrdata;
			if (kempe_val == 0) printf("**");
			printf("It=%5d: ",nbiter);
			printinfo(member);
			printpurity(member);

			i = getrusage(RUSAGE_SELF,&tmp);
			tsecs = tmp.ru_utime.tv_sec - seconds;
			tmicros = tmp.ru_utime.tv_usec - microsecs;
			printf("COLOR %3d IT= %5d CPU= %5.2f\n",minclr,nbiter, 
				tsecs+(tmicros/1000000.0));
			fflush(stdout);
			/* if the global has improved use it next */
		}
		else if (mintotal > member->clrdata.total) {
			mintotal = member->clrdata.total;
			bestiter = nbiter;
			retry = 0;
			for(i=0;i<order;i++) bestmem.vc[i] = member->vc[i];
			bestmem.clrdata = member->clrdata;
			if (kempe_val == 0) printf("**");
			printf("It=%5d: ",nbiter);
			printinfo(member);
			fflush(stdout);
		}


		if ( ( nbiter % out_freq) == 0) {
			if( ( nbiter % (out_freq*10)) == 0) printpurity(member);
			printf("ITR %5d ",nbiter);
			printinfo(member);
			fflush(stdout);
		}

		if ( (switchback < (nbiter - bestiter))
		  && (retry < retrylimit) ) {
			/* revert to best again */
			for(i=0;i<order;i++) member->vc[i] = bestmem.vc[i];
			member->clrdata = bestmem.clrdata;
			retry++;
		}

#ifdef DEBUG
		if (DEBUG & 2) {
			printcoloring(member);
			printf("\n");
		}
		if(DEBUG & 1)
			printinfo(member);
#endif
	}
#ifdef DEBUG
	printf("---\n");
	fflush(stdout);
#endif
	printf("Final iterations = %d\n",nbiter);
	printpurity(member);
}


