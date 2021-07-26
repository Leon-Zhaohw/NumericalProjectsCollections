/* COLOR ROUTINES */

/*
        Title: Color Routines Source File.
        file: colorrtns.c
	does:
        A Collection of routines for manipulating vectors
        of colors, vertices, degree sequences etc.
        used here and there throughout the system

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
#include "graph.h"
#include "colorrtns.h"

/*
	A Collection of routines for manipulating vectors
	of colors, vertices, degree sequences etc.
	used here and there throughout the system
*/


void printinfo( popmembertype *member)
{
	printf("CLRS =%d\tCLRSUM = %d\n",member->clrdata.numcolors,
		member->clrdata.total);
}

/* maximum colors printed per line */
#define MAXPERLINE  20

void printcoloring( popmembertype *member)
/*
	Prints the sequence of colors, in order of vertices 
	first to last 
*/
{
	colortype c[MAXVERTEX];
	int i;
	int fmt;

	for(i=0;i<order;i++) {
		c[member->vc[i].vertex] = member->vc[i].color;

	}
	fmt = 0;
	for(i=0;i<order;i++) {
		printf("%3d ",c[i]);
		fmt++;
		if ( fmt >= MAXPERLINE ) {
			printf("\n");
			fmt = 0;
		}
	}
}

void printpurity( popmembertype *member)
{
        int i,j, purity;
	int purityset[MAXVERTEX];
        vertextype v,w;
        colortype cv,cw;

	if (partitionflag) {
        	purity = 0;
		for(i=0;i<=partitionnumber;i++) purityset[i] =0;
        	for (i=0;i<order;i++) {
                	v= member->vc[i].vertex;
                	cv=member->vc[i].color;
                	for(j=0;j<i;j++) {
                        	w= member->vc[j].vertex;
                        	cw=member->vc[j].color;
                        	if ((partset[v] == partset[w]) && (cv == cw)){
					purity++;
					purityset[partset[v]]++;
                        	}
                	}
        	}
        	printf("P= %4d:",purity);
		for(i=0;i<=partitionnumber;i++)
			printf("%3d ",purityset[i]);
		printf("\n");

	}
}


void getcolorinfo( popmembertype *member)
{
	vertextype i;
	struct vrtxandclr *v;

	/* get color info */
	member->clrdata.numcolors = 1;
	member->clrdata.total = 0;
	v = member->vc;

	for(i=0;i<order;i++) {
		if (v[i].color > member->clrdata.numcolors) 
			member->clrdata.numcolors = v[i].color;
		member->clrdata.total += v[i].color;
	}
	member->clrdata.total += member->clrdata.numcolors * order;
}

/* permute a subrange of a member 			*/
/* permutes the elements in the range [first..last-1] 	*/
/* all permutations equally likely 			*/
void permute(
popmembertype *member,
vertextype first, vertextype last)
{
	vertextype i,j;
	struct vrtxandclr k;
	
	for(i=first; i<last-1; i++) {
		j = i + ((vertextype) random() % (last - i));
		k = member->vc[i];
		member->vc[i] = member->vc[j];
		member->vc[j] = k;
	}
}

void trivial_color( popmembertype *m)
{
        vertextype i;
        for(i=0;i<order;i++)
                m->vc[i].color = i+1;
}

void verifycolor(popmembertype *m)
{
	char verifyset[MAXVERTEX];
	int clrflagerror;
	int i,j;

	for(i=0;i<order;i++) verifyset[i] = 0;
	clrflagerror = 0;
	for(i=0;i<order;i++) {
		verifyset[m->vc[i].vertex] = 1;
		for(j=i+1;j<order;j++) {
			if ((edge(m->vc[i].vertex,m->vc[j].vertex)) && 
				(m->vc[i].color == m->vc[j].color))
					clrflagerror++;
		}
	}
	
	if (clrflagerror > 0) printf("COLORING ERRORS %d\n",clrflagerror);
	else {
		clrflagerror = 0;
		for(i=0;i<order;i++) if (verifyset[i] == 0) clrflagerror++;
		if (clrflagerror > 0) 
			printf("ERROR: %d missing vertices\n",clrflagerror);
		else
			printf("Coloring Verified\n");
	}
			
}

void getacoloring(popmembertype *m, char *name, int *which)
{
	int i,cnt, colnum;

	FILE *fp;

	char fname[100], oneline[100],*s;

	
	for(i=0;i<order;i++)
		m->vc[i].vertex = i;

	strncpy(fname,name,80);
	strncat(fname,".res",5);

	if ( (fp = fopen(fname,"r") ) == NULL) {
		printf("WARNING: getacoloring - cannot open file %s \n",
			fname);
		printf("\tUSING TRIVIAL COLORING INSTEAD\n");
		trivial_color(m);
	} else {
		cnt = 0;
		while ( !feof(fp) ) {
			s = fgets( oneline, 81 , fp);
			if (0 == strncmp(oneline, "CLRS",4)) {
				cnt++;
				printf("[ %d] %s\n",cnt, oneline);
			}
		}
		
		if (cnt > 0) {
			printf("There are %d colorings. Which do you want? ",
				cnt);
			scanf("%d",&colnum);

			/* default take the last */
			if ( (colnum > cnt) || (colnum < 1) ) colnum = cnt;
			printf("%d\n",colnum);

			*which = colnum;

			rewind(fp);
			cnt = 0;
			while (cnt < colnum) {
				s = fgets( oneline, 81 , fp);
				if (0 == strncmp(oneline, "CLRS",4)) {
					cnt++;
				}
			}
			for(i=0;i<order;i++) {
				if ( 0 == fscanf(fp,"%hu",&(m->vc[i].color)) ) {
					printf("BAD FILE FORMAT\n");
					exit(1);
				}
			}
		} else {
			printf("NO COLORINGS PRESENT IN FILE\n");
			printf("Using Trivial Color\n");
			trivial_color(m);
		}
		fclose(fp);
	}
			
}

void printperm(popmembertype *m)
{
	int i;
	for(i=0;i<order;i++)
		printf("%d %d\n",m->vc[i].vertex, m->vc[i].color);
}
		

int degseq[MAXVERTEX+1];

int decdeg( struct vrtxandclr *a, struct vrtxandclr *b)
/*
        Comparison routine for sorting by degree downwards
*/
{
        if (degseq[a->vertex] < degseq[b->vertex]) return(1);
        else if (degseq[a->vertex] == degseq[b->vertex]) return(0);
        else return(-1);
}

int incdeg(struct vrtxandclr *a, struct vrtxandclr *b)
/*
        Comparison routine for sorting by degree upwards
*/
{
        if (degseq[a->vertex] > degseq[b->vertex]) return(1);
        else if (degseq[a->vertex] == degseq[b->vertex]) return(0);
        else return(-1);
}

void computedeg()
/*
         Compute degree sequence.
*/

{
        int i,j;
        adjacencytype *x;
        for(i=0;i<order;i++)
                degseq[i] =0;
        for(i=1;i<order;i++) {
                initnbr(x,i);
                for(j=0;j<i;j++)
                        if (isnbr(x,j)) {
                                degseq[i]++;
                                degseq[j]++;
                        }
        }
}

void fileres(char *name, popmembertype *m, char *prgm)
{
	char newname[100];
	FILE *fp;

	int i;
	colortype c[MAXVERTEX];
	int fmt;

	strncpy(newname,name,95);
	strncat(newname,".res",5);

	if ( (fp = fopen(newname,"a")) == NULL) {
		printf("unable to open file %s for update\n",newname);
		printf("dumping coloing information to standard output");
		printinfo(m);
		printcoloring(m);
		exit(1);
	}
	fprintf(fp,"CLRS %d FROM %s\n",m->clrdata.numcolors, prgm);
	
	for(i=0;i<order;i++) {
                c[m->vc[i].vertex] = m->vc[i].color;
	}

	fmt = 0;
        for(i=0;i<order;i++) {
                fprintf(fp,"%3d ",c[i]);
                fmt++;
                if ( fmt >= MAXPERLINE ) {
                        fprintf(fp,"\n");
                        fmt = 0;
                }
        }
	if (fmt != 0) fprintf(fp,"\n");

	fclose(fp);
}
	
	
void about(char *pgrmname)
{
	printf("J. Culberson's Implemenation of\n");
	printf("\t\t%s\n",pgrmname);
	printf("A program for coloring graphs.\n");

	printf("For more information visit the webpages at: \n\n");

 	printf("\thttp://www.cs.ualberta.ca/~joe/Coloring/index.html\n");
	printf("\nThis program is available for research ");
	printf("and educational purposes only.\n");
	printf("There is no warranty of any kind.\n");
/*
	printf("This research is supported in part by a grant from the\n");
	printf("Natural Sciences and Engineering Research Council.\n");
*/
	printf("\n\tEnjoy!\n\n");
}
