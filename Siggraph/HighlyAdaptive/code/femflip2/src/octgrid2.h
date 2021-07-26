/*
 *	octgrid2.h
 *
 *	Created by Ryoichi Ando on Dec 25 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "octree2.h"

#ifndef _OCTGRID2_H
#define _OCTGRID2_H

template <class T> class octgrid2 {
public:
	octgrid2() { octree = NULL; }
	virtual void setOctree( const octree2 &octree ) {
		this->octree = &octree;
	}
	virtual void setValue( const std::vector<T> &q ) {
		this->q = q;
	}
	virtual T getValue(vec2d p) const {
		for( uint dim=0; dim<DIM; dim++ ) p[dim] = fmin(1.0-1e-10,fmax(1e-10,p[dim]));
		int hit = octree->hitTest(p);
		if( hit >= 0 ) {
			octree2::leaf2 *leaf = octree->terminals[hit];
			FLOAT64 dx = leaf->dx/(FLOAT64)octree->resolution;
			vec2d origin = octree->nodes[leaf->corners[0][0]];
			FLOAT64 x = (p[0]-origin[0])/dx;
			FLOAT64 y = (p[1]-origin[1])/dx;
			T value = (1.0-y)*((1.0-x)*q[leaf->corners[0][0]]+x*q[leaf->corners[1][0]])+y*((1.0-x)*q[leaf->corners[0][1]]+x*q[leaf->corners[1][1]]);
			return value;
		}
		printf( "octlevelset2: hit test failed (%f,%f)\n", p[0],p[1]);
		return T();
	}
	const octree2 *octree;
	std::vector<T> q;
};

#endif
