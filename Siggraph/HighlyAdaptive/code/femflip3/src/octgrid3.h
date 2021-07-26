/*
 *	octgrid3.h
 *
 *	Created by Ryoichi Ando on Dec 25 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "octree3.h"

#ifndef _OCTGRID3_H
#define _OCTGRID3_H

template <class T> class octgrid3 {
public:
	octgrid3() { octree = NULL; }
	virtual void setOctree( const octree3 &octree ) {
		this->octree = &octree;
	}
	virtual void setValue( const std::vector<T> &q ) {
		this->q = q;
	}
	virtual T getValue(vec3d p) const {
		for( uint dim=0; dim<DIM; dim++ ) p[dim] = fmin(1.0-1e-10,fmax(1e-10,p[dim]));
		int hit = octree->hitTest(p);
		if( hit >= 0 ) {
			octree3::leaf3 *leaf = octree->terminals[hit];
			FLOAT64 dx = leaf->dx/(FLOAT64)octree->resolution;
			vec3d origin = octree->nodes[leaf->corners[0][0][0]];
			FLOAT64 x = (p[0]-origin[0])/dx;
			FLOAT64 y = (p[1]-origin[1])/dx;
			FLOAT64 z = (p[2]-origin[2])/dx;
			T value = (1.0-z)*((1.0-y)*((1.0-x)*q[leaf->corners[0][0][0]]+x*q[leaf->corners[1][0][0]])+y*((1.0-x)*q[leaf->corners[0][1][0]]+x*q[leaf->corners[1][1][0]]))+
			(z)*((1.0-y)*((1.0-x)*q[leaf->corners[0][0][1]]+x*q[leaf->corners[1][0][1]])+y*((1.0-x)*q[leaf->corners[0][1][1]]+x*q[leaf->corners[1][1][1]]));
			return value;
		}
		printf( "octlevelset3: hit test failed (%f,%f,%f)\n", p[0],p[1],p[2]);
		return T();
	}
	const octree3 *octree;
	std::vector<T> q;
};

#endif
