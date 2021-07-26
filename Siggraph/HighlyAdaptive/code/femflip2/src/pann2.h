/*
 *	pann2.h
 *
 *	Created by Ryoichi Ando on 8 Nov 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <stdio.h>
#include <vector>
#include <math.h>
#include "particle2.h"
#include "ANN.h"

#ifndef _PANN2_H
#define _PANN2_H

class pann2 {
public:
	pann2 () {
		max_level = 0;
	}
	void buildKDTree( const std::vector<particle2 *> particles ) {
		// Remove all the existing kdtrees
		for( uint n=0; n<max_level; n++ ) {
			if( numbers[n] ) delete kdtrees[n];
		}
		kdtrees.resize(0);
		// Remove all the exisiting dataPoints
		for( uint n=0; n<max_level; n++ ) {
			if( numbers[n] ) {
				delete [] dataPts[n][0];
				delete [] dataPts[n];
			}
			numbers[n] = 0;
		}
		max_level = 0;
		if( particles.empty() ) return;
		
		// Count the maximum particle radius
		std::vector<uint> indices;
		
		for( uint n=0; n<particles.size(); n++ ) {
			max_level = fmax(particles[n]->r,max_level);
		}
	
		// Allocate the KDtree
		kdtrees.resize(max_level);
		numbers.resize(max_level);
		indices.resize(max_level);
		
		for( uint level=0; level<max_level; level++ ) {
			kdtrees[level] = NULL;
			numbers[level] = 0;
			indices[level] = 0;
		}
		
		// Count numbers for each level
		for( uint n=0; n<particles.size(); n++ ) {
			uint r = particles[n]->r-1;
			numbers[r] ++;
		}
		
		// Create data point set
		dataPts.resize(max_level);
		for( uint n=0; n<max_level; n++ ) {
			if( numbers[n] ) dataPts[n] = annAllocPts(numbers[n],DIM);
		}
		
		// Assign point information
		for( uint n=0; n<particles.size(); n++ ) {
			uint r = particles[n]->r-1;
			for( uint dim=0; dim<DIM; dim++) dataPts[r][indices[r]][dim] = particles[n]->p[dim];
			indices[r] ++;
		}
		
		// Build KDTree
		for( uint n=0; n<max_level; n++ ) {
			if( numbers[n] ) kdtrees[n] = new ANNkd_tree(dataPts[n],numbers[n],DIM);
		}
	}
	FLOAT64 evalAbsLevelset( vec2d pos, FLOAT64 dpx ) const {
		FLOAT64 min_phi = 1e6;
		if( ! max_level ) return min_phi;
		
		for( uint n=0; n<max_level; n++ ) {
			if( numbers[n] ) {
				FLOAT64 radius = dpx*(n+1);
				ANNkd_tree *kdtree = kdtrees[n];
				ANNpoint queryPt = annAllocPt(DIM);
				for( uint dim=0; dim<DIM; dim++ ) queryPt[dim] = pos[dim];
				ANNidx nnIdx;
				ANNdist dist;
				kdtree->annkSearch(queryPt,1,&nnIdx,&dist,0.0);
				delete [] queryPt;
				min_phi = fmin(min_phi,sqrt(dist)-0.5*radius);
			}
		}
		return min_phi;
	}
	uint max_level;
	std::vector<uint> numbers;
	std::vector<ANNpointArray> dataPts;
	std::vector<ANNkd_tree *> kdtrees;
};

#endif