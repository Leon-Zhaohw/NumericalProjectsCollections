/*
 *	ann2.cpp
 *
 *	Created by Ryoichi Ando on 9 Nov 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "ANN.h"
#include "vec2.h"
#include "util2.h"

#ifndef _ANN2_H
#define _ANN2_H

class ann2 {
public:
	ann2() {
		numbers = 0;
		kdtree = NULL;
	}
	virtual ~ann2() {
		clear();
	}
	void clear() {
		if( numbers ) {
			delete kdtree;
			delete [] dataPts[0];
			delete [] dataPts;
		}
		numbers = 0;
	}
	void sort( const std::vector<vec2d> &array ) {
		clear();
		numbers = array.size();
		if( ! numbers ) return;
	
		// Create data point set
		dataPts = annAllocPts(numbers,DIM);
	
		// Assign point information
		for( uint n=0; n<numbers; n++ ) {
			for( uint dim=0; dim<DIM; dim++) dataPts[n][dim] = array[n][dim];
		}
	
		// Build KD Tree
		kdtree = new ANNkd_tree(dataPts,numbers,DIM);
	}

	std::vector<ANNidx> getNeighbors( vec2d p, FLOAT64 r ) const {
		if( ! numbers ) return std::vector<ANNidx>();
		ANNpoint queryPt = annAllocPt(DIM);
		for( uint dim=0; dim<DIM; dim++ ) queryPt[dim] = p[dim];
		std::vector<ANNidx> neighbors_idx = kdtree->annkFRSearch(queryPt,r*r,0.0);
		delete [] queryPt;
		return neighbors_idx;
	}

	std::vector<ANNidx> getkNeighbors( vec2d p, uint n ) const {
		std::vector<ANNidx> res(n);
		if( ! numbers ) return std::vector<ANNidx>();
		
		ANNpoint queryPt = annAllocPt(DIM);
		for( uint dim=0; dim<DIM; dim++ ) queryPt[dim] = p[dim];
		ANNidx *nnIdx = new ANNidx[n];
		ANNdist *dists = new ANNdist[n];
		kdtree->annkSearch(queryPt,n,nnIdx,dists,0.0);
		for( uint k=0; k<n; k++ ) {
			res[k] = nnIdx[k];
		}
		delete [] queryPt;
		delete [] nnIdx;
		delete [] dists;
		return res;
	}
	uint numbers;
	ANNpointArray dataPts;
	ANNkd_tree * kdtree;
};

#endif