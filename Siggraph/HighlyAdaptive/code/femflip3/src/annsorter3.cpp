/*
 *	annsorter2.cpp
 *
 *	Created by Ryoichi Ando on 9 Nov 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "annsorter3.h"
#include "particle3.h"
#include "util3.h"

static void warnAndDie() {
	email::print("Sorter operation was attemped while particles are dirty !\n" );
	email::print("Please sort particles.\n");
	email::send();
	exit(0);
}

annsorter3::annsorter3() {
	numbers = 0;
	kdtree = NULL;
	dirty = true;
}

void annsorter3::sortParticles( std::vector<particle3 *> &particles ) {
	if( ! dirty ) return;
	if( numbers ) {
		delete kdtree;
		delete [] dataPts[0];
		delete [] dataPts;
	}
	numbers = particles.size();
	this->particles = particles;
	if( ! numbers ) return;
	
	// Create data point set
	dataPts = annAllocPts(numbers,DIM);
	
	// Assign point information
	for( uint n=0; n<numbers; n++ ) {
		for( uint dim=0; dim<DIM; dim++) dataPts[n][dim] = particles[n]->p[dim];
	}
	
	// Build KD Tree
	kdtree = new ANNkd_tree(dataPts,numbers,DIM);
	dirty = false;
}

std::vector<particle3 *> annsorter3::getNeighbors( vec3d p, FLOAT64 r ) const {
	if( dirty ) warnAndDie();
	ANNpoint queryPt = annAllocPt(DIM);
	for( uint dim=0; dim<DIM; dim++ ) queryPt[dim] = p[dim];
	std::vector<ANNidx> neighbors_idx = kdtree->annkFRSearch(queryPt,r*r,0.0);
	std::vector<particle3 *> neighbors(neighbors_idx.size());
	for( uint n=0; n<neighbors_idx.size(); n++ ) {
		neighbors[n] = particles[neighbors_idx[n]];
	}
	delete [] queryPt;
	return neighbors;
}

std::vector<particle3 *> annsorter3::getkNeighbors( vec3d p, uint n ) const {
	std::vector<particle3 *> res(n);
	ANNpoint queryPt = annAllocPt(DIM);
	for( uint dim=0; dim<DIM; dim++ ) queryPt[dim] = p[dim];
	ANNidx *nnIdx = new ANNidx[n];
	ANNdist *dists = new ANNdist[n];
	kdtree->annkSearch(queryPt,n,nnIdx,dists,0.0);
	for( uint k=0; k<n; k++ ) {
		res[k] = particles[nnIdx[k]];
	}
	delete [] queryPt;
	delete [] nnIdx;
	delete [] dists;
	return res;
}

void annsorter3::setDirty( bool isDirty ) {
	dirty = isDirty;
}