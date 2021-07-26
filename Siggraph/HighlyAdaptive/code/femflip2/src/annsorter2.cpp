/*
 *	annsorter2.cpp
 *
 *	Created by Ryoichi Ando on 9 Nov 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "annsorter2.h"
#include "particle2.h"
#include "util2.h"

static void warnAndDie() {
	email::print("Sorter operation was attemped while particles are dirty !\n" );
	email::print("Please sort particles.\n");
	email::send();
	exit(0);
}

annsorter2::annsorter2() {
	dirty = true;
}

void annsorter2::sortParticles( std::vector<particle2 *> &particles ) {
	if( ! dirty ) return;
	this->particles = particles;
	std::vector<vec2d> array(particles.size());
	for( uint n=0; n<particles.size(); n++ ) {
		array[n] = particles[n]->p;
	}
	ann.sort(array);
	dirty = false;
}

std::vector<particle2 *> annsorter2::getNeighbors( vec2d p, FLOAT64 r ) const {
	if( dirty ) warnAndDie();
	std::vector<ANNidx> neighbors_idx = ann.getNeighbors(p,r);
	std::vector<particle2 *> neighbors(neighbors_idx.size());
	for( uint n=0; n<neighbors_idx.size(); n++ ) {
		neighbors[n] = particles[neighbors_idx[n]];
	}
	return neighbors;
}

std::vector<particle2 *> annsorter2::getkNeighbors( vec2d p, uint n ) const {
	std::vector<particle2 *> res(n);
	std::vector<ANNidx> idx = ann.getkNeighbors(p,n);
	for( uint k=0; k<n; k++ ) {
		res[k] = particles[idx[k]];
	}
	return res;
}

void annsorter2::setDirty( bool isDirty ) {
	dirty = isDirty;
}