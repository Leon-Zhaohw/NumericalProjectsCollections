/*
 *	particle2.cpp
 *	
 *	Created by Ryoichi Ando on 11/9/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "particle2.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

std::list<particle2 *> particle2::removeList;

particle2::particle2 () {
	n = 1;
	levelset = -1.0;
	rm = false;
	r = 1.0;
	remesh_r = 0.0;
    ref = NULL;
	index = 0;
	merge = false;
	split = false;
	isolated = false;
	curvature[0] = 0.0;
	curvature[1] = 0.0;
}

// Remove particles which "remove" flag is turned on
bool particle2::cleanParticles( std::vector<particle2 *> &particles ) {
	// Now find particles to remove
    bool didRemove = false;
	for( uint n=0; n<particles.size(); n++ ) {
		if( particles[n]->rm ) {
			// Don't delete here, first just collect into the list
			removeList.push_back(particles[n]);
			particles[n] = NULL;
            didRemove = true;
		}
	}
	if( didRemove ) {
		uint cnt=0;
		for( uint n=0; n<particles.size(); n++ ) {
			if( particles[n] ) particles[cnt++] = particles[n];
		}
		particles.resize(cnt);
	}
	
	// Then, finally remove those particles
	list<particle2 *>::iterator it;
	for( it=removeList.begin(); it!=removeList.end(); it++ ) {
		delete *it;
	}
	removeList.clear();
	std::list<particle2 *>().swap(removeList);
    return didRemove;
}

void particle2::setRemovable( bool remove ) {
	rm = remove;
}

void particle2::computeRadius() {
	r = pow(n,1.0/DIM);
}