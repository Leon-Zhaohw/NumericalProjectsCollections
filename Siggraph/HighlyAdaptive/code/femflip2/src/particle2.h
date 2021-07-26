/*
 *	particle2.h
 *	
 *	Created by Ryoichi Ando on 11/3/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "vec2.h"
#include "svd2.h"
#include <vector>
#include <list>

#ifndef _PARTICLE2_H
#define _PARTICLE2_H

class particle2;
class particle2 {
public:
	particle2 ();
	static bool cleanParticles( std::vector<particle2 *> &particles );	// Call this to release removed particles from RAM when possible
	void setRemovable( bool remove );	// Call this function to set particle removable when cleanParticles() is called
										// When called, the particle is immediately pop out from the sim particle list, but stays 
										// in RAM if the reference counter is not zero. 
										// Don't worry, just call this to remove the particle from the sim.
										// Also don't forget to occasionally call cleanParticles() to release memory.
	
	void computeRadius();				// Recompute particle radius based on the number of contained particles
	vec2d p;							// Position
	vec2d u;							// Velocity
	FLOAT64 levelset;					// Fluid levelset
	vec2d gradient;						// Gradient of levelset
	FLOAT64 curvature[2];				// Curvature info 0: current 1: previous timestep
	FLOAT64 r;							// Particle radius (particle radius of smallest particle should be multiplied for actual use)
	FLOAT64 remesh_r;					// Previous remesh radius
	uint n;								// Number of merged particles
	svd2 svd;							// SVD information
	uint index;							// Index
	
	// Variables below are used as temporary space
    particle2 *ref;						// Reference to the closest particle
	bool rm;							// Remove flag
	bool merge;							// Merge flag
	bool split;							// Split flag
	bool isolated;						// Splash flag
	static std::list<particle2 *> removeList; // Temporary storage to track particle to remove
};

#endif