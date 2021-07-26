/*
 *	particle3.h
 *	
 *	Created by Ryoichi Ando on 1/8/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "vec3.h"
#include "svd3.h"
#include <vector>
#include <list>

#ifndef _PARTICLE3_H
#define _PARTICLE3_H

class particle3;
class particle3 {
public:
	particle3 ();
	static bool cleanParticles( std::vector<particle3 *> &particles );	// Call this to release removed particles from RAM when possible
	void setRemovable( bool remove );	// Call this function to set particle removable when cleanParticles() is called
										// When called, the particle is immediately pop out from the sim particle list, but stays
										// in RAM if the reference counter is not zero.
										// Don't worry, just call this to remove the particle from the sim.
										// Also don't forget to occasionally call cleanParticles() to release memory.
	
	void computeRadius();				// Recompute particle radius based on the number of contained particles
	vec3d p;							// Position
	vec3d u;							// Velocity
	FLOAT64 levelset;					// Fluid levelset
	vec3d gradient;						// Fluid levelset gradient
	FLOAT64 curvature[2];				// Curvature info 0: current 1: previous timestep
	FLOAT64 r;							// Particle radius (particle radius of smallest particle should be multiplied for actual use)
	FLOAT64 remesh_r;					// Previous remeshing radius
	uint n;								// Number of merged particles
	uint index;							// Index
	
	// Variables below are used as temporary space
    particle3 *ref;						// Reference to the closest particle
	bool rm;							// Remove flag
	bool merge;							// Merge flag
	bool split;							// Split flag
	bool isolated;						// Splash flag
	static std::list<particle3 *> removeList; // Temporary storage to track particle to remove
};

#endif