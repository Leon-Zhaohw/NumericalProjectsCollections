/*
 *	corrector3.cpp
 *	
 *	Created by Ryoichi Ando on 1/10/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "particle3.h"
#include "corrector3.h"
#include "levelset3.h"
#include "sorter3.h"
#include "vec3.h"
#include "kernel.h"
#include "util3.h"
using namespace std;

corrector3::corrector3() {
}

static bool sampleVelocity( sorter3& sorter, vec3d pos, FLOAT64 sample_visc, FLOAT64 dpx, const std::vector<particle3 *> &neighbors, vec3d &vel ) {
	FLOAT64 wsum = 0.0;
	vec3d out_velocity;
	FOR_EACH_PARTICLE(neighbors) {
		FLOAT64 w = p.n*kernel::sharp_kernel( (pos-p.p).len2(), sample_visc*p.r*dpx);
		out_velocity += w*p.u;
		wsum += w;
	} END_FOR
	if( wsum ) {
		vel = out_velocity/wsum;
	} else {
		FOR_EACH_PARTICLE(neighbors) {
			FLOAT64 w = 1.0;
			out_velocity += w*p.u;
			wsum += w;
		} END_FOR
		if( wsum ) {
			vel = out_velocity/wsum;
		}
	}
	return wsum > 0.0;
}

void corrector3::correct( sorter3 &sorter, int num, const levelset3 *level, std::vector<particle3 *> &particles, FLOAT64 stiffness, FLOAT64 sample_visc,
						 FLOAT64 dt, uint gsize, FLOAT64 dpx, const levelset3* fluid, const levelset3* solid, uint numSample ) {
	tick(); dump(">>> Correcting particle positions\n");
	
	// Sort particles
	sorter.sortParticles(particles);
    
	// If stiffness is zero, probably we don't want to correct particle positions at all...
	if( ! stiffness ) return;
	
	tick(); dump("Computing new particle positions (%d times)...", num);
	
	// Copy particle positions
	std::vector<vec3d> old_positions(particles.size());
	std::vector<vec3d> new_positions(particles.size());
	PARALLEL_FOR for( int n=0; n<particles.size(); n++ ) {
		old_positions[n] = new_positions[n] = particles[n]->p;
		particles[n]->index = n;
	}
	
	std::vector<std::vector<particle3 *> > neighbors_array(particles.size());
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		neighbors_array[pcount] = sorter.getNeighbors(p.p,1.01*sample_visc*dpx*powf(2,level->evalLevelset(p.p)-1));
	} END_FOR
	
	for( uint k=0; k<num; k++ ) {	
		// Now push with weak spring forces
		FLOAT64 st = 1.0/sqrt(DIM) * powf(2,DIM);
		PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
			vec3d f;
			
			// Search for nearby particles
			vec3d pos = old_positions[pcount];
			const std::vector<particle3 *> &neighbors = neighbors_array[pcount];
			FOR_EACH_NEIGHBOR_PARTICLE(neighbors) { // np <- neighbors
				if( &p == &np ) continue;
				FLOAT64 r = 0.5*(p.r+np.r)*dpx;
				vec3d rp = old_positions[p.index]-old_positions[np.index];
				FLOAT64 len = rp.len();
				if( len < r && len > 0) {
					vec3d force = r * stiffness * rp / len * kernel::smooth_kernel(len*len,r);
					if( force.len() > r ) force = force.normal() * r;
					f += force;
				}
			} END_FOR
			
			// Don't forget the virtual solid particles
			FLOAT64 phi = solid->evalLevelset(p.p);
			vec3d normal = solid->evalGradient(p.p).normal();
			FLOAT64 len = fabs(phi);
			vec3d rp = normal*len;
			FLOAT64 r = 0.5*dpx*p.r;
			if( len < r && len > 0 ) {
				vec3d force = st * r * stiffness * rp / len * kernel::smooth_kernel(len*len,r);
				if( force.len() > r ) force = force.normal() * r;
				f += force;
			}
			if( p.levelset > -p.r*dpx ) {
				vec3d grad = p.gradient;
				if( grad * f > 0.0 ) f -= grad*(f*grad);
			}
			new_positions[pcount] = old_positions[pcount] + f;			
#if 1
			// Move to surface as much as possible
			if( k==0 && p.levelset < -p.r*dpx && p.levelset > -3.0*p.r*dpx && solid->evalLevelset(p.p) > p.r*dpx ) {
				FLOAT64 min_neighbor = 1e9;
				FLOAT64 dist_surface = fabs(p.levelset)-0.5*dpx*p.r;
				vec3d avg_gradient;
				FOR_EACH_NEIGHBOR_PARTICLE(neighbors) {
					if( &np != &p ) avg_gradient += np.gradient;
				} END_FOR
				if( avg_gradient * p.gradient > 0.0 ) {
					FOR_EACH_NEIGHBOR_PARTICLE(neighbors) {
						if( &np != &p && (np.p-p.p) * p.gradient > 0.0 ) {
							min_neighbor = fmin(min_neighbor,(np.p-p.p).len()-0.5*dpx*(p.r+np.r));
						}
					} END_FOR
					if( min_neighbor > 0 ) {
						FLOAT64 displace = fmin(min_neighbor,0.5*dist_surface);
						if( displace > 0.1*dpx*p.r ) new_positions[pcount] += fmin(0.5*dpx*p.r,stiffness*displace)*p.gradient;
					}
				}
			}
#endif
		} END_FOR
		
		// Assign new particle positions
		PARALLEL_FOR for( int n=0; n<particles.size(); n++ ) {
			old_positions[n] = new_positions[n];
			for( uint dim=0; dim<DIM; dim++ ) new_positions[n][dim] = fmin(1.0,fmax(0.0,new_positions[n][dim]));
			particles[n]->levelset = fluid->evalLevelset(new_positions[n]);
			particles[n]->gradient = fluid->evalGradient(new_positions[n]);
		}
	}
	
	dump("Done. Took %s.\n", stock("correct_pos"));
	tick(); dump("Computing new particle velocities...");
	
	// Copy particle velocities
	std::vector<vec3d> velocities(particles.size());
	PARALLEL_FOR for( int n=0; n<particles.size(); n++ ) {
		velocities[n] = particles[n]->u;
	}
	
	// Resample new particle velocity
	if( sample_visc > 0 ) {
		PARALLEL_FOR for( int n=0; n<particles.size(); n++ ) {
			sampleVelocity(sorter,new_positions[n],sample_visc,dpx,neighbors_array[n],velocities[n]);
		}
	}
	
	// Assign new particle velocities and positions
	PARALLEL_FOR for( int n=0; n<particles.size(); n++ ) {
		particles[n]->p = new_positions[n];
		particles[n]->u = velocities[n];
	}
	dump("Done. Took %s.\n", stock("correct_resample"));
	
	// Don't forget to sort particles !
	if( sample_visc ) sorter.setDirty();
	
	dump("<<< Done. Took %s.\n", stock("correct"));
}