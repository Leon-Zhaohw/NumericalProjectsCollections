/*
 *	corrector2.cpp
 *	
 *	Created by Ryoichi Ando on 11/6/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "particle2.h"
#include "corrector2.h"
#include "levelset2.h"
#include "sorter2.h"
#include "surf2.h"
#include "vec2.h"
#include "kernel.h"
#include "util2.h"
using namespace std;

corrector2::corrector2() {
}

static bool sampleVelocity( sorter2& sorter, vec2d pos, FLOAT64 sample_visc, FLOAT64 dpx, const std::vector<particle2 *> &neighbors, vec2d &vel ) {
	vec2d out_velocity;
	FLOAT64 wsum = 0.0;
	FOR_EACH_PARTICLE(neighbors) {
		FLOAT64 w = p.n*kernel::sharp_kernel( (pos-p.p).len2(), sample_visc*p.r*dpx);
		out_velocity += w*p.u;
		wsum += w;
	} END_FOR
	if( wsum ) {
		vel = out_velocity/wsum;
	}
	return wsum > 0.0;
}

void corrector2::correct( sorter2 &sorter, int num, const levelset2 *level, std::vector<particle2 *> &particles, FLOAT64 stiffness, FLOAT64 sample_visc,
						  FLOAT64 dt, uint gsize, FLOAT64 dpx, const levelset2* fluid, const levelset2* solid, uint numSample ) {
    // Sort particles
	sorter.sortParticles(particles);
    
	// If stiffness is zero, probably we don't want to correct particle positoins at all...
	if( ! stiffness ) return;

	// Copy particle positions
	std::vector<vec2d> old_positions(particles.size());
	std::vector<vec2d> new_positions(particles.size());
	PARALLEL_FOR for( int n=0; n<particles.size(); n++ ) {
		old_positions[n] = new_positions[n] = particles[n]->p;
		particles[n]->index = n;
	}
	
	std::vector<std::vector<particle2 *> > neighbors_array(particles.size());
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		neighbors_array[pcount] = sorter.getNeighbors(p.p,1.01*sample_visc*dpx*powf(2,level->evalLevelset(p.p)-1));
	} END_FOR
	
	for( uint k=0; k<num; k++ ) {
		// Now push with weak spring forces
		FLOAT64 st = 1.0/sqrt(DIM) * powf(2,DIM);
		PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
			vec2d f;
			
			// Search for nearby particles
			vec2d pos = old_positions[pcount];
			const std::vector<particle2 *> &neighbors = neighbors_array[pcount];
			FOR_EACH_NEIGHBOR_PARTICLE(neighbors) { // np <- neighbors
				if( &p == &np ) continue;
				FLOAT64 r = 0.5*(p.r+np.r)*dpx;
				vec2d rp = old_positions[p.index]-old_positions[np.index];
				FLOAT64 len = rp.len();
				if( len < r && len > 0) {
					vec2d force = r * stiffness * rp / len * kernel::smooth_kernel(len*len,r);
					if( force.len() > r ) force = force.normal() * r;
					f += force;
				}
			} END_FOR
			
			// Don't forget the virtual solid particles
			FLOAT64 phi = solid->evalLevelset(p.p);
			vec2d normal = solid->evalGradient(p.p).normal();
			FLOAT64 len = fabs(phi);
			vec2d rp = normal*len;
			FLOAT64 r = 0.5*dpx*p.r;
			if( len < r && len > 0 ) {
				vec2d force = st * r * stiffness * rp / len * kernel::smooth_kernel(len*len,r);
				if( force.len() > r ) force = force.normal() * r;
				f += force;
			}
			if( p.levelset > -p.r*dpx ) {
				vec2d grad = p.gradient;
				if( grad * f > 0.0 ) f -= grad*(f*grad);
			}
			new_positions[pcount] = old_positions[pcount] + f;
			
			// Move to surface as much as possible
#if 1
			if( k==0 && p.levelset < -p.r*dpx && p.levelset > -3.0*p.r*dpx && solid->evalLevelset(p.p) > p.r*dpx ) {
				FLOAT64 min_neighbor = 1e9;
				FLOAT64 dist_surface = fabs(p.levelset)-0.5*dpx*p.r;
				vec2d avg_gradient;
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
	
	// Copy particle velocities
	std::vector<vec2d> velocities(particles.size());
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
	
	// Don't forget to sort particles !
	if( sample_visc ) sorter.setDirty();
}