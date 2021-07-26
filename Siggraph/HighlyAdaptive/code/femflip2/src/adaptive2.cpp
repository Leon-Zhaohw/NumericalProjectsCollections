/*
 *	adaptive2.cpp
 *	
 *	Created by Ryoichi Ando on 11/9/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "adaptive2.h"
#include "particle2.h"
#include "levelset2.h"
#include "sorter2.h"
#include "kernel.h"
#include "pcgsolver/util.h"
#include "util2.h"
#include "ann2.h"

#include <stdlib.h>
using namespace std;

particle2 *adaptive2::findClosestParticle( const sorter2& sorter, const particle2 &p, FLOAT64 r ) {
    particle2 *nearest = NULL;
    FLOAT64 mind2 = 1.0;
    FLOAT64 r2 = r*r;
	std::vector<particle2 *> particles = sorter.getkNeighbors(p.p,3);
    FOR_EACH_NEIGHBOR_PARTICLE(particles) {
		if( np.rm ) continue;
		FLOAT64 d2 = (np.p-p.p).len2();
		if( &p != &np ) {
			if( d2 < mind2 && d2 < r2 ) {
				nearest = &np;
				mind2 = d2;
			}
		}
    } END_FOR
    return nearest;
}

static void fitParticles( particle2 &p, const levelset2 *levelset, FLOAT64 dpx ) {
	for( uint n=0; n<3; n++ ) {
		FLOAT64 lv = levelset->evalLevelset(p.p);
		FLOAT64 r = 0.5*dpx*p.r;
		FLOAT64 dist = fabs(lv+r);
		if( dist < 1.25*r && ! p.isolated ) {
			vec2d grad = levelset->evalGradient(p.p).normal();
			p.p -= (lv+r)*grad;
		}
	}
	p.levelset = levelset->evalLevelset(p.p);
	p.gradient = levelset->evalGradient(p.p).normal();
}

static bool sampleVelocity( sorter2& sorter, vec2d pos, FLOAT64 dpx, const std::vector<particle2 *> &neighbors, vec2d &vel ) {
	FLOAT64 wsum = 0.0;
	vec2d out_velocity;
	FOR_EACH_PARTICLE(neighbors) {
		FLOAT64 w = p.n*kernel::sharp_kernel( (pos-p.p).len2(), p.r*dpx);
		out_velocity += w*p.u;
		wsum += w;
	} END_FOR
	if( wsum ) {
		vel = out_velocity/wsum;
	}
	return wsum > 0.0;
}

void adaptive2::merge(sorter2& sorter, const levelset2 *level, const levelset2 *surface, std::vector<particle2 *> &particles, FLOAT64 dpx, FLOAT64 dx ) {
    // Do until no particles are left to be merged
	sorter.sortParticles(particles);
	
	for( uint k=0; k<5; k++ ) {
		uint sum_before = 0;
		FOR_EACH_PARTICLE(particles) {
			if( ! p.rm ) sum_before += p.n;
		} END_FOR
		
        // Find unique closest particle pairs
        PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
			p.ref = NULL;
			if( p.isolated || p.levelset > 0.0 ) continue;
			if( ! p.rm ) {
				// This particle lies in the deeper zone, let's search for a nearest particle
				uint depth = level->evalLevelset(p.p);
				if( powf(powf(dx/dpx,depth),DIM) > powf(dx/dpx,DIM)*p.n ) p.ref = findClosestParticle(sorter,p,1.5*dpx*powf(2,depth-1));
			}
        } END_FOR
    
        // Now merge !
		bool merged = false;
		FOR_EACH_PARTICLE(particles) {
			if( p.isolated || p.levelset > 0.0 ) continue;
            // Make sure the pair is referring each other
            if( p.ref ) {
				particle2 &np = *p.ref;
				if( p.ref == &np && np.ref == &p ) {
					if( &p > &np ) {
						vec2d center = (np.n*np.p+p.n*p.p)/(np.n+p.n);
						FLOAT64 depth = level->evalLevelset(center);
						vec2d vel;
						vel = (np.n*np.u+p.n*p.u)/(np.n+p.n);
						FLOAT64 lv = (np.n*np.levelset+p.n*p.levelset)/(np.n+p.n);
						vec2d grad = (np.n*np.gradient+p.n*p.gradient)/(np.n+p.n);
						if( powf(powf(dx/dpx,depth),DIM) > powf(dx/dpx,DIM)*(np.n+p.n) ) {
							if( &p > &np ) {
								p.p = center;
								p.u = vel;
								p.n = np.n+p.n;
								p.gradient = grad;
								p.computeRadius();
								p.levelset = fmin(lv,surface->evalLevelset(p.p));
								p.curvature[0] = 0.5*(np.curvature[0]+p.curvature[0]);
								p.curvature[1] = 0.5*(np.curvature[1]+p.curvature[1]);
								p.remesh_r = 0.5*(np.remesh_r+p.remesh_r);
								np.setRemovable(true);
								np.ref = NULL;
								merged = true;
								fitParticles(p,surface,dpx);
							}
						}
					}
				}
            }
        } END_FOR
		
		uint sum_after = 0;
		FOR_EACH_PARTICLE(particles) {
			if( ! p.rm ) sum_after += p.n;
		} END_FOR
		if( sum_before != sum_after ) {
			email::print("Merge mass not conserved !\n");
			email::send();
			exit(0);
		}
		
		if( ! merged ) break;
    };
	// Cleanup
	if( particle2::cleanParticles(particles)) sorter.setDirty();
}

void adaptive2::split( sorter2& sorter, const levelset2 *level, const levelset2 *surface, std::vector<particle2 *> &particles, FLOAT64 dpx, FLOAT64 dx, bool splitAll ) {
	uint sum_before = 0;
	FOR_EACH_PARTICLE(particles) {
		sum_before += p.n;
	} END_FOR

    // Now perfrorm split
	std::vector<vec2d> thin_array;
	std::vector<FLOAT64> radius_array;
	FOR_EACH_PARTICLE(particles) {
		if( ! p.isolated && ! p.rm && p.levelset > -p.r*dpx ) {
			thin_array.push_back(p.p);
			radius_array.push_back(p.r);
		}
	} END_FOR

	for( uint k=0; k<5; k++ ) {
		std::vector<particle2 *> splitted;
		PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
			p.split = false;
			if( p.isolated ) continue;
			if( ! p.rm && p.n > 1 ) {				
				// If this particle lies in the shallow zone, let's split into the two one-level small particles
				FLOAT64 depth = level->evalLevelset(p.p)-1;
				if( powf(powf(dx/dpx,depth),DIM) < p.n || splitAll ) {
					p.split = true;
				}
			}
		} END_FOR
		
		FOR_EACH_PARTICLE(particles) {
			if( p.split ) {
				// Now split !
				if( fmod(p.n,powf(2,DIM)) == 0.0 ) {
					FOR_EACH(2,2) {
						particle2 *newp = new particle2();
						FLOAT64 r = 0.5*dpx*p.r;
						vec2d pos( p.p[0]+(i==0?-1:1)*0.5*r, p.p[1]+(j==0?-1:1)*0.5*r );
						newp->p = pos;
						newp->u = p.u;
						newp->n = p.n/powf(2,DIM);
						newp->computeRadius();
						splitted.push_back(newp);
						for( uint dim=0; dim<DIM; dim++ ) newp->p[dim] = fmin(1.0,fmax(0.0,newp->p[dim]));
						FLOAT64 surface_lv = surface->evalLevelset(newp->p);
						newp->levelset = surface_lv < 0.0 ? surface_lv : p.levelset;
						newp->gradient = p.gradient;
						newp->curvature[0] = p.curvature[0];
						newp->curvature[1] = p.curvature[1];
						newp->remesh_r = p.remesh_r;
						fitParticles(*newp,surface,dpx);
					} END_FOR
				} else {
					int pn = p.n;
					vec2d normal = level->evalGradient(p.p).normal();
					vec2d direction = vec2d((rand()%201)/100.0-1.0,(rand()%201)/100.0-1.0);
					if( p.levelset > -p.r*dpx ) direction = (direction-(direction*normal)*normal).normal();
					for( uint k=0; k<2; k++ ) {
						particle2 *newp = new particle2();
						// Set random value that ranges -1 ~ 1
						newp->p = p.p+dpx*p.r*direction*(FLOAT64)(k-0.5);
						newp->u = p.u;
						newp->n = k == 0 ? pn/2 : pn;
						newp->computeRadius();
						splitted.push_back(newp);
						pn -= newp->n;
						for( uint dim=0; dim<DIM; dim++ ) newp->p[dim] = fmin(1.0,fmax(0.0,newp->p[dim]));
						FLOAT64 surface_lv = surface->evalLevelset(newp->p);
						newp->levelset = surface_lv < 0.0 ? surface_lv : p.levelset;
						newp->gradient = p.gradient;
						newp->curvature[0] = p.curvature[0];
						newp->curvature[1] = p.curvature[1];
						newp->remesh_r = p.remesh_r;
						fitParticles(*newp,surface,dpx);
					}
				}
				p.setRemovable(true);
			}
		} END_FOR

		if( splitted.size() ) {
			ann2 thin_ann;
			thin_ann.sort(thin_array);
			
			std::vector<vec2d> new_position(splitted.size());
			PARALLEL_FOR FOR_EACH_PARTICLE(splitted) {
				new_position[pcount] = p.p;
				if( p.levelset > -p.r*dpx ) {
					std::vector<ANNidx> neighbors = thin_ann.getkNeighbors(p.p,6);
					vec2d pos;
					FLOAT64 wsum = 0.0;
					for( uint n=0; n<neighbors.size(); n++ ) {
						vec2d op_p = thin_array[neighbors[n]];
						if( (op_p-p.p).len2() > 0 ) {
							vec2d mid = (op_p+p.p)*0.5;
							if( surface->evalLevelset(mid) < 0.0 ) {
								std::vector<ANNidx> neighbors2 = thin_ann.getkNeighbors(mid,1);
								if( neighbors2.size() ) {
									ANNidx min_idx = neighbors2[0];
									FLOAT64 op_r = radius_array[min_idx];
									FLOAT64 w = sqr((thin_array[min_idx]-mid).len()-0.5*dpx*(op_r+p.r));
									wsum += w;
									pos += w*mid;
								}
							}
						}
					}
					if( wsum ) {
						new_position[pcount] = pos / wsum;
					}
				}
			} END_FOR
			PARALLEL_FOR for( uint n=0; n<splitted.size(); n++ ) {
				splitted[n]->p = new_position[n];
				fitParticles(*splitted[n],surface,dpx);
			}
			FOR_EACH_PARTICLE(splitted) {
				if( ! p.isolated && ! p.rm && p.levelset > -p.r*dpx ) {
					thin_array.push_back(p.p);
					radius_array.push_back(p.r);
				}
			} END_FOR
		}
		particles.insert( particles.end(), splitted.begin(), splitted.end() );
		if( splitted.empty()) break;
	};
	
	// Clean removed particles
    particle2::cleanParticles(particles);
	
	uint sum_after = 0;
	FOR_EACH_PARTICLE(particles) {
		sum_after += p.n;
	} END_FOR
	if( sum_before != sum_after ) {
		email::print("Split mass not conserved !\n");
		email::send();
		exit(0);
	}
	
    // Don't forget to sort particles
	sorter.setDirty();
}