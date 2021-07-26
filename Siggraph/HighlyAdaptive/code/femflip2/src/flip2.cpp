/*
 *  flip2.cpp
 *	
 *	Created by Ryoichi Ando on 11/3/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "particle2.h"
#include "flip2.h"
#include "kernel.h"
#include "util2.h"
#include "util.h"
#include "testcase.h"
#include "sizefunc2.h"
#include <vector>
#include <math.h>
#include <stdlib.h>
using namespace std;

// Level strategy
class discreteDepthLevelset2 : public levelset2 {
public:
	discreteDepthLevelset2( const flip2 &sim ) : sim(sim) {
	}
	virtual FLOAT64 evalLevelset(vec2d pos) const {
		for( uint dim=0; dim<DIM; dim++ ) pos[dim] = fmin(1.0,fmax(0.0,pos[dim]));
		const octree2 *octree = sim.fluidSolver->getOctree();
		if( octree ) {
			FLOAT64 depth = log2(octree->terminals[octree->hitTest(pos)]->dx/(FLOAT64)octree->resolution / sim.dx) + 1;
			return depth;
		} else {
			return 1;
		}
	}
	virtual vec2d evalGradient(vec2d p) const {
		return sim.surf.evalGradient(p);
	}
	const flip2 &sim;
};

// Remeshing strategy
class remeshLevelset2 : public levelset2 {
public:
	remeshLevelset2( const flip2 &sim ) : sim(sim) {
	}
	void removePrecomputedSpheres() {
		for( uint n=0; n<spheres.size(); n++ ) {
			delete spheres[n];
		}
		spheres.resize(0);
		std::vector<octree2::sphere2 *>().swap(spheres);
	}
	virtual bool precomputeRemeshAreaBySpheres( const std::vector<particle2 *> &particles, int step=1 ) {
		removePrecomputedSpheres();
		spheres.resize(particles.size());
		PARALLEL_FOR FOR_EACH_PARTICLE(particles) {		
			octree2::sphere2 *sphere = new octree2::sphere2;
			if( step == 0 ) {
				sphere->r = p.remesh_r ? p.remesh_r : evalLevelset(p.p);
			} else {
				sphere->r = evalLevelset(p.p);
				p.remesh_r = sphere->r;
			}
			sphere->p = p.p;
			spheres[pcount] = sphere;
		} END_FOR
		return true;
	}
	virtual const std::vector<octree2::sphere2 *> *getSpheres () const {
		if( spheres.empty()) return NULL;
		return &spheres;
	}
	virtual FLOAT64 evalLevelset(vec2d pos) const {
		for( uint dim=0; dim<DIM; dim++ ) pos[dim] = fmin(1.0-dx,fmax(dx,pos[dim]));
		FLOAT64 levelset = sim.fluid->evalLevelset(pos);
		FLOAT64 curvature = sim.fluid->evalCurvature(pos);
		FLOAT64 diff_curvature = sim.getDiffCurvature(pos);
		FLOAT64 strain = sim.fluidSolver->getStrain(pos);
		if( sim.fluidSolver->getOctree() ) {
			// Compute distance from surface
			levelset = copysign(fmax(fabs(levelset),-sim.solid->evalLevelset(pos)),levelset);
			FLOAT64 dx = fmax(sim.dx,size_func(pos,levelset,curvature,diff_curvature,strain,
										 sim.camera.getProjectedCoord(pos),
										 sim.solid->evalLevelset(pos),
										 sim.solid->evalCurvature(pos),
										 sim.dx,1.0/sim.min_resolution,sim.step==0));
			return dx;
		} else {
			return fabs(levelset) < sim.dx ? -1 : 1;
		}
	}
	const flip2 &sim;
	std::vector<octree2::sphere2 *> spheres;
};

flip2::flip2() {
	surf_order = 2;
	remeshRate = 10;
	doubleRemesh = false;
	numSample = 6;
	totalTime = 0.0;
	variation = true;
	correction = true;
    adaptiveSampling = true;
	doRemesh = false;
	extrapolation_dist = 1.0;
	solid = NULL;
	fluid = NULL;
	fluidSolver = &femSolver;
	sorter = &annsorter;
	surf_sorter = &surf_annsorter;
	remeshLevelset = new remeshLevelset2(*this);
	discreteDepthLevelset = new discreteDepthLevelset2(*this);
	
	// Set up camera
	camera.lookAt(vec2d(0.3,0.5), vec2d(0.85,0.75), vec2d(0.0,1.0), 28.0);
}

void flip2::parameterChanged() {
	fluidSolver->setParameter(fluid2::SURF_ORDER,surf_order);
	fluidSolver->setParameter(fluid2::VARIATION,variation);
}

void flip2::setBoundaryAccuracy( uint surf_order ) {
	this->surf_order = surf_order;
	parameterChanged();
}

void flip2::setVariation( bool variation ) {
	this->variation = variation;
	parameterChanged();
}

void flip2::setCorrection( bool correction ) {
	this->correction = correction;
	parameterChanged();
}

void flip2::setAdaptiveSampling( bool enabled ) {
    this->adaptiveSampling = enabled;
	parameterChanged();
}

void flip2::setRemesh( bool enabled ) {
	doRemesh = enabled;
	parameterChanged();
}

void flip2::setFluidSolver( uint name ) {
	fluid2 *old = fluidSolver;
	if( name == 0 ) {
		fluidSolver = &macSolver;
	} else if( name == 1 ) {
		fluidSolver = &femSolver;
	} else if( name == 2 ) {
		fluidSolver = &fvmSolver;
	} else if( name == 3 ) {
		fluidSolver = &octSolver;
	}
	
	if( old != fluidSolver ) {
		// Make sure the simulation is initialized by taking a look of "solid" pointer
		if( solid ) {
			fluidSolver->init(gn,doRemesh ? remeshLevelset : NULL);
			fluidSolver->setupSolidLevelset(solid);
			
			// Save old levelset
			surf0 = surf;
			// Save old levelset
			if( fluidSolver->getOctree() ) {
				resampleLevelsetOnOctree(fluid,octLevelset,*fluidSolver->getOctree());
			}
			
			// Set levelset
			fluidSolver->setupFluidLevelset(this);
			
			// Set sorter
			sorter->sortParticles(particles);

			// Build surface
			sortSurfaceParticles(*surf_sorter,particles,surf_particles);
			surf.setResolution(gm);
			surf.buildSurface(this,solid,fluidSolver);
			surf_sorter->setDirty();
		}
	}
}

void flip2::init( uint gsize, levelset2* fluid, levelset2* solid, FLOAT64 gscale ) {
	// Set simulation parameters
	srand(110);
	step = 0;
	gn = gsize;
	gm = 2*gn;
	dx = 1.0/gn;
	remesh_dt = 0.0;
	N0 = 4;
	dpx = 0.5*dx;
	anisotropyRad = 5.0;
	CFL = 1.0;
	totalTime = 0.0;
	gravity = 9.8*gscale;
	stiffness = 0.2;
	sample_visc = 1.0;
	flip_visc = 0.0;
	min_stretch = 0.5*sqrt(2.0);
	svd_ratio[0] = 0.2;
	svd_ratio[1] = 0.4;
	grid_tol = 0.0;
	msec = 0.0;
	min_resolution = 16;
	this->solid = solid;
	this->fluid = fluid;
	
	const octree2 *octree = fluidSolver->getOctree();	
    sorter->setDirty();
	solid->resize(gn);
	fluid->resize(gn);
	remeshLevelset->removePrecomputedSpheres();
	
	// Set solver levelset for fluid solver
	fluidSolver->init(gn,doRemesh ? remeshLevelset : NULL);
	fluidSolver->setupSolidLevelset(solid);
	fluidSolver->setupFluidLevelset(fluid);
	parameterChanged();
	
	// Now seed particles
	if( ! adaptiveSampling ) octree = NULL;
	removeAllParticles(particles);
	seedParticle(particles,fluid,solid,octree);
	fitParticles(particles,fluid);
	sorter->sortParticles(particles);
	
	// Compute levelsets on particles
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		p.levelset = fluid->evalLevelset(p.p);
		p.gradient = fluid->evalGradient(p.p);
		p.curvature[0] = fluid->evalCurvature(p.p);
		p.curvature[1] = p.curvature[0];
	} END_FOR
	
	// Build surface
	sortSurfaceParticles(*surf_sorter,particles,surf_particles);
	surf.setResolution(gn);
	surf.buildSurface(fluid,solid,fluidSolver);
	
	// Save old levelset
	if( fluidSolver->getOctree() ) {
		resampleLevelsetOnOctree(fluid,octLevelset,*fluidSolver->getOctree());
	} else {
		surf0 = surf;
	}
	
	// Precompute particle remeshing size
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		p.remesh_r = remeshLevelset->evalLevelset(p.p);
	} END_FOR
	
	// Change fluid levelset
	this->fluid = &surf;
	
	// Clean up
	sorter->sortParticles(particles);
}

void flip2::fitParticles( std::vector<particle2 *> &particles, const levelset2 *levelset ) {
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		FLOAT64 r = 0.5*dpx*p.r;
		FLOAT64 lv = levelset->evalLevelset(p.p);
		if( fabs(lv+r) < r ) {
			for( uint k=0; k<1; k++ ) {
				lv = levelset->evalLevelset(p.p);
				vec2d grad = levelset->evalGradient(p.p).normal();
				p.p += -(r+lv)*grad;
				p.levelset = levelset->evalLevelset(p.p);
				p.gradient = grad;
			}
		}
	} END_FOR
}

void flip2::computeParameters(FLOAT64 tframe, FLOAT64 maxdt, FLOAT64 mindt) {
	// Compute timstep
	maxdt = fmin(0.5*sqrt(dx/gravity),maxdt);
	FLOAT64 max_u = 0.0;
	FOR_EACH_PARTICLE(particles) {
		if( p.isolated ) continue;
		max_u = fmax(max_u,p.u.len());
	} END_FOR
	dt = fmin(maxdt,dx/max_u);
	if( dt < mindt ) {
		dt = mindt;
	}
	CFL = max_u*dt/dx;
	
	// Compute timestep related parameters
	flip_visc = 0.3*dt*gn;
	
	// Set extrapolation degree
	extrapolation_dist = 2.0*fmax(1,CFL)*dx;
	fluidSolver->setParameter(fluid2::EXTRAPOLATE_DIST,extrapolation_dist);
}

class ballLevelset2 : public levelset2 {
public:
	virtual FLOAT64 evalLevelset(vec2d pos) const {
		return (pos-vec2d(0.5,0.4)).len() - 0.025;
	}
};

FLOAT64 flip2::simStep( FLOAT64 tframe, FLOAT64 maxdt, FLOAT64 mindt ) {
	// Start recording simTime
	FLOAT64 msec0 = util2::getMilliseconds();
	
	// Compute timestep and parameters
	computeParameters(tframe,maxdt,mindt);
	
    // Advect particles through grid with 2nd order Runge=Kutta samplnig.
	advectParticles(*sorter,solid,particles,*fluidSolver,dt);
	sorter->sortParticles(particles);
	
	// Save old levelset
	if( fluidSolver->getOctree() ) {
		resampleLevelsetOnOctree(fluid,octLevelset,*fluidSolver->getOctree());
	} else {
		surf0 = surf;
	}
	
	// Collect splash particles into the sorter
	collectSplashParticles(bullet_particles,clustered_particles,ann,*sorter,particles,fluid);
	
	// Resample splashing particle velocity
	if( bullet_particles.size()) {
		std::vector<vec2d> new_velocity(bullet_particles.size());
		std::vector<bool> resampled(bullet_particles.size());
		PARALLEL_FOR for( uint n=0; n<bullet_particles.size(); n++ ) {
			FLOAT64 r = bullet_particles[n]->r;
			new_velocity[n] = bullet_particles[n]->u;
			resampled[n] = sampleVelocity(*sorter,bullet_particles[n]->p,1.5*r*sample_visc,dpx,new_velocity[n]);
		}
		for( uint n=0; n<bullet_particles.size(); n++ ) {
			if( resampled[n] ) bullet_particles[n]->u = new_velocity[n];
		}
	}
	
	// Compute new surface
	sorter->sortParticles(particles);
	sortSurfaceParticles(*surf_sorter,particles,surf_particles);
	surf.buildSurface(this,solid,fluidSolver);
	surf_sorter->setDirty();
	
	// Compute levelsets on particles
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		p.levelset = fluid->evalLevelset(p.p);
		p.gradient = fluid->evalGradient(p.p);
	} END_FOR
    
	if( doRemesh && (step+1) % remeshRate == 0 ) {
		if( doubleRemesh ) {
			// Prepare remesh with spheres
			remeshLevelset->precomputeRemeshAreaBySpheres(particles,0);
			
			// Remesh
			fluidSolver->init(gn,remeshLevelset);
			fluidSolver->setupSolidLevelset(solid);
			
			sortSurfaceParticles(*surf_sorter,particles,surf_particles);
			surf.buildSurface(this,solid,fluidSolver);
			surf_sorter->setDirty();
			
			// Prepare remesh with spheres
			remeshLevelset->precomputeRemeshAreaBySpheres(particles,1);
		} else {
			// Prepare remesh with spheres
			remeshLevelset->precomputeRemeshAreaBySpheres(particles);
		}
		// Remesh
		fluidSolver->init(gn,remeshLevelset);
		fluidSolver->setupSolidLevelset(solid);
		
		// Merge particles for adaptive sampling
		if( adaptiveSampling ) adaptive.merge(*sorter,discreteDepthLevelset,fluid,particles,dpx,dx);
		
		// Sort particles
		sorter->sortParticles(particles);
		
		// Split particles for adaptive sampling
		adaptive.split(*sorter,discreteDepthLevelset,fluid,particles,dpx,dx,adaptiveSampling==0);
		
		remesh_dt = dt;
		FOR_EACH_PARTICLE(particles) {
			p.curvature[1] = p.curvature[0];
			p.curvature[0] = fluid->evalCurvature(p.p);
		} END_FOR
		
		// Remove prepared remesh spheres
		remeshLevelset->removePrecomputedSpheres();
	} else {
		// Accumulate curvature on particles
		PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
			p.curvature[0] += fluid->evalCurvature(p.p);
		} END_FOR
		remesh_dt += dt;
	}
	
	// Correct particle position
	if( correction ) corrector.correct(*sorter,floor(fmax(0.1,fmin(3,CFL))+1.0-1e-8),discreteDepthLevelset,particles,
									    stiffness,sample_visc,dt,gn,dpx,fluid,solid,numSample);
	
	// Add force
	addForce(*sorter,particles,extForcePos,extForces,gravity,dt);
	
	// Solve FLIP fluid
	projectFLIP(*fluidSolver,*sorter,particles,solid,fluid,flip_visc,sample_visc,dpx,doRemesh,dt);
	
	// Finish recording the simtime for this step
	msec = util2::getMilliseconds()-msec0;
	
	// Increment simulation step
	step++;
	
	// Accumulate elappsed time
	totalTime += dt;
	
	// Return a total time
	return totalTime;
}

void flip2::projectFLIP(fluid2 &fluidSolver, sorter2& sorter, std::vector<particle2 *> &particles,
						const levelset2* solid, const levelset2* fluid, FLOAT64 flip_visc, FLOAT64 sample_visc, FLOAT64 dpx, bool doRemesh, FLOAT64 dt ) {
	// Sort particles
	sorter.sortParticles(particles);
	
	// Collect grid positions
	std::vector<vec2d> pos;
	std::vector<vec2d> vel;
	std::vector<bool> mapped;
	fluidSolver.getSamplePoints(pos);
	vel.resize(pos.size());
	mapped.resize(pos.size());
	
	// Set velocity
	PARALLEL_FOR for( uint n=0; n<pos.size(); n++ ) {
		FLOAT64 r = dx*powf(2,discreteDepthLevelset->evalLevelset(pos[n])-1);
		if( fluid->evalLevelset(pos[n]) <= r ) {
			vel[n] = vec2d();
			mapped[n] = sampleVelocity(sorter,pos[n],sample_visc,dpx,vel[n]);
		} else {
			mapped[n] = false;
			vel[n] = vec2d();
		}
	}
	fluidSolver.setupVelocity(vel,mapped);
	fluidSolver.setupFluidLevelset(fluid);
	
	// Now project on grid
	fluidSolver.setTimestep(dt);
	fluidSolver.project();
	
	// Compute particle new FLIP velocity
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		if( p.levelset < 0.0 ) {
			vec2d pos = p.p;
			vec2d gu = fluidSolver.getVelocity(pos);
			vec2d pu = p.u;
			vec2d gradp = fluidSolver.getPressureGradient(pos);
			FLOAT64 visc = fmax(0.01,fmin(1.0,flip_visc/sqr(p.r)));
			p.u = visc*gu+(1.0-visc)*pu-dt*gradp;
		}
	} END_FOR
}

void flip2::addForce(sorter2& sorter, std::vector<particle2 *> &particles,
					 std::vector<vec2d> &extForcePos, std::vector<vec2d> &extForces, FLOAT64 gravity, FLOAT64 dt ) {
	// Sort particles
	sorter.sortParticles(particles);
	
	// Add gravity force
	FOR_EACH_PARTICLE(particles) {
		p.u[1] -= dt*gravity;
	} END_FOR
	
	// Apply user-defined external forces...
	for( uint m=0; m<extForcePos.size(); m++ ) {
		vec2d position = extForcePos[m];
		vec2d f = extForces[m];
		FLOAT64 r = dx;
		FOR_EACH_PARTICLE(particles) {
			if( (p.p-position).len2() < r*r ) p.u += dt*f;
		} END_FOR
	}
	extForcePos.clear();
	extForces.clear();
}

bool flip2::getClosestSurfacePos(vec2d &pos) const {
	FLOAT64 min_phi = 1e9;
	std::vector<particle2 *> neighbors = surf_sorter->getkNeighbors(pos,numSample/2);
	FLOAT64 r = 0.0;
	FLOAT64 rsum = 0.0;
	for( uint n=0; n<neighbors.size(); n++ ) {
		r += neighbors[n]->r;
		rsum += 1.0;
	}
	if( rsum ) r = r / rsum;
	if( solid->evalLevelset(pos) <= r*dpx) return false;
	vec2d p = pos;
	min_phi = fmin(min_phi,levelset2::getParticleLevelset(neighbors,dpx,p,pos,solid));
	return fabs(min_phi) < 0.5*dpx*r;
}

void flip2::sortSurfaceParticles( sorter2 &sorter, const std::vector<particle2 *> &particles, std::vector<particle2 *> &surf_particles ) {
	surf_particles.clear();
	FOR_EACH_PARTICLE(particles) {
		if( p.levelset > -2.0*p.r*dpx ) {
			surf_particles.push_back(&p);
		}
	} END_FOR
	sorter.sortParticles(surf_particles);
}

FLOAT64 flip2::evalLevelset(vec2d pos) const {
	// Compute old levelset
	FLOAT64 oldLevelset;
	if( fluidSolver->getOctree() ) oldLevelset = octLevelset.evalLevelset(pos);
	else oldLevelset = surf0.evalLevelset(pos);
		
	// Compute base grid size
	FLOAT64 base_dx = 1.0;
	const octree2 *octree = fluidSolver->getOctree();
	if( octree ) {
		std::vector<uint> array;
		if( octree->hitTest( pos, 0.5*dx, array ) ) {
			for( uint n=0; n<array.size(); n++ ) {
				int hit = array[n];
				base_dx = fmin(base_dx,octree->terminals[hit]->dx/(FLOAT64)octree->resolution);
			}
		} else {
			dump( "flip2::evalLevelset base grid size failed to fetch at (%f,%f)\n", pos[0], pos[1] );
		}
	} else {
		base_dx = dx;
	}
	
	// Compute solid levelset there
	FLOAT64 solid_phi = solid->evalLevelset(pos);
	
	// Shift position
	int max_atempt = 5;
	FLOAT64 min_dist = 0.01*base_dx;
	while( solid_phi < min_dist && max_atempt >= 0 ) {
		pos += (min_dist-solid_phi)*solid->evalGradient(pos);
		solid_phi = solid->evalLevelset(pos);
		max_atempt --;
	}
	
	// First let's see the closest particle position
	FLOAT64 LCFL = 1.0;
	FLOAT64 min_phi = 1e9;
	FLOAT64 shrink = 0.75;
	std::vector<particle2 *> neighbors = sorter->getkNeighbors(pos,numSample);

	FOR_EACH_PARTICLE(neighbors) {
		vec2d out = pos;
		min_phi = fmin(min_phi,levelset2::getParticleConvextHullLevelset(out,p,dpx,shrink,false));
		LCFL = fmax(LCFL,p.u.len()*dt/base_dx);
	} END_FOR
	
	if( fabs(oldLevelset) > fmax(1.0,LCFL)*base_dx ) return fmin(min_phi,oldLevelset);
	
	// Evaluate exact distance
	if( fabs(min_phi) < fmax(1.0,LCFL)*base_dx ) {
		vec2d outpos;
		neighbors = surf_sorter->getkNeighbors(pos,numSample);
		min_phi = fmin(min_phi,levelset2::getParticleLevelset(neighbors,dpx,pos,outpos,solid,shrink));
	}
	
	// Compute splash particle levelset
	if( min_phi > 0 ) {
		min_phi = fmin(min_phi,ann.evalAbsLevelset(pos,dpx));
	}
	return min_phi;
}

void flip2::advectParticles( sorter2& sorter, const levelset2* solid, 
							 std::vector<particle2 *> &particles, const fluid2 &fluidSolver, FLOAT64 dt ) {
	// Advect using 2nd order RK sampling
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		FLOAT64 r = dx*powf(2,discreteDepthLevelset->evalLevelset(p.p)-1);
		if( ! p.isolated && p.levelset <= 0.0 && dpx*p.r > grid_tol*r ) {
			vec2d v1 = fluidSolver.getVelocity(p.p)-dt*fluidSolver.getPressureGradient(p.p);
			vec2d p2 = p.p-dt*v1;
			for( uint dim=0; dim<DIM; dim++ ) p2[dim] = fmin(1.0,fmax(0.0,p2[dim]));
			vec2d v2 = fluidSolver.getVelocity(p2)-dt*fluidSolver.getPressureGradient(p2);
			p.p += 0.5*dt*(v1+v2);
		} else {
			p.p += dt*p.u;
		}
	} END_FOR

	// Push out particles out of objects, easy because we know the levelset value
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		FLOAT64 phi;
		if( (phi = solid->evalLevelset(p.p)-0.5*p.r*dpx) < 0.0 ) {
			vec2d normal = solid->evalGradient(p.p).normal();
			p.p += -1.01*phi*normal;
			p.u -= (normal*p.u)*normal;
		}
		for( uint dim=0; dim<DIM; dim++ ) p.p[dim] = fmin(1.0,fmax(0.0,p.p[dim]));
	} END_FOR
	
	// Set sorter as dirty
	sorter.setDirty();
}

void flip2::removeAllParticles( std::vector<particle2 *> &particles ) {
	for( uint n=0; n<particles.size(); n++ ) {
		particles[n]->setRemovable(true);
	}
	particle2::cleanParticles(particles);
	particles.clear();
}

void flip2::seedParticle(std::vector<particle2 *> &particles, const levelset2* fluid, const levelset2* solid, const octree2 *octree ) {
	FLOAT64 e = 1e-4*dpx;
	if( ! octree ) {
		FOR_EACH(gn,gn) {
			for( FLOAT64 x = 0.5*dpx; x < dx; x += dpx ) {
				for( FLOAT64 y = 0.5*dpx; y < dx; y += dpx ) {
					// Place particles if levelset is negative there
					vec2d p( dx*i+x+e*nrand(), dx*j+y+e*nrand() );
					if( solid->evalLevelset(p) > 0 && fluid->evalLevelset(p) < 0.0 ) {
						particle2 *newp = new particle2;
						newp->p = p;
						newp->u = vec2d();
						newp->levelset = fluid->evalLevelset(p);
						particles.push_back(newp);
					}
				}
			}
		} END_FOR
	} else {
		// Seed particles per terminal cell
		for( uint n=0; n<octree->terminals.size(); n++ ) {
			const octree2::leaf2 *leaf = octree->terminals[n];
			FLOAT64 dx = leaf->dx/(FLOAT64)octree->resolution;
			vec2d corner = vec2d(leaf->center)/octree->resolution-0.5*dx*vec2d(1,1);
			int r = dx > dpx ? 2 : 1;
			FOR_EACH(r,r) {
				FLOAT64 padding = 0.5*dx/r;
				vec2d p(corner[0]+padding+i*dx/r+e*nrand(),corner[1]+padding+j*dx/r+e*nrand());
				if( solid->evalLevelset(p) > 0.0 && fluid->evalLevelset(p) < 0.0 ) {
					particle2 *newp = new particle2;
					newp->p = p;
					newp->u = vec2d();
					newp->n = fmax(1.0,powf(0.5*dx/dpx,DIM));
					newp->computeRadius();
					particles.push_back(newp);
				}
			} END_FOR
		}
	}
}

bool flip2::sampleVelocity( sorter2& sorter, vec2d pos, FLOAT64 sample_visc, FLOAT64 dpx, vec2d &vel ) {
	FLOAT64 scale = fluidSolver->getMesh() == NULL ? 2.0 : 1.0;
	FLOAT64 base_r = dpx*powf(2,discreteDepthLevelset->evalLevelset(pos)-1);
	FLOAT64 r = scale*sample_visc*base_r;
	std::vector<particle2 *> neighbors = sorter.getNeighbors(pos,r);
	FLOAT64 wsum = 0.0;
	vec2d out_velocity;
	FOR_EACH_PARTICLE(neighbors) {
		if( ! p.isolated && dpx*p.r > grid_tol*base_r ) {
			FLOAT64 w = p.n*kernel::sharp_kernel( (pos-p.p).len2(), sample_visc*p.r*dpx );
			out_velocity += w*p.u;
			wsum += w;
		}
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

void flip2::resampleLevelsetOnOctree(const levelset2 *levelset ,octlevelset2 &octLevelset, const octree2 &octree) {
	octLevelset.setOctree(octree);
	std::vector<FLOAT64> q(octree.nodes.size());
	PARALLEL_FOR for( uint n=0; n<octree.nodes.size(); n++ ) {
		q[n] = levelset->evalLevelset(octree.nodes[n]);
	}
	octLevelset.setNodalLevelset(q);
}

void flip2::collectSplashParticles( std::vector<particle2 *> &bullet_particles, std::vector<particle2 *> &clustered_particles,
								    pann2 &ann, const sorter2 &sorter, const std::vector<particle2 *> &particles, const levelset2 *fluid ) {
	// Collect isolated particles
	bullet_particles.clear();
	clustered_particles.clear();
	
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		FLOAT64 levelset = fluid->evalLevelset(p.p);
		FLOAT64 r = dpx*p.r;
		if( levelset < 0.0 ) r = fmax(r,dpx*powf(2,discreteDepthLevelset->evalLevelset(p.p)-1));
		std::vector<particle2 *> neighbors = sorter.getkNeighbors(p.p,2);
		p.isolated = true;
		for( uint n=0; n<neighbors.size(); n++ ) {
			if( neighbors[n] != &p && (p.p-neighbors[n]->p).len2() < sqr(3.0*r) ) {
				p.isolated = false;
				break;
			}
		}
	} END_FOR
	
	FOR_EACH_PARTICLE(particles) {
		if( p.isolated ) bullet_particles.push_back(&p);
	} END_FOR
	
	// Collect cluster particles
	FOR_EACH_PARTICLE(particles) {
		if( p.levelset > 0.5*p.r*dpx && ! p.isolated ) {
			bullet_particles.push_back(&p);
			clustered_particles.push_back(&p);
		}
	} END_FOR
	
	// Put into ANN
	ann.buildKDTree(clustered_particles);
}

FLOAT64 flip2::getDiffCurvature(vec2d pos) const {
	if( surf_particles.empty() ) return 0.0;
	std::vector<particle2 *> neighbors = surf_sorter->getkNeighbors(pos,numSample);
	FLOAT64 curvature[2] = { 0.0, 0.0 };
	FLOAT64 wsum = 0.0;
	FLOAT64 max_r = 0.0;
	FOR_EACH_PARTICLE(neighbors) {
		max_r = fmax(max_r,(p.p-pos).len());
	} END_FOR
	FOR_EACH_PARTICLE(neighbors) {
		FLOAT64 w = p.n*kernel::smooth_kernel((p.p-pos).len2(),4*max_r);
		curvature[0] = w*p.curvature[0];
		curvature[1] = w*p.curvature[1];
		wsum += w;
	} END_FOR
	if( wsum ) {
		curvature[0] /= wsum;
		curvature[1] /= wsum;
		FLOAT64 diff = curvature[0]-curvature[1];
		if( diff ) {
			return  diff / fmax(remesh_dt,1e-18);
		}
	}
	return 0.0;
}

void flip2::addExternalForce( vec2d p, vec2d f ) {
	extForcePos.push_back(p);
	extForces.push_back(f/dt);
}