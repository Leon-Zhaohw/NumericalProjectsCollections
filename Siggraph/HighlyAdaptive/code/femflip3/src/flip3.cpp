/*
 *	flip3.cpp
 *	
 *	Created by Ryoichi Ando on 1/9/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "particle3.h"
#include "flip3.h"
#include "kernel.h"
#include "util3.h"
#include "util.h"
#include "matutil.h"
#include "sizefunc3.h"
#include "livemem.h"
#include "exporter3.h"
#include <vector>
#include <math.h>
#include <stdlib.h>
using namespace std;

// Level strategy
class discreteDepthLevelset3 : public levelset3 {
public:
	discreteDepthLevelset3( const flip3 &sim ) : sim(sim) {
	}
	virtual FLOAT64 evalLevelset(vec3d pos) const {
		for( uint dim=0; dim<DIM; dim++ ) pos[dim] = fmin(1.0,fmax(0.0,pos[dim]));
		const octree3 *octree = sim.fluidSolver->getOctree();
		if( octree ) {
			FLOAT64 depth = log2(octree->terminals[octree->hitTest(pos)]->dx/(FLOAT64)octree->resolution / sim.dx) + 1;
			return depth;
		} else {
			return 1;
		}
	}
	virtual vec3d evalGradient(vec3d p) const {
		return sim.surf.evalGradient(p);
	}
	const flip3 &sim;
};

// Remeshing strategy
class remeshLevelset3 : public levelset3 {
public:
	remeshLevelset3( const flip3 &sim ) : sim(sim) {
	}
	void removePrecomputedSpheres() {
		for( uint n=0; n<spheres.size(); n++ ) {
			delete spheres[n];
		}
		spheres.resize(0);
		std::vector<octree3::sphere3 *>().swap(spheres);
	}
	virtual bool precomputeRemeshAreaBySpheres( const std::vector<particle3 *> &particles, int step=1 ) {
		removePrecomputedSpheres();
		spheres.resize(particles.size());
		PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
			octree3::sphere3 *sphere = new octree3::sphere3;
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
	virtual const std::vector<octree3::sphere3 *> *getSpheres () const {
		if( spheres.empty()) return NULL;
		return &spheres;
	}
	virtual FLOAT64 evalLevelset(vec3d pos) const {
		for( uint dim=0; dim<DIM; dim++ ) pos[dim] = fmin(1.0,fmax(0.0,pos[dim]));
		FLOAT64 levelset = sim.fluid->evalLevelset(pos);
		FLOAT64 curvature = sim.fluid->evalCurvature(pos);
		FLOAT64 diff_curvature = sim.getDiffCurvature(pos);
		FLOAT64 strain = sim.fluidSolver->getStrain(pos);
		bool visible = true;
		for( uint n=0; n<sim.objects.size(); n++ ) {
			const object3 &object = *sim.objects[n];
			if( object.type == object3::SOLID ) {
				vec3d camera_pos = sim.camera->origin;
				if( object.intersect(pos, camera_pos)) {
					visible = false;
					break;
				}
			}
		}
		if( sim.fluidSolver->getOctree() ) {
			// Compute distance from surface
			levelset = copysign(fmax(fabs(levelset),-sim.solid->evalLevelset(pos)),levelset);
			FLOAT64 dx = fmax(sim.dx,size_func(pos,levelset,curvature,diff_curvature,strain,
											   sim.camera->origin, sim.camera->target,
											   sim.camera->getProjectedCoord(pos),visible,
											   sim.solid->evalLevelset(pos),
											   sim.solid->evalCurvature(pos),
											   sim.dx,1.0/sim.min_resolution,sim.step==0));
			return dx;
		} else {
			return 0;
		}
	}
	const flip3 &sim;
	std::vector<octree3::sphere3 *> spheres;
};

flip3::flip3() {
	solid = NULL;
	fluid = NULL;
	remeshRate = 10;
	doubleRemesh = false;
	correction = true;
	adaptiveSampling = true;
	doRemesh = true;
	extrapolation_dist = 1.0;
	grid_tol = 0.0;
	step = 0;
	msec = 0.0;
	elapsedTime = 0.0;
	frameNumber = 0;
	avg_msec = 0.0;
	totalTime = 0.0;
	fluidSolver = &femSolver;
	remeshLevelset = new remeshLevelset3(*this);
	discreteDepthLevelset = new discreteDepthLevelset3(*this);
}

void flip3::setElapsedTime( FLOAT64 elapsedTime ) {
	this->elapsedTime = elapsedTime;
}

void flip3::setVideoFrameNumber( uint frame ) {
	this->frameNumber = frame;
}

void flip3::init( uint gsize, std::vector<object3 *> objects, FLOAT64 gravity, const camera3 *camera ) {
	tick(); dump( "--------- Initializing... --------\n" );
	setTimestepNumber(0);
	
	// Set simulation parameters
	gn = gsize;
	gm = gn;
	numSample = 18;
	dx = 1.0/gn;
	min_resolution = 16;
	N0 = 8;
	dpx = 0.5*dx;
	anisotropyRad = 5.0;
	CFL = 1.0;
	totalTime = 0.0;
	sitffness = 0.2;
	dt = 0.16/gn;
	remesh_dt = 0.0;
	sample_visc = 1.0;
	flip_visc = 0.0;
	min_stretch = 0.5*sqrt(DIM);
	svd_ratio[0] = 0.2;
	svd_ratio[1] = 0.4;
	this->solid = new objLevelset3(object3::SOLID,objects);
	this->fluid = new objLevelset3(object3::FLUID,objects);
	this->camera = camera;
	this->gravity = gravity;
	this->tension = tension;
	
	// Set sorter
	const octree3 *octree = fluidSolver->getOctree();
	sorter = &annsorter;
	surf_sorter = &surf_annsorter;
	
	sorter->setDirty();
	solid->resize(gn);
	fluid->resize(gn);
	remeshLevelset->removePrecomputedSpheres();
	
	// Dump out setup configuration
	dump( ">>> Configuration:\n" );
	dump( "Adaptive maximal resolution: %dx%dx%d\n", gn, gn, gn );
	dump( "Adaptive minimal resolution: %dx%dx%d (surface)\n", min_resolution, min_resolution, min_resolution );
	dump( "Grid solver name: %s\n", fluidSolver->getName() );
	dump( "Positive-definite matrix solver: %s\n", "ICCG" );
	dump( "Adaptive sampling enabled: %s\n", adaptiveSampling ? "YES" : "NO" );
	dump( "Position correction enabled: %s\n", correction ? "YES" : "NO" );
	dump( "Double remeshing enabled: %s\n", doubleRemesh ? "YES" : "NO" );
	dump( "<<<\n" );
	
	// Allocate domain
	this->objects = objects;
    sorter->setDirty();
	solid->resize(gn);
	fluid->resize(gn);
	
	// Set solver levelset for fluid solver
	sorter->sortParticles(particles);
	fluidSolver->init(gn,doRemesh ? remeshLevelset : NULL);
	fluidSolver->setupSolidLevelset(solid);
	fluidSolver->setupFluidLevelset(fluid);
	
	// Now seed particles
	tick(); dump("Seeding particles...");
	if( ! adaptiveSampling ) octree = NULL;
	removeAllParticles(particles);
	seedParticle(particles,fluid,solid,octree);
	fitParticles(particles,fluid);
	sorter->sortParticles(particles);
	dump("Done. Seeded %d particles. Took %s.\n", particles.size(), stock("main_seeding"));
	
	tick(); dump("Computing levelset at particles...");
	// Compute levelsets on particles
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		p.levelset = fluid->evalLevelset(p.p);
		p.gradient = fluid->evalGradient(p.p);
	} END_FOR
	dump("Done. Took %s.\n", stock());
	
	// Build surface
	sortSurfaceParticles(*surf_sorter,particles,surf_particles);
	surf.setResolution(gm);
	sorter->sortParticles(particles);
	surf.buildSurface(fluidSolver,fluid,solid);
	surf_sorter->setDirty();
	
	// Save old levelset
	if( fluidSolver->getOctree() ) {
		resampleLevelsetOnOctree(fluid,octLevelset,*fluidSolver->getOctree());
	} else {
		surf0 = surf;
	}
	
	// Precompute particle remeshing size
	tick(); dump("Precomputing remesh radius at particles...");
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		p.remesh_r = remeshLevelset->evalLevelset(p.p);
	} END_FOR
	dump("Done. Took %s.\n", stock());
	
	// Change fluid levelset
	this->fluid = &surf;

	// Clean up
	sorter->sortParticles(particles);
	dump( "Initialization done. Took %s.\n", stock("time") );
}

void flip3::fitParticles( std::vector<particle3 *> &particles, const levelset3 *levelset ) {
	FOR_EACH_PARTICLE(particles) {
		FLOAT64 r = 0.5*dpx*p.r;
		FLOAT64 lv = levelset->evalLevelset(p.p);
		if( fabs(lv+r) < r ) {
			for( uint k=0; k<1; k++ ) {
				lv = levelset->evalLevelset(p.p);
				vec3d grad = levelset->evalGradient(p.p).normal();
				p.p += -(r+lv)*grad;
			}
		}
	} END_FOR
}

void flip3::computeParameters( FLOAT64 tframe, FLOAT64 maxdt, FLOAT64 mindt) {
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
	CFL = max_u*dt/dpx;
	
	// Compute timestep dependent parameters
	flip_visc = 0.3*dt*gn;
	
	// Set extrapolation degree
	extrapolation_dist = 2.0*fmax(1,CFL)*dx;
	fluidSolver->setParameter(fluid3::EXTRAPOLATE_DIST,extrapolation_dist);
	
	dump( ">>> Timestep dependent parameters:\n" );
	dump( "dt = %f (approx %.2f step/frame, %.2f steps until next frame)\n", dt, tframe/dt, maxdt/dt );
	dump( "time/frame = %f, maxdt = %f, mindt = %f, max_u = %f\n", tframe, maxdt, mindt, max_u );
	dump( "CFL = %f\n", CFL );
	dump( "flip visc = %f\n", flip_visc );
	dump( "<<<\n" );
	
	// Make sure the CFL number is acceptible
	if( CFL > 100 ) {
		uint extreme_num = 0;
		FLOAT64 min_levelset = 1e8;
		FLOAT64 max_height = 0.0;
		FOR_EACH_PARTICLE(particles) {
			if( p.isolated ) continue;
			if( p.u.len()*dt/dpx > 100 ) {
				extreme_num ++;
				min_levelset = fmin(p.levelset,min_levelset);
				max_height = fmax(p.p[1],max_height);
			}
		} END_FOR
		dump( "CFL number is too large ! something is wrong. Number of extremely fast moving particle is %d. min(levelset) = %e. max(height)=%f\n", extreme_num, min_levelset, max_height );
	}
}

void flip3::writeStatus() const {
	uint num_isolated = 0;
	FOR_EACH_PARTICLE(particles) {
		if( p.isolated ) num_isolated++;
	} END_FOR
	writeNumber("status_particle",particles.size());
	writeNumber("status_bullet", bullet_particles.size());
	writeNumber("status_clustered", clustered_particles.size());
	writeNumber("status_isolated", num_isolated);
	writeNumber("status_CFL", CFL);
	writeNumber("status_dt", dt);
	
	// Get system memory
	double vm, rss;
	double gbscale = 1024*1024;
	process_mem_usage(vm,rss);
	writeNumber("status_memory_vm", vm/gbscale);
	writeNumber("status_memory_rss", rss/gbscale);
	writeNumber("status_memory_sum", (vm+rss)/gbscale);
}

FLOAT64 flip3::simStep( FLOAT64 tframe, FLOAT64 maxdt, FLOAT64 mindt ) {
	// Start recording simTime
	writeStatus();
	setTimestepNumber(step+1);
	tick(); dump( "------------ %s step -------------\n", nth(step+1) );
	
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
	
	// Resample splashing particle velocit
	tick(); dump("Resampling %d bullet particles velocities...", bullet_particles.size());
	if( bullet_particles.size()) {
		std::vector<vec3d> new_velocity(bullet_particles.size());
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
	dump("Done. Took %s.\n", stock("main_resample_bullet_particle_velocity"));
	
	// Generate a mesh
	sortSurfaceParticles(*surf_sorter,particles,surf_particles);
	sorter->sortParticles(particles);
	surf.buildSurface(fluidSolver,this,solid);
	surf_sorter->setDirty();
	
	// Compute levelsets on particles
	tick(); dump("Computing levelset at particles...");
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		p.levelset = fluid->evalLevelset(p.p);
		p.gradient = fluid->evalGradient(p.p);
	} END_FOR
	dump("Done. Took %s.\n", stock("main_computing_levelset_particle"));
	
	uint num_particles = particles.size();
	if( doRemesh && (step+1) % remeshRate == 0 ) {
		tick(); dump( ">>> Remeshing started... (Remeshrate = %d Doubled = %s)\n", remeshRate, doubleRemesh ? "YES" : "NO" );
		if( doubleRemesh ) {
			// Prepare remesh with spheres
			tick(); dump("Precomputing remeshing sizing functions at particles (1st pass)...");
			remeshLevelset->precomputeRemeshAreaBySpheres(particles,0);
			dump("Done. Took %s.\n", stock("main_computing_remesh_prepare_1st"));
			
			// Remesh
			fluidSolver->init(gn,remeshLevelset);
			fluidSolver->setupSolidLevelset(solid);
			
			// Prepare remesh with spheres again			
			sortSurfaceParticles(*surf_sorter,particles,surf_particles);
			surf.buildSurface(fluidSolver,this,solid);
			surf_sorter->setDirty();
			
			// Prepare remesh with spheres
			tick(); dump("Precomputing remeshing sizing functions at particles (2nd pass)...");
			remeshLevelset->precomputeRemeshAreaBySpheres(particles,1);
			dump("Done. Took %s.\n", stock("main_computing_remesh_prepare_2nd"));
		} else {
			// Prepare remesh with spheres
			tick(); dump("Precomputing remeshing sizing functions at particles...");
			remeshLevelset->precomputeRemeshAreaBySpheres(particles);
			dump("Done. Took %s.\n", stock("main_computing_remesh_prepare"));
		}
		fluidSolver->init(gn,remeshLevelset);
		fluidSolver->setupSolidLevelset(solid);
		sorter->sortParticles(particles);
		dump( "<<< Done. Took %s.\n", stock("main_remesh"));
		
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
	if( correction ) corrector.correct(*sorter,floor(fmax(0.1,fmin(6,CFL))+1.0-1e-8),discreteDepthLevelset,particles,
									    sitffness,sample_visc,dt,gn,dpx,fluid,solid,numSample);
    
	// Add force
	addForce(particles,gravity,dt);
	
	// Solve FLIP fluid
	projectFLIP(*fluidSolver,*sorter,particles,solid,fluid,flip_visc,sample_visc,dpx,doRemesh,dt);
	
	// Finish recording the simtime for this step
	msec = tock();
	
	avg_msec += msec;
	FLOAT64 time_per_particle = 0.0;
	if( (step+1) % remeshRate == 0 ) {
		writeNumber("avg_time",avg_msec/(FLOAT64)remeshRate);
		time_per_particle = avg_msec/(FLOAT64)remeshRate/num_particles;
		writeNumber("avg_time_per_particle", time_per_particle);
		avg_msec = 0.0;
	}
	
	dump( "%s step done. Took %s.", nth(step+1), tstr(msec) );
	if( elapsedTime ) dump( " Elapsed time: %s.", tstr(elapsedTime+msec) );
	if( frameNumber ) dump( " Video frame: %d.", frameNumber );
	dump( "\n" );
	if( time_per_particle ) dump("Time per particle = %s\n", tstr(time_per_particle) );
		
	writeNumber("time",msec);
	if( elapsedTime ) writeNumber("time_total",elapsedTime+msec);
	
	// Advance simulation step
	step++;
	
	// Accumulate total time
	totalTime += dt;
	
	// Return a total time
	return totalTime;
}

void flip3::projectFLIP( fluid3 &fluidSolver, sorter3& sorter, std::vector<particle3 *> &particles,
						 const levelset3* solid, const levelset3* fluid, FLOAT64 flip_visc, FLOAT64 sample_visc, FLOAT64 dpx, bool doRemesh, FLOAT64 dt ) {
	// Sort particles
	sorter.sortParticles(particles);
	
	tick(); dump( ">>> FLIP projection started...\n" );
	tick(); dump( "Mapping velocity..." );
	
	// Collect grid positions
	vector<vec3d> pos;
	vector<vec3d> vel;
	vector<bool> mapped;
	fluidSolver.getSamplePoints(pos);
	vel.resize(pos.size());
	mapped.resize(pos.size());
	
    // Set velocity
	FLOAT64 max_u = 0.0;
	PARALLEL_FOR for( uint n=0; n<pos.size(); n++ ) {
		FLOAT64 r = dx*powf(2,discreteDepthLevelset->evalLevelset(pos[n])-1);
		if( fluid->evalLevelset(pos[n]) <= r ) {
			vel[n] = vec3d();
			mapped[n] = sampleVelocity(sorter,pos[n],sample_visc,dpx,vel[n]);
			if( mapped[n] ) max_u = fmax(max_u,vel[n].len());
		} else {
			mapped[n] = false;
			vel[n] = vec3d();
		}
	}
    fluidSolver.setupVelocity(vel,mapped);
    
	dump( "Done. Took %s.\n", stock("main_velocity_map"));
	tick(); dump( ">>> Grid projection (%s).\n", fluidSolver.getName());
	
	// Set variables
	fluidSolver.setTimestep(dt);
	fluidSolver.setupFluidLevelset(fluid);
	
	// Now project on grid
	fluidSolver.project();
	
	dump( "<<< Grid projection done. Took %s.\n", stock("main_grid_projection"));
	tick(); dump( "Fetching gradient and computing new FLIP velocity..." );
	// Compute particle new FLIP velocity
	FLOAT64 max_grad = 0.0;
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		if( p.levelset < 0.0 ) {
			vec3d pos = p.p;
			vec3d gu = fluidSolver.getVelocity(pos);
			vec3d pu = p.u;
			vec3d gradp = fluid->evalLevelset(pos) < 0.0 ? fluidSolver.getPressureGradient(pos) : vec3d();
			max_grad = fmax(max_grad,gradp.len());
			FLOAT64 visc = fmax(0.1,fmin(1.0,flip_visc/sqr(p.r)));
			p.u = visc*gu+(1.0-visc)*pu-dt*gradp;
		}
	} END_FOR
	dump( "Done. Took %s.\n", stock("main_pressure_gradient_fetch"));
	dump( "<<< FLIP projection done. Took %s.\n", stock("main_projection"));
}

void flip3::removeAllParticles( std::vector<particle3 *> &particles ) {
	for( uint n=0; n<particles.size(); n++ ) {
		particles[n]->setRemovable(true);
	}
	particle3::cleanParticles(particles);
	particles.clear();
}

void flip3::seedParticle( std::vector<particle3 *> &particles, const levelset3* fluid, const levelset3* solid, const octree3 *octree ) {
	FLOAT64 e = 1e-3*dpx;
	if( ! octree ) {
		FOR_EACH(gn,gn,gn) {
			for( FLOAT64 x = 0.5*dpx; x < dx; x += dpx ) {
				for( FLOAT64 y = 0.5*dpx; y < dx; y += dpx ) {
					for( FLOAT64 z = 0.5*dpx; z < dx; z += dpx ) {
						// Place particles if levelset is negative there
						vec3d p( dx*i+x, dx*j+y, dx*k+z );
						if( solid->evalLevelset(p) > 0 && fluid->evalLevelset(p) < 0.0 ) {
							particle3 *newp = new particle3;
							newp->p = p;
							newp->u = vec3d();
							particles.push_back(newp);
						}
					}
				}
			}
		} END_FOR
	} else {
		// Seed particles per terminal cell
		for( uint n=0; n<octree->terminals.size(); n++ ) {
			const octree3::leaf3 *leaf = octree->terminals[n];
			FLOAT64 dx = leaf->dx/(FLOAT64)octree->resolution;
			vec3d corner = vec3d(leaf->center)/octree->resolution-0.5*dx*vec3d(1,1,1);
			int r = dx > dpx ? 2 : 1;
			FOR_EACH(r,r,r) {
				FLOAT64 padding = 0.5*dx/r;
				vec3d p(corner[0]+padding+i*dx/r,corner[1]+padding+j*dx/r,corner[2]+padding+k*dx/r);
				p += e*vec3d(nrand(),nrand(),nrand());
				if( solid->evalLevelset(p) > 0.0 && fluid->evalLevelset(p) < 0.0 ) {
					particle3 *newp = new particle3;
					newp->p = p;
					newp->u = vec3d();
					newp->n = fmax(1.0,powf(0.5*dx/dpx,DIM));
					newp->computeRadius();
					particles.push_back(newp);
				}
			} END_FOR
		}
	}
}

void flip3::addForce( std::vector<particle3 *> &particles, FLOAT64 gravity, FLOAT64 dt ) {
	// Add gravity force
	FOR_EACH_PARTICLE(particles) {
		p.u[1] -= dt*gravity;
	} END_FOR
}

void flip3::advectParticles( sorter3& sorter, const levelset3* solid,
							 std::vector<particle3 *> &particles, const fluid3 &fluidSolver, FLOAT64 dt ) {
	tick(); dump("Advecting %d particles...", particles.size() );
	
	// Advect using 2nd order RK sampling
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		FLOAT64 r = dpx*powf(2,discreteDepthLevelset->evalLevelset(p.p)-1);
		if( ! p.isolated && p.levelset <= 0.0 && dpx*p.r > grid_tol*r ) {
			vec3d v1 = fluidSolver.getVelocity(p.p)-dt*fluidSolver.getPressureGradient(p.p);
			vec3d p2 = p.p-dt*v1;
			for( uint dim=0; dim<DIM; dim++ ) p2[dim] = fmin(1.0,fmax(0.0,p2[dim]));
			vec3d v2 = fluidSolver.getVelocity(p2)-dt*fluidSolver.getPressureGradient(p2);
			// Why backward sampling ? technically speaking, this should be forward sampling.
			// But in practice, forward sampling weiredly changes the behavior, and it seems that the backward sampling
			// acts more naturally ( I guess because this is a family of upwind sampling ? ).
			// This may be because we extrapolate the velocity by nearest neighbor fashion.
			// Accurately solving Eikonal equation may improve the issue.
			p.p += 0.5*dt*(v1+v2);
		} else {
			p.p += dt*p.u;
		}
	} END_FOR
	
	// Push out particles out of objects, easy because we know the levelset value
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		FLOAT64 phi;
		if( (phi = solid->evalLevelset(p.p)-0.5*dpx) < 0.0 ) {
			vec3d normal = solid->evalGradient(p.p).normal();
			p.p += -phi*normal;
			p.u -= (normal*p.u)*normal;
		}
		for( uint dim=0; dim<DIM; dim++ ) p.p[dim] = fmin(1.0,fmax(0.0,p.p[dim]));
	} END_FOR
	
	// Set sorter as dirty
	sorter.setDirty();
	dump("Done. Took %s.\n",stock("main_advection"));
}

bool flip3::getClosestSurfacePos(vec3d &pos) const {
	FLOAT64 min_phi = 1e9;
	std::vector<particle3 *> neighbors = surf_sorter->getkNeighbors(pos,numSample/2);
	FLOAT64 r = 0.0;
	FLOAT64 rsum = 0.0;
	for( uint n=0; n<neighbors.size(); n++ ) {
		r += neighbors[n]->r;
		rsum += 1.0;
	}
	if( rsum ) r = r / rsum;
	vec3d p = pos;
	min_phi = fmin(min_phi,levelset3::getParticleLevelset(neighbors,dpx,p,pos,solid));
	return fabs(min_phi) < 0.5*dpx*r;
}

void flip3::sortSurfaceParticles( sorter3 &sorter, const std::vector<particle3 *> &particles, std::vector<particle3 *> &surf_particles ) {
	surf_particles.clear();
	FOR_EACH_PARTICLE(particles) {
		if( p.levelset > -2.0*p.r*dpx ) {
			surf_particles.push_back(&p);
		}
	} END_FOR
	sorter.sortParticles(surf_particles);
}

FLOAT64 flip3::evalLevelset( vec3d pos ) const {	
	// Compute old levelset
	FLOAT64 oldLevelset;
	if( fluidSolver->getOctree() ) oldLevelset = octLevelset.evalLevelset(pos);
	else oldLevelset = surf0.evalLevelset(pos);
	
	// Compute base grid size
	FLOAT64 base_dx = 1.0;
	const octree3 *octree = fluidSolver->getOctree();
	if( octree ) {
		std::vector<uint> array;
		if( octree->hitTest( pos, 0.5*dx, array ) ) {
			for( uint n=0; n<array.size(); n++ ) {
				int hit = array[n];
				base_dx = fmin(base_dx,octree->terminals[hit]->dx/(FLOAT64)octree->resolution);
			}
		} else {
			dump( "flip3::evalLevelset base grid size failed to fetch at (%f,%f,%f)\n", pos[0], pos[1], pos[2] );
		}
	} else {
		base_dx = dx;
	}
	
	// Compute solid levelset there
	FLOAT64 solid_phi = solid->evalLevelset(pos);
	
	// Slide position
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
	std::vector<particle3 *> neighbors = sorter->getkNeighbors(pos,numSample);
	
	vec3d gradient;
	FOR_EACH_PARTICLE(neighbors) {
		vec3d out = pos;
		min_phi = fmin(min_phi,levelset3::getParticleConvextHullLevelset(out,p,dpx,shrink,false));
		LCFL = fmax(LCFL,p.u.len()*dt/base_dx);
		gradient += p.gradient;
	} END_FOR
	
	if( fabs(oldLevelset) > fmax(1.0,LCFL)*base_dx ) return fmin(min_phi,oldLevelset);
	
	// Evaluate exact distance
	if( fabs(min_phi) < fmax(1.0,LCFL)*base_dx ) {
		vec3d outpos;
		neighbors = surf_sorter->getkNeighbors(pos,numSample);
		min_phi = fmin(min_phi,levelset3::getParticleLevelset(neighbors,dpx,pos,outpos,solid,shrink));
	}
	
	// Compute splash particle levelset
	if( min_phi > 0 ) {
		min_phi = fmin(min_phi,ann.evalAbsLevelset(pos,dpx));
	}
	return min_phi;
}

FLOAT64 flip3::getDiffCurvature(vec3d pos) const {
	if( surf_particles.empty() ) return 0.0;
	std::vector<particle3 *> neighbors = surf_sorter->getkNeighbors(pos,numSample);
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

bool flip3::sampleVelocity( sorter3& sorter, vec3d pos, FLOAT64 sample_visc, FLOAT64 dpx, vec3d &vel ) {
	FLOAT64 scale = fluidSolver->getMesh() == NULL ? 2.0 : 1.0;
	FLOAT64 base_r = dpx*powf(2,discreteDepthLevelset->evalLevelset(pos)-1);
	FLOAT64 r = scale*sample_visc*base_r;
	std::vector<particle3 *> neighbors = sorter.getNeighbors(pos,r);
	FLOAT64 wsum = 0.0;
	vec3d out_velocity;
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

void flip3::resampleLevelsetOnOctree(const levelset3 *levelset, octlevelset3 &octLevelset, const octree3 &octree) {
	tick(); dump("Sampling levelset on octree...");
	octLevelset.setOctree(octree);
	std::vector<FLOAT64> q(octree.nodes.size());
	PARALLEL_FOR for( uint n=0; n<octree.nodes.size(); n++ ) {
		q[n] = levelset->evalLevelset(octree.nodes[n]);
	}
	octLevelset.setNodalLevelset(q);
	dump("Done. Took %s.\n",stock());
}

void flip3::collectSplashParticles( std::vector<particle3 *> &bullet_particles, std::vector<particle3 *> &clustered_particles,
								   pann3 &ann, const sorter3 &sorter, const std::vector<particle3 *> &particles, const levelset3 *fluid ) {
	tick(); dump("Collecting splash particles...");
	// Collect isolated particles
	bullet_particles.clear();
	clustered_particles.clear();
	
	PARALLEL_FOR FOR_EACH_PARTICLE(particles) {
		FLOAT64 levelset = fluid->evalLevelset(p.p);
		FLOAT64 r = dpx*p.r;
		if( levelset < 0.0 ) r = fmax(r,dpx*powf(2,discreteDepthLevelset->evalLevelset(p.p)-1));
		std::vector<particle3 *> neighbors = sorter.getkNeighbors(p.p,2);
		p.isolated = true;
		for( uint n=0; n<neighbors.size(); n++ ) {
			if( neighbors[n] != &p && (p.p-neighbors[n]->p).len2() < sqr(3.0*r) ) {
				p.isolated = false;
				break;
			}
		}
	} END_FOR
	
	uint isolated_num = 0;
	FOR_EACH_PARTICLE(particles) {
		if( p.isolated ) {
			bullet_particles.push_back(&p);
			isolated_num ++;
		}
	} END_FOR
	
	// Collect cluster particles
	FOR_EACH_PARTICLE(particles) {
		if( p.levelset > 0.5*p.r*dpx && ! p.isolated ) {
			bullet_particles.push_back(&p);
			clustered_particles.push_back(&p);
		}
	} END_FOR
	dump("Done. Collected %d bullet particles %d isolated particles and %d clustered particles. Took %s.\n",
			bullet_particles.size(), isolated_num, clustered_particles.size(), stock("main_splash_collect"));
	
	// Put into ANN
	ann.buildKDTree(clustered_particles);
}

