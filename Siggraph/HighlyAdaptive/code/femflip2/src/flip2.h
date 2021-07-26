/*
 *	flip2.h
 *	
 *	Created by Ryoichi Ando on 11/3/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "vec2.h"
#include "array2.h"
#include "annsorter2.h"
#include "corrector2.h"
#include "adaptive2.h"
#include "macfluid2.h"
#include "femfluid2.h"
#include "fvmfluid2.h"
#include "octfluid2.h"
#include "levelset2.h"
#include "surf2.h"
#include "camera2.h"
#include "octlevelset2.h"
#include "pann2.h"
#include <vector>

#ifndef _FLIP2_H_
#define _FLIP2_H_

// Don't want to include them here...
class particle2;
class glviewer;
class glhelpviewer;
class remeshLevelset2;
class discreteDepthLevelset2;
class flip2 : public levelset2 {
public:
	flip2();							// Don't forget to call init() !
	void init( uint gsize, levelset2* fluid, levelset2* solid, FLOAT64 gscale ); // Setup simulation domain (Can be called anytime)
	FLOAT64 simStep( FLOAT64 tframe, FLOAT64 maxdt, FLOAT64 mindt );			 // Advance a timestep
	// Set functions
	void setBoundaryAccuracy( uint surf_order );	// Accuracy of boundary condition. 2nd order default
	void setVariation( bool variation );			// Whether to use variational framework for pressure solve
	void setCorrection( bool enabled ); // Whether to use weak spring correction to improve particle spacing
    void setAdaptiveSampling( bool enabled );       // Whether to use adaptive sampling
	void setRemesh( bool enabled );					// Whether to remesh if available
	void setFluidSolver( uint name );				// Which solver to solve fluid flow

	uint	gn;							// Underlying grid resolution ( gn x gn )
	uint	gm;							// Surface resolution
	uint	step;						// Simulation step
	uint	min_resolution;				// Minimum resolution
	uint	remeshRate;					// Remesh rate
	bool	doubleRemesh;				// Double remeshing
	uint	numSample;					// Sample number of levelset
	FLOAT64	dx;							// Grid cell width
	FLOAT64 dpx;						// Particle space
	uint	N0;							// Particle number in a cell at initial placement
	FLOAT64 dt;							// Timestep size
	FLOAT64 remesh_dt;					// Timestep size from previously remeshed time
	FLOAT64 CFL;						// CFL Number
	FLOAT64 extrapolation_dist;			// Extrapolation distance
	FLOAT64 totalTime;					// Accumulated time
	FLOAT64 gravity;					// Gravity acceleration.
	FLOAT64	flip_visc;					// FLIP viscosity
	FLOAT64	sample_visc;				// Sampling viscosity
	FLOAT64 anisotropyRad;				// Anisotropy Radius
	FLOAT64 msec;						// Simulation time per timestep in milliseconds
	uint	surf_order;					// Accuracy of boundary condition
	FLOAT64 densityRad;					// Density sampling radius
	FLOAT64	stiffness;					// Position correction stiffness
	FLOAT64 svd_ratio[2];				// Stretch threshold to split and collapse particles
	FLOAT64 min_stretch;				// Minimum stretch of particles
	FLOAT64 grid_tol;					// Minimal particle size to be sampled on grid
	
	bool	variation;					// Whether to use variational frame work or not for pressure solve
	bool	correction;					// Whether to correct particle position using weak spring
	bool	anisotropic;				// Whether to use anisotropic spring correction
    bool    adaptiveSampling;           // Whether to do adaptive sampling
	bool	doRemesh;					// Whether to do remeshing
	
	levelset2* solid;					// Solid levelset evaluator
	levelset2* fluid;					// Fluid (initial) levelset evaluator
	remeshLevelset2* remeshLevelset;	// Lelvelset for remesh
	levelset2 *discreteDepthLevelset;
	surf2 surf;							// Surface levelset
	surf2 surf0;						// Old surface levelset
	octlevelset2 octLevelset;			// Octree levelset
	std::vector<particle2 *> particles;	// Particle array
	std::vector<particle2 *> surf_particles;			// Surface particles
	std::vector<particle2 *> bullet_particles;			// Bullet particle array
	std::vector<particle2 *> clustered_particles;		// Clustered particle array
	sorter2 *sorter;					// Particle sorter
	annsorter2 annsorter;				// ANN-based particle sorter
	sorter2 *surf_sorter;				// Surface particle sorter
	annsorter2 surf_annsorter;			// ANN-based surface particle sorter
	adaptive2 adaptive;					// Adaptive merging and splitting operator
	corrector2 corrector;				// Particle position corrector
	fluid2 *fluidSolver;				// Abstracted fluid solver
	macfluid2 macSolver;				// MAC fluid solver
    femfluid2 femSolver;				// FEM fluid solver
	fvmfluid2 fvmSolver;				// FVM fluid solver
	octfluid2 octSolver;				// OCT fluid solver
	pann2 ann;							// ANN particle searcher
	camera2	camera;						// 2D virtual camera
	
	std::vector<vec2d> extForcePos;		// External force poistion list
	std::vector<vec2d> extForces;		// External forces list

	FLOAT64 evalLevelset(vec2d pos) const;
	virtual bool getClosestSurfacePos(vec2d &pos) const;
	FLOAT64 getDiffCurvature(vec2d pos) const;
	// Utility functions
	void addExternalForce( vec2d p, vec2d f );
private:
	// Private functions
	void fitParticles( std::vector<particle2 *> &particles, const levelset2 *levelset );
	void computeParameters( FLOAT64 tframe, FLOAT64 maxdt, FLOAT64 mindt );
	void parameterChanged();
	void removeAllParticles( std::vector<particle2 *> &particles );
	void seedParticle( std::vector<particle2 *> &particles, const levelset2* fluid, const levelset2* solid, const octree2 *octree=NULL );
	void advectParticles( sorter2& sorter, const levelset2* solid, std::vector<particle2 *> &particles, const fluid2 &fluidSolver, FLOAT64 dt );
	void addForce( sorter2& sorter, std::vector<particle2 *> &particles, std::vector<vec2d> &extForcePos, std::vector<vec2d> &extForces, FLOAT64 gravity, FLOAT64 dt );
	void projectFLIP( fluid2 &fluidSolver, sorter2& sorter, std::vector<particle2 *> &particles, const levelset2* solid, const levelset2* fluid, FLOAT64 flip_visc, FLOAT64 sample_visc, FLOAT64 dpx, bool doRemesh, FLOAT64 dt);
	bool sampleVelocity( sorter2& sorter, vec2d pos, FLOAT64 sample_visc, FLOAT64 dpx, vec2d &vel );
	void resampleLevelsetOnOctree(const levelset2 *levelset,octlevelset2 &octLevelset, const octree2 &octree);
	void collectSplashParticles( std::vector<particle2 *> &bullet_particles, std::vector<particle2 *> &clustered_particles,
								 pann2 &ann, const sorter2 &sorter, const std::vector<particle2 *> &particles, const levelset2 *fluid );
	void sortSurfaceParticles( sorter2 &sorter, const std::vector<particle2 *> &particles, std::vector<particle2 *> &surf_particles );
};

#endif