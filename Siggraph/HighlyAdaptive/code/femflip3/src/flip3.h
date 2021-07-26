/*
 *	flip3.h
 *	
 *	Created by Ryoichi Ando on 1/9/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "vec3.h"
#include "array3.h"
#include "annsorter3.h"
#include "cellsurf3.h"
#include "corrector3.h"
#include "adaptive3.h"
#include "levelset3.h"
#include "surf3.h"
#include "camera3.h"
#include "octlevelset3.h"
#include "object3.h"
#include "macfluid3.h"
#include "femfluid3.h"
#include "fvmfluid3.h"
#include "octfluid3.h"
#include "pann3.h"
#include <vector>

#ifndef _FLIP3_H_
#define _FLIP3_H_

class particle3;
class glviewer;
class exporter3;
class remeshLevelset3;
class discreteDepthLevelset3;
class flip3 : public levelset3 {
public:
	flip3();																	// Don't forget to call init() !
	void init( uint gsize, std::vector<object3 *> objects, FLOAT64 gravity, const camera3 *camera );	// Setup simulation
	FLOAT64 simStep( FLOAT64 tframe, FLOAT64 maxdt, FLOAT64 mindt );			// Advance a timestep
	void setElapsedTime( FLOAT64 elapsedTime );									// Set elapsed time (optional)
	void setVideoFrameNumber( uint frame );										// Set video frame number (optional)

	uint	gn;							// Underlying grid resolution ( gn^3 )
	uint	gm;							// Surface mesh resolution
	uint	numSample;					// Number of samples to evaluate levelset
	FLOAT64	dx;							// Grid cell width
	uint	min_resolution;				// Minimum resolution
	uint	remeshRate;					// Remesh rate
	FLOAT64 gravity;					// Gravity
	FLOAT64 tension;					// Surface tension
	FLOAT64 dpx;						// Particle space
	uint	N0;							// Particle number in a cell at initial placement
	FLOAT64 dt;							// Timestep size
	FLOAT64 remesh_dt;					// Timestep size from previously remeshed time
	FLOAT64 CFL;						// CFL Number
	FLOAT64 extrapolation_dist;			// Extrapolation distance
	FLOAT64	flip_visc;					// FLIP viscosity
	FLOAT64	sample_visc;				// Sampling viscosity
	FLOAT64 tension_force;				// Surface tension force
	FLOAT64 anisotropyRad;				// Anisotropy Radius
	FLOAT64 msec;						// Simulation time per timestep in milliseconds
	FLOAT64 avg_msec;					// Simulation time per (remeshRate) timestep in milliseconds
	FLOAT64	sitffness;					// Position correction stiffness
	FLOAT64 svd_ratio[2];				// Stretch threshold to split particles
	FLOAT64 min_stretch;				// Minimum stretch of particles
	FLOAT64 grid_tol;					// Minimal particle size to be sampled on grid
	bool	anisotropic;				// Whether to use anisotropic spring correction
	bool	correction;					// Whether to correct particle position using weak spring
	bool    adaptiveSampling;           // Whether to do adaptive sampling
	bool	doRemesh;					// Whether to remesh
	bool	doubleRemesh;				// Whether to use precise remeshing

	std::vector<particle3 *> particles;				// Particle array
	std::vector<particle3 *> surf_particles;		// Surface particles
	std::vector<particle3 *> bullet_particles;		// Bullet particle array
	std::vector<particle3 *> clustered_particles;	// Clustered particle array
	sorter3 *sorter;					// Particle sorter
	annsorter3 annsorter;				// ANN-based particle sorter
	sorter3 *surf_sorter;				// Surface particle sorter
	annsorter3 surf_annsorter;			// ANN-based surface particle sorter
	adaptive3 adaptive;					// Adaptive merging and splitting operator
	corrector3 corrector;				// Particle position corrector
	std::vector<object3 *> objects;		// Objects array
	
	levelset3* solid;					// Solid levelset evaluator
	levelset3* fluid;					// Fluid (initial) levelset evaluator
	remeshLevelset3* remeshLevelset;	// Lelvelset for remesh
	discreteDepthLevelset3 *discreteDepthLevelset;
	surf3 surf;							// Surface levelset
	surf3 surf0;						// Old surface levelset
	octlevelset3 octLevelset;			// Octree levelset
	fluid3 *fluidSolver;				// Abstracted fluid solver
	macfluid3 macSolver;				// MAC fluid solver
	femfluid3 femSolver;				// FEM fluid solver
	fvmfluid3 fvmSolver;				// FVM fluid solver
	octfluid3 octSolver;				// OCT fluid solver
	pann3	ann;						// ANN particle searcher
	const camera3 *camera;				// 3D virtual camera
	uint	step;						// Simulation step
	FLOAT64 totalTime;					// Total time
	uint	frameNumber;				// Video frame number
	FLOAT64 elapsedTime;				// Elapsed time
	FLOAT64 evalLevelset( vec3d pos ) const;
	virtual bool getClosestSurfacePos(vec3d &pos) const;
	FLOAT64 getDiffCurvature(vec3d pos) const;
private:
	void fitParticles( std::vector<particle3 *> &particles, const levelset3 *levelset );
	void computeParameters( FLOAT64 tframe, FLOAT64 maxdt, FLOAT64 mindt );
	void removeAllParticles( std::vector<particle3 *> &particles );
	void seedParticle( std::vector<particle3 *> &particles, const levelset3* fluid, const levelset3* solid, const octree3 *octree=NULL );
	void advectParticles( sorter3& sorter, const levelset3* solid, std::vector<particle3 *> &particles, const fluid3 &fluidSolver, FLOAT64 dt );
	void addForce( std::vector<particle3 *> &particles, FLOAT64 gravity, FLOAT64 dt );
	void projectFLIP( fluid3 &fluidSolver, sorter3& sorter, std::vector<particle3 *> &particles, const levelset3* solid, const levelset3* fluid, FLOAT64 flip_visc, FLOAT64 sample_visc, FLOAT64 dpx, bool doRemesh, FLOAT64 dt );
	bool sampleVelocity( sorter3& sorter, vec3d pos, FLOAT64 sample_visc, FLOAT64 dpx, vec3d &vel );
	void resampleLevelsetOnOctree(const levelset3 *levelset,octlevelset3 &octLevelset, const octree3 &octree);
	void collectSplashParticles( std::vector<particle3 *> &bullet_particles, std::vector<particle3 *> &clustered_particles, pann3 &ann, const sorter3 &sorter, const std::vector<particle3 *> &particles, const levelset3 *fluid );
	void writeStatus() const;
	void sortSurfaceParticles( sorter3 &sorter, const std::vector<particle3 *> &particles, std::vector<particle3 *> &surf_particles );
};

#endif
