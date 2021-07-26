/*
 *	fluid3.h
 *	
 *	Created by Ryoichi Ando on 2012/05/03
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <stdlib.h>
#include "macros.h"
#include "vec3.h"
#include <vector>

#ifndef _FLUID3_H
#define _FLUID3_H

class levelset3;
class bccmesher3;
class octree3;
class mesher3;
class fluid3 {
public:
	// Parameter list
	enum { EXTRAPOLATE_DIST, SURF_ORDER, VARIATION };
	// Render list
	enum { PRESSURE, DIVERGENCE, LEVELSET, VELOCITY, MESH, NAME };
	
	fluid3() {}
	virtual ~fluid3() {}
	
	// 1st step: setup simulation configuration. You have to call them in defined order.
	virtual void init( uint gn, const levelset3 *hint ) {}
	virtual void setParameter( int name, FLOAT64 value ) {}
	virtual void setTimestep( FLOAT64 dt ) {}
	virtual void setupSolidLevelset( const levelset3 *solid ) {}
	virtual void setupFluidLevelset( const levelset3 *fluid ) {}
	virtual void getSamplePoints( std::vector<vec3d> &pos ) {}
	virtual void setupVelocity( const std::vector<vec3d> &vel, const std::vector<bool> &mapped ) {}
	
	// 2nd step: project
	virtual void project() {}
	
	// 3rd step: retrieve information
	virtual vec3d getPressureGradient( vec3d p ) const { return vec3d(); }
	virtual vec3d getVelocity( vec3d p ) const { return vec3d(); }
	virtual FLOAT64 getDivergence( vec3d p ) const { return 0.0; }
	virtual FLOAT64 getStrain( vec3d p ) const { return 0.0; }
	
	// Render variables
	virtual void render( int name ) const {}
	virtual const char *getName() const { return ""; }
	
	// Retrive some solver info
	virtual const mesher3 *getMesh() const { return NULL; }
	virtual const bccmesher3 *getBCCMesh() const { return NULL; }
	virtual const octree3 *getOctree() const { return NULL; }
};

#endif
