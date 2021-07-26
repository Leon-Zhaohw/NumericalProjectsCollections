/*
 *	fluid2.h
 *	
 *	Created by Ryoichi Ando on 12/9/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "vec2.h"
#include <vector>

#ifndef _FLUID2_H
#define _FLUID2_H

class levelset2;
class bccmesher2;
class octree2;
class mesher2;
class fluid2 {
public:
	// Parameter list
	enum { EXTRAPOLATE_DIST, SURF_ORDER, VARIATION };
	// Render list
	enum { PRESSURE, CURVATURE, DIVERGENCE, FLUID, SOLID, VELOCITY, MESH, MATRIX_CONNECTION, NAME };
	
	fluid2() {}
	virtual ~fluid2() {}
	
	// 1st step: setup simulation configuration. You have to call them in defined order.
	virtual void init( uint gn, const levelset2 *hint ) {}
	virtual void setParameter( int name, FLOAT64 value ) {}
	virtual void setTimestep( FLOAT64 dt ) {}
	virtual void setupSolidLevelset( const levelset2 *solid ) {}
	virtual void setupFluidLevelset( const levelset2 *fluid ) {}
	virtual void getSamplePoints( std::vector<vec2d> &pos ) {}
	virtual void setupVelocity( const std::vector<vec2d> &vel, const std::vector<bool> &mapped ) {}
	
	// 2nd step: project
	virtual void project() {}
	
	// 3rd step: retrieve information
	virtual vec2d getPressureGradient( vec2d p ) const { return vec2d(); }
	virtual vec2d getVelocity( vec2d p ) const { return vec2d(); }
	virtual FLOAT64 getDivergence( vec2d p ) const { return 0.0; }
	virtual FLOAT64 getStrain( vec2d p ) const { return 0.0; }
	
	// Render variables
	virtual void render( int name, vec2d mousePos=vec2d() ) const {}
	
	// Retrive some solver info
	virtual const mesher2 *getMesh() const { return NULL; }
	virtual const bccmesher2 *getBCCMesh() const { return NULL; }
	virtual const octree2 *getOctree() const { return NULL; }
};

#endif