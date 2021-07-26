/*
 *	macfluid2.h
 *	
 *	Created by Ryoichi Ando on 12/7/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "fluid2.h"
#include "array2.h"

#ifndef _MACFLUID2_H
#define _MACFLUID2_H

class macfluid2 : public fluid2 {
public:
	macfluid2();
	
	// 1st step: setup simulation configuration
	virtual void init( uint gn, const levelset2 *hint );
	virtual void setParameter( int name, FLOAT64 value );
	virtual void setTimestep( FLOAT64 dt );
	virtual void setupSolidLevelset( const levelset2 *solid );
	virtual void setupFluidLevelset( const levelset2 *fluid );
	virtual void getSamplePoints( std::vector<vec2d> &pos );
	virtual void setupVelocity( const std::vector<vec2d> &vel, const std::vector<bool> &mapped );
	
	// 2nd step: project
	virtual void project();
	
	// 3rd step: retrieve information
	virtual vec2d getPressureGradient( vec2d p ) const;
	virtual vec2d getVelocity( vec2d p ) const;
	virtual FLOAT64 getDivergence( vec2d p ) const;
	
	// Render variables
	virtual void render( int name, vec2d mousePos=vec2d() ) const;
protected:
	uint gn;
	FLOAT64 dx;
	FLOAT64 dt;
	const levelset2 *solid;				// Solid levelset
	const levelset2 *fluid;				// Fluid levelset
	array2<FLOAT64> u[DIM];				// Cellface fluid velocity
	array2<bool> flag[DIM];				// Cellface map flag
	array2<FLOAT64> gradp[DIM];			// Cellface pressure gradient
	array2<FLOAT64> area[DIM];			// Area fraction
	array2<FLOAT64> volume;				// Volume fraction
	array2<FLOAT64> cellfluidLS;		// Cell centered levelset
	array2<FLOAT64> facetfluidLS[DIM];	// Facet centered fluid levelset
	array2<FLOAT64> facetsolidLS[DIM];	// Facet centered solid levelset
	array2<FLOAT64> cellCurvature;
	
	FLOAT64 extrapolate_dist;			// Extrapolate distance
	uint surf_order;
	uint variation;
	FLOAT64 tension;
	FLOAT64 frac_eps;
private:
	array2<FLOAT64> p;                   // Pressure in 2D
	array2<FLOAT64> div;                 // Divergence in 2D
    
	FLOAT64 computeVolume( const array2<FLOAT64> &levelset );
	FLOAT64 computeDivergenceControl( FLOAT64 volume0, FLOAT64 accmVolumeErr, FLOAT64 volume );
	void computeFractions( array2<FLOAT64> area[], array2<FLOAT64> &volume, const levelset2 *solid );
	void extrapolate( array2<FLOAT64> u[], const array2<bool> flag[] );
	vec2d getVelocity( vec2d p, const array2<FLOAT64> u[] ) const;
	void copyGridVelocity( array2<FLOAT64> dest[], const array2<FLOAT64> src[] );
	void project(array2<FLOAT64> u[], const array2<FLOAT64> &cellfluidLS, const levelset2* solid,
				 bool variation, uint surf_order, FLOAT64 dt );
};

#endif