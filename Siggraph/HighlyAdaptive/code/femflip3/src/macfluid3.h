/*
 *	macfluid3.h
 *	
 *	Created by Ryoichi Ando on 2012/05/03
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "fluid3.h"
#include "array3.h"

#ifndef _MACFLUID3_H
#define _MACFLUID3_H

class macfluid3 : public fluid3 {
public:
	macfluid3();
	
	// 1st step: setup simulation configuration
	virtual void init( uint gn, const levelset3 *hint );
	virtual void setParameter( int name, FLOAT64 value );
	virtual void setTimestep( FLOAT64 dt );
	virtual void setupSolidLevelset( const levelset3 *solid );
	virtual void setupFluidLevelset( const levelset3 *fluid );
	virtual void getSamplePoints( std::vector<vec3d> &pos );
	virtual void setupVelocity( const std::vector<vec3d> &vel, const std::vector<bool> &mapped );
	
	// 2nd step: project
	virtual void project();
	
	// 3rd step: retrieve information
	virtual vec3d getPressureGradient( vec3d p ) const;
	virtual vec3d getVelocity( vec3d p ) const;
	virtual FLOAT64 getDivergence( vec3d p ) const;
	
	// Render variables
	virtual void render( int name ) const;
	virtual const char *getName() const { return "MAC"; }
protected:
	uint gn;
	FLOAT64 dx;
	FLOAT64 dt;
	const levelset3 *solid;				// Solid levelset
	const levelset3 *fluid;				// Fluid levelset
	array3<FLOAT64> u[DIM];				// Cellface fluid velocity
	array3<bool> flag[DIM];				// Cellface map flag
	array3<FLOAT64> gradp[DIM];			// Cellface pressure gradient
	array3<FLOAT64> area[DIM];			// Area fraction
	array3<FLOAT64> volume;				// Volume fraction
	array3<FLOAT64> cellfluidLS;		// Cell centered levelset
	array3<FLOAT64> facetfluidLS[DIM];	// Facet centered fluid levelset
	array3<FLOAT64> facetsolidLS[DIM];	// Facet centered solid levelset
	
	FLOAT64 extrapolate_dist;			// Extrapolate distance
	bool wall_separate;
	uint surf_order;
	uint variation;
private:
	array3<FLOAT64> p;                   // Pressure in 3D
	array3<FLOAT64> div;                 // Divergence in 3D
    
	void computeFractions( array3<FLOAT64> area[], array3<FLOAT64> &volume, const levelset3 *solid );
	void extrapolate( array3<FLOAT64> u[], const array3<bool> flag[] );
	vec3d getVelocity( vec3d p, const array3<FLOAT64> u[] ) const;
	void copyGridVelocity( array3<FLOAT64> dest[], const array3<FLOAT64> src[] );
	void project(array3<FLOAT64> u[], const array3<FLOAT64> &cellfluidLS, const levelset3* solid,
				 bool variation, uint surf_order, bool wallSepration, FLOAT64 dt );
};

#endif
