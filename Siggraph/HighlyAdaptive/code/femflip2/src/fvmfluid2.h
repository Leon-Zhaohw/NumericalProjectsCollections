/*
 *	fvmfluid2.h
 *	
 *	Created by Ryoichi Ando on 12/27/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "pcgsolver/sparse_matrix.h"
#include "fluid2.h"
#include "umesher2.h"
#include "bccmesher2.h"

#ifndef _FVMFLUID2_H
#define _FVMFLUID2_H

class fvmfluid2 : public fluid2 {
public:
	fvmfluid2();
	
	// 1st step: setup simulation configuration. You have to call them in defined order.
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
	virtual const mesher2 *getMesh() const { return &g; }
protected:
	typedef struct _facet2 {
		vec2d center;
		vec2d normal;
		FLOAT64 velocity;
		FLOAT64 gradient;
		FLOAT64 fraction;
		bool mapped;
	} facet2;
	std::vector<facet2> facets;
	
	typedef struct _cell2 {
		FLOAT64 pressure;
		FLOAT64 divergence;
		FLOAT64 solidVolume;
		FLOAT64 solidLS;
		FLOAT64 fluidLS;
		vec2d velocity;
		vec2d gradient;
		bool mapped;
	} cell2;
	std::vector<cell2> cells;
	
	typedef struct _node2 {
		vec2d velocity;
		vec2d gradient;
		bool mapped;
		FLOAT64 fluidLS;
		FLOAT64 solidLS;
		FLOAT64 volume;
	} node2;
	std::vector<node2> nodes;
	
	umesher2 g;		// Mesh geometry
	uint gn;		// Base resolution
	FLOAT64 dx;
	FLOAT64 dt;
	FLOAT64 extrapolate_dist;
	const levelset2 *solid;
	const levelset2 *fluid;
	bool variation;
	uint surf_order;

	SparseMatrix<FLOAT64> LhsMatrix;
	SparseMatrix<FLOAT64> RhsMatrix;

	void initMesh(const levelset2 *hint);
protected:
	void buildMatrix();
	void computeElementVelocity( const mesher2 &g, const std::vector<facet2> &facets, std::vector<cell2> &cells, int kind );
	void extrapolate(const mesher2 &g, std::vector<node2> &nodes, int kind);
	vec2d getVector( vec2d p, uint kind ) const;
	void resample( uint type, uint kind );
};

#endif