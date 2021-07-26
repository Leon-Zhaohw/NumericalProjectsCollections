/*
 *	fvmfluid3.h
 *	
 *	Created by Ryoichi Ando on 2012/05/03
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "pcgsolver/sparse_matrix.h"
#include "fluid3.h"
#include "bccmesher3.h"

#ifndef _FVMFLUID3_H
#define _FVMFLUID3_H

class fvmfluid3 : public fluid3 {
public:
	fvmfluid3();
	
	// 1st step: setup simulation configuration. You have to call them in defined order.
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
	virtual const bccmesher3 *getBCCMesh() const { return &g; }
	virtual const mesher3 *getMesh() const { return &g; }
	virtual const octree3 *getOctree() const { return &g.octree; }
	virtual const char *getName() const { return "FVM"; }
	
protected:
	typedef struct _facet2 {
		vec3d center;
		vec3d normal;
		FLOAT64 velocity;
		FLOAT64 gradient;
		FLOAT64 fraction;
		bool mapped;
	} facet3;
	std::vector<facet3> facets;
	
	typedef struct _cell3 {
		FLOAT64 pressure;
		FLOAT64 divergence;
		FLOAT64 solidVolume;
		FLOAT64 solidLS;
		FLOAT64 fluidLS;
		vec3d velocity;
		vec3d gradient;
		bool mapped;
	} cell3;
	std::vector<cell3> cells;
	
	typedef struct _node3 {
		vec3d velocity;
		vec3d gradient;
		bool mapped;
		FLOAT64 fluidLS;
		FLOAT64 solidLS;
		FLOAT64 volume;
	} node3;
	std::vector<node3> nodes;
	
	bccmesher3 g;	// Mesh geometry
	uint gn;	// Base resolution
	FLOAT64 dx;
	FLOAT64 dt;
	FLOAT64 extrapolate_dist;
	const levelset3 *fluid;
	const levelset3 *solid;
	bool variation;
	uint surf_order;
	FLOAT64 tension;
	
	SparseMatrix<FLOAT64> LhsMatrix;
	SparseMatrix<FLOAT64> RhsMatrix;
	
	void initMesh(const levelset3 *hint);
protected:
	void buildMatrix();
	void computeElementVelocity( const mesher3 &g, const std::vector<facet3> &facets, std::vector<cell3> &cells, int kind );
	void extrapolate(const mesher3 &g, std::vector<node3> &nodes, int kind);
	vec3d getVector( vec3d p, uint kind ) const;
	void resample( uint type, uint kind );
};

#endif
