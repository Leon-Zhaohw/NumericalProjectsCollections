/*
 *	femfluid3.h
 *	
 *	Created by Ryoichi Ando on 2012/05/03
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "fluid3.h"
#include "mesher3.h"
#include "bccmesher3.h"
#include "pcgsolver/sparse_matrix.h"
#include <vector>

#ifndef _FEMFLUID3_H
#define _FEMFLUID3_H

class bcc3;
class femfluid3 : public fluid3 {
public:
	femfluid3();
    
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
	
    // Render variables
	virtual void render( int name ) const;
	virtual const char *getName() const { return "FEM"; }
	
	// 3rd step: retrieve information
	virtual vec3d getPressureGradient( vec3d p ) const;
	virtual vec3d getVelocity( vec3d p ) const;
	virtual FLOAT64 getDivergence( vec3d p ) const;
	virtual FLOAT64 getStrain( vec3d p ) const;
	
	// Retrieve some info
	virtual const bccmesher3 *getBCCMesh() const { return &g; }
	virtual const mesher3 *getMesh() const { return &g; }
	virtual const octree3 *getOctree() const { return &g.octree; }
	
protected:
	typedef struct _node3 {
		vec3d velocity;			// Velocity
		vec3d gradient;			// Gradient
		bool mapped;
		FLOAT64 pressure;		// Pressure value
		FLOAT64 divergence;		// Divergence
		FLOAT64 solidLS;		// Solid levelset
		FLOAT64 fluidLS;		// Fluid levelset
		FLOAT64 volume;			// One ring volume
		FLOAT64 strain;			// Strain
	} node3;
	std::vector<node3> nodes;
	
	typedef struct _element3 {
		vec3d velocity;			 // Velocity
		vec3d gradient;			 // Gradient
		FLOAT64 solidLS;		 // Solid levelset
		FLOAT64 fluidLS;		 // Fluid levelset
		bool mapped;									// Whether the velocity is mapped
		FLOAT64 VBtB[NUM_VERT][NUM_VERT];				// Laplace matrix
		FLOAT64 rho;									// Fluid density  (0 ~ 1)
		FLOAT64 frac;									// Solid fraction (0 ~ 1)
	} element3;
	std::vector<element3> elements;

	bccmesher3 g;
	uint gn;
    FLOAT64 dt;
    const levelset3 *solid;
    const levelset3 *fluid;
	SparseMatrix<FLOAT64>		LhsMatrix;			// Left hand side matrix (Num nodes) x (Num nodes)
	SparseMatrix<FLOAT64>		RhsMatrix;			// Right hand side matrix (Num nodes) x (Num nodes x dim)
	
	uint surf_order;
	FLOAT64 tension;
	FLOAT64 extrapolate_dist;
	bool variation;
private:
	FLOAT64 dx;
	void initMesh(const levelset3 *hint);
	void resample( uint type, uint kind ); //type 0: e->v 1: v->e
	void extrapolate( const mesher3 &g, std::vector<node3> &nodes, int kind );
	int hitElements( const std::vector<vec3d> &nodes, const std::vector<std::vector<uint> > &elements, const array3<std::vector<uint> > &hash, const vec3d p ) const;
	bool hitElement( const std::vector<vec3d> &nodes, const std::vector<uint> &element, vec3d p ) const;
	bool computeGhostPressure( uint elm_idx, uint index, FLOAT64 coef[NUM_VERT]	) const;
	void embedGhostPressure( uint elm_idx, FLOAT64 press[NUM_VERT] ) const;
	void buildMatrix();
	void printDivergence( int kind, int step );
	vec3d getVector( vec3d p, uint kind ) const;
	void computeElementLiquidFraction();
	void computeStrain();
};

#endif
