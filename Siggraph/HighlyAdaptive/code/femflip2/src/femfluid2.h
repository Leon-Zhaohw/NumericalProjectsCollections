/*
 *	femfluid2.h
 *	
 *	Created by Ryoichi Ando on 12/9/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "fluid2.h"
#include "umesher2.h"
#include "bccmesher2.h"
#include "pcgsolver/sparse_matrix.h"
#include <vector>

#ifndef _FEMFLUID2_H
#define _FEMFLUID2_H

class femfluid2 : public fluid2 {
public:
	femfluid2();
    
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
	
    // Render variables
	virtual void render( int name, vec2d mousePos=vec2d() ) const;
	
	// 3rd step: retrieve information
	virtual vec2d getPressureGradient( vec2d p ) const;
	virtual vec2d getVelocity( vec2d p ) const;
	virtual FLOAT64 getDivergence( vec2d p ) const;
	virtual FLOAT64 getStrain( vec2d p ) const;
	
	// Retrieve some info
	virtual const bccmesher2 *getBCCMesh() const { return &g; }
	virtual const mesher2 *getMesh() const { return &g; }
	virtual const octree2 *getOctree() const { return &g.octree; }
	
protected:	
	typedef struct _node2 {
		vec2d velocity;			// Velocity
		vec2d gradient;			// Gradient
		bool mapped;
		FLOAT64 pressure;		// Pressure value
		FLOAT64 divergence;		// Divergence value
		FLOAT64 solidLS;		// Solid levelset
		FLOAT64 fluidLS;		// Fluid levelset
		FLOAT64 volume;			// One ring volume
		FLOAT64 strain;			// Strain tensor
	} node2;
	std::vector<node2> nodes;
	
	typedef struct _element2 {
		vec2d velocity;			 // Velocity
		vec2d gradient;			 // Gradient
		FLOAT64 solidLS;		 // Solid levelset
		FLOAT64 fluidLS;		 // Fluid levelset
		bool mapped;									// Whether the velocity is mapped
		FLOAT64 VBtB[NUM_VERT][NUM_VERT];				// Laplace matrix
		FLOAT64 rho;									// Fluid density  (0 ~ 1)
		FLOAT64 frac;									// Solid fraction (0 ~ 1)
	} element2;
	std::vector<element2> elements;

public:
	bccmesher2 g;
	uint gn;
    FLOAT64 dt;
    const levelset2 *solid;
    const levelset2 *fluid;
	SparseMatrix<FLOAT64>		LhsMatrix;			// Left hand side matrix (Num nodes) x (Num nodes)
	SparseMatrix<FLOAT64>		RhsMatrix;			// Right hand side matrix (Num nodes) x (Num nodes x dim)
	
	uint surf_order;
	FLOAT64 extrapolate_dist;
	bool variation;
private:
	FLOAT64 dx;
	void initMesh(const levelset2 *hint);
	void resample( uint type, uint kind ); //type 0: e->v 1: v->e
	void extrapolate( const mesher2 &g, std::vector<node2> &nodes, int kind ); // 0 for velocity 1 for gradient
	int hitElements( const std::vector<vec2d> &nodes, const std::vector<std::vector<uint> > &elements, const array2<std::vector<uint> > &hash, const vec2d p ) const;
	bool hitElement( const std::vector<vec2d> &nodes, const std::vector<uint> &element, vec2d p ) const;
	bool computeGhostPressure( uint elm_idx, uint index, FLOAT64 coef[NUM_VERT] ) const;
	void embedGhostPressure( uint elm_idx, FLOAT64 press[NUM_VERT] ) const;
	void buildMatrix();
	void printDivergence( int kind, int step );
	void computeElementLiquidFraction();
	void computeStrain();
	vec2d getVector( vec2d p, uint kind ) const;
};

#endif