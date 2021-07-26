/*
 *	octfluid2.h
 *
 *	Created by Ryoichi Ando on 6/9/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "fluid2.h"
#include "octree2.h"
#include <list>
#include <map>
#include <vector>
#include "pcgsolver/sparse_matrix.h"

#ifndef OCTFLUID2_H
#define OCTFLUID2_H

class levelset2;
class octfluid2 : public fluid2 {
public:
	octfluid2();
	
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
	virtual const octree2 *getOctree() const { return &octree; }
protected:
	// Cell structure
	struct _facet2;
	struct _node2;
	typedef struct _cell2 {
		vec2i p;
		uint dx;
		FLOAT64 pressure;
		FLOAT64 divergence;
		FLOAT64 fluidLS;
		std::list<_facet2 *> facets;
		_node2 *nodes[2][2];
		// Used to extrapolate
		bool mapped;
		vec2d velocity;
		vec2d gradient;
		uint index;
	} cell2;
	
	// Facet structure
	typedef struct _facet2 {
		char dim; // X, Y, Z
		vec2i p;
		FLOAT64 velocity;
		FLOAT64 gradient;
		FLOAT64 fraction;
		bool mapped;
		bool edge;
		uint dx;
		uint index;
		std::list<cell2 *> cells;
		std::list<_node2 *> nodes;
		// Used for pressure solve
		FLOAT64 ds;
		cell2 *Lcell;
	} facet2;
	
	typedef struct _sample2 {
		vec2i p;
		facet2 *facets[2];
		FLOAT64 weights[2];
		char num;
	} sample2;

	// Node structure
	typedef struct _node2 {
		vec2i p;			// Position
		uint dx;			// Minimum dx
		FLOAT64 solidLS;	// Solid levelset
		FLOAT64 fluidLS;	// Fluid levelset
		vec2d velocity;
		vec2d gradient;
		// Terms for the extrapolation
		std::list<_node2 *> nodes;
		std::vector<sample2 *> samples[DIM];
		std::vector<FLOAT64> weights[DIM];
		bool mapped;
		bool edge[DIM];
		uint index;
	} node2;
	
	// Cell and facet array
	std::vector<cell2 *> cells;
	std::vector<facet2 *> facets;
	std::map<uint64,node2 *> nodes;
	std::map<uint64,sample2 *> samples[DIM];
	
	// Solver info
	FLOAT64 dt;
	FLOAT64 extrapolate_dist;
	FLOAT64 free_e;
	int surf_order;
	bool variation;
	
	// Matrix
	SparseMatrix<FLOAT64> LhsMatrix;	// Left hand side matrix
	SparseMatrix<FLOAT64> RhsMatrix;	// Right hand side matrix
	
	const levelset2 *solid;
	const levelset2 *fluid;
public:
	// Octree info
	octree2 octree;
	uint maxdepth;
private:
	uint64 computeIndex( vec2i p );
	void drawLevelset( int kind ) const;
	void resample( int kind, int type ); // type 0: node -> facet, 1: facet -> node
	vec2d sampleVector( vec2d p, char kind ) const;
	void extrapolate( char kind );
	void buildMatrix();
};

#endif