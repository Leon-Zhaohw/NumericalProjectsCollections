/*
 *	octfluid3.h
 *
 *	Created by Ryoichi Ando on 6/29/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "fluid3.h"
#include "octree3.h"
#include <list>
#include <map>
#include <vector>
#include "pcgsolver/sparse_matrix.h"

#ifndef OCTFLUID3_H
#define OCTFLUID3_H

class levelset3;
class octfluid3 : public fluid3 {
public:
	octfluid3();
	
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
	virtual const char *getName() const { return "OCT"; }
	virtual const octree3 *getOctree() const { return &octree; }
protected:
	// Cell structure
	struct _facet3;
	struct _node3;
	typedef struct _cell3 {
		vec3i p;
		uint dx;
		FLOAT64 pressure;
		FLOAT64 divergence;
		FLOAT64 fluidLS;
		std::list<_facet3 *> facets;
		_node3 *nodes[2][2][2];
		// Used to extrapolate
		bool mapped;
		vec3d velocity;
		vec3d gradient;
		uint index;
	} cell3;
	
	// Facet structure
	typedef struct _facet3 {
		char dim; // X, Y, Z
		vec3i p;
		FLOAT64 velocity;
		FLOAT64 gradient;
		FLOAT64 fraction;
		bool mapped;
		bool edge;
		uint dx;
		uint index;
		std::list<cell3 *> cells;
		std::list<_node3 *> nodes;
		// Used for pressure solve
		FLOAT64 ds;
		cell3 *Lcell;
	} facet3;
	
	typedef struct _sample3 {
		vec3i p;
		facet3 *facets[4];
		FLOAT64 weights[4];
		char num;
	} sample3;
	
	// Node structure
	typedef struct _node3 {
		vec3i p;			// Position
		uint dx;			// Minimum dx
		FLOAT64 solidLS;	// Solid levelset
		FLOAT64 fluidLS;	// Fluid levelset
		vec3d velocity;
		vec3d gradient;
		// Terms for the extrapolation
		std::list<_node3 *> nodes;
		std::vector<sample3 *> samples[DIM];
		std::vector<FLOAT64> weights[DIM];
		bool mapped;
		bool edge[DIM];
		uint index;
	} node3;
	
	// Cell and facet array
	std::vector<cell3 *> cells;
	std::vector<facet3 *> facets;
	std::map<uint64,node3 *> nodes;
	std::map<uint64,sample3 *> samples[DIM];
	
	// Solver info
	FLOAT64 dt;
	FLOAT64 extrapolate_dist;
	FLOAT64 free_e;
	int surf_order;
	bool variation;
	FLOAT64 tension;
	
	// Matrix
	SparseMatrix<FLOAT64> LhsMatrix;	// Left hand side matrix
	SparseMatrix<FLOAT64> RhsMatrix;	// Right hand side matrix
	
	const levelset3 *solid;
	const levelset3 *fluid;
public:
	// Octree info
	octree3 octree;
	uint maxdepth;
private:
	uint64 computeIndex( vec3i p );
	void drawLevelset( int kind ) const;
	void resample( int kind, int type ); // type 0: node -> facet, 1: facet -> node
	vec3d sampleVector( vec3d p, char kind ) const;
	void extrapolate( char kind );
	void buildMatrix();
	void getTangentVectors( int dim, vec3i &vec0, vec3i &vec1 );
	bool computeWeights( const std::vector<sample3 *> &samples, std::vector<FLOAT64> &weights, int dim, vec3i p );
	void writeConnections( const char *path, uint dim, node3 *target_node );
};

#endif