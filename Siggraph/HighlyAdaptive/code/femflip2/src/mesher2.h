/*
 *	mesher2.h
 *	
 *	Created by Ryoichi Ando on 12/26/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "vec2.h"
#include "array2.h"
#include <vector>

#ifndef _MESHER2_H
#define _MESHER2_H

class levelset2;
class mesher2 {
public:
	enum { BARYCENTRIC = 0, CIRCUMCENTRIC = 1 };
	enum { NODAL = 0, ELEMENT = 1 };
	
	mesher2();
	void setGenHash( bool state );
	void setGenFacet( bool state );
	void setCenterType( int type );						// 0 for barycentric, 1 for circumcenter
	virtual void generateMesh( uint gn, const levelset2 *hint=NULL );
	virtual int	 hitElements(vec2d p) const;		// Returns hit element index (-1 if not hit)

	uint gn;											// Base resolution
	std::vector<vec2d> nodes;							// Grid node array
	std::vector<std::vector<uint> > elements;			// Elements array
	std::vector<vec2d> centers;							// Element center positions
	
	typedef struct {
		FLOAT64 m[NUM_VERT][NUM_VERT];
	} shapeM;
	std::vector<shapeM> matrix;							// Element shape function matrix
	
	// Information below may not be available
	std::vector<std::vector<uint> > node_elements;		// Node to element array. This can be interpreted as Voronoi elements
	std::vector<std::vector<uint> > facets;				// Elements facets
	std::vector<FLOAT64>			facetArea;			// Facet area
	std::vector<std::vector<uint> > facet_elements;		// Facet to element array
	std::vector<std::vector<uint> > edges;				// Edges
	std::vector<std::vector<uint> > element_facets;		// Element to facet array
	std::vector<std::vector<uint> > element_edges;		// Element to edges
	std::vector<std::vector<uint> > element_elements;	// Adjacent element - element array
	std::vector<std::vector<uint> > node2node;			// Node to node connection array
    std::vector<FLOAT64> volumes;						// Element volumes
	array2<std::vector<uint> > hash;					// Element spatial hash
	
	virtual void drawMesh() const;
	virtual void drawCenters() const;
	virtual void drawElement2Facet() const;
	virtual void drawLevelset( const std::vector<FLOAT64> &levelsets ) const;
	virtual void drawScalar( const std::vector<FLOAT64> &scalar, uint type=NODAL, const std::vector<bool> &mask=std::vector<bool>(), FLOAT64 minv=0.0, FLOAT64 maxv=0.0 ) const;
	virtual void drawVector( const std::vector<vec2d> &vecs, uint type=NODAL, FLOAT64 scale=1.0 ) const;
	virtual void updateConnection();
	virtual void write_matlab( const char *name, const std::vector<FLOAT64> &scalar );
protected:
	int centerType;
	bool genHash;
	bool genFacet;
	bool hitElement( uint index, vec2d p ) const;
	void computeCenterPosition( const std::vector<vec2d> &nodes, const std::vector<std::vector<uint> > &elements, std::vector<vec2d> &centers, int centerType );
	virtual void generateElements( std::vector<vec2d> &nodes, std::vector<std::vector<uint> > &elements, const levelset2 *hint, uint gn );
	void cleanup();
};

#endif