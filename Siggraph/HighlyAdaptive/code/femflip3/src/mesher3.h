/*
 *	mesher3.h
 *	
 *	Created by Ryoichi Ando on 2/1/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "vec3.h"
#include "array3.h"
#include <vector>

#ifndef _MESHER3_H
#define _MESHER3_H

class levelset3;
class mesher3 {
public:
	enum { BARYCENTRIC = 0, CIRCUMCENTRIC = 1 };
	enum { NODAL = 0, ELEMENT = 1 };
	
	mesher3();
	void setGenHash( bool state );
	void setGenFacet( bool state );
	void setCenterType( int type );						// 0 for barycentric, 1 for circumcenter
	virtual void generateMesh( uint gn, const levelset3 *hint=NULL );
	virtual int	 hitElements(vec3d p) const;		// Returns hit element index (-1 if not hit)
	virtual void updateConnection();
	
	uint gn;											// Base resolution
	std::vector<vec3d> nodes;							// Node array
	std::vector<std::vector<uint> > elements;			// Elements array
	std::vector<vec3d> centers;							// Element centers
	
	typedef struct {
		FLOAT64 m[NUM_VERT][NUM_VERT];
	} shapeM;
	std::vector<shapeM> matrix;							// Element shape function matrix
	
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
	array3<std::vector<uint> > hash;					// Element spatial hash
	
	// utility function
	FLOAT64 getLevelsetVolume( uint elmid, FLOAT64 isoval[4] );
	virtual void drawMesh() const;
	void writeObj( const char *path, int type=0 ) const;
protected:
	int centerType;
	bool genHash;
	bool genFacet;
	bool hitElement( uint index, vec3d p ) const;
	void computeCenterPosition( const std::vector<vec3d> &nodes, const std::vector<std::vector<uint> > &elements, std::vector<vec3d> &centers, int centerType );
	virtual void generateElements( std::vector<vec3d> &nodes, std::vector<std::vector<uint> > &elements, const levelset3 *hint, uint gn );
	void cleanup();
};

#endif