/*
 *	octree2.h
 *
 *	Created by Ryoichi Ando on 5/30/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "vec2.h"
#include <vector>
#ifndef _OCTREE2_H
#define _OCTREE2_H

class levelset2;
class octree2 {
public:
	// Octree cell structure
	typedef struct _leaf2 {
		_leaf2 *children[2][2];
		vec2i position;
		bool subdivided;
		vec2i center;
		uint depth;
		int dx;
		uint64 corners[2][2];	// index reference to the corner indices
		uint index;
	} leaf2;
	
	typedef struct {
		vec2d p;
		FLOAT64 r;
	} sphere2;
	
	octree2();
	octree2( const octree2 &octree );
	void operator=( const octree2 &octree );
	~octree2();
	
	bool buildOctree( const levelset2 *hint, uint maxdepth );
	bool buildOctree( const std::vector<sphere2 *> &spheres, uint maxdepth );
	bool buildOctree( const levelset2 *hint, const std::vector<sphere2 *> &spheres, uint maxdepth );
	
	int hitTest( vec2d p ) const;
	bool hitTest( vec2d p, FLOAT64 r, std::vector<uint> &array, leaf2 *leaf=NULL ) const;
	void drawOctree() const;
	
	leaf2 *root;
	std::vector<leaf2 *> terminals; // Octree terminal list
	std::vector<vec2d> nodes;	// Octree corner list
	
	uint maxdepth;
	uint resolution;
	void subdivide( leaf2 *leaf, const levelset2 *hint, const std::vector<sphere2 *> &spheres, uint maxdepth );
	void enforceWeakBalance();
private:
	void countNumTerminal( leaf2 *leaf, uint &index );
	void buildArray( leaf2 *leaf, uint &index );
	void buildNodes();
	bool clearData();
	bool releaseChilren( leaf2 *leaf );
	leaf2* allocLeaf( vec2i center, uint depth, vec2i position );
	void drawOctree( const leaf2 *leaf ) const;
	// Box levelset
	FLOAT64 box( vec2d p, vec2d p0, vec2d p1 ) const;
	void copy( leaf2 *src, leaf2 *dest );
	uint64 computeCornerIndex( vec2i p ) const;
	bool checkSubdivision( vec2i pos, uint dx, const levelset2 *hint, int threshold, int depth, uint max_nest ) const;
};

#endif