/*
 *	octree3.h
 *
 *	Created by Ryoichi Ando on 5/30/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "vec3.h"
#include <vector>
#include <list>
#ifndef _OCTREE3_H
#define _OCTREE3_H

class levelset3;
class octree3 {
public:
	// Octree cell structure
	typedef struct _leaf3 {
		_leaf3 *children[2][2][2];
		vec3i position;
		bool subdivided;
		vec3i center;
		uint depth;
		int dx;
		uint64 corners[2][2][2];	// index reference to the corner indices
		uint index;
	} leaf3;
	
	typedef struct {
		vec3d p;
		FLOAT64 r;
	} sphere3;
	
	octree3();
	octree3( const octree3 &octree );
	void operator=( const octree3 &octree );
	~octree3();
	
	bool buildOctree( const levelset3 *hint, uint maxdepth );
	bool buildOctree( const std::vector<sphere3 *> &spheres, uint maxdepth );
	bool buildOctree( const levelset3 *hint, const std::vector<sphere3 *> &spheres, uint maxdepth );
	void dontEnforceWeakBalance();
	
	int hitTest( vec3d p ) const;
	bool hitTest( vec3d p, FLOAT64 r, std::vector<uint> &array, leaf3 *leaf=NULL ) const;
	void drawOctree() const;
	void writeObj( const char *path ) const;
	
	leaf3 *root;
	std::vector<leaf3 *> terminals; // Octree terminal list
	std::vector<vec3d> nodes;	// Octree corner list
	
	uint maxdepth;
	uint resolution;
	bool dontEnfortceWB;
	void subdivide( leaf3 *leaf, const levelset3 *hint, const std::vector<sphere3 *> &spheres, uint maxdepth );
	void enforceWeakBalance();
private:
	void countNumTerminal( leaf3 *leaf, uint &index );
	void buildArray( leaf3 *leaf, uint &index );
	void buildNodes();
	bool clearData();
	bool releaseChilren( leaf3 *leaf );
	leaf3* allocLeaf( vec3i center, uint depth, vec3i position );
	void drawOctree( const leaf3 *leaf ) const;
	FLOAT64 box( vec3d p, vec3d p0, vec3d p1 ) const;
	void copy( leaf3 *src, leaf3 *dest );
	uint64 computeCornerIndex( vec3i p ) const;
	bool checkSubdivision( vec3i pos, uint dx, const levelset3 *hint, int threshold, int depth, uint max_nest ) const;
};

#endif