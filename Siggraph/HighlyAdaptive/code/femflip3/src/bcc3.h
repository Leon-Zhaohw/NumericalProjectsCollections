/*
 *	bcc3.h
 *
 *	Created by Ryoichi Ando on 5/28/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "octree3.h"
#include <tr1/unordered_map>
#ifndef _BCC3_H
#define _BCC3_H

#define CENTERS_DEF		0		// Whether or not to mark the BCC center positions

class levelset3;
class bcc3 {
public:
	bool buildBCC( const octree3 &octree );
	const std::vector<vec3d> &getNodes() const { return cleanNodes; }
	const std::vector<std::vector<uint> > &getElements() const { return cleanElements; }
	void freeTets();
	
	virtual void drawElements();
	int hitTest(vec3d p) const;
	virtual void writeObj( const char *path ) const;
	octree3 octree;
	std::tr1::unordered_map<uint64,vec3i> nodes;
	std::tr1::unordered_map<uint64,bool> centers;
	std::vector<std::vector<uint64> > elements;
	std::vector<vec3d> cleanNodes;
#if CENTERS_DEF
	std::vector<bool> cleanCenters;
#endif
	std::vector<std::vector<uint> > cleanElements;
	std::vector<std::vector<uint> > hash;
	
	typedef struct {
		FLOAT64 m[NUM_VERT][NUM_VERT];
	} shapeM;
	std::vector<shapeM> matrix;

protected:
	bool clearData();
	void buildCornerArray( const std::vector<octree3::leaf3 *> &terminals );
	void buildElements( const std::vector<octree3::leaf3 *> &terminals, const octree3 &octree, std::vector<std::vector<uint> > &hash );
	void cleanupNodeAndElements();
	void precomputeShapeMatrix();
	void buildTets( octree3::leaf3 *leaf, const octree3 &octree, std::vector<std::vector<uint> > &hash );
private:
	uint resolution;
	uint64 computeCornerIndex( vec3i p ) const;
};

#endif