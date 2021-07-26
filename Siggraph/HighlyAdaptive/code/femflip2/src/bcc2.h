/*
 *	bcc2.h
 *
 *	Created by Ryoichi Ando on 5/26/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "octree2.h"
#include <vector>
#include <tr1/unordered_map>
#ifndef _BCC2_H
#define _BCC2_H

#define CENTERS_DEF		0		// Whether or not to mark the BCC center positions

class bcc2 {
public:
	bool buildBCC( const octree2 &octree );
	const std::vector<vec2d> &getNodes() const { return cleanNodes; }
	const std::vector<std::vector<uint> > &getElements() const { return cleanElements; }
	void freeTets();
	void drawElements();
	int hitTest(vec2d p) const;
	octree2 octree;
	std::tr1::unordered_map<uint64,vec2i> nodes;
	std::tr1::unordered_map<uint64,bool> centers;
	std::vector<std::vector<uint64> > elements;
	std::vector<vec2d> cleanNodes;
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
	void buildCornerArray( const std::vector<octree2::leaf2 *> &terminals );
	void buildElements( const std::vector<octree2::leaf2 *> &terminals, const octree2 &octree, std::vector<std::vector<uint> > &hash );
	void cleanupNodeAndElements();
	void precomputeShapeMatrix();
private:
	uint resolution;
	uint64 computeCornerIndex( vec2i p ) const;
};

#endif