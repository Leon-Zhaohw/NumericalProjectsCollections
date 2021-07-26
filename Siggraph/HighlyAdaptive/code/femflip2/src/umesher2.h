/*
 *	umesher2.h
 *	
 *	Created by Ryoichi Ando on 1/16/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "mesher2.h"

#ifndef _UMESHER2_H
#define _UMESHER2_H

class umesher2 : public mesher2 {
public:
	enum { BIG, SUBDIV, DOWN, NEIGHBOR, UP };
	virtual void generateElements( std::vector<vec2d> &nodes, std::vector<std::vector<uint> > &elements, const levelset2 *hint, uint gn );
	virtual void drawMesh() const;
	virtual void updateElements( std::vector<uint> &changedElements );
	
	std::vector<uint>	elementMap;							// Access like elements_bank[elementMap[elements[n]]]
	std::vector<vec2d>  nodes_bank;							// Grid node array
	std::vector<std::vector<uint> > elements_bank;			// Elements array
protected:
	uint w;
	uint h;
	uint numBig;
	
	std::vector<std::vector<int> > local_indices;
	std::vector<std::vector<int> > local_neighbors;
	std::vector<std::vector<int> > local_elements;
	std::vector<std::vector<int> > local_status;
	std::vector<int>  status;		// numBig
	std::vector<char> directions;	// numBig
	std::vector<bool> subdivided;	// numBig
	std::vector<bool> enabled;		// num elements_bank
	std::vector<bool> nodes_enabled;	// per node
	void drawBigMesh(uint n) const;
};

#endif