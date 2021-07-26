/*
 *	bccmesher2.cpp
 *
 *	Created by Ryoichi Ando on 5/28/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <math.h>
#include "bccmesher2.h"
#include "levelset2.h"

void bccmesher2::generateElements( std::vector<vec2d> &nodes, std::vector<std::vector<uint> > &elements, const levelset2 *hint, uint gn ) {
	// Now, let's consider 'gn' as a maximum resolution. Therefore, maxdepth will be approx log(2,gn)
	uint maxdepth = log2(gn)+1;

	// Now, do it
	if( hint && hint->getSpheres()) {
		octree.buildOctree(*hint->getSpheres(),maxdepth);
	} else {
		octree.buildOctree(hint,maxdepth);
	}
	bcc.buildBCC(octree);
	
	// Set nodes
	this->nodes = bcc.getNodes();
	
	// Set elements
	this->elements = bcc.getElements();
	
	// Free bcc tets data
	bcc.freeTets();
	
	// Set no hash generation
	setGenHash(false);
}

int bccmesher2::hitElements(const vec2d p) const {
	return bcc.hitTest(p);
}