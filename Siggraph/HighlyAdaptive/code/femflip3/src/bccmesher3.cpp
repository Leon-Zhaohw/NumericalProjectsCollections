/*
 *	bccmesher3.cpp
 *
 *	Created by Ryoichi Ando on 5/29/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <math.h>
#include <stdlib.h>
#include "bccmesher3.h"
#include "levelset3.h"

void bccmesher3::generateElements( std::vector<vec3d> &nodes, std::vector<std::vector<uint> > &elements, const levelset3 *hint, uint gn ) {
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
	
	// Deallocate BCC tets
	bcc.freeTets();
	
	// Set no hash generation
	setGenHash(false);
}

int bccmesher3::hitElements(const vec3d p) const {
	return bcc.hitTest(p);
}