/*
 *	octree2.cpp
 *
 *	Created by Ryoichi Ando on 5/30/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <tr1/unordered_map>
#include <algorithm>
#include "octree2.h"
#include "util.h"
#include "levelset2.h"
#include "opengl.h"
using namespace std;

class defaultLevelset2 : public levelset2 {
public:
	defaultLevelset2( uint maxdepth ) : maxdepth(maxdepth) {}
	virtual FLOAT64 evalLevelset(vec2d p) const {
		return 1.0/powf(2.0,maxdepth);
	}
	uint maxdepth;
};

octree2::octree2() {
	root = NULL;
	clearData();
}

octree2::octree2( const octree2 &octree ) {
	root = NULL;
	*this = octree;
}

octree2::~octree2() {
	clearData();
}

void octree2::operator=( const octree2 &octree ) {
	clearData();
	if( octree.root ) {
		maxdepth = octree.maxdepth;
		resolution = octree.resolution;
		root = new leaf2;
		copy(octree.root,root);
		terminals.resize(octree.terminals.size());
		
		// Build terminal array
		uint index = 0;
		countNumTerminal(root,index);
		terminals.resize(index);
		index = 0;
		buildArray(root,index);
		
		// Copy nodes
		nodes = octree.nodes;
	}
}

void octree2::copy( leaf2 *src, leaf2 *dest ) {
	*dest = *src;
	FOR_EACH(2,2) {
		if( src->children[i][j] ) {
			dest->children[i][j] = new leaf2;
			copy(src->children[i][j],dest->children[i][j]);
		}
	} END_FOR
}

bool octree2::buildOctree( const levelset2 *hint, const std::vector<sphere2 *> &spheres, uint maxdepth ) {
	// Make default levelset
	defaultLevelset2 defaultLevelset(maxdepth-1);
	if( ! hint && spheres.empty() ) hint = &defaultLevelset;
	
	// Clear the octree first
	if( ! clearData()) return false;
	
	// Set the maximal depth
	this->maxdepth = maxdepth;
	resolution = powf(2,maxdepth+2);
	
	// Allocate root leaf
	root = allocLeaf(vec2i(resolution/2,resolution/2),0,vec2i(0,0));
	
	// Subdivide this root leaf...
	subdivide(root,hint,spheres,maxdepth);
	
	// Enforce Weak Balance
	enforceWeakBalance();
	
	// Build terminal array
	uint index = 0;
	countNumTerminal(root,index);
	terminals.resize(index);
	index = 0;
	buildArray(root,index);
	
	// Build corner nodes and its references
	buildNodes();
	return true;
}

bool octree2::buildOctree( const levelset2 *hint, uint maxdepth ) {
	const std::vector<sphere2 *> spheres;
	return buildOctree(hint,spheres,maxdepth);
}

bool octree2::buildOctree( const std::vector<sphere2 *> &spheres, uint maxdepth ) {
	return buildOctree(NULL,spheres,maxdepth);
}

int octree2::hitTest( vec2d p ) const {
	for( uint dim=0; dim<DIM; dim++ ) p[dim] = fmin(1.0,fmax(0.0,p[dim]));
	p = p * resolution;
	std::vector<uint> array;
	hitTest(p,0.0,array,root);
	if( array.empty() ) return -1;
	return *array.begin();
}

bool octree2::hitTest( vec2d p, FLOAT64 r, std::vector<uint> &array, leaf2 *leaf ) const {
	if( ! leaf ) {
		p = p * resolution;
		leaf = root;
	}
	vec2d center = vec2d(leaf->center[0],leaf->center[1]);
	vec2d p0 = center - leaf->dx*vec2d(0.5,0.5);
	vec2d p1 = center + leaf->dx*vec2d(0.5,0.5);
	bool hit = box(p,p0,p1) <= resolution*r;
	if( hit ) {
		if( leaf->subdivided ) {
			FOR_EACH(2,2) {
				hitTest(p,r,array,leaf->children[i][j]);
			} END_FOR
		} else {
			array.push_back(leaf->index);
		}
	}
	return array.size();
}

octree2::leaf2* octree2::allocLeaf( vec2i center, uint depth, vec2i position ) {
	leaf2 *leaf = new leaf2;
	leaf->center = center;
	leaf->position = position;
	leaf->depth = depth;
	leaf->dx = resolution/powf(2,depth);
	leaf->subdivided = false;
	FOR_EACH(2,2) {
		leaf->children[i][j] = NULL;
		leaf->corners[i][j] = 0;
	} END_FOR
	return leaf;
}

static vec2i computeCenterPos( vec2i center, uint dx, uint i, uint j ) {
	return center-0.25*dx*vec2i(1,1)+0.5*dx*vec2i(i,j);
}

static vec2i computeCornerPos( vec2i center, uint dx, uint i, uint j ) {
	return center-0.5*dx*vec2i(1,1)+dx*vec2i(i,j);
}

bool octree2::checkSubdivision( vec2i pos, uint dx, const levelset2 *hint, int threshold, int depth, uint max_nest ) const {
	if( ! max_nest ) return false;
	if( hint->evalLevelset(vec2d(pos)/(FLOAT64)resolution) < powf(0.5,depth) ) {
		return true;
	}
	if( dx > threshold ) {
		// Compute the center position of this children
		FOR_EACH(2,2) {
			// Compute the levelset at this position
			vec2i sub_pos = computeCenterPos(pos,dx,i,j);
			if( checkSubdivision(sub_pos,dx/2,hint,threshold,depth,max_nest-1)) return true;
		} END_FOR
	}
	return false;
}

void octree2::subdivide( leaf2 *leaf, const levelset2 *hint, const std::vector<sphere2 *> &spheres, uint maxdepth ) {
	bool doSubdivision = false;
	
	// Compute the center position of this children
	FLOAT64 dx = leaf->dx/(FLOAT64)resolution;
	
	// See this octree contains large small particles
	for( uint n=0; n<spheres.size(); n++ ) {
		const sphere2 &sphere = *spheres[n];
		if( sphere.r <= 0.5*dx ) {
			doSubdivision = true;
			break;
		}
	}

	// See hint indicates a subdivision
	if( hint && ! doSubdivision ) {
		int threshold = 8;
		doSubdivision = checkSubdivision(leaf->center,leaf->dx,hint,threshold,leaf->depth,5);
	}
	
	// If to subdivide, do it
	if( doSubdivision ) {
		uint depth = leaf->depth+1;
		if( depth <= maxdepth ) {
			leaf->subdivided = true;
			FOR_EACH(2,2) {
				// Compute the center position for this children
				vec2i center = computeCenterPos(leaf->center,leaf->dx,i,j);
				leaf2 *child = allocLeaf(center,depth,vec2i(i,j));
				
				// Make a new sphere array for this child
				std::vector<sphere2 *> child_spheres;
				for( uint n=0; n<spheres.size(); n++ ) {
					const sphere2 &sphere = *spheres[n];
					vec2d child_pos = vec2d(center)/(FLOAT64)resolution;
					FLOAT64 child_dx = child->dx/(FLOAT64)resolution;
					if( (child_pos-sphere.p).len2() < sqr(child_dx) ) {
						child_spheres.push_back(spheres[n]);
					}
				}

				leaf->children[i][j] = child;
				subdivide(child,hint,child_spheres,maxdepth);
			} END_FOR
		} 
	}
}

void octree2::enforceWeakBalance() {
	// Repeat while more than 2-level T-junction exists
	while(true) {
		// Collect terminals
		std::tr1::unordered_map<uint,leaf2 *> terminals_collapse;
		
		// Build terminal array
		uint index = 0;
		countNumTerminal(root,index);
		terminals.resize(index);
		index = 0;
		buildArray(root,index);
		
		// For each terminal
		for( uint n=0; n<terminals.size(); n++ ) {
			// Look for neighbors
			FOR_EACH(2,2) {
				vec2d p = vec2d(terminals[n]->center+terminals[n]->dx*vec2d(2*i-1,2*j-1))/resolution;
				int neigh = hitTest(p);
				if( neigh >= 0 ) {
					if( terminals[neigh]->dx > 2*terminals[n]->dx && terminals_collapse.find(terminals[neigh]->index) == terminals_collapse.end()) {
						terminals_collapse[terminals[neigh]->index] = terminals[neigh];
					}
				}
			} END_FOR
		}
		if( terminals_collapse.empty()) break;
		
		// Collapse
		std::tr1::unordered_map<uint,leaf2 *>::iterator it;
		for( it=terminals_collapse.begin(); it!=terminals_collapse.end(); it++ ) {
			leaf2 *leaf = it->second;
			uint depth = leaf->depth+1;
			if( depth <= maxdepth ) {
				leaf->subdivided = true;
				FOR_EACH(2,2) {
					// Compute the center position for this children
					vec2i center = computeCenterPos(leaf->center,leaf->dx,i,j);
					leaf2 *child = allocLeaf(center,depth,vec2i(i,j));
					leaf->children[i][j] = child;
				} END_FOR
			}
		}
	}
}

void octree2::countNumTerminal( leaf2 *leaf, uint &count ) {
	if( leaf->subdivided ) {
		FOR_EACH(2,2) {
			countNumTerminal(leaf->children[i][j],count);
		} END_FOR
	} else {
		count++;
	}
}

void octree2::buildArray( leaf2 *leaf, uint &index ) {
	if( leaf->subdivided ) {
		FOR_EACH(2,2) {
			buildArray(leaf->children[i][j],index);
		} END_FOR
	} else {
		terminals[index] = leaf;
		leaf->index = index;
		index ++;
	}
}

void octree2::buildNodes() {
	std::tr1::unordered_map<uint64,uint64> nodeDictionary;
	uint64 index = 0;
	for( uint n=0; n<terminals.size(); n++ ) {
		leaf2 *leaf = terminals[n];
		FOR_EACH(2,2) {
			// Compute the center position for this children
			vec2i corner = computeCornerPos(leaf->center,leaf->dx,i,j);
			uint64 idx = computeCornerIndex(corner);
			if( nodeDictionary.find(idx) == nodeDictionary.end() ) {
				nodeDictionary[idx] = index;
				leaf->corners[i][j] = index;
				index ++;
			} else {
				leaf->corners[i][j] = nodeDictionary[idx];
			}
		} END_FOR
	}
	nodes.resize(index);
	for( uint n=0; n<terminals.size(); n++ ) {
		leaf2 *leaf = terminals[n];
		FOR_EACH(2,2) {
			uint index = leaf->corners[i][j];
			nodes[index] = vec2d(computeCornerPos(leaf->center,leaf->dx,i,j))/resolution;
		} END_FOR
	}
}

uint64 octree2::computeCornerIndex( vec2i p ) const {
	uint64 R = resolution;
	return p[0]+p[1]*R;
}

bool octree2::clearData() {
	if( root ) {
		releaseChilren(root);
		root = NULL;
	}
	maxdepth = 1;
	terminals.clear();
	nodes.clear();
	return true;
}

bool octree2::releaseChilren( leaf2 *leaf ) {
	if( ! leaf ) return false;
	// Make sure we release all the chilren first
	FOR_EACH(2,2) {
		if( leaf->children[i][j] ) {
			releaseChilren(leaf->children[i][j]);
			leaf->children[i][j] = NULL;
		}
	} END_FOR
	
	// After that we deallocate this structure
	delete leaf;
	return true;
}

void octree2::drawOctree() const {
	glPushMatrix();
	glScaled(1.0/resolution,1.0/resolution,1.0);
	if( root ) drawOctree(root);
	glPopMatrix();
}

void octree2::drawOctree( const leaf2 *leaf ) const {
	FOR_EACH(2,2) {
		if( leaf->subdivided ) {
			drawOctree(leaf->children[i][j]);
		} else {
			int dx = leaf->dx;
			vec2i center = leaf->center;
			glBegin(GL_LINE_LOOP);
			glVertex2d(center[0]-0.5*dx,center[1]-0.5*dx);
			glVertex2d(center[0]-0.5*dx,center[1]+0.5*dx);
			glVertex2d(center[0]+0.5*dx,center[1]+0.5*dx);
			glVertex2d(center[0]+0.5*dx,center[1]-0.5*dx);
			glEnd();
		}
	} END_FOR
}

// Box levelset
FLOAT64 octree2::box( vec2d p, vec2d p0, vec2d p1 ) const {
	FLOAT64 sd = -9999.0;
	sd = fmax(sd,p0[0]-p[0]);
	sd = fmax(sd,p0[1]-p[1]);
	sd = fmax(sd,p[0]-p1[0]);
	sd = fmax(sd,p[1]-p1[1]);
	return sd;
}