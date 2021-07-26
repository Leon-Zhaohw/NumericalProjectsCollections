/*
 *	octree3.cpp
 *
 *	Created by Ryoichi Ando on 5/30/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <tr1/unordered_map>
#include <algorithm>
#include "macros.h"
#include "util3.h"
#include "util.h"
#include "octree3.h"
#include "levelset3.h"
#include "opengl.h"
using namespace std;

class defaultLevelset3 : public levelset3 {
public:
	defaultLevelset3( uint maxdepth ) : maxdepth(maxdepth) {}
	virtual FLOAT64 evalLevelset(vec3d p) const {
		return 1.0/powf(2.0,maxdepth);
	}
	uint maxdepth;
};

octree3::octree3() {
	root = NULL;
	dontEnfortceWB = false;
	clearData();
}

octree3::octree3( const octree3 &octree ) {
	root = NULL;
	*this = octree;
}

octree3::~octree3() {
	clearData();
}

void octree3::operator=( const octree3 &octree ) {
	clearData();
	if( octree.root ) {
		maxdepth = octree.maxdepth;
		resolution = octree.resolution;
		root = new leaf3;
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

void octree3::copy( leaf3 *src, leaf3 *dest ) {
	*dest = *src;
	FOR_EACH(2,2,2) {
		if( src->children[i][j][k] ) {
			dest->children[i][j][k] = new leaf3;
			copy(src->children[i][j][k],dest->children[i][j][k]);
		}
	} END_FOR
}

bool octree3::buildOctree( const levelset3 *hint, const std::vector<sphere3 *> &spheres, uint maxdepth ) {
	// Make default levelset
	defaultLevelset3 defaultLevelset(maxdepth-1);
	if( ! hint && spheres.empty()) hint = &defaultLevelset;
	
	// Clear the octree first
	if( ! clearData()) return false;
	
	tick(); dump( ">>> Building octree started...\n" );
	
	// Set the maximal depth
	this->maxdepth = maxdepth;
	resolution = powf(2,maxdepth+2);
	
	// Allocate root leaf
	root = allocLeaf(vec3i(resolution/2,resolution/2,resolution/2),0,vec3i(0,0,0));
	
	// Subdivide this root leaf...
	tick(); dump( "Subdividing octree..." );
	subdivide(root,hint,spheres,maxdepth);
	dump( "Done. Took %s.\n", stock("octree_subdivision_hinted"));
	
	// Enforce weak balance
	enforceWeakBalance();
	
	// Build terminal array
	tick(); dump( "Building terminal list..." );
	uint index = 0;
	countNumTerminal(root,index);
	terminals.resize(index);
	index = 0;
	buildArray(root,index);
	dump( "Done. Found %d terminals. Took %s.\n", index, stock());
	writeNumber("octree_terminal_num", index);
	
	// Build corner nodes and its references
	tick(); dump( "Building nodes list..." );
	buildNodes();
	dump( "Done. Took %s.\n", stock("octree_node_list_hinted"));
	writeNumber("octree_node_num", nodes.size());
	
	dump( "<<< Octree done. Took %s.\n", stock("octree"));
	return true;
}

bool octree3::buildOctree( const levelset3 *hint, uint maxdepth ) {
	const std::vector<sphere3 *> spheres;
	return buildOctree(hint,spheres,maxdepth);
}

bool octree3::buildOctree( const std::vector<sphere3 *> &spheres, uint maxdepth ) {
	return buildOctree(NULL,spheres,maxdepth);
}

void octree3::dontEnforceWeakBalance() {
	dontEnfortceWB = true;
}

int octree3::hitTest( vec3d p ) const {
	for( uint dim=0; dim<DIM; dim++ ) p[dim] = fmin(1.0,fmax(0.0,p[dim]));
	p = p * resolution;
	std::vector<uint> array;
	hitTest(p,0.0,array,root);
	if( array.empty() ) return -1;
	return *array.begin();
}

bool octree3::hitTest( vec3d p, FLOAT64 r, std::vector<uint> &array, leaf3 *leaf ) const {
	if( ! leaf ) {
		p = p * resolution;
		leaf = root;
	}
	vec3d center = vec3d(leaf->center[0],leaf->center[1],leaf->center[2]);
	vec3d p0 = center - leaf->dx*vec3d(0.5,0.5,0.5);
	vec3d p1 = center + leaf->dx*vec3d(0.5,0.5,0.5);
	bool hit = box(p,p0,p1) <= resolution*r;
	if( hit ) {
		if( leaf->subdivided ) {
			FOR_EACH(2,2,2) {
				hitTest(p,r,array,leaf->children[i][j][k]);
			} END_FOR
		} else {
			array.push_back(leaf->index);
		}
	}
	return array.size();
}

octree3::leaf3* octree3::allocLeaf( vec3i center, uint depth, vec3i position ) {
	leaf3 *leaf = new leaf3;
	leaf->center = center;
	leaf->position = position;
	leaf->depth = depth;
	leaf->dx = resolution/powf(2,depth);
	leaf->subdivided = false;
	FOR_EACH(2,2,2) {
		leaf->children[i][j][k] = NULL;
		leaf->corners[i][j][k] = 0;
	} END_FOR
	return leaf;
}

static vec3i computeCenterPos( vec3i center, uint dx, uint i, uint j, uint k ) {
	return center-0.25*dx*vec3i(1,1,1)+0.5*dx*vec3i(i,j,k);
}

static vec3i computeCornerPos( vec3i center, uint dx, uint i, uint j, uint k ) {
	return center-0.5*dx*vec3i(1,1,1)+dx*vec3i(i,j,k);
}

bool octree3::checkSubdivision( vec3i pos, uint dx, const levelset3 *hint, int threshold, int depth, uint max_nest ) const {
	if( ! max_nest ) return false;
	if( hint->evalLevelset(vec3d(pos)/(FLOAT64)resolution) < powf(0.5,depth) ) {
		return true;
	}
	if( dx > threshold ) {
		// Compute the center position of this children
		FOR_EACH(2,2,2) {
			// Compute the levelset at this position
			vec3i sub_pos = computeCenterPos(pos,dx,i,j,k);
			if( checkSubdivision(sub_pos,dx/2,hint,threshold,depth,max_nest-1)) return true;
		} END_FOR
	}
	return false;
}

void octree3::subdivide( leaf3 *leaf, const levelset3 *hint, const std::vector<sphere3 *> &spheres, uint maxdepth ) {
	bool doSubdivision = false;
	
	// Compute the center position of this children
	FLOAT64 dx = leaf->dx/(FLOAT64)resolution;
	
	// See this octree contains large small particles
	for( uint n=0; n<spheres.size(); n++ ) {
		const sphere3 &sphere = *spheres[n];
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
			PARALLEL_FOR FOR_EACH(2,2,2) {
				// Compute the center position for this children
				vec3i center = computeCenterPos(leaf->center,leaf->dx,i,j,k);
				leaf3 *child = allocLeaf(center,depth,vec3i(i,j,k));
				// Make a new sphere array for this child
				std::vector<sphere3 *> child_spheres;
				for( uint n=0; n<spheres.size(); n++ ) {
					const sphere3 &sphere = *spheres[n];
					vec3d child_pos = vec3d(center)/(FLOAT64)resolution;
					FLOAT64 child_dx = child->dx/(FLOAT64)resolution;
					if( (child_pos-sphere.p).len2() < sqr(child_dx) ) {
						child_spheres.push_back(spheres[n]);
					}
				}				
				leaf->children[i][j][k] = child;
				subdivide(child,hint,child_spheres,maxdepth);
			} END_FOR
		}
	}
}

void octree3::enforceWeakBalance() {
	if( dontEnfortceWB ) return;
	
	tick(); dump( "Enforcing Weak Balance condition..." );
	uint itnum = 0;
	uint first_num = 0;
	uint last_num = 0;
	uint subdiv_num = 0;
	// Repeat while more than 2-level T-junction exists
	while(true) {
		// Collect terminals
		std::tr1::unordered_map<uint,leaf3 *> terminals_collapse;
		
		// Build terminal array
		uint index = 0;
		countNumTerminal(root,index);
		terminals.resize(index);
		if( ! first_num ) first_num = index;
		last_num = index;
		index = 0;
		buildArray(root,index);
		
		// For each terminal
		for( uint n=0; n<terminals.size(); n++ ) {
			// Look for neighbors
			FOR_EACH(2,2,2) {
				vec3d p = vec3d(terminals[n]->center+terminals[n]->dx*vec3d(2*i-1,2*j-1,2*k-1))/resolution;
				int neigh = hitTest(p);
				if( neigh >= 0 ) {
					if( terminals[neigh]->dx > 2*terminals[n]->dx && terminals_collapse.find(terminals[neigh]->index) == terminals_collapse.end()) {
						terminals_collapse[terminals[neigh]->index] = terminals[neigh];
						subdiv_num ++;
					}
				}
			} END_FOR
		}
		if( terminals_collapse.empty() ) break;
		
		// Collapse
		std::tr1::unordered_map<uint,leaf3 *>::iterator it;
		for( it=terminals_collapse.begin(); it!=terminals_collapse.end(); it++ ) {
			leaf3 *leaf = it->second;
			uint depth = leaf->depth+1;
			if( depth <= maxdepth ) {
				leaf->subdivided = true;
				FOR_EACH(2,2,2) {
					// Compute the center position for this children
					vec3i center = computeCenterPos(leaf->center,leaf->dx,i,j,k);
					leaf3 *child = allocLeaf(center,depth,vec3i(i,j,k));
					leaf->children[i][j][k] = child;
				} END_FOR
			}
		}
		itnum ++;
	}
	dump( "Done. Looped %d times. %d terminals are subdivided and extra %d terminals are generated. Took %s.\n", itnum, subdiv_num, last_num-first_num, stock("octree_weakbalance"));
	writeNumber("octree_weakbalace_generated", last_num-first_num);
}

void octree3::countNumTerminal( leaf3 *leaf, uint &count ) {
	if( leaf->subdivided ) {
		FOR_EACH(2,2,2) {
			countNumTerminal(leaf->children[i][j][k],count);
		} END_FOR
	} else {
		count++;
	}
}

void octree3::buildArray( leaf3 *leaf, uint &index ) {
	if( leaf->subdivided ) {
		FOR_EACH(2,2,2) {
			buildArray(leaf->children[i][j][k],index);
		} END_FOR
	} else {
		terminals[index] = leaf;
		leaf->index = index;
		index ++;
	}
}

void octree3::buildNodes() {
	std::tr1::unordered_map<uint64,uint64> nodeDictionary;
	uint64 index = 0;
	for( uint n=0; n<terminals.size(); n++ ) {
		leaf3 *leaf = terminals[n];
		FOR_EACH(2,2,2) {
			// Compute the center position for this children
			vec3i corner = computeCornerPos(leaf->center,leaf->dx,i,j,k);
			uint64 idx = computeCornerIndex(corner);
			if( nodeDictionary.find(idx) == nodeDictionary.end() ) {
				nodeDictionary[idx] = index;
				leaf->corners[i][j][k] = index;
				index ++;
			} else {
				leaf->corners[i][j][k] = nodeDictionary[idx];
			}
		} END_FOR
	}
	nodes.resize(index);
	for( uint n=0; n<terminals.size(); n++ ) {
		leaf3 *leaf = terminals[n];
		FOR_EACH(2,2,2) {
			uint64 index = leaf->corners[i][j][k];
			nodes[index] = vec3d(computeCornerPos(leaf->center,leaf->dx,i,j,k))/resolution;
		} END_FOR
	}
}

uint64 octree3::computeCornerIndex( vec3i p ) const {
	uint64 R = resolution;
	return p[0]+p[1]*R+p[2]*R*R;
}

bool octree3::clearData() {
	if( root ) {
		releaseChilren(root);
		root = NULL;
	}
	maxdepth = 1;
	terminals.clear();
	nodes.clear();
	return true;
}

bool octree3::releaseChilren( leaf3 *leaf ) {
	if( ! leaf ) return false;
	// Make sure we release all the chilren first
	FOR_EACH(2,2,2) {
		if( leaf->children[i][j][k] ) {
			releaseChilren(leaf->children[i][j][k]);
			leaf->children[i][j][k] = NULL;
		}
	} END_FOR
	
	// After that we deallocate this structure
	delete leaf;
	return true;
}

FLOAT64 octree3::box( vec3d p, vec3d p0, vec3d p1 ) const {
	FLOAT64 sd = -9999.0;
	sd = fmax(sd,p0[0]-p[0]);
	sd = fmax(sd,p0[1]-p[1]);
	sd = fmax(sd,p0[2]-p[2]);
	sd = fmax(sd,p[0]-p1[0]);
	sd = fmax(sd,p[1]-p1[1]);
	sd = fmax(sd,p[2]-p1[2]);
	return sd;
}

void octree3::drawOctree() const {
	glPushMatrix();
	glScaled(1.0/resolution,1.0/resolution,1.0/resolution);
	if( root ) drawOctree(root);
	glPopMatrix();
}

void octree3::drawOctree( const leaf3 *leaf ) const {
	glColor4d(0.5,0.5,0.5,0.5);
	FOR_EACH(2,2,2) {
		if( leaf->subdivided ) {
			drawOctree(leaf->children[i][j][k]);
		} else {
			int dx = leaf->dx;
			vec3i center = leaf->center;
			// GL box wire box here...
			glPushMatrix();
			glTranslatef(center[0],center[1],center[2]);
			glutWireCube(dx);
			glPopMatrix();
		}
	} END_FOR
}

void octree3::writeObj( const char *path ) const {
	FILE *fp = fopen(path,"w");
	if( ! fp ) return;
	uint index=0;
	for( uint n=0; n<terminals.size(); n++ ) {
		leaf3* leaf = terminals[n];
		vec3d p = vec3d(terminals[n]->center)/resolution;
		FLOAT64 dx = leaf->dx/(FLOAT64)resolution;
		for( uint dim=0; dim<DIM; dim++ ) {
			vec3d vec0, vec1;
			if( dim == 0 ) {
				vec0 = vec3d(0,1,0);
				vec1 = vec3d(0,0,1);
			} else if( dim == 1 ) {
				vec0 = vec3d(1,0,0);
				vec1 = vec3d(0,0,1);
			} else if( dim == 2 ) {
				vec0 = vec3d(1,0,0);
				vec1 = vec3d(0,1,0);
			}
			for( int dir=-1; dir<=1; dir+=2 ) {
				vec3d q[4] = { -0.5*vec0-0.5*vec1, 0.5*vec0-0.5*vec1, 0.5*vec0+0.5*vec1, -0.5*vec0+0.5*vec1 };
				for( uint m=0; m<4; m++ ) {
					vec3d corner = p+dir*0.5*dx*vec3d(dim==0,dim==1,dim==2)+dx*q[m];
					fprintf(fp,"v %f %f %f\n", corner[0], corner[1], corner[2] );
					index ++;
				}
				fprintf(fp,"f %d %d %d %d\n", index-3, index-2, index-1, index );
			}
		}
	}
	fclose(fp);
}

