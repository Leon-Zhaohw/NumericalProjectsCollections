/*
 *	bcc2.cpp
 *
 *	Created by Ryoichi Ando on 5/26/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "util2.h"
#include "bcc2.h"
#include "matutil.h"
#include "levelset2.h"
#include "opengl.h"
using namespace std;

bool bcc2::buildBCC( const octree2 &octree ) {
	// Copy octree
	this->octree = octree;
	
	// Get resolution
	resolution = octree.resolution;
	
	// Clear data
	if( ! clearData()) return false;
	
	// Build an array of corner nodes
	buildCornerArray(octree.terminals);
	
	// Build elements
	hash.resize(octree.terminals.size());
	buildElements(octree.terminals,octree,hash);
	
	// Cleanup node indices
	cleanupNodeAndElements();
	
	// Precompute shape matrix
	precomputeShapeMatrix();
	return true;
}

void bcc2::buildCornerArray( const std::vector<octree2::leaf2 *> &terminals ) {
	for( uint n=0; n<terminals.size(); n++ ) {
		int dx = terminals[n]->dx;
		vec2i center = terminals[n]->center;
		// Insert corners
		FOR_EACH(2,2) {
			vec2i p = center-0.5*dx*vec2i(1,1)+dx*vec2i(i,j);
			uint64 index = computeCornerIndex(p);
			if( nodes.find(index) == nodes.end() ) {
				nodes[index] = p;
#if CENTERS_DEF
				centers[index] = false;
#endif
			}
		} END_FOR

		// Insert a center
		uint64 index = computeCornerIndex(center);
		if( nodes.find(index) == nodes.end() ) {
			nodes[index] = center;
#if CENTERS_DEF
			centers[index] = true;
#endif
		}
	}
}

void bcc2::buildElements( const std::vector<octree2::leaf2 *> &terminals, const octree2 &octree, std::vector<std::vector<uint> > &hash ) {
	elements.reserve(terminals.size()*4);
	for( uint t=0; t<terminals.size(); t++ ) {
		int dx = terminals[t]->dx;
		vec2i center = terminals[t]->center;
		int idx = octree.hitTest(vec2d(center)/octree.resolution);
		
		// Build meshes here !
		uint64 v0 = computeCornerIndex(center);
		uint64 indices[] = {  
			computeCornerIndex(center+0.5*dx*vec2i(-1,-1)),
			computeCornerIndex(center+0.5*dx*vec2i(-1,1)),
			computeCornerIndex(center+0.5*dx*vec2i(1,1)),
			computeCornerIndex(center+0.5*dx*vec2i(1,-1)) };
		for( uint n=0; n<4; n++ ) {
			uint64 v1 = indices[n];
			uint64 v2 = indices[(n+1)%4];
			uint64 vt = computeCornerIndex(0.5*(nodes[v1]+nodes[v2]));
			if( nodes.find(vt) != nodes.end() ) {
				// Transition start
				vector<uint64> element1(3);
				vector<uint64> element2(3);
				element1[0] = v0;
				element1[1] = vt;
				element1[2] = v1;
				elements.push_back(element1);
				hash[idx].push_back(elements.size()-1);
				element2[0] = v0;
				element2[1] = vt;
				element2[2] = v2;
				elements.push_back(element2);
				hash[idx].push_back(elements.size()-1);
			} else {
				// Regular BCC tet
				vector<uint64> element(3);
				element[0] = v0;
				element[1] = v1;
				element[2] = v2;
				elements.push_back(element);
				hash[idx].push_back(elements.size()-1);
			}
		}
	}
}

void bcc2::cleanupNodeAndElements() {
	std::tr1::unordered_map<uint64,uint> remap;
	
	// Compute remap array
	uint size = nodes.size();
	cleanNodes.resize(size);
#if CENTERS_DEF
	cleanCenters.resize(size);
#endif
	uint index=0;
	std::tr1::unordered_map<uint64,vec2i>::iterator it = nodes.begin();
	while( it != nodes.end() ) {
		remap[(*it).first] = index;
		vec2i p = (*it).second;
		cleanNodes[index] = vec2d(p[0],p[1]) / (FLOAT64)resolution;
		it ++;
		index ++;
	}
#if CENTERS_DEF
	std::tr1::unordered_map<uint64,bool>::iterator it2 = centers.begin();
	while( it2 != centers.end() ) {
		cleanCenters[remap[(*it2).first]] = (*it2).second;
		it2 ++;
	}
#endif
	// Deallocate random nodes
	std::tr1::unordered_map<uint64,vec2i>().swap(nodes);
#if CENTERS_DEF
	std::tr1::unordered_map<uint64,bool>().swap(centers);
#endif
	// Cleanup elements reference
	cleanElements.resize(elements.size());
	for( uint n=0; n<elements.size(); n++ ) {
		cleanElements[n].resize(elements[n].size());
		for( uint m=0; m<elements[n].size(); m++ ) {
			cleanElements[n][m] = remap[elements[n][m]];
		}
	}	
	// Deallocate random elements
	std::vector<std::vector<uint64> >().swap(elements);
}

void bcc2::freeTets() {
	std::vector<vec2d>().swap(cleanNodes);
#if CENTERS_DEF
	std::vector<bool>().swap(cleanCenters);
#endif
	std::vector<std::vector<uint> >().swap(cleanElements);
}

void bcc2::precomputeShapeMatrix() {
	// Build shape function
	matrix.resize(cleanElements.size());
	for( uint n=0; n<cleanElements.size(); n++ ) {
		FLOAT64 A[NUM_VERT][NUM_VERT];
		for( uint i=0; i<NUM_VERT; i++ ) {
			for( uint dim=0; dim<DIM; dim++ ) {
				A[dim][i] = cleanNodes[cleanElements[n][i]][dim];
			}
			A[DIM][i] = 1.0;
		}
		if( ! invert3x3( A, matrix[n].m )) {
			printf( "Failed inverse 3x3!\n" );
		}
	}
}

bool bcc2::clearData() {
	for( uint n=0; n<hash.size(); n++ ) hash[n].clear();
	hash.clear();
	nodes.clear();
#if CENTERS_DEF
	centers.clear();
#endif
	cleanNodes.clear();
	for( uint n=0; n<elements.size(); n++ ) elements[n].clear();
	elements.clear();
	for( uint n=0; n<cleanElements.size(); n++ ) cleanElements[n].clear();
	cleanElements.clear();
	return true;
}

uint64 bcc2::computeCornerIndex( vec2i p ) const {
	uint64 R = resolution;
	for( uint dim=0; dim<DIM; dim++ ) {
		if( p[dim] < 0 || p[dim] > R ) return 0;
	}
	return p[0]+p[1]*R+1;
}

int bcc2::hitTest(vec2d p) const {
	FLOAT64 e = 1e-8;
	for( uint dim=0; dim<DIM; dim++) p[dim] = fmin(1.0-e,fmax(e,p[dim]));
	int t = octree.hitTest(p);
	if( t >= 0 ) {
		std::vector<uint>::const_iterator it;
		for( it=hash[t].begin(); it!=hash[t].end(); it++ ) {
			int index = *it;
			bool skip = false;
			FLOAT64 x[NUM_VERT] = { p[0], p[1], 1.0 };
			for( uint i=0; i<NUM_VERT; i++ ) {
				FLOAT64 t = 0.0;
				for( uint j=0;j<NUM_VERT;j++) {
					t += matrix[index].m[i][j]*x[j];
				}
				if( t < -1e-8 || t > 1.0+1e-8 ) {
					skip = true;
					break;
				}
			}
			if( ! skip ) return index;
		}
	}
	email::print("BCC hit test failed (%f,%f)\n", p[0], p[1] );
	email::send();
	exit(0);
	return -1;
}

void bcc2::drawElements() {
	for( uint n=0; n<cleanElements.size(); n++ ) {
		glBegin(GL_LINE_LOOP);
		for( uint m=0; m<NUM_VERT; m++ ) {
			vec2d p = cleanNodes[cleanElements[n][m]];
			glVertex2dv(p.v);
		}
		glEnd();
	}
}