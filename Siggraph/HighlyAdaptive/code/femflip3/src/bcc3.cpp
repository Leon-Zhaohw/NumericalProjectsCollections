/*
 *	bcc3.cpp
 *
 *	Created by Ryoichi Ando on 5/28/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "util3.h"
#include "bcc3.h"
#include "matutil.h"
#include "levelset3.h"
#include "opengl.h"
using namespace std;

bool bcc3::buildBCC( const octree3 &octree ) {
	// Copy octree
	this->octree = octree;
	
	// Get resolution
	resolution = octree.resolution;
	
	// Clear the octree first
	if( ! clearData()) return false;
	
	tick(); dump( ">>> Building BCC tets...\n" );
	
	// Get the resolution
	resolution = octree.resolution;
	
	// Build an array of corner nodes
	tick(); dump( "Building corner nodes..." );
	buildCornerArray(octree.terminals);
	dump( "Done. Took %s.\n", stock());
	
	// Build elements
	tick(); dump( "Building elements..." );
	hash.resize(octree.terminals.size());
	buildElements(octree.terminals,octree,hash);
	dump( "Done. Took %s.\n", stock());
	
	// Cleanup node indices
	tick(); dump( "Cleaning up array..." );
	cleanupNodeAndElements();
	dump( "Done. Took %s.\n", stock());
	
	// Precompute shape matrix
	tick(); dump( "Precomputing BCC shape matrix..." );
	precomputeShapeMatrix();
	dump( "Done. Took %s.\n", stock());
	dump( "<<< Building BCC tets done. Took %s.\n", stock());
	
	// Write info
	writeNumber("BCC_tet_num", cleanElements.size());
	writeNumber("BCC_node_num", cleanNodes.size());
	
	return true;
}

void bcc3::buildCornerArray( const std::vector<octree3::leaf3 *> &terminals ) {
	for( uint n=0; n<terminals.size(); n++ ) {
		int dx = terminals[n]->dx;
		vec3i center = terminals[n]->center;
		// Insert corners
		FOR_EACH(2,2,2) {
			vec3i p = center-0.5*dx*vec3i(1,1,1)+dx*vec3i(i,j,k);
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

void bcc3::buildElements( const std::vector<octree3::leaf3 *> &terminals, const octree3 &octree, std::vector<std::vector<uint> > &hash ) {
	for( uint n=0; n<terminals.size(); n++ ) {
		// Build meshes here !
		buildTets(terminals[n],octree,hash);
	}
}

void bcc3::buildTets( octree3::leaf3 *leaf, const octree3 &octree, std::vector<std::vector<uint> > &hash ) {
	vec3i center = leaf->center;
	uint idx = leaf->index;
	int dx = leaf->dx;
	uint64 center_index = computeCornerIndex(center);
	// 1. Pick up a direction
	for( uint dim=0; dim<DIM; dim++ ) {
		// Compute nullspace vector
		vec3i vec1, vec2, dirvec;
		uint axis[2];
		if( dim==0 ) {
			vec1 = vec3i(0,1,0); axis[0] = 1;
			vec2 = vec3i(0,0,1); axis[1] = 2;
		} else if( dim==1 ) {
			vec1 = vec3i(1,0,0); axis[0] = 0;
			vec2 = vec3i(0,0,1); axis[1] = 2;
		} else if( dim==2 ) {
			vec1 = vec3i(1,0,0); axis[0] = 0;
			vec2 = vec3i(0,1,0); axis[1] = 1;
		}
		// For both side
		for( int dir=-1; dir<=1; dir+=2 ) {
			// 2. Move by the half of the cell size, see if a node exist there
			vec3i facepos = center+0.5*dx*dir*vec3i(dim==0,dim==1,dim==2);
			uint64 face_index = computeCornerIndex(facepos);
			if( nodes.find(face_index) != nodes.end() ) {
				// 3. If found, grade it
				int dirset[][2] = { {-1,-1}, {1,-1}, {1,1}, {-1,1} };
				uint64 indices[4];
				vec3i pos[4];
				// Fetch surrounding corner nodes
				for( uint n=0; n<4; n++ ) {
					pos[n] = facepos+0.5*dx*vec1*dirset[n][0]+0.5*dx*vec2*dirset[n][1];
					indices[n] = computeCornerIndex(pos[n]);
				}
				// Build grading mesh
				std::vector<uint64> element(4);
				element[0] = center_index;
				element[1] = face_index;
				for( uint n=0; n<4; n++ ) {
					vec3i cutpos = (pos[n]+pos[(n+1)%4])/2;
					uint64 cutindex = computeCornerIndex(cutpos);
					if( nodes.find(cutindex) == nodes.end() ) {
						element[2] = indices[n];
						element[3] = indices[(n+1)%4];
						elements.push_back(element);
						hash[idx].push_back(elements.size()-1);
					} else {
						element[2] = indices[n];
						element[3] = cutindex;
						elements.push_back(element);
						hash[idx].push_back(elements.size()-1);
						element[2] = indices[(n+1)%4];
						element[3] = cutindex;
						elements.push_back(element);
						hash[idx].push_back(elements.size()-1);
					}
				}
			} else {
				// 4. If not found, move further by the half of the cell size and try find a node there
				vec3i adjacent_pos = center+dx*dir*vec3i(dim==0,dim==1,dim==2);
				uint64 adjacent_index = computeCornerIndex(adjacent_pos);
				if( nodes.find(adjacent_index) == nodes.end() ) {
					// 5. If not found, generate grading tets there
					std::vector<uint64> element(4);
					element[0] = center_index;
					uint position[2] = { leaf->position[axis[0]], leaf->position[axis[1]] };
					vec3i diagonal_vec;
					vec3i non_diagonal_vec;
					if( position[0] == position[1] ) {
						diagonal_vec = 0.5*dx*vec1+0.5*dx*vec2;
						non_diagonal_vec = 0.5*dx*vec1-0.5*dx*vec2;
					} else {
						diagonal_vec = 0.5*dx*vec1-0.5*dx*vec2;
						non_diagonal_vec = 0.5*dx*vec1+0.5*dx*vec2;
					}
					uint64 indices[4];
					indices[0] = center_index;
					indices[1] = computeCornerIndex(facepos+diagonal_vec);
					indices[2] = computeCornerIndex(facepos-diagonal_vec);
					for( int non_diagonal_dir=-1; non_diagonal_dir<=1; non_diagonal_dir+=2 ) {
						indices[3] = computeCornerIndex(facepos+non_diagonal_dir*non_diagonal_vec);
						// Look cut point, and decide whether we should split there
						char cuttable[2] = { 0, 0 };
						uint64 cutindex[2] = { 0, 0 };
						for( int cutdir=-1; cutdir<=2; cutdir+=2 ) {
							vec3i cutpos = facepos+(non_diagonal_dir*non_diagonal_vec+cutdir*diagonal_vec)/2;
							uint n = cutdir==-1 ? 0 : 1;
							cutindex[n] = computeCornerIndex(cutpos);
							if( nodes.find(cutindex[n]) != nodes.end() ) {
								cuttable[n] = 1;
							}
						}
						std::vector<uint64> element(4);
						if( cuttable[0] == 0 && cuttable[1] == 0 ) {
							for( uint n=0; n<4; n++ ) {
								element[n] = indices[n];
							}
							elements.push_back(element);
							hash[idx].push_back(elements.size()-1);
						} else if( cuttable[0] == 1 && cuttable[1] == 0) {
							element[0] = indices[0];
							element[1] = indices[1];
							element[2] = indices[2];
							element[3] = cutindex[0];
							elements.push_back(element);
							hash[idx].push_back(elements.size()-1);
							element[0] = indices[0];
							element[1] = indices[1];
							element[2] = indices[3];
							element[3] = cutindex[0];
							elements.push_back(element);
							hash[idx].push_back(elements.size()-1);
						} else if( cuttable[0] == 0 && cuttable[1] == 1 ) {
							element[0] = indices[0];
							element[1] = indices[1];
							element[2] = indices[2];
							element[3] = cutindex[1];
							elements.push_back(element);
							hash[idx].push_back(elements.size()-1);
							element[0] = indices[0];
							element[1] = indices[2];
							element[2] = indices[3];
							element[3] = cutindex[1];
							elements.push_back(element);
							hash[idx].push_back(elements.size()-1);
						} else if( cuttable[0] == 1 && cuttable[1] == 1 ) {
							// What a worst case
							element[0] = indices[0];
							element[1] = indices[2];
							element[2] = cutindex[1];
							element[3] = indices[1];
							elements.push_back(element);
							hash[idx].push_back(elements.size()-1);
							element[0] = indices[0];
							element[1] = cutindex[1];
							element[2] = indices[2];
							element[3] = cutindex[0];
							elements.push_back(element);
							hash[idx].push_back(elements.size()-1);
							element[0] = indices[0];
							element[1] = indices[3];
							element[2] = cutindex[0];
							element[3] = cutindex[1];
							elements.push_back(element);
							hash[idx].push_back(elements.size()-1);
						}
					}
				} else if( adjacent_index < center_index ) {
					vec3i adjacent_pos = center+dx*dir*vec3i(dim==0,dim==1,dim==2);
					int nidx = octree.hitTest(vec3d(adjacent_pos)/octree.resolution);
					
					// 6. If found, pick up another direction and move by the half of the grid cell, and see if a node exists there
					vec3i cutpos[4] = { facepos+0.5*dx*vec1, facepos+0.5*dx*vec2, facepos-0.5*dx*vec1, facepos-0.5*dx*vec2 };
					vec3i slide_vec[4] = { 0.5*dx*vec2, 0.5*dx*vec1, 0.5*dx*vec2, 0.5*dx*vec1 };
					for( uint n=0; n<4; n++ ) {
						uint64 cut_index = computeCornerIndex(cutpos[n]);
						if( nodes.find(cut_index) != nodes.end() ) {
							// If found, build two tets there
							std::vector<uint64> element(4);
							element[0] = center_index;
							element[1] = adjacent_index;
							element[2] = cut_index;
							for( int slide_dir=-1; slide_dir<=1; slide_dir+=2 ) {
								element[3] = computeCornerIndex(cutpos[n]+slide_dir*slide_vec[n]);
								elements.push_back(element);
								hash[idx].push_back(elements.size()-1);
								if( nidx >= 0 ) hash[nidx].push_back(elements.size()-1);
							}
						} else {
							// If not found, build a regular BCC tet there
							std::vector<uint64> element(4);
							element[0] = center_index;
							element[1] = adjacent_index;
							element[2] = computeCornerIndex(cutpos[n]+slide_vec[n]);
							element[3] = computeCornerIndex(cutpos[n]-slide_vec[n]);
							elements.push_back(element);
							hash[idx].push_back(elements.size()-1);
							if( nidx >= 0 ) hash[nidx].push_back(elements.size()-1);
						}
					}
				}
			}
		}
	}
}

void bcc3::cleanupNodeAndElements() {
	std::tr1::unordered_map<uint64,uint> remap;
	
	// Compute remap array
	uint size = nodes.size();
	cleanNodes.resize(size);
#if CENTERS_DEF
	cleanCenters.resize(size);
#endif
	uint index=0;
	std::tr1::unordered_map<uint64,vec3i>::iterator it = nodes.begin();
	while( it != nodes.end() ) {
		remap[(*it).first] = index;
		vec3i p = (*it).second;
		cleanNodes[index] = vec3d(p[0],p[1],p[2]) / (FLOAT64)resolution;
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
	std::tr1::unordered_map<uint64,vec3i>().swap(nodes);
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

void bcc3::freeTets() {
	std::vector<vec3d>().swap(cleanNodes);
#if CENTERS_DEF
	std::vector<bool>().swap(cleanCenters);
#endif
	std::vector<std::vector<uint> >().swap(cleanElements);
}

static vec3i computeCenterPos( vec3i center, uint dx, uint i, uint j, uint k ) {
	return center-0.25*dx*vec3i(1,1,1)+0.5*dx*vec3i(i,j,k);
}

void bcc3::precomputeShapeMatrix() {
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
		if( ! invert4x4( A, matrix[n].m )) {
			printf( "Failed inverse 4x4!\n" );
		}
	}
}

bool bcc3::clearData() {
	for( uint n=0; n<hash.size(); n++ ) hash[n].clear();
	hash.clear();
	std::vector<std::vector<uint> >().swap(hash);
	matrix.clear();
	std::vector<shapeM>().swap(matrix);
	nodes.clear();
#if CENTERS_DEF
	centers.clear();
#endif
	cleanNodes.clear();
	for( uint n=0; n<elements.size(); n++ ) elements[n].clear();
	elements.clear();
	for( uint n=0; n<cleanElements.size(); n++ ) cleanElements[n].clear();
	cleanElements.clear();
	freeTets();
	return true;
}

uint64 bcc3::computeCornerIndex( vec3i p ) const {
	uint64 R = resolution;
	for( uint dim=0; dim<DIM; dim++ ) {
		if( p[dim] < 0 || p[dim] > R ) return 0;
	}
	return p[0]+p[1]*R+p[2]*R*R+1;
}

int bcc3::hitTest(vec3d p) const {
	FLOAT64 e = 1e-8;
	for( uint dim=0; dim<DIM; dim++) p[dim] = fmin(1.0-e,fmax(e,p[dim]));
	int t = octree.hitTest(p);
	if( t >= 0 ) {
		std::vector<uint>::const_iterator it;
		for( it=hash[t].begin(); it!=hash[t].end(); it++ ) {
			int index = *it;
			bool skip = false;
			FLOAT64 x[NUM_VERT] = { p[0], p[1], p[2], 1.0 };
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
	email::print("BCC hit test failed (%f,%f,%f)\n", p[0], p[1], p[2] );
	email::send();
	exit(0);
	return -1;
}

// Rest of lines are for the visualization...

void bcc3::drawElements() {
	glPushMatrix();
	glScaled(1.0/resolution,1.0/resolution,1.0/resolution);
	glColor4d(0.5,0.5,0.5,1.0);
	for( uint n=0; n<elements.size(); n++ ) {
		glBegin(GL_LINES);
		for( uint m0=0; m0<elements[n].size(); m0++ ) {
			for( uint m1=m0+1; m1<elements[n].size(); m1++ ) {
				vec3i p0 = nodes[elements[n][m0]];
				vec3i p1 = nodes[elements[n][m1]];
				glVertex3d(p0[0],p0[1],p0[2]);
				glVertex3d(p1[0],p1[1],p1[2]);
			}
		}
		glEnd();
	}
	glPopMatrix();
}

void bcc3::writeObj( const char *path ) const {
	FILE *fp = fopen(path,"w");
	if( ! fp ) return;
	// Write vertices
	for( uint n=0; n<cleanNodes.size(); n++ ) {
		fprintf(fp, "v %f %f %f\n", cleanNodes[n][0], cleanNodes[n][1], cleanNodes[n][2] );
	}
	// Write elements
	for( uint n=0; n<elements.size(); n++ ) {
		// Compute center position
		vec3d center = (cleanNodes[cleanElements[n][0]]+
						cleanNodes[cleanElements[n][1]]+
						cleanNodes[cleanElements[n][2]]+
						cleanNodes[cleanElements[n][3]])/4.0;
		if( center[2] < 0.5 ) {
			for( uint v1=0; v1<4; v1++ ) for( uint v2=v1+1; v2<4; v2++ ) for( uint v3=v2+1; v3<4; v3++ ) {
				fprintf(fp, "f %d %d %d\n", cleanElements[n][v1]+1, cleanElements[n][v2]+1, cleanElements[n][v3]+1 );
			}
		}
	}
	fclose(fp);
}
