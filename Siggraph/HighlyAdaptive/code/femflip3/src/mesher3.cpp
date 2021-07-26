/*
 *	mesher3.cpp
 *	
 *	Created by Ryoichi Ando on 2/1/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "tetgen.h"
#include "mesher3.h"
#include "levelset3.h"
#include "util3.h"
#include "pcgsolver/matutil.h"
#include "opengl.h"
#include <stdio.h>
#include <stdlib.h>
using namespace std;

mesher3::mesher3() {
	centerType = BARYCENTRIC;
	genHash = true;
	genFacet = true;
}

void mesher3::setCenterType( int type ) {
	centerType = type;
}

void mesher3::setGenHash( bool state ) {
	genHash = state;
}

void mesher3::setGenFacet( bool state ) {
	genFacet = state;
}

static void buildTetgenBoxInput( tetgenio &in, uint gn ) {
	FLOAT64 dx = 1.0/(gn-1);
	vector<vec3d> grids;
	// Create tet nodes
	for( uint k=0; k<=gn; k++ ) {
		for( uint j=0; j<=gn; j++ ) {
			for( uint i=0; i<=gn; i++ ) {
				// Corner nodes
				grids.push_back(vec3d(dx*i,dx*j,dx*k));
			}
		}
	}
	// Insert center nodes
	for( uint j=0; j<gn; j++ ) {
		for( uint i=0; i<gn; i++ ) {
			for( uint k=0; k<gn; k++ ) {
				grids.push_back(vec3d(dx*(i+0.5),dx*(j+0.5),dx*(k+0.5)));
			}
		}
	}
	
	// Generate grid points
	vec3d maxp(0.0,0.0,0.0);
	for( uint n=0; n<grids.size(); n++ ) {
		maxp[0] = fmax(grids[n][0],maxp[0]);
		maxp[1] = fmax(grids[n][1],maxp[1]);
		maxp[2] = fmax(grids[n][2],maxp[2]);
	}
	vec3d shift(maxp[0]/2-0.5,maxp[1]/2-0.5,maxp[2]/2-0.5);
	for( uint n=0; n<grids.size(); n++ ) {
		grids[n] -= shift;
	}
	
	in.firstnumber = 1;
	in.numberofpoints = grids.size();
	in.pointlist = new REAL[in.numberofpoints * 3];
	for( uint n=0; n<grids.size(); n++ ) {
		in.pointlist[3*n+0] = grids[n][0];
		in.pointlist[3*n+1] = grids[n][1];
		in.pointlist[3*n+2] = grids[n][2];
	}
}

static bool lu_decmp(FLOAT64 lu[4][4], int n, int* ps, FLOAT64* d, int N)
{
	FLOAT64 scales[4];
	FLOAT64 pivot, biggest, mult, tempf;
	int pivotindex = 0;
	int i, j, k;
	
	*d = 1.0;                                      // No row interchanges yet.
	
	for (i = N; i < n + N; i++) {                             // For each row.
		// Find the largest element in each row for row equilibration
		biggest = 0.0;
		for (j = N; j < n + N; j++)
			if (biggest < (tempf = fabs(lu[i][j])))
				biggest  = tempf;
		if (biggest != 0.0)
			scales[i] = 1.0 / biggest;
		else {
			scales[i] = 0.0;
			return false;                            // Zero row: singular matrix.
		}
		ps[i] = i;                                 // Initialize pivot sequence.
	}
	
	for (k = N; k < n + N - 1; k++) {                      // For each column.
		// Find the largest element in each column to pivot around.
		biggest = 0.0;
		for (i = k; i < n + N; i++) {
			if (biggest < (tempf = fabs(lu[ps[i]][k]) * scales[ps[i]])) {
				biggest = tempf;
				pivotindex = i;
			}
		}
		if (biggest == 0.0) {
			return false;                         // Zero column: singular matrix.
		}
		if (pivotindex != k) {                         // Update pivot sequence.
			j = ps[k];
			ps[k] = ps[pivotindex];
			ps[pivotindex] = j;
			*d = -(*d);                          // ...and change the parity of d.
		}
		
		// Pivot, eliminating an extra variable  each time
		pivot = lu[ps[k]][k];
		for (i = k + 1; i < n + N; i++) {
			lu[ps[i]][k] = mult = lu[ps[i]][k] / pivot;
			if (mult != 0.0) {
				for (j = k + 1; j < n + N; j++)
					lu[ps[i]][j] -= mult * lu[ps[k]][j];
			}
		}
	}
	
	// (lu[ps[n + N - 1]][n + N - 1] == 0.0) ==> A is singular.
	return lu[ps[n + N - 1]][n + N - 1] != 0.0;
}

static FLOAT64 tetVolume( vec3d p0, vec3d p1, vec3d p2, vec3d p3 ) {
	// Set the edge vectors: V[0], ..., V[5]
	FLOAT64 V[6][3];
	FLOAT64 A[4][4];
	FLOAT64 D;
	int indx[4];
	int i, j;
    for (i = 0; i < 3; i++) V[0][i] = p0[i] - p3[i]; // V[0]: p3->p0.
    for (i = 0; i < 3; i++) V[1][i] = p1[i] - p3[i]; // V[1]: p3->p1.
    for (i = 0; i < 3; i++) V[2][i] = p2[i] - p3[i]; // V[2]: p3->p2.
    for (i = 0; i < 3; i++) V[3][i] = p1[i] - p0[i]; // V[3]: p0->p1.
    for (i = 0; i < 3; i++) V[4][i] = p2[i] - p1[i]; // V[4]: p1->p2.
    for (i = 0; i < 3; i++) V[5][i] = p0[i] - p2[i]; // V[5]: p2->p0.
	// Set the matrix A = [V[0], V[1], V[2]]^T.
    for (j = 0; j < 3; j++) {
		for (i = 0; i < 3; i++) A[j][i] = V[j][i];
    }
	// Decompose A just once.
    if( lu_decmp(A, 3, indx, &D, 0)) return fabs(A[indx[0]][0] * A[indx[1]][1] * A[indx[2]][2]) / 6.0;
	else return 0.0;
}

static void putHash( vec3d left_bottom, vec3d right_top, array3<std::vector<uint> > &hash, uint id ) {
	FLOAT64 e = 1.0e-8;
	uint gn[3] = { hash.size().w, hash.size().h, hash.size().d };
	uint min_i = gn[0]*fmin(1.0-e,fmax(0.0,left_bottom[0]+e));
	uint min_j = gn[1]*fmin(1.0-e,fmax(0.0,left_bottom[1]+e));
	uint min_k = gn[2]*fmin(1.0-e,fmax(0.0,left_bottom[2]+e));
	uint max_i = gn[0]*fmin(1.0-e,fmax(0.0,right_top[0]-e));
	uint max_j = gn[1]*fmin(1.0-e,fmax(0.0,right_top[1]-e));
	uint max_k = gn[2]*fmin(1.0-e,fmax(0.0,right_top[2]-e));
	
	max_i = imin(max_i+1,gn[0]);
	max_j = imin(max_j+1,gn[1]);
	max_k = imin(max_k+1,gn[2]);

	for( uint i=min_i; i<max_i; i++ ) for( uint j=min_j; j<max_j; j++ ) for( uint k=min_k; k<max_k; k++ ) {
		hash[i][j][k].push_back(id);
	}
}

static void getBoundingBox( vector<vec3d> &points, vec3d &left_bottom, vec3d &right_top ) {
	FLOAT64 r1x = 9999.0;
	FLOAT64 r1y = 9999.0;
	FLOAT64 r1z = 9999.0;
	FLOAT64 r2x = -999.0;
	FLOAT64 r2y = -999.0;
	FLOAT64 r2z = -999.0;
	for( uint n=0; n<points.size(); n++ ) {
		FLOAT64 x = points[n][0];
		FLOAT64 y = points[n][1];
		FLOAT64 z = points[n][2];
		if( x < r1x ) r1x = x;
		if( y < r1y ) r1y = y;
		if( z < r1z ) r1z = z;
		if( x > r2x ) r2x = x;
		if( y > r2y ) r2y = y;
		if( z > r2z ) r2z = z;
	}
	left_bottom[0] = r1x;
	left_bottom[1] = r1y;
	left_bottom[2] = r1z;
	right_top[0] = r2x;
	right_top[1] = r2y;
	right_top[2] = r2z;
}

static bool checkElementQuality( vec3d vertices[] ) {
	bool passed = true;
	FLOAT64 A[DIM+1][DIM+1];
	for( uint i=0; i<DIM+1; i++ ) for( uint j=0; j<DIM+1; j++ ) {
		if( i < DIM ) A[i][j] = vertices[j][i];
		else A[i][j] = 1.0;
	}
	FLOAT64 Ainv[DIM+1][DIM+1];
	if( invert4x4(A,Ainv) ) {
		FLOAT64 B[DIM][DIM+1];
		for( uint i=0; i<DIM+1; i++ ) for( uint j=0; j<DIM; j++ ) {
			B[i][j] = Ainv[j][i];
		}
		// BtB = B^t * B
		FLOAT64 BtB[DIM+1][DIM+1];
		for( uint i=0; i<DIM+1; i++ ) for( uint j=0; j<DIM+1; j++ ) {
			BtB[i][j] = 0.0;
			for( uint k=0; k<DIM; k++ ) BtB[i][j] += B[k][i]*B[k][j];
		}
		// Check that the BtB matrix does not have any positive non-diagonal term
		passed = true;
		for( uint i=0; i<DIM+1; i++ ) for( uint j=0; j<DIM+1; j++ ) {
			if( i < j && BtB[i][j] > 1e-8 ) {
				passed = false;
			}
		}
#if 0
		if( ! passed ) {
			dump( "WARNING: Badly shaped element detected !\n" );
			dump( "A = \n" );
			dump( "[%.2f %.2f %.2f %.2f;\n", A[0][0], A[0][1], A[0][2], A[0][3] );
			dump( " %.2f %.2f %.2f %.2f;\n", A[1][0], A[1][1], A[1][2], A[1][3] );
			dump( " %.2f %.2f %.2f %.2f;\n", A[2][0], A[2][1], A[2][2], A[2][3] );
			dump( " %.2f %.2f %.2f %.2f]\n", A[3][0], A[3][1], A[3][2], A[3][3] );
			dump( "BtB = \n" );
			dump( "[%.2f %.2f %.2f %.2f;\n", BtB[0][0], BtB[0][1], BtB[0][2], BtB[0][3] );
			dump( " %.2f %.2f %.2f %.2f;\n", BtB[1][0], BtB[1][1], BtB[1][2], BtB[1][3] );
			dump( " %.2f %.2f %.2f %.2f;\n", BtB[2][0], BtB[2][1], BtB[2][2], BtB[2][3] );
			dump( " %.2f %.2f %.2f %.2f]\n", BtB[3][0], BtB[3][1], BtB[3][2], BtB[3][3] );
		}
#endif
	} else {
		passed = false;
	}
	return passed;
}

void mesher3::cleanup() {
	std::vector<vec3d>().swap(nodes);							// Grid node array
	std::vector<std::vector<uint> >().swap(elements);			// Elements array
	std::vector<vec3d>().swap(centers);							// Element center positions
	std::vector<shapeM>().swap(matrix);							// Element shape function matrix
	
	// Information below may not be available
	std::vector<std::vector<uint> >().swap(node_elements);		// Node to element array. This can be interpreted as Voronoi elements
	std::vector<std::vector<uint> >().swap(facets);				// Elements facets
	std::vector<FLOAT64>().swap(facetArea);						// Facet area
	std::vector<std::vector<uint> >().swap(facet_elements);		// Facet to element array
	std::vector<std::vector<uint> >().swap(edges);				// Edges
	std::vector<std::vector<uint> >().swap(element_facets);		// Element to facet array
	std::vector<std::vector<uint> >().swap(element_edges);		// Element to edges
	std::vector<std::vector<uint> >().swap(element_elements);	// Adjacent element - element array
	std::vector<std::vector<uint> >().swap(node2node);			// Node to node connection array
    std::vector<FLOAT64>().swap(volumes);						// Element volumes
	hash.clear();												// Element spatial hash
}

void mesher3::generateElements( std::vector<vec3d> &nodes, std::vector<std::vector<uint> > &elements, const levelset3 *hint, uint gn ) {
	tetgenio in, out;
	buildTetgenBoxInput(in,gn);
	tetrahedralize((char *)"Qnnfq1.41",&in,&out);
#if 0
	// Output mesh to files ’barout.node’, ’barout.ele’ and ’barout.face’.
	out.save_nodes("box.node");
	out.save_elements("box.ele");
	out.save_faces("box.face");
#endif
	// Build node informaion
	nodes.clear();
	nodes.resize(out.numberofpoints);
	for( uint n=0; n<out.numberofpoints; n++ ) {
		nodes[n] = vec3d(out.pointlist[3*n+0],out.pointlist[3*n+1],out.pointlist[3*n+2]);
	}
	
	// Build elements information
	elements.clear();
	elements.resize(out.numberoftetrahedra);
	for( uint n=0; n<out.numberoftetrahedra; n++ ) {
		elements[n].resize(out.numberofcorners);
		for( uint m=0; m<out.numberofcorners; m++ ) {
			elements[n][m] = out.tetrahedronlist[n*out.numberofcorners+m]-1;
		}
	}
}

void mesher3::updateConnection() {
	// Quality check for debug
#if 0
	for( uint n=0; n<elements.size(); n++ ) {
		vec3d vertices[DIM+1];
		for( uint m=0; m<elements[n].size(); m++ ) vertices[m] = nodes[elements[n][m]];
		vec3d center = (vertices[0]+vertices[1]+vertices[2]+vertices[3])/4.0;
		if( ! checkElementQuality(vertices)) {
			const char *name = "bad_element.m";
			FILE *fp = fopen(name,"w");
			if( fp ) {
				fprintf(fp,"T=[1 2 3 4];\n");
				fprintf(fp,"V=[");
				for( uint i=0; i<4; i++ ) {
					fprintf(fp,"%f %f %f;\n",vertices[i][0],vertices[i][1],vertices[i][2]);
				}
				fprintf(fp,"];\n");
				for( uint i=0; i<4; i++ ) {
					fprintf(fp,"p%d = [%f %f %f];\n",i+1,vertices[i][0],vertices[i][1],vertices[i][2]);
				}
				fprintf(fp,"tetinfo(p1,p2,p3,p4);\n");
				fclose(fp);
			}
			dump("MATLAB file \"%s\" exported.\n", name );
			exit(0);
		}
	}
#endif
	
	// Build node to elements matrices that tells which elements a node share
	node_elements.clear();
	node_elements.resize(nodes.size());
	for( uint n=0; n<node_elements.size(); n++ ) node_elements[n].clear();
	for( uint n=0; n<elements.size(); n++ ) {
		for( uint m=0; m<elements[n].size(); m++ ) {
			node_elements[elements[n][m]].push_back(n);
		}
	}
	
	// Build shape function
	matrix.resize(elements.size());
	for( uint n=0; n<elements.size(); n++ ) {
		FLOAT64 A[NUM_VERT][NUM_VERT];
		for( uint i=0; i<NUM_VERT; i++ ) {
			for( uint dim=0; dim<DIM; dim++ ) {
				A[dim][i] = nodes[elements[n][i]][dim];
			}
			A[DIM][i] = 1.0;
		}
		if( ! invert4x4( A, matrix[n].m )) {
			printf( "Failed inverse 4x4!\n" );
		}
	}
	
	// Build node to node matrix based on the elements
	node2node.clear();
	node2node.resize(nodes.size());
	for( uint n=0; n<node2node.size(); n++ ) node2node[n].clear();
	for( uint n=0; n<elements.size(); n++ ) {
		for( uint m1=0; m1<elements[n].size(); m1++ ) {
			uint idx1 = elements[n][m1];
			for( uint m2=0; m2<elements[n].size(); m2++ ) {
				uint idx2 = elements[n][m2];
				if( idx1 < idx2 ) {
					bool duplicated = false;
					for( uint k=0; k<node2node[idx1].size(); k++ ) {
						if( node2node[idx1][k] == idx2 ) {
							duplicated = true;
							break;
						}
					}
					if( ! duplicated ) {
						node2node[idx1].push_back(idx2);
						node2node[idx2].push_back(idx1);
					}
				}
			}
		}
	}
			
	//Build facet information
	if( genFacet ) {
		facets.clear();
		facet_elements.clear();
		for( uint n=0; n<facets.size(); n++ ) node2node[n].clear();
		for( uint n=0; n<facet_elements.size(); n++ ) node2node[n].clear();
		for( uint n=0; n<node_elements.size(); n++ ) {
			for( uint m0=0; m0<node_elements[n].size(); m0++ ) for( uint m1=m0+1; m1<node_elements[n].size(); m1++ ) {
				uint element0 = node_elements[n][m0];
				uint element1 = node_elements[n][m1];
				// Check that these two elements are adjacent each other
				uint count = 0;
				uint node_indices[NUM_VERT];
				for( uint v0=0; v0<NUM_VERT; v0++ ) for( uint v1=0; v1<NUM_VERT; v1++ ) {
					if( elements[element0][v0] == elements[element1][v1] ) {
						node_indices[count] = elements[element0][v0];
						count ++;
					}
				}
				// See if adjacent
				if( count == DIM ) {
					std::vector<uint> face(DIM);
					std::vector<uint> adjacent(2);
					for( uint k=0; k<DIM; k++ ) face[k] = node_indices[k];
					facets.push_back(face);
					adjacent[0] = element0;
					adjacent[1] = element1;
					facet_elements.push_back(adjacent);
				}
			}
		}
		
		// Compute facet area
		facetArea.clear();
		facetArea.resize(facets.size());
		for( uint n=0; n<facets.size(); n++ ) {
			vec3d face[3];
			for( uint m=0; m<3; m++ ) face[m] = nodes[facets[n][m]];
			facetArea[n] = 0.5*((face[1]-face[0])^(face[2]-face[0])).len();
		}
		
		// Build element_elements and element_facets
		element_elements.clear();
		element_elements.resize(elements.size());
		for( uint n=0; n<element_elements.size(); n++ ) element_elements[n].clear();
		element_facets.clear();
		element_facets.resize(elements.size());
		for( uint n=0; n<element_facets.size(); n++ ) element_facets[n].clear();
		for( uint n=0; n<facet_elements.size(); n++ ) {
			for( uint m1=0; m1<facet_elements[n].size(); m1++ ) {
				for( uint m2=0; m2<facet_elements[n].size(); m2++ ) {
					if( m1 == m2 ) continue;
					element_elements[facet_elements[n][m1]].push_back(facet_elements[n][m2]);
					element_facets[facet_elements[n][m1]].push_back(n);
				}
			}
		}
	}
	
	// Build edge array
	edges.clear();
	for( uint n=0; n<edges.size(); n++ ) edges[n].clear();
	element_edges.resize(elements.size());
	for( uint n=0; n<element_edges.size(); n++ ) element_edges[n].clear();
	for( uint n=0; n<node2node.size(); n++ ) {
		for( uint m=0; m<node2node[n].size(); m++ ) {
			if( n < node2node[n][m] ) {
				// New edge
				std::vector<uint> edge(2);
				edge[0] = n;
				edge[1] = node2node[n][m];
				edges.push_back(edge);
				for( uint i=0; i<node_elements[n].size(); i++) {
					uint eidx = node_elements[n][i];
					for( uint j=0; j<elements[eidx].size(); j++ ) {
						if( elements[eidx][j] == node2node[n][m] ) {
							element_edges[eidx].push_back(edges.size()-1);
							if( element_edges[eidx].size() > 6 ) {
								printf( "Something wrong !! %d\n", element_edges[eidx].size());
								exit(0);
							}
						}
					}
				}
			}
		}
	}
	
	// Compute volume
	volumes.clear();
	volumes.resize(elements.size());
	for( uint n=0; n<elements.size(); n++ ) {
		volumes[n] = tetVolume(nodes[elements[n][0]],nodes[elements[n][1]],nodes[elements[n][2]],nodes[elements[n][3]]);
	}
	
	// Build an uniform hash
	if( genHash ) {
		hash.resize(gn,gn,gn);
		hash.clear();
		for( uint n=0; n<elements.size(); n++ ) {
			vector<vec3d> posarray(elements[n].size());
			for( uint m=0; m<elements[n].size(); m++ ) {
				FLOAT64 x = nodes[elements[n][m]][0];
				FLOAT64 y = nodes[elements[n][m]][1];
				FLOAT64 z = nodes[elements[n][m]][2];
				posarray[m] = vec3d(x,y,z);
			}
			vec3d leftbottom;
			vec3d righttop;
			getBoundingBox(posarray,leftbottom,righttop);
			putHash(leftbottom,righttop,hash,n);
		}
	}
	
	// Compute center positions
	computeCenterPosition(nodes,elements,centers,centerType);
}

void mesher3::generateMesh( uint gn, const levelset3 *hint ) {
	this->gn = gn;
	
	// Generate nodes and matrix
	generateElements(nodes,elements,hint,gn);
	
	// Build element connections
	tick(); dump("Organizing connectivity...");
	updateConnection();
	dump("Done. Took %s.\n",stock("mesher_organize"));
}

void mesher3::computeCenterPosition( const std::vector<vec3d> &nodes, const std::vector<std::vector<uint> > &elements, std::vector<vec3d> &centers, int centerType ) {
	// Build element centers
	centers.clear();
	centers.resize(elements.size());
	for( uint n=0; n<elements.size(); n++ ) {
		if( centerType == BARYCENTRIC ) {
			vec3d cntr;
			for( uint m=0; m<elements[n].size(); m++ ) {
				cntr += nodes[elements[n][m]];
			}
			centers[n] = cntr / elements[n].size();
		} else {
			REAL radius;
			tetgenmesh mesh;
			vec3<double> nodepos[NUM_VERT];
			for( uint k=0; k<NUM_VERT; k++ ) for( uint dim=0; dim<DIM; dim++ ) {
				nodepos[k][dim] = nodes[elements[n][k]][dim];
			}
			
			vec3<double> center;
			mesh.circumsphere(nodepos[0].v,nodepos[1].v,nodepos[2].v,nodepos[3].v,center.v,&radius);
			for( uint dim=0; dim<DIM; dim++ ) {
				centers[n][dim] = center[dim];
			}
		}
	}
}

bool mesher3::hitElement( uint index, vec3d p ) const {
	FLOAT64 x[NUM_VERT] = { p[0], p[1], p[2], 1.0 };
	for( uint i=0; i<NUM_VERT; i++ ) {
		FLOAT64 t = 0.0;
		for( uint j=0;j<NUM_VERT;j++) {
			t += matrix[index].m[i][j]*x[j];
		}
		if( t < -1e-8 || t > 1.0+1e-8 ) {
			return false;
		}
	}
	return true;
}

int	mesher3::hitElements(vec3d p) const {
	int ref = -1;
	if( genHash ) {
		FLOAT64 e = 1e-8;
		uint gn = hash.size().w;
		uint i = gn*fmin(1.0-e,fmax(0.0,p[0]));
		uint j = gn*fmin(1.0-e,fmax(0.0,p[1]));
		uint k = gn*fmin(1.0-e,fmax(0.0,p[2]));
		for( uint n=0; n<hash[i][j][k].size(); n++ ) {
			if( hitElement(hash[i][j][k][n],p) ) {
				ref = hash[i][j][k][n];
				break;
			}
		}
	} else {
		printf( "Spatical hash generation is disabled.\n" );
	}
	return ref;
}

static FLOAT64 fraction( FLOAT64 phi0, FLOAT64 phi1 ) {
	return phi0*phi1 > 0 ? 1 - (phi0 > 0) : -fmin(phi0,phi1)/fmax(fabs(phi1-phi0),1.0e-6);
}

FLOAT64 mesher3::getLevelsetVolume( uint elmid, FLOAT64 isoval[4] ) {
	uint cnt = 0;
	for( uint n=0; n<4; n++ ) {
		if( isoval[n] < 0 ) cnt++;
	}
	switch(cnt) {
		case 0:
			return 0.0;
		case 1: {
			vec3d v[4];
			uint idx=0;
			for( uint n=0; n<4; n++ ) {
				if( isoval[n] < 0 ) {
					if( idx >= 4 ) {
						dump( "Illegal index number at getLevelsetVolume:1.1\n");
						exit(0);
					}
					v[idx++] = nodes[elements[elmid][n]];
					for( uint m=0; m<4; m++ ) {
						if( n==m ) continue;
						if( isoval[m] >= 0 ) {
							FLOAT64 t = fraction(isoval[n],isoval[m]);
							if( idx >= 4 ) {
								dump( "Illegal index number at getLevelsetVolume:1.2\n");
								exit(0);
							}
							v[idx++] = (1.0-t)*nodes[elements[elmid][n]] + t*nodes[elements[elmid][m]];
						}
					}
				}
			}
			return tetVolume(v[0],v[1],v[2],v[3]);
			break;
		}
		case 2: {
			// Gather floor and ceiling position
			vec3d v[6];
			uint idx=0;
			for( uint n=0; n<4; n++ ) {
				if( isoval[n] < 0 ) {
					if( idx >= 6 ) {
						dump( "Illegal index number at getLevelsetVolume:2.1\n");
						exit(0);
					}
					v[idx++] = nodes[elements[elmid][n]];
					for( uint m=0; m<4; m++ ) {
						if( n==m ) continue;
						if( isoval[m] >= 0 ) {
							FLOAT64 t = fraction(isoval[n],isoval[m]);
							if( idx >= 6 ) {
								dump( "Illegal index number at getLevelsetVolume:2.2\n");
								exit(0);
							}
							v[idx++] = (1.0-t)*nodes[elements[elmid][n]] + t*nodes[elements[elmid][m]];
						}
					}
				}
			}
			// Compute volume
			return tetVolume(v[0],v[1],v[2],v[3])+tetVolume(v[1],v[2],v[3],v[4])+tetVolume(v[2],v[3],v[4],v[5]);
			break;
		}
		case 3: {
			vec3d v[4];
			uint idx=0;
			for( uint n=0; n<4; n++ ) {
				if( isoval[n] >= 0 ) {
					if( idx >= 4 ) {
						dump( "Illegal index number at getLevelsetVolume:3.1\n");
						exit(0);
					}
					v[idx++] = nodes[elements[elmid][n]];
					for( uint m=0; m<4; m++ ) {
						if( n==m ) continue;
						if( isoval[m] < 0 ) {
							FLOAT64 t = fraction(isoval[n],isoval[m]);
							if( idx >= 4 ) {
								dump( "Illegal index number at getLevelsetVolume:3.2\n");
								exit(0);
							}
							v[idx++] = t*nodes[elements[elmid][n]] + (1.0-t)*nodes[elements[elmid][m]];
						}
					}
				}
			}
			return volumes[elmid]-tetVolume(v[0],v[1],v[2],v[3]);
			break;
		}
	}
	return volumes[elmid];
}

static void setLighting( int enabled ) {
	if( enabled ) {
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);
		
		GLfloat light0pos[] = { 0.0, 3.0, 5.0, 1.0 };
		GLfloat light1pos[] = { 5.0, 3.0, 0.0, 1.0 };
		
		glLightfv(GL_LIGHT0, GL_POSITION, light0pos);
		glLightfv(GL_LIGHT1, GL_POSITION, light1pos);
	} else {
		glDisable(GL_LIGHTING);
		glDisable(GL_LIGHT0);
		glDisable(GL_LIGHT1);
	}
}

static void setMaterial( int enabled ) {
	if( enabled ) {
		GLfloat water[] = { 0.4, 0.5, 0.8, 1.0 };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, water);
		glShadeModel(GL_FLAT);
	}
}

void mesher3::drawMesh() const {
	// Draw edge
#if 0
	glColor4f(1.0,1.0,1.0,1.0);
	glPointSize(2);
	glBegin(GL_POINTS);
	for( uint n=0; n<nodes.size(); n++ ) {
		glVertex3dv(nodes[n].v);
	}
	glEnd();
#endif
	
	// Draw facet
#if 0
	glColor4f(1.0,1.0,1.0,0.3);
	glLineWidth(1.0);
	for( uint n=0; n<facets.size(); n++ ) {
		glBegin(GL_LINE_LOOP);
		for( uint m=0; m<3; m++ ) {
			glVertex3dv(nodes[facets[n][m]].v);
		}
		glEnd();
	}
	glLineWidth(1.0);
#endif
	
	// Draw center points
#if 0
	glPointSize(2);
	glBegin(GL_POINTS);
	for( uint n=0; n<centers.size(); n++ ) {
		glColor3dv(centers[n].v);
		glVertex3dv(centers[n].v);
	}
	glEnd();
#endif
	
#if 0
	// Draw adjacent elements
	glBegin(GL_LINES);
	for( uint n=0; n<elements.size(); n++ ) {
		for( uint m=0; m<element_elements[n].size(); m++ ) {
			glColor3dv(centers[n].v);
			glVertex3dv(centers[n].v);
			glColor3dv(centers[element_elements[n][m]].v);
			glVertex3dv(centers[element_elements[n][m]].v);
		}
	}
	glEnd();
#endif
	
#if 0
	// Draw illegal circumceters
	glPointSize(2);
	glBegin(GL_POINTS);
	for( uint n=0; n<elements.size(); n++ ) {
		for( uint m=0; m<element_elements[n].size(); m++ ) {
			if( ! (centers[n]-centers[element_elements[n][m]]).len() ) {
				glColor3dv(centers[n].v);
				glVertex3dv(centers[n].v);
			}
		}
	}
	glEnd();
#endif
	
#if 0
	// Draw adjacent facets
	glBegin(GL_LINES);
	for( uint n=0; n<element_facets.size(); n++ ) {
		for( uint m=0; m<element_facets[n].size(); m++ ) {
			// Compute facet center
			vec3d ctnr;
			uint idx = element_facets[n][m];
			for( uint k=0; k<3; k++ ) ctnr += nodes[facets[idx][k]]/3.0;
			glColor3dv(centers[n].v);
			glVertex3dv(centers[n].v);
			glColor3dv(ctnr.v);
			glVertex3dv(ctnr.v);
		}
	}
	glEnd();
#endif
	
#if 0
	// Draw adjacent elements from facets side
	glBegin(GL_LINES);
	for( uint n=0; n<facet_elements.size(); n++ ) {
		// Compute facet center
		vec3d ctnr;
		for( uint k=0; k<3; k++ ) ctnr += nodes[facets[n][k]]/3.0;
		for( uint m=0; m<facet_elements[n].size(); m++ ) {
			uint idx = facet_elements[n][m];
			glColor3dv(ctnr.v);
			glVertex3dv(ctnr.v);
			glColor3dv(centers[idx].v);
			glVertex3dv(centers[idx].v);
		}
	}
	glEnd();
#endif
}

void mesher3::writeObj( const char *path, int type ) const {
	FILE *fp = fopen(path,"w");
	if( ! fp ) return;
	if( type == 0 ) {
		// Write vertices
		for( uint n=0; n<nodes.size(); n++ ) {
			fprintf(fp, "v %f %f %f\n", nodes[n][0], nodes[n][1], nodes[n][2] );
		}
		// Write elements
		for( uint n=0; n<elements.size(); n++ ) {
			
			vec3d vertices[DIM+1];
			for( uint m=0; m<elements[n].size(); m++ ) vertices[m] = nodes[elements[n][m]];
			vec3d center = (vertices[0]+vertices[1]+vertices[2]+vertices[3])/4.0;
			// Compute center position
			if( center[2] < 0.6 ) {
				for( uint v1=0; v1<4; v1++ ) for( uint v2=v1+1; v2<4; v2++ ) for( uint v3=v2+1; v3<4; v3++ ) {
					fprintf(fp, "f %d %d %d\n", elements[n][v1]+1, elements[n][v2]+1, elements[n][v3]+1 );
				}
			}
		}
	} else {
		util3::writeContainer(fp);
		// Write vertices
		for( uint n=0; n<elements.size(); n++ ) {
			fprintf(fp, "v %f %f %f\n", centers[n][0], centers[n][1], centers[n][2] );
		}
	}
	fclose(fp);
}
