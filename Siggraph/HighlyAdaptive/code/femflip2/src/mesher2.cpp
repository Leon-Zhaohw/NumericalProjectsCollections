/*
 *	mesher2.cpp
 *	
 *	Created by Ryoichi Ando on 12/26/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "mesher2.h"
#include "matutil.h"
#include "util2.h"
#include "opengl.h"
#include "pcgsolver/matutil.h"
#include "levelset2.h"
#include <stdio.h>
#include <stdlib.h>
using namespace std;

static FLOAT64 computeElementVolume( const std::vector<vec2d> &nodes, const std::vector<uint> &element ) {
	uint elmnum = element.size();
	std::vector<vec2d> points(elmnum);
	for( uint i=0; i<elmnum; i++ ) {
		points[i] = nodes[element[i]];
	}
	return util2::computeVolume(points);
}

static std::vector<uint> getNeighbors( const array2<vector<uint> > &sorter, const vec2d &p, uint gn ) {
	int i = gn*p[0]; i = imin(gn-1,imax(0,i));
	int j = gn*p[1]; j = imin(gn-1,imax(0,j));
	std::vector<uint> res;
	int r = 2;
	for( int w=-r; w<=r; w++ ) for( int h=-r; h<=r; h++ ) {
		int fi = i+w;
		int fj = j+h;
		if( fi<0 || fi>gn-1 || fj<0 || fj>gn-1 ) continue;
		res.insert(res.end(),sorter[fi][fj].begin(),sorter[fi][fj].end());
	}
	return res;
}

static std::vector<std::vector<uint> > makeDelaunayTriangles( const std::vector<vec2d> &points, uint gn ) {
	array2<vector<uint> > sorter(gn,gn);
	for( uint n=0; n<points.size(); n++ ) {
		int i = gn*points[n][0]; i = imin(gn-1,imax(0,i));
		int j = gn*points[n][1]; j = imin(gn-1,imax(0,j));
		sorter[i][j].push_back(n);
	}
	
	vector<vector<uint> > triangles;
	if( points.size() < 3 ) return triangles;
	for( uint i=0; i<points.size(); i++ ) {
		vector<uint> around_i = getNeighbors( sorter, points[i], gn );
		for( uint ni=0; ni<around_i.size(); ni++ ) {
			int j = around_i[ni];
			if( j <= i ) continue;
			for( uint nj=0; nj<around_i.size(); nj++ ) {
				int k = around_i[nj];
				if( k <= j ) continue;
				bool outside = true;
				for( uint nk=0; nk<around_i.size(); nk++ ) {
					int n = around_i[nk];
					if( n==i || n==j || n==k ) continue;
					FLOAT64 det = util2::detDelaunay( points[i], points[j], points[k], points[n] );
					if( det <= 0.0 ) {
						outside = false;
						break;
					}
				}
				if( outside ) {
					vector<uint> tri(3);
					tri[0] = i;
					tri[1] = j;
					tri[2] = k;
					triangles.push_back(tri);
				}
			}
		}
	}
	return triangles;
}

mesher2::mesher2() {
	centerType = BARYCENTRIC;
	genHash = true;
	genFacet = true;
}

void mesher2::setCenterType( int type ) {
	centerType = type;
}

void mesher2::setGenHash( bool state ) {
	genHash = state;
}

void mesher2::setGenFacet( bool state ) {
	genFacet = state;
}

static void putHash( vec2d left_bottom, vec2d right_top, array2<std::vector<uint> > &hash, uint id ) {
	FLOAT64 e = 1.0e-8;
	uint gn[2] = { hash.size().w, hash.size().h };
	uint min_i = gn[0]*fmin(1.0-e,fmax(0.0,left_bottom[0]+e));
	uint min_j = gn[1]*fmin(1.0-e,fmax(0.0,left_bottom[1]+e));
	uint max_i = gn[0]*fmin(1.0-e,fmax(0.0,right_top[0]-e));
	uint max_j = gn[1]*fmin(1.0-e,fmax(0.0,right_top[1]-e));
	
	max_i = imin(max_i+1,gn[0]);
	max_j = imin(max_j+1,gn[1]);
	
	for( uint i=min_i; i<max_i; i++ ) for( uint j=min_j; j<max_j; j++ ) {
		hash[i][j].push_back(id);
	}
}

static void getBoundingBox( vector<vec2d> &points, vec2d &left_bottom, vec2d &right_top ) {
	FLOAT64 r1x = 9999.0;
	FLOAT64 r1y = 9999.0;
	FLOAT64 r2x = -999.0;
	FLOAT64 r2y = -999.0;
	for( uint n=0; n<points.size(); n++ ) {
		FLOAT64 x = points[n][0];
		FLOAT64 y = points[n][1];
		if( x < r1x ) r1x = x;
		if( y < r1y ) r1y = y;
		if( x > r2x ) r2x = x;
		if( y > r2y ) r2y = y;
	}
	left_bottom[0] = r1x;
	left_bottom[1] = r1y;
	right_top[0] = r2x;
	right_top[1] = r2y;
}

void mesher2::cleanup() {
	std::vector<vec2d>().swap(nodes);							// Grid node array
	std::vector<std::vector<uint> >().swap(elements);			// Elements array
	std::vector<vec2d>().swap(centers);							// Element center positions
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

void mesher2::generateElements( std::vector<vec2d> &nodes, std::vector<std::vector<uint> > &elements, const levelset2 *hint, uint gn ) {
	FLOAT64 dx = 1.0/(gn-1);
	
	// Generate nodal points
	cleanup();
	
	// Create triangle nodes
	FLOAT64 w = sqrt(3.0)/2.0*dx;
	FLOAT64 rs = 0.35;
	for( uint j=0; j<gn; j++ ) for( uint i=0; i<gn*sqrt(3.0); i++ ) {
		if( i % 2 == 0 ) {
			if( j!=gn-1 ) {
				FLOAT64 x = i*w;
				FLOAT64 y = (j+0.5)*dx;
				x += (x>dx&&x<1.0-2*dx)*rs*nrand()*dx;
				y += (x>dx&&y<1.0-2*dx)*rs*nrand()*dx;
				if(x<1.0+dx&&y<1.0+dx) nodes.push_back(vec2d(fmin(1.0+(gn>65)*dx,x),fmin(1.0,y)));
			}
		} else {
			FLOAT64 x = i*w;
			FLOAT64 y = j*dx;
			x += (x>dx&&j!=0&&j!=gn-1&&x<1.0-2*dx)*rs*nrand()*dx;;
			y += (x>dx&&j!=0&&y<1.0-2*dx)*rs*nrand()*dx;
			if(x<1.0+dx&&y<1.0+dx) nodes.push_back(vec2d(fmin(1.0+(gn>65)*dx,x),fmin(1.0,y)));
		}
	}
	nodes.push_back(vec2d(0.0,0.0));
	nodes.push_back(vec2d(0.0,1.0));
	
	// Create triangle elements
	vector<vec2d> triPoints;
	vector<uint> triIndices;
	for( uint i=0; i<nodes.size(); i++ ) {
		triPoints.push_back(nodes[i]);
		triIndices.push_back(i);
	}
	vector<vector<uint> > triangles = makeDelaunayTriangles(triPoints,gn);
	for( uint n=0; n<triangles.size(); n++ ) {
		vector<uint> elm(3);
		for( uint k=0; k<3; k++ ) elm[k] = triIndices[triangles[n][k]];
		elements.push_back(elm);
	}
}

static bool checkElementQuality( vec2d vertices[] ) {
	bool passed = true;
	FLOAT64 A[DIM+1][DIM+1];
	for( uint i=0; i<DIM+1; i++ ) for( uint j=0; j<DIM+1; j++ ) {
		if( i < DIM ) A[i][j] = vertices[j][i];
		else A[i][j] = 1.0;
	}
	FLOAT64 Ainv[DIM+1][DIM+1];
	if( invert3x3(A,Ainv) ) {
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
		if( ! passed ) {
			printf( "WARNING: Badly shaped element detected !\n" );
			printf( "BtB = \n" );
			printf( "[%.2f %.2f %.2f;\n", BtB[0][0], BtB[0][1], BtB[0][2] );
			printf( " %.2f %.2f %.2f;\n", BtB[1][0], BtB[1][1], BtB[1][2] );
			printf( " %.2f %.2f %.2f;\n", BtB[2][0], BtB[2][1], BtB[2][2] );
			printf( "-----------\n" );
		}
	} else {
		passed = false;
	}
	return passed;
}

void mesher2::updateConnection() {
	// Quality check for debug
#if 0
	for( uint n=0; n<elements.size(); n++ ) {
		vec2d vertices[DIM+1];
		for( uint m=0; m<elements[n].size(); m++ ) vertices[m] = nodes[elements[n][m]];
		if( ! checkElementQuality(vertices)) exit(0);
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
		if( ! invert3x3( A, matrix[n].m )) {
			printf( "Failed inverse 3x3!\n" );
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
	
	// Build facet marix based on the node - node matrix
	if( genFacet ) {
		facets.clear();
		facet_elements.clear();
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
		
		// Compute facet volume
		facetArea.resize(facets.size());
		for( uint n=0; n<facets.size(); n++ ) {
			FLOAT64 len = (nodes[facets[n][0]]-nodes[facets[n][1]]).len();
			facetArea[n] = len;
		}
		
		// Build element_elements
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
	element_edges.clear();
	edges.reserve(nodes.size());
	for( uint n=0; n<edges.size(); n++ ) edges[n].clear();
	element_edges.clear();
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
						}
					}
				}
			}
		}
	}
	
	// Compute the volume of elements
    volumes.resize(elements.size());
    for( uint n=0; n<elements.size(); n++ ) {
		volumes[n] = computeElementVolume(nodes,elements[n]);
	}
	
	// Build an uniform hash
	if( genHash ) {
		hash.clear();
		hash.resize(gn-1,gn-1);
		for( uint n=0; n<elements.size(); n++ ) {
			vector<vec2d> posarray(elements[n].size());
			for( uint m=0; m<elements[n].size(); m++ ) {
				FLOAT64 x = nodes[elements[n][m]][0];
				FLOAT64 y = nodes[elements[n][m]][1];
				posarray[m] = vec2d(x,y);
			}
			vec2d leftbottom;
			vec2d righttop;
			getBoundingBox(posarray,leftbottom,righttop);
			putHash(leftbottom,righttop,hash,n);
		}
	}
	
	// Compute center positions
	computeCenterPosition(nodes,elements,centers,centerType);
}

void mesher2::generateMesh( uint gn, const levelset2 *hint ) {
	this->gn = gn;
	
	// Generate nodes and matrix
	generateElements(nodes,elements,hint,gn);
	
	// Build element connections
	updateConnection();
}

void mesher2::computeCenterPosition( const std::vector<vec2d> &nodes, const std::vector<std::vector<uint> > &elements, std::vector<vec2d> &centers, int centerType ) {
	// Compute element center position
	centers.clear();
	centers.resize(elements.size());
	for( uint n=0; n<elements.size(); n++ ) {
		// Circumcenters
		if( centerType == BARYCENTRIC ) {
			// Barycenters
			vec2d barypos;
			for( uint k=0; k<elements[n].size(); k++ ) barypos += nodes[elements[n][k]]/elements[n].size();
			centers[n] = barypos;
		} else if( centerType == CIRCUMCENTRIC ) {
			vec2d pos;
			// Compute circumcenters directly
			if( util2::centerCircumcircle(nodes[elements[n][0]],nodes[elements[n][1]],nodes[elements[n][2]],pos)) {
				centers[n] = pos;
			}
		}
	}
}

bool mesher2::hitElement( uint index, vec2d p ) const {
	FLOAT64 x[NUM_VERT] = { p[0], p[1], 1.0 };
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

int mesher2::hitElements(vec2d p) const {
	int ref = -1;
	if( genHash ) {
		FLOAT64 e = 1e-8;
		uint gn = hash.size().w;
		uint i = gn*fmin(1.0-e,fmax(0.0,p[0]));
		uint j = gn*fmin(1.0-e,fmax(0.0,p[1]));
		for( uint n=0; n<hash[i][j].size(); n++ ) {
			if( hitElement(hash[i][j][n],p) ) {
				ref = hash[i][j][n];
				break;
			}
		}
	} else {
		printf( "Spatical hash generation is disabled.\n" );
	}
	return ref;
}

void mesher2::drawMesh() const {
	glBegin(GL_LINES);
	for( uint n=0; n<node2node.size(); n++ ) {
		for( uint m=0; m<node2node[n].size(); m++ ) {
			if( n < node2node[n][m] ) {
				vec2d p0 = nodes[n];
				vec2d p1 = nodes[node2node[n][m]];
				glVertex2f(p0[0],p0[1]);
				glVertex2f(p1[0],p1[1]);
			}
		}
	}
	glEnd();
}

void mesher2::drawElement2Facet() const {
	glBegin(GL_LINES);
	for( uint n=0; n<element_facets.size(); n++ ) {
		for( uint m=0; m<element_facets[n].size(); m++ ) {
			vec2d p0 = centers[n];
			vec2d p1 = 0.5*(nodes[facets[element_facets[n][m]][0]]+nodes[facets[element_facets[n][m]][1]]);
			glVertex2f(p0[0],p0[1]);
			glVertex2f(p1[0],p1[1]);
		}
	}
	glEnd();
}

void mesher2::drawCenters() const {
	glBegin(GL_POINTS);
	for( uint n=0; n<centers.size(); n++ ) {
		glVertex2f(centers[n][0],centers[n][1]);
	}
	glEnd();
}

void mesher2::drawLevelset( const std::vector<FLOAT64> &levelsets ) const {
	if( levelsets.size() != nodes.size() ) return;
	for( uint n=0; n<elements.size(); n++ ) {
		// Render levelset
		vector<vec2d> points(elements[n].size());
		vector<FLOAT64> LS(elements[n].size());
		for( uint i=0; i<elements[n].size(); i++ ) {
			points[i] = nodes[elements[n][i]];
			LS[i] = levelsets[elements[n][i]];
		}
		vector<vec2d> contourLines = util2::marchPoints(points,LS);					
		glBegin(GL_POLYGON);
		for( uint n=0; n<contourLines.size(); n++ ) {
			glVertex2f(contourLines[n][0],contourLines[n][1]);
		}
		glEnd();
	}
}

void mesher2::write_matlab( const char *name, const std::vector<FLOAT64> &scalar ) {
	FILE *fp = fopen(name,"w");
	if(fp) {
		// Write mesh config
		fprintf( fp, "tri = [\n" );
		for( uint n=0; n<elements.size(); n++ ) {
			for( uint m=0; m<elements[n].size(); m++ ) {
				fprintf(fp, "%d ", elements[n][m]+1 );
			}
			fprintf( fp, ";\n" );
		}
		fprintf( fp, "];\n" );
		fprintf( fp, "x = [\n" );
		for( uint n=0; n<nodes.size(); n++ ) {
			fprintf( fp, "%f;\n", nodes[n][0] );
		}
		fprintf( fp, "];\n" );
		fprintf( fp, "y = [\n" );
		for( uint n=0; n<nodes.size(); n++ ) {
			fprintf( fp, "%f;\n", nodes[n][1] );
		}
		fprintf( fp, "];\n" );
		fprintf( fp, "z = [\n" );
		for( uint n=0; n<nodes.size(); n++ ) {
			fprintf( fp, "%f;\n", scalar[n] );
		}
		fprintf( fp, "];\n" );
		fprintf( fp, "trisurf(tri,x,y,z);\n" );
	}
}

void mesher2::drawScalar( const std::vector<FLOAT64> &scalar, uint type, const std::vector<bool> &mask, FLOAT64 minv, FLOAT64 maxv ) const {
	FLOAT64 alpha = 1.0;
	if( type == NODAL ) {
		if( scalar.size() != nodes.size() ) return;
		if( maxv-minv == 0.0 ) {
			minv = 99999.0;
			maxv = -999999.0;
			for( uint n=0; n<nodes.size(); n++ ) {
				if( (mask.empty()||mask[n]) && scalar[n] < minv ) minv = scalar[n];
				if( (mask.empty()||mask[n]) && scalar[n] > maxv ) maxv = scalar[n];
			}
		}
		FLOAT64 det = maxv - minv;
		for( uint n=0; n<elements.size(); n++ ) {
			bool drawable = true;
#if 0
			for( uint m=0; m<elements[n].size(); m++ ) {
				if( ! mask[elements[n][m]] ) {
					drawable = false;
					break;
				}
			}
#endif
			if( drawable ) {
				glBegin(GL_POLYGON);
				for( uint m=0; m<elements[n].size(); m++ ) {
					FLOAT64 normp = 2.0*(scalar[elements[n][m]]-minv)/det-1.0;
					glColor4d(normp>0,0.3,normp<=0,alpha*fabs(normp)*(mask.empty()||mask[elements[n][m]]));
					vec2d p = nodes[elements[n][m]];
					glVertex2f(p[0],p[1]);
				}
				glEnd();
			}
		}
	} else if( type == ELEMENT ) {
		if( scalar.size() != elements.size() ) return;
		if( maxv-minv == 0.0 ) {
			minv = 99999.0;
			maxv = -999999.0;
			for( uint n=0; n<elements.size(); n++ ) {
				if( (mask.empty()||mask[n]) && scalar[n] < minv ) minv = scalar[n];
				if( (mask.empty()||mask[n]) && scalar[n] > maxv ) maxv = scalar[n];
			}
		}
		FLOAT64 det = maxv - minv;
		for( uint n=0; n<elements.size(); n++ ) {
			if(mask.empty()||mask[n]) {
				glBegin(GL_POLYGON);
				FLOAT64 normp = 2.0*(scalar[n]-minv)/det-1.0;
				glColor4d(normp>0,0.3,normp<=0,alpha*fabs(normp));
				for( uint m=0; m<elements[n].size(); m++ ) {
					vec2d p = nodes[elements[n][m]];
					glVertex2f(p[0],p[1]);
				}
				glEnd();
			}
		}
	}
}

void mesher2::drawVector( const std::vector<vec2d> &vecs, uint type, FLOAT64 scale ) const {
	if( type == NODAL ) {
		if( vecs.size() == nodes.size() ) {
			glBegin(GL_LINES);
			for( uint n=0; n<nodes.size(); n++ ) {
				vec2d p0 = nodes[n];
				vec2d p1 = nodes[n]+scale*vecs[n];
				glVertex2f(p0[0],p0[1]);
				glVertex2f(p1[0],p1[1]);
			}
			glEnd();
		}
	} else if( type == ELEMENT ) {
		if( vecs.size() == elements.size() ) {
			glBegin(GL_LINES);
			for( uint n=0; n<elements.size(); n++ ) {
				vec2d p0 = centers[n];
				vec2d p1 = centers[n]+scale*vecs[n];
				glVertex2dv(p0.v);
				glVertex2dv(p1.v);
			}
			glEnd();
		}
	}
}
