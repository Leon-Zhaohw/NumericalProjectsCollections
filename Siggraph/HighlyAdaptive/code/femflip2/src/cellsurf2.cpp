/*
 *	cellsurf2.cpp
 *	
 *	Created by Ryoichi Ando on 2/11/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "cellsurf2.h"
#include "levelset2.h"
#include "opengl.h"
#include "util2.h"
#include "kernel.h"
#include "fastmarch2.h"
using namespace std;

cellsurf2::cellsurf2() {
	gn = 0;
	dx = 0.0;
	extrapolate_dist = 1.0;
}

cellsurf2::cellsurf2( const cellsurf2 &cellsurf ) {
	*this = cellsurf;
}

void cellsurf2::setResolution( uint gn ) {
	this->gn = gn;
}

void cellsurf2::setExtrapolationDist( FLOAT64 dist ) {
	extrapolate_dist = dist;
}

static void putHash( vec2d left_bottom, vec2d right_top, array2<std::vector<uint> > &hash, uint id ) {
	FLOAT64 e = 1.0e-12;
	uint gn[2] = { hash.size().w, hash.size().h };
	uint min_i = gn[0]*fmin(1.0-e,fmax(0.0,left_bottom[0]+e));
	uint min_j = gn[1]*fmin(1.0-e,fmax(0.0,left_bottom[1]+e));
	uint max_i = gn[0]*fmin(1.0-e,fmax(0.0,right_top[0]-e));
	uint max_j = gn[1]*fmin(1.0-e,fmax(0.0,right_top[1]-e));
	
	max_i = imin(max_i+1,gn[0]);
	max_j = imin(max_j+1,gn[1]);
	
#if 0
	if( max_i-min_i != 1 ) {
		printf( "cellsurf2.cpp: x %d %d\n", min_i, max_i );
	}
	if( max_j-min_j != 1 ) {
		printf( "cellsurf2.cpp: y %d %d\n", min_j, max_j );
	}
#endif
	
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

static FLOAT64 distance( const vec2d &p0, const vec2d &p1, vec2d &p ) {
	vec2d a = p1-p0;
	vec2d b = p-p0;
	FLOAT64 alen = a.len();
	if( alen ) {
		FLOAT64 dot = a/alen * b;
		dot = fmin(alen,fmax(0.0,dot));
		vec2d o = p;
		p = p0 + (a/alen) * dot;
		return (p-o).len();
	} else {
		p = p0;
		return a.len();
	}
}

void cellsurf2::buildSurface( std::vector<vec2d> &vertices, std::vector<vec2d> &normals, std::vector<vec2i> &faces, const levelset2 *fluid, const levelset2 *solid, bool enclose ) {
	dx = 1.0/gn;
	array2<int> vindex[2];
	vindex[0].resize(gn+1,gn);
	vindex[1].resize(gn,gn+1);
	
	// Compute nodal levelset
	nodalLevelset.resize(gn+1,gn+1);
	nodalLevelset.clear();
	FOR_EACH(gn+1,gn+1) {
		vec2d p(dx*i,dx*j);
		FLOAT64 phi = fluid->evalLevelset(p);
		if( fabs(phi) < 1e-2*fluid->dx ) phi = 0.0;
		nodalLevelset[i][j] = phi;
	} END_FOR
	
	// Compute X face vertices
	vertices.clear();
	FOR_EACH(gn+1,gn) {
		FLOAT64 LS[2] = { nodalLevelset[i][j], nodalLevelset[i][j+1] };
		vindex[0][i][j] = -1;
		if( copysign(1.0,LS[0]) * copysign(1.0,LS[1]) <= 0 ) {
			FLOAT64 a = LS[0]/(LS[0]-LS[1]);
			vec2d pos = (1.0-a)*vec2d(i*dx,j*dx)+a*vec2d(i*dx,(j+1)*dx);
			vindex[0][i][j] = vertices.size();
			vertices.push_back(pos);
		}
	} END_FOR
	
	// Compute Y face vertices
	FOR_EACH(gn,gn+1) {
		FLOAT64 LS[2] = { nodalLevelset[i][j], nodalLevelset[i+1][j] };
		vindex[1][i][j] = -1;
		if( copysign(1.0,LS[0]) * copysign(1.0,LS[1]) <= 0 ) {
			FLOAT64 a = LS[0]/(LS[0]-LS[1]);
			vec2d pos = (1.0-a)*vec2d(i*dx,j*dx)+a*vec2d((i+1)*dx,j*dx);
			vindex[1][i][j] = vertices.size();
			vertices.push_back(pos);
		}
	} END_FOR
	
	// Now connect them with edges
	faces.clear();
	FOR_EACH(gn,gn) {
		int f[4] = {vindex[0][i][j],vindex[1][i][j+1],vindex[0][i+1][j],vindex[1][i][j]};
		int num_v = 0;
		for( uint n=0; n<4; n++ ) if(f[n]>=0) num_v++;
		if( num_v == 2 ) {
			vector<uint> edge(2);
			int idx=0;
			for( uint n=0; n<4; n++ ) {
				if(f[n]>=0) {
					edge[idx++] = f[n];
				}
			}
			faces.push_back(vec2i(edge[0],edge[1]));
		} else if( num_v == 4 ) {
			vector<uint> edge(2);
			int idx=0;
			int s = nodalLevelset[i][j] < 0 ? 3 : 0;
			for( uint n=0; n<4; n++ ) {
				if(f[(s+n)%4]>=0) {
					edge[idx++] = f[(s+n)%4];
					if( idx == 2 ) {
						faces.push_back(vec2i(edge[0],edge[1]));
						idx=0;
					}
				}
			}
		}
	} END_FOR
	
	// Smooth surface
	smoothMesh(vertices,normals,faces,1);
	
	// Build levelset structure
	buildLevelsetInfo(vertices,normals,faces);
}

static bool is_nan( vec2d p ) {
	return is_nan(p[0]) || is_nan(p[1]) || is_nan(p[2]);
}

void cellsurf2::smoothMesh( std::vector<vec2d> &vertices, std::vector<vec2d> &normals, std::vector<vec2i> &faces, uint iterations ) {
	std::vector<vec2d> old_vertices = vertices;
	std::vector<vec2d> old_normals = normals;
	std::vector<std::vector<uint> > connections;
	std::vector<FLOAT64> areas;
	connections.resize(vertices.size());
	areas.resize(vertices.size());
	for( uint n=0; n<vertices.size(); n++ ) {
		areas[n] = 0.0;
	}
	for( uint n=0; n<faces.size(); n++ ) {
		connections[faces[n][0]].push_back(faces[n][1]);
		connections[faces[n][1]].push_back(faces[n][0]);
		FLOAT64 area = (vertices[faces[n][1]]-vertices[faces[n][0]]).len();
		for( uint j=0; j<2; j++ ) {
			areas[faces[n][j]] += area / 2.0;
		}
	}
	
	// Smooth out
	for( uint k=0; k<iterations; k++ ) {
		for( uint n=0; n<vertices.size(); n++ ) {
			bool edge = false;
			for( uint dim=0; dim<DIM; dim++) {
				if( vertices[n][dim] <= 1e-4 || vertices[n][dim] >= 1.0-1e-4 ) edge = true;
			}
			if( edge ) continue;
			
			FLOAT64 wsum = 0.0;
			vec2d pos;
			// Add self
			FLOAT64 w = areas[n];
			if( ! is_nan(old_vertices[n]) ) {
				pos += w*old_vertices[n];
				wsum += w;
			}
			// Add neighbors
			for( uint m=0; m<connections[n].size(); m++ ) {
				uint idx = connections[n][m];
				if( ! is_nan(old_vertices[idx]) ) {
					FLOAT64 w = areas[idx];
					pos += w*old_vertices[idx];
					wsum += w;
				}
			}
			if( wsum ) vertices[n] = pos / wsum;
		}
		old_vertices = vertices;
#if 0
		for( uint n=0; n<normals.size(); n++ ) {
			bool edge = false;
			for( uint dim=0; dim<DIM; dim++) {
				if( vertices[n][dim] <= 1e-4 || vertices[n][dim] >= 1.0-1e-4 ) edge = true;
			}
			if( edge ) continue;
			
			FLOAT64 wsum = 0.0;
			vec2d normal;
			// Add self
			FLOAT64 w = areas[n];
			if( ! is_nan(old_normals[n]) ) {
				normal += w*old_normals[n];
				wsum += w;
			}
			// Add neighbors
			for( uint m=0; m<connections[n].size(); m++ ) {
				uint idx = connections[n][m];
				if( ! is_nan(old_normals[idx]) ) {
					FLOAT64 w = areas[idx];
					normal += w*old_normals[idx];
					wsum += w;
				}
			}
			if( wsum ) normals[n] = normal.normal();
		}
		old_normals = normals;
#endif
	}
}

void cellsurf2::buildLevelsetInfo( std::vector<vec2d> &vertices, std::vector<vec2d> &normals, std::vector<vec2i> &faces ) {
	// Build spatial hash
	hash.resize(gn,gn);
	hash.clear();
	for( uint n=0; n<faces.size(); n++ ) {
		vec2d left_bottom;
		vec2d right_top;
		vector<vec2d> points(2);
		points[0] = vertices[faces[n][0]];
		points[1] = vertices[faces[n][1]];
		getBoundingBox(points,left_bottom,right_top);
		putHash(left_bottom,right_top,hash,n);
	}
	
	vhash.resize(gn,gn);
	vhash.clear();
	for( uint n=0; n<vertices.size(); n++ ) {
		int i = imax(0,imin(gn-1,vertices[n][0]*gn));
		int j = imax(0,imin(gn-1,vertices[n][1]*gn));
		vhash[i][j].push_back(n);
	}
	
	// Build vertex to vertex info
	v2v.clear();
	v2v.resize(vertices.size());
	for( uint n=0; n<faces.size(); n++ ) {
		for( uint m1=0; m1<DIM; m1++ ) for( uint m2=0; m2<DIM; m2++ ) {
			if( m1 != m2 ) v2v[faces[n][m1]].push_back(faces[n][m2]);
		}
	}
	// Gather narrow bands
	vector<vec2i> narrowBands;
	FOR_EACH(gn+1,gn+1) {
		bool found = false;
		for( int ki=0; ki<2; ki++ ) for( int kj=0; kj<2; kj++ ) {
			if( found ) continue;
			int ni = i-1+ki;
			int nj = j-1+kj;
			if( ni<0 || nj<0 || ni>=gn || nj>=gn ) continue;
			if( hash[ni][nj].size() ) {
				narrowBands.push_back(vec2i(i,j));
				found = true;
			}
		}
		nodalLevelset[i][j] = (nodalLevelset[i][j] >= 0.0 ? 1.0 : -1.0)*1e8;
	} END_FOR
	
	// Fix the narrow bands
	referencePos.resize(gn+1,gn+1);
	referencePos.clear();
	PARALLEL_FOR for( uint n=0; n<narrowBands.size(); n++ ) {
		uint i = narrowBands[n][0];
		uint j = narrowBands[n][1];
		vec2d p(dx*i,dx*j);
		FLOAT64 dist = findClosestSurfPosition(vertices,faces,i,j,p);
		FLOAT64 sgn = nodalLevelset[i][j] >= 0.0 ? 1.0 : -1.0;
		nodalLevelset[i][j] = dist*sgn;
		referencePos[i][j] = p;
	}
	
	if( extrapolate_dist ) {
		// Setup a network
		FLOAT64 dx = 1.0/gn;
		array2<fastmarch2<FLOAT64>::node2> nodeArray;
		nodeArray.resize(gn+1,gn+1);
		PARALLEL_FOR FOR_EACH(gn+1,gn+1) {
			vec2d p = vec2d(i*dx,j*dx);
			nodeArray[i][j].p = p;
			nodeArray[i][j].fixed = fabs(nodalLevelset[i][j])<1.0;
			nodeArray[i][j].levelset = nodalLevelset[i][j];
			nodeArray[i][j].value = 0.0;
			int q[][DIM] = {{i-1,j},{i+1,j},{i,j-1},{i,j+1}};
			for( uint n=0; n<4; n++ ) {
				int ni = q[n][0];
				int nj = q[n][1];
				if( ni>=0 && ni<gn+1 && nj>=0 && nj<gn+1 ) {
					nodeArray[i][j].p2p.push_back(&nodeArray[ni][nj]);
				}
			}
		} END_FOR
		std::vector<fastmarch2<FLOAT64>::node2 *> nodes((gn+1)*(gn+1));
		uint index = 0;
		FOR_EACH(gn+1,gn+1) {
			nodes[index++] = &nodeArray[i][j];
		} END_FOR
		
		// Perform fast march
		fastmarch2<FLOAT64>::fastMarch(nodes,extrapolate_dist,-extrapolate_dist,1);
		
		// Pickup extrapolated levelsets
		PARALLEL_FOR FOR_EACH(gn+1,gn+1) {
			nodalLevelset[i][j] = nodeArray[i][j].levelset;
		} END_FOR
	}
	
	// Now evaluate normals
	normals.resize(vertices.size());
	for( uint n=0; n<vertices.size(); n++ ) {
		normals[n] = -1.0 * evalGradient(vertices[n]);
	}
	
	// Build vertex to vertex info
	std::vector<std::vector<uint> > v2v(vertices.size());
	for( uint n=0; n<faces.size(); n++ ) {
		for( uint m1=0; m1<DIM; m1++ ) for( uint m2=0; m2<DIM; m2++ ) {
			if( m1 != m2 ) v2v[faces[n][m1]].push_back(faces[n][m2]);
		}
	}
}


FLOAT64 cellsurf2::evalLevelset(vec2d p) const {
	return util2::interp<FLOAT64>(gn*p,nodalLevelset);
}

vec2d cellsurf2::evalGradient(vec2d p) const {
	uint w, h;
	w = nodalLevelset.size().w;
	h = nodalLevelset.size().h;
	const array2<FLOAT64> &q = nodalLevelset;
	FLOAT64 x = fmax(0.0,fmin(w,gn*p[0]));
	FLOAT64 y = fmax(0.0,fmin(h,gn*p[1]));
	int i = imin(x,w-2);
	int j = imin(y,h-2);
	FLOAT64 gx = (q[i+1][j]-q[i][j])*(j+1-y)+(q[i+1][j+1]-q[i][j+1])*(y-j);
	FLOAT64 gy = (q[i][j+1]-q[i][j])*(i+1-x)+(q[i+1][j+1]-q[i+1][j])*(x-i);
	return vec2d(gx,gy).normal();
}

FLOAT64 cellsurf2::findClosestSurfPosition( std::vector<vec2d> &vertices, std::vector<vec2i> &faces, int nodal_i, int nodal_j, vec2d &p ) {
	FLOAT64 dist = 1e8;
	vec2d dest = p;
	FOR_EACH(2,2) {
		int ni = nodal_i-1+i;
		int nj = nodal_j-1+j;
		if( ni<0 || nj<0 || ni>=gn || nj>=gn ) continue;
		int size = hash[ni][nj].size();
		for( int n=0; n<size; n++ ) {
			uint id = hash[ni][nj][n];
			const vec2d &p0 = vertices[faces[id][0]];
			const vec2d &p1 = vertices[faces[id][1]];
			vec2d crossp = p;
			FLOAT64 d = distance(p0,p1,crossp);
			if( d < dist ) {
				dist = d;
				dest = crossp;
			}
		}
	} END_FOR
	p = dest;
	return dist;
}

void cellsurf2::drawGL( const levelset2 *solid, uint kind ) const {
	// Fill this levelset
	if( kind == 0 ) {
		glColor4d(0.3,0.3,0.8,0.6);
		FOR_EACH(gn-1,gn-1) {
			uint q[4][2] = {{i,j},{i+1,j},{i+1,j+1},{i,j+1}};
#if 1
			vector<vec2d> nodes(4);
			vector<FLOAT64> levelsets(4);
			for( uint k=0; k<4; k++ ) {
				nodes[k] = vec2d(dx*q[k][0],dx*q[k][1]);
				levelsets[k] = fmax(nodalLevelset[q[k][0]][q[k][1]],solid ? -solid->evalLevelset(nodes[k]) : -1e9 );
			}
			vector<vec2d> points = util2::marchPoints(nodes,levelsets);
			glBegin(GL_POLYGON);
			for( int m=0; m<points.size(); m++ ) glVertex2f(points[m][0],points[m][1]);
			glEnd();
#else
			glBegin(GL_QUADS);
			for( uint k=0; k<4; k++ ) {
				FLOAT64 levelset = nodalLevelset[q[k][0]][q[k][1]];
				glColor4d(levelset>0,0.0,levelset<0.0,7.0*fabs(levelset));
				glVertex2d(dx*q[k][0],dx*q[k][1]);
			}
			glEnd();
#endif
		} END_FOR
	}
#if 0
	// Draw points
	glColor4d(1.0,1.0,1.0,1.0);
	glPointSize(3.0);
	glBegin(GL_POINTS);
	for( uint n=0; n<vertices.size(); n++ ) {
		glVertex2dv(vertices[n].v);
	}
	glEnd();
	glPointSize(1.0);
#endif
	
#if 0
	// Turn on here to see closest positions
	glBegin(GL_faces);
	FOR_EACH(gn+1,gn+1) {
		FLOAT64 levelset = nodalLevelset[i][j];
		if( referencePos[i][j][0] ) {
			glColor4d(levelset>0,0.5,levelset<=0.0,1.0);
			vec2d p(dx*i,dx*j);
			glVertex2dv(p.v);
			glVertex2dv(referencePos[i][j].v);
		}
	} END_FOR
	glEnd();
#endif
}
