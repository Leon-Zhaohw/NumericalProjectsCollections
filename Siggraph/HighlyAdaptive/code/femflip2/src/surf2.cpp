/*
 *	surf2.h
 *
 *	Created by Ryoichi Ando on 8/19/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "surf2.h"
#include "fluid2.h"
#include "opengl.h"
#include <stdlib.h>

surf2::surf2() {
	surf = NULL;
}

surf2::surf2( const surf2 &surf ) {
	*this = surf;
}

void surf2::setResolution( uint gm ) {
	this->gm = gm;
	cellsurf.setResolution(gm);
	bccsurf.dx = 1.0/gm;
	bccsurf.meshsurf.dx = 1.0/gm;
	meshsurf.dx = 1.0/gm;
}

void surf2::setExtrapolationDist( FLOAT64 dist ) {
	bccsurf.setExtrapolationDist(dist);
	meshsurf.setExtrapolationDist(dist);
	cellsurf.setExtrapolationDist(dist);
}

void surf2::buildSurface(const levelset2 *fluid, const levelset2 *solid, fluid2 *fluidSolver, bool enclose ) {
	if( fluidSolver->getBCCMesh() ) {
		surf = &bccsurf;
		bccsurf.setBCC(*fluidSolver->getBCCMesh());
		bccsurf.buildSurface(vertices,normals,faces,fluid,solid,enclose);
	} else if( fluidSolver->getMesh() ) {
		surf = &meshsurf;
		mesh = *fluidSolver->getMesh();
		meshsurf.setReference(&mesh);
		meshsurf.buildSurface(vertices,normals,faces,fluid,solid,enclose);
	} else {
		surf = &cellsurf;
		cellsurf.buildSurface(vertices,normals,faces,fluid,solid,enclose);
	}
	
	// Build ANN
	ann.sort(vertices);
	
	// Compute vertices to face info
	v2f.clear();
	v2f.resize(vertices.size());
	for( uint n=0; n<faces.size(); n++ ) {
		for( uint m=0; m<DIM; m++ ) {
			v2f[faces[n][m]].push_back(n);
		}
	}
	
	// Compute curvature
	computeCurvature(solid);
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

FLOAT64 surf2::evalLevelset(vec2d p) const {
	if( ! surf ) {
		printf( "Surface not initialized !\n");
		exit(0);
	}
	return surf->evalLevelset(p);
}

vec2d surf2::evalGradient(vec2d p) const {
	if( ! surf ) {
		printf( "Surface not initialized !\n");
		exit(0);
	}
	return surf->evalGradient(p);
}

FLOAT64 surf2::evalCurvature(vec2d p) const {
	if( ! surf ) {
		printf( "Surface not initialized !\n");
		exit(0);
	}
	// Pick up from surfaces by smooth kernel balls
	std::vector<ANNidx> neighbors = ann.getkNeighbors(p,3);
	FLOAT64 sum = 0.0;
	FLOAT64 wsum = 0.0;
	for( uint n=0; n<neighbors.size(); n++ ) {
		FLOAT64 w = areas[neighbors[n]];
		sum += w*curvature[neighbors[n]];
		wsum += w;
	}
	if( wsum ) return sum/wsum;
	return 0.0;
}

void surf2::drawGL( const levelset2 *solid, uint kind ) const {
	if( surf == &bccsurf ) bccsurf.drawGL(solid,kind);
	else if( surf == &meshsurf ) meshsurf.drawGL(solid,kind);
	else if( surf == &cellsurf ) cellsurf.drawGL(solid,kind);
	
	if( kind == 1 ) {
		// Draw outlines
		glColor4d(1.0,1.0,1.0,1.0);
		glBegin(GL_LINES);
		for( uint n=0; n<faces.size(); n++ ) {
			FLOAT64 len = (vertices[faces[n][0]]-vertices[faces[n][1]]).len();
			if( ! solid || (solid->evalLevelset(vertices[faces[n][0]]) > -len &&
				solid->evalLevelset(vertices[faces[n][1]]) > -len ))  {
				   glVertex2dv(vertices[faces[n][0]].v);
				   glVertex2dv(vertices[faces[n][1]].v);
			}
		}
		glEnd();
	}
}

void surf2::operator=( const surf2 &surf ) {
	gm = surf.gm;
	if( surf.surf == &surf.bccsurf ) {
		this->surf = &bccsurf;
		bccsurf = surf.bccsurf;
	} else if( surf.surf == &surf.meshsurf ) {
		this->mesh = surf.mesh;
		this->surf = &meshsurf;
		meshsurf = surf.meshsurf;
		meshsurf.setReference(&mesh);
	} else if( surf.surf == &surf.cellsurf ) {
		this->surf = &cellsurf;
		cellsurf = surf.cellsurf;
	}
}

void surf2::computeCurvature( const levelset2 *solid ) {
	// Build vertex to vertex info
	std::vector<std::vector<uint> > v2v(vertices.size());
	for( uint n=0; n<faces.size(); n++ ) {
		for( uint m1=0; m1<DIM; m1++ ) for( uint m2=0; m2<DIM; m2++ ) {
			if( m1 != m2 ) v2v[faces[n][m1]].push_back(faces[n][m2]);
		}
	}
	
	// Compute curvature at every vertex
	curvature.resize(vertices.size());
	areas.resize(vertices.size());
	for( uint n=0; n<vertices.size(); n++ ) {
		curvature[n] = 0.0;
		areas[n] = 0.0;
		if( solid->evalLevelset(vertices[n]) > 0.0 ) {
			// Get two connected opposite vertex
			vec2d v0 = vertices[n];
			vec2d v1 = vertices[v2v[n][0]];
			vec2d v2 = vertices[v2v[n][1]];
			FLOAT64 len0 = (v1-v0).len()+1e-18;
			FLOAT64 len1 = (v2-v0).len()+1e-18;
			FLOAT64 area = 0.5*(len0+len1);
			// Stretch
			v1 = v0+(v1-v0)/len0;
			v2 = v0+(v2-v0)/len1;
			vec2d p = v0;
			FLOAT64 d = distance(v1,v2,p);
			curvature[n] = d/area*(evalGradient(v0)*(p-v0)>=0.0?1.0:-1.0);
			areas[n] = area;
		}
	}
	
#if 1
	// Smoothing part
	for( uint k=0; k<2; k++ ) {
		std::vector<FLOAT64> old_curvature = curvature;
		PARALLEL_FOR for( uint n=0; n<vertices.size(); n++ ) {
			FLOAT64 new_curvature = 0.0;
			FLOAT64 wsum = 0.0;
			for( uint m=0; m<v2v[n].size(); m++ ) {
				FLOAT64 w = areas[v2v[n][m]];
				wsum += w;
				new_curvature += w*old_curvature[v2v[n][m]];
			}
			if( wsum ) {
				curvature[n] = new_curvature / wsum;
			}
		}
	}
#endif
}
