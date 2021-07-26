/*
 *	surf3.cpp
 *
 *	Created by Ryoichi Ando on 2012/08/04
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */
 
#include "util3.h"
#include "surf3.h"
#include "fluid3.h"
#include "opengl.h"
#include <stdlib.h>

surf3::surf3() {
	surf = NULL;
}

surf3::surf3( const surf3 &surf ) {
	*this = surf;
}

void surf3::setResolution( uint gm ) {
	this->gm = gm;
	cellsurf.setResolution(gm);
}

void surf3::setExtrapolationDist( FLOAT64 dist ) {
	bccsurf.setExtrapolationDist(dist);
	meshsurf.setExtrapolationDist(dist);
	cellsurf.setExtrapolationDist(dist);
}

void surf3::buildSurface(const fluid3 *fluidSolver, const levelset3 *fluid, const levelset3 *solid, bool enclose ) {
	tick(); dump(">>> Building surfaces started...\n");
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
	tick(); dump("Building ANN structure for mesh vertices...");
	ann.sort(vertices);
	dump("Done. Took %s.\n", stock("surf_ANN"));
	
	// Compute vertices to face info
	v2f.clear();
	v2f.resize(vertices.size());
	for( uint n=0; n<faces.size(); n++ ) {
		for( uint m=0; m<DIM; m++ ) {
			v2f[faces[n][m]].push_back(n);
		}
	}
	
	dump("<<< Done. Took %s.\n",stock("surf_build"));
	
	tick(); dump("Computing mesh-based curvature (%d vertices %d faces)...", vertices.size(), faces.size());
	computeCurvature(solid);
	dump("Done. Took %s.\n", stock("surf_curvature"));
	
	// Write some info
	writeNumber("surf_vertex_num", vertices.size());
	writeNumber("surf_face_num", faces.size());
}

// Thank you Christopher Batty !!
// Code picked up from:
// https://github.com/christopherbatty/SDFGen
//

// find distance x0 is from segment x1-x2
static FLOAT64 point_segment_distance(const vec3d &x0, const vec3d &x1, const vec3d &x2, vec3d &out ) {
	vec3d dx(x2-x1);
	FLOAT64 m2=dx.len2();
	if( ! m2 ) {
		out = x1;
		return (x0-out).len();
	}
	// find parameter value of closest point on segment
	FLOAT64 s12=(x2-x0)*dx/m2;
	if(s12<0){
		s12=0;
	}else if(s12>1){
		s12=1;
	}
	// and find the distance
	out = s12*x1+(1-s12)*x2;
	return (x0-out).len();
}

// find distance x0 is from triangle x1-x2-x3
static FLOAT64 point_triangle_distance(const vec3d &x0, const vec3d &x1, const vec3d &x2, const vec3d &x3, vec3d &out ) {
	// first find barycentric coordinates of closest point on infinite plane
	vec3d x13(x1-x3), x23(x2-x3), x03(x0-x3);
	FLOAT64 m13=x13.len2(), m23=x23.len2(), d=x13*x23;
	FLOAT64 invdet=1.0/fmax(m13*m23-d*d,1e-30);
	FLOAT64 a=x13*x03, b=x23*x03;
	// the barycentric coordinates themselves
	FLOAT64 w23=invdet*(m23*a-d*b);
	FLOAT64 w31=invdet*(m13*b-d*a);
	FLOAT64 w12=1-w23-w31;
	if(w23>=0 && w31>=0 && w12>=0){ // if we're inside the triangle
		out = w23*x1+w31*x2+w12*x3;
		return (x0-out).len();
	}else{ // we have to clamp to one of the edges
		if(w23>0) // this rules out edge 2-3 for us
			return fmin(point_segment_distance(x0,x1,x2,out), point_segment_distance(x0,x1,x3,out));
		else if(w31>0) // this rules out edge 1-3
			return fmin(point_segment_distance(x0,x1,x2,out), point_segment_distance(x0,x2,x3,out));
		else // w12 must be >0, ruling out edge 1-2
			return fmin(point_segment_distance(x0,x1,x3,out), point_segment_distance(x0,x2,x3,out));
	}
}

FLOAT64 surf3::evalLevelset(vec3d p) const {
	if( ! surf ) {
		dump( "Surface not initialized !\n");
		exit(0);
	}
	return surf->evalLevelset(p);
}

vec3d surf3::evalGradient(vec3d p) const {
	if( ! surf ) {
		dump( "Surface not initialized !\n");
		exit(0);
	}
	return surf->evalGradient(p);
}

void surf3::drawGL( const levelset3 *solid, uint kind ) const {
	for( int i=0; i<faces.size(); i++ ) {
		if( kind == 0 ) glBegin(GL_LINES);
		else glBegin(GL_POLYGON);
		bool clear = true;
		for( int idx=0; idx<3; idx++ ) {
			int face = faces[i][idx];
			if( solid->evalLevelset(vertices[face]) < -0.5*dx ) {
				clear = false;
				break;
			}
		}
		if( clear ) {
			for( int idx=0; idx<3; idx++ ) {
				int face = faces[i][idx];
				glVertex3dv(vertices[face].v);
			}
		}
		glEnd();
	}
}

void surf3::operator=( const surf3 &surf ) {
	tick(); dump(">>> Copying surface mesh...\n");
	gm = surf.gm;
	if( surf.surf == &surf.bccsurf ) {
		this->surf = &bccsurf;
		bccsurf = surf.bccsurf;
	} else if( surf.surf == &surf.meshsurf ) {
		this->mesh = surf.mesh;
		this->surf = &meshsurf;
		tick(); dump("Copying mesh surface...");
		meshsurf = surf.meshsurf;
		meshsurf.setReference(&mesh);
		dump("Done. Took %s.\n",stock());
	} else if( surf.surf == &surf.cellsurf ) {
		this->surf = &cellsurf;
		cellsurf = surf.cellsurf;
	}
	dump("<<< Done. Took %s.\n",stock());
}

FLOAT64 surf3::evalCurvature(vec3d p) const {
	if( ! surf ) {
		printf( "Surface not initialized !\n");
		exit(0);
	}
	// Pick up from surfaces by smooth kernel balls
	std::vector<ANNidx> neighbors = ann.getkNeighbors(p,4);
	FLOAT64 sum = 0.0;
	FLOAT64 wsum = 0.0;
	for( uint n=0; n<neighbors.size(); n++ ) {
		FLOAT64 w = ringArea[neighbors[n]];
		sum += w*curvature[neighbors[n]].len();
		wsum += w;
	}
	if( wsum ) return sum/wsum;
	return 0.0;
}

static FLOAT64 vec_angle( vec3d v1, vec3d v2 ) {
	FLOAT64 len1 = v1.len();
	FLOAT64 len2 = v2.len();
	FLOAT64 dot = len1 * len2;
	if( dot ) {
		FLOAT64 res = acos(fmin(1.0-1e-8,fmax(-1.0+1e-8,v1*v2/dot)));
		return res;
	}
	return 0.0;
}

void surf3::computeCurvature( const levelset3 *solid ) {
	// Build vertex to vertex info
	std::vector<std::vector<uint> > v2v(vertices.size());
	for( uint n=0; n<faces.size(); n++ ) {
		for( uint m1=0; m1<DIM; m1++ ) for( uint m2=0; m2<DIM; m2++ ) {
			if( m1 != m2 ) {
				bool duplicated = false;
				for( uint i=0; i<v2v[faces[n][m1]].size(); i++ ) {
					if( v2v[faces[n][m1]][i] == faces[n][m2] ) {
						duplicated = true;
						break;
					}
				}
				if( ! duplicated ) v2v[faces[n][m1]].push_back(faces[n][m2]);
			}
		}
	}
	
	// Build vertex to face info
	std::vector<std::vector<uint> > v2faces(vertices.size());
	for( uint n=0; n<faces.size(); n++ ) {
		for( uint m=0; m<DIM; m++ ) {
			v2faces[faces[n][m]].push_back(n);
		}
	}
	
	// Build segments and vertex to segment info
	std::vector<std::vector<uint> > segments;
	std::vector<std::vector<uint> > v2segs(vertices.size());
	for( uint n=0; n<v2v.size(); n++ ) {
		for( uint m=0; m<v2v[n].size(); m++ ) {
			if( n < v2v[n][m] ) {
				std::vector<uint> seg(2);
				seg[0] = n;
				seg[1] = v2v[n][m];
				segments.push_back(seg);
				v2segs[seg[0]].push_back(segments.size()-1);
				v2segs[seg[1]].push_back(segments.size()-1);
			}
		}
	}
	
	// Build seg2face matrix
	std::vector<std::vector<uint> > seg2faces(segments.size());
	for( uint n=0; n<segments.size(); n++ ) {
		uint idx0 = segments[n][0];
		uint idx1 = segments[n][1];
		for( uint m=0; m<v2faces[idx0].size(); m++ ) {
			for( int k=DIM-1; k>=0; k-- ) {
				if( faces[v2faces[idx0][m]][k] == idx1 ) {
					seg2faces[n].push_back(v2faces[idx0][m]);
					break;
				}
			}
		}
	}
	
	// Now begin computing the curvature
	curvature.resize(vertices.size());
	PARALLEL_FOR for( uint n=0; n<vertices.size(); n++ ) {
		vec3d kn;
		for( uint m=0; m<v2segs[n].size(); m++ ) {
			vec3d x0, x1;
			uint seg = v2segs[n][m];
			x0 = vertices[n];
			for( uint i=0; i<2; i++ ) if( segments[seg][i] != n ) {
				x1 = vertices[segments[seg][i]];
			}
			for( uint f=0; f<seg2faces[seg].size(); f++ ) {
				uint face = seg2faces[seg][f];
				// Find opposite vertex
				int opv = -1;
				for( uint i=0; i<3; i++ ) {
					if( faces[face][i] != segments[seg][0] && faces[face][i] != segments[seg][1] ) {
						opv = i;
						break;
					}
				}
				if( opv >= 0 ) {
					vec3d v1 = vertices[faces[face][(opv+1)%3]]-vertices[faces[face][opv]];
					vec3d v2 = vertices[faces[face][(opv+2)%3]]-vertices[faces[face][opv]];
					FLOAT64 angle = vec_angle(v1,v2);
					FLOAT64 k = cos(angle)/(1e-16+sin(angle));
					kn += k*(x1-x0);
				}
			}
		}
		curvature[n] = kn;
	}
	
	ringArea.resize(vertices.size());
	PARALLEL_FOR for( uint n=0; n<vertices.size(); n++ ) {
		FLOAT64 area_sum = 0.0;
		for( uint m=0; m<v2faces[n].size(); m++ ) {
			uint face = v2faces[n][m];
			vec3d v1 = vertices[faces[face][1]]-vertices[faces[face][0]];
			vec3d v2 = vertices[faces[face][2]]-vertices[faces[face][0]];
			FLOAT64 area = 0.5 * (v1 ^ v2).len();
			if( ! area ) vertices[faces[face][0]] += 1e-8*vec3d(nrand(),nrand(),nrand());
			area_sum += area;
		}
		area_sum = fmax(1e-18,area_sum);
		ringArea[n] = area_sum/3.0;
		curvature[n] = curvature[n] / ringArea[n];
	}
	
	// Set curvature zero inside the solid
	if( solid ) {
		for( uint n=0; n<vertices.size(); n++ ) {
			if( solid->evalLevelset(vertices[n]) < 0.0 ) {
				curvature[n] = vec3d();
				ringArea[n] = 0.0;
			}
		}
	}
	
#if 1
	// Smoothing part
	for( uint k=0; k<2; k++ ) {
		std::vector<vec3d> old_curvature = curvature;
		PARALLEL_FOR for( uint n=0; n<vertices.size(); n++ ) {
			vec3d new_curvature;
			FLOAT64 wsum = 0.0;
			for( uint m=0; m<v2v[n].size(); m++ ) {
				FLOAT64 w = ringArea[v2v[n][m]];
				wsum += w;
				new_curvature += w*old_curvature[v2v[n][m]];
			}
			FLOAT64 w = ringArea[n];
			new_curvature += w*old_curvature[n];
			wsum += w;
			if( wsum ) {
				curvature[n] = new_curvature / wsum;
			}
		}
	}
#endif
}