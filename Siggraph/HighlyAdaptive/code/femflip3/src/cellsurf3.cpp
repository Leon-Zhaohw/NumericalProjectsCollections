/*
 *	cellsurf3.cpp
 *	
 *	Created by Ryoichi Ando on 1/9/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "util3.h"
#include "cellsurf3.h"
#include "array3.h"
#include "CIsoSurface.h"
#include "kernel.h"
#include "fastmarch3.h"
#include "opengl.h"
#include <vector>
using namespace std;

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

static void putHash( vec3d left_bottom, vec3d right_top, array3<std::vector<uint> > &hash, uint id ) {
	FLOAT64 e = 1.0e-12;
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
	
#if 0
	if( max_i-min_i != 1 ) {
		printf( "cellsurf3.cpp: x %d %d\n", min_i, max_i );
		exit(0);
	}
	if( max_j-min_j != 1 ) {
		printf( "cellsurf3.cpp: y %d %d\n", min_j, max_j );
		exit(0);
	}
	if( max_k-min_k != 1 ) {
		printf( "cellsurf3.cpp: z %d %d\n", min_k, max_k );
		exit(0);
	}
#endif
	
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

cellsurf3::cellsurf3() {
	extrapolate_dist = 1.0;
}

cellsurf3::cellsurf3( const cellsurf3& cellsurf ) {
	tick(); dump("Copying cell surf...");
	*this = cellsurf;
	dump("Done. Took %s.\n",stock());
}

void cellsurf3::setResolution( uint gn ) {
	this->gn = gn;
}

void cellsurf3::setExtrapolationDist( FLOAT64 dist ) {
	extrapolate_dist = dist;
}

FLOAT64 cellsurf3::evalLevelset(vec3d p) const {
	return util3::interp<FLOAT64>(gn*p,nodalLevelset);
}

vec3d cellsurf3::evalGradient(vec3d p) const {
	uint w, h, d;
	w = nodalLevelset.size().w;
	h = nodalLevelset.size().h;
	d = nodalLevelset.size().d;
	const array3<FLOAT64> &q = nodalLevelset;
	FLOAT64 x = fmax(0.0,fmin(w,gn*p[0]));
	FLOAT64 y = fmax(0.0,fmin(h,gn*p[1]));
	FLOAT64 z = fmax(0.0,fmin(d,gn*p[2]));
	int i = imin(x,w-2);
	int j = imin(y,h-2);
	int k = imin(z,d-2);
	FLOAT64 gx = (k+1-z)*((q[i+1][j][k]-q[i][j][k])*(j+1-y)+(q[i+1][j+1][k]-q[i][j+1][k])*(y-j))+
				 (z-k)*((q[i+1][j][k+1]-q[i][j][k+1])*(j+1-y)+(q[i+1][j+1][k+1]-q[i][j+1][k+1])*(y-j));
	FLOAT64 gy = (k+1-z)*((q[i][j+1][k]-q[i][j][k])*(i+1-x)+(q[i+1][j+1][k]-q[i+1][j][k])*(x-i))+
				 (z-k)*((q[i][j+1][k+1]-q[i][j][k+1])*(i+1-x)+(q[i+1][j+1][k+1]-q[i+1][j][k+1])*(x-i));
	FLOAT64 gz = (i+1-x)*(j+1-y)*(q[i][j][k+1]-q[i][j][k])+(x-i)*(j+1-y)*(q[i+1][j][k+1]-q[i+1][j][k])+
				 (i+1-x)*(y-j)*(q[i][j+1][k+1]-q[i][j+1][k])+(x-i)*(y-j)*(q[i+1][j+1][k+1]-q[i+1][j+1][k]);
	return vec3d(gx,gy,gz).normal();
}

FLOAT64 cellsurf3::findClosestSurfPosition( std::vector<vec3d> &vertices, std::vector<vec3i> &faces, int nodal_i, int nodal_j, int nodal_k, vec3d &p ) {
	FLOAT64 dist = 1e8;
	vec3d dest = p;
	FOR_EACH(2,2,2) {
		int ni = nodal_i-1+i;
		int nj = nodal_j-1+j;
		int nk = nodal_k-1+k;
		if( ni<0 || nj<0 || nk<0 || ni>=gn || nj>=gn || nk>=gn ) continue;
		int size = hash[ni][nj][nk].size();
		for( int n=0; n<size; n++ ) {
			uint id = hash[ni][nj][nk][n];
			const vec3d &p0 = vertices[faces[id][0]];
			const vec3d &p1 = vertices[faces[id][1]];
			const vec3d &p2 = vertices[faces[id][2]];
			vec3d out;
			FLOAT64 d = point_triangle_distance(p,p0,p1,p2,out);
			if( d < dist ) {
				dist = d;
				dest = out;
			}
		}
	} END_FOR
	p = dest;
	return dist;
}

void cellsurf3::buildSurface( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, const levelset3 *fluid, const levelset3 *solid, bool enclose ) {
	tick(); dump(">>> Building surface mesh...\n");
	dx = 1.0/gn;

	// Compute nodal levelset
	tick();	dump("Evaluating cubic levelset %dx%dx%d...", gn,gn,gn );
	
	nodalLevelset.resize(gn+1,gn+1,gn+1);
	nodalLevelset.clear();
	FOR_EACH(gn+1,gn+1,gn+1) {
		vec3d p(dx*i,dx*j,dx*k);
		FLOAT64 phi = fluid->evalLevelset(p);
		if( enclose ) {
			dump( "Enclose mode is still not implemeted.\n");
			exit(0);
		}
		nodalLevelset[i][j][k] = phi;
	} END_FOR
	dump("Done. Took %s.\n",stock("cellsurf_levelset"));
	tick(); dump("Marching cube mesher %dx%dx%d...", gn,gn,gn );
	
	// Density Voxels
    vector<Voxel<FLOAT64> > Voxels;
    static array3<FLOAT64> phi = array3<FLOAT64>(gn+1,gn+1,gn+1);
	
	// Create Density Field
    PARALLEL_FOR FOR_EACH(gn+1,gn+1,gn+1) {
        FLOAT64 value = nodalLevelset[i][j][k];
        phi[i][j][k] = value;	
    } END_FOR
	
	FOR_EACH(gn,gn,gn) {
        Voxel<FLOAT64> aVoxel;
        aVoxel.x[0] = i;
        aVoxel.x[1] = j;
        aVoxel.x[2] = k;
        bool FoundNegative = false;
        bool FoundPositive = false;
        for( int ix=0; ix<2; ix++ ) {
            for( int iy=0; iy<2; iy++ ) {
                for( int iz=0; iz<2; iz++ ) {
                    aVoxel.v[ix][iy][iz] = phi[i+ix][j+iy][k+iz];
                    if( aVoxel.v[ix][iy][iz] < 0.0 ) {
                        FoundNegative = true;
                    } else {
                        FoundPositive = true;
                    }
                }
            }
        } 
        if( FoundNegative && FoundPositive ) Voxels.push_back(aVoxel);
    } END_FOR
	
	// Extract mesh
	vertices.clear();
	normals.clear();
	faces.clear();
	CIsoSurface<FLOAT64> surface;
	surface.GenerateSurface(Voxels,0,gn,gn,gn,dx,dx,dx);
	
	for( int i=0; i<surface.m_nVertices; i++ ) {
		vertices.push_back(vec3d(surface.m_ppt3dVertices[i][0],
								 surface.m_ppt3dVertices[i][1],
								 surface.m_ppt3dVertices[i][2]));
	}
	
	for( int i=0; i<surface.m_nNormals; i++ ) {
		normals.push_back(vec3d(surface.m_pvec3dNormals[i][0],
								surface.m_pvec3dNormals[i][1],
								surface.m_pvec3dNormals[i][2]));
	}

	for (uint i = 0; i < surface.m_nTriangles; i++) {
		uint id0, id1, id2;
		id0 = surface.m_piTriangleIndices[i*3];
		id1 = surface.m_piTriangleIndices[i*3+1];
		id2 = surface.m_piTriangleIndices[i*3+2];
        faces.push_back(vec3i(id0,id1,id2));
	}
	
	smoothMesh(vertices,normals,faces,solid,1);
	
	dump("Done. Took %s.\n",stock());
	buildLevelsetInfo(vertices,normals,faces);
	dump("<<< Building surface mesh done. Took %s.\n",stock("cellsurf_surfacemesh"));
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

void cellsurf3::buildLevelsetInfo( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces ) {
	tick();	dump("Building hash table...");
	
	// Build spatial hash
	hash.resize(gn,gn,gn);
	hash.clear();
	for( uint n=0; n<faces.size(); n++ ) {
		vec3d left_bottom;
		vec3d right_top;
		vector<vec3d> points(3);
		points[0] = vertices[faces[n][0]];
		points[1] = vertices[faces[n][1]];
		points[2] = vertices[faces[n][2]];
		getBoundingBox(points,left_bottom,right_top);
		putHash(left_bottom,right_top,hash,n);
	}
	
	vhash.resize(gn,gn,gn);
	vhash.clear();
	for( uint n=0; n<vertices.size(); n++ ) {
		int i = imax(0,imin(gn-1,vertices[n][0]*gn));
		int j = imax(0,imin(gn-1,vertices[n][1]*gn));
		int k = imax(0,imin(gn-1,vertices[n][2]*gn));
		vhash[i][j][k].push_back(n);
	}
	
	dump("Done. Took %s.\n",stock("cellsurf_hash"));
	tick(); dump("Computing SLC levelset at narrow bands...");
	
	// Allocate nodal grids
	referencePos.resize(gn+1,gn+1,gn+1);
	referencePos.clear();
	
	// Gather narrow bands
	vector<vec3i> narrowBands;
	FOR_EACH(gn+1,gn+1,gn+1) {
		bool found = false;
		for( int ki=0; ki<2; ki++ ) for( int kj=0; kj<2; kj++ ) for( int kk=0; kk<2; kk++ ) {
			if( found ) continue;
			int ni = i-1+ki;
			int nj = j-1+kj;
			int nk = k-1+kk;
			if( ni<0 || nj<0 || nk<0 || ni>=gn || nj>=gn || nk>=gn ) continue;
			if( hash[ni][nj][nk].size() ) {
				narrowBands.push_back(vec3i(i,j,k));
				found = true;
			}
		}
		vec3d p(dx*i,dx*j,dx*k);
		nodalLevelset[i][j][k] = 1e8*(util3::interp<FLOAT64>(gn*p,nodalLevelset)>0.0?1.0:-1.0);
	} END_FOR
	
	// Fix the narrow bands
	referencePos.resize(gn+1,gn+1,gn+1);
	referencePos.clear();
	PARALLEL_FOR for( uint n=0; n<narrowBands.size(); n++ ) {
		uint i = narrowBands[n][0];
		uint j = narrowBands[n][1];
		uint k = narrowBands[n][2];
		vec3d p(dx*i,dx*j,dx*k);
		FLOAT64 dist = findClosestSurfPosition(vertices,faces,i,j,k,p);
		FLOAT64 sgn = nodalLevelset[i][j][k] >= 0.0 ? 1.0 : -1.0;
		nodalLevelset[i][j][k] = dist*sgn;
		referencePos[i][j][k] = p;
	}
	
	dump("Done. Took %s.\n",stock("cellsurf_SLC"));
	
	if( extrapolate_dist ) {
		// Fast march
		tick(); dump("Fast march...");
		
		// Setup a network
		FLOAT64 dx = 1.0/gn;
		array3<fastmarch3<FLOAT64>::node3> nodeArray;
		nodeArray.resize(gn+1,gn+1,gn+1);
		PARALLEL_FOR FOR_EACH(gn+1,gn+1,gn+1) {
			vec3d p = vec3d(i*dx,j*dx,k*dx);
			nodeArray[i][j][k].p = p;
			nodeArray[i][j][k].fixed = fabs(nodalLevelset[i][j][k])<1.0;
			nodeArray[i][j][k].levelset = nodalLevelset[i][j][k];
			nodeArray[i][j][k].value = 0.0;
			int q[][DIM] = {{i-1,j,k},{i+1,j,k},{i,j-1,k},{i,j+1,k},{i,j,k-1},{i,j,k+1}};
			for( uint n=0; n<6; n++ ) {
				int ni = q[n][0];
				int nj = q[n][1];
				int nk = q[n][2];
				if( ni>=0 && ni<gn+1 && nj>=0 && nj<gn+1 && nk>=0 && nk<gn+1 ) {
					nodeArray[i][j][k].p2p.push_back(&nodeArray[ni][nj][nk]);
				}
			}
		} END_FOR
		std::vector<fastmarch3<FLOAT64>::node3 *> nodes((gn+1)*(gn+1)*(gn+1));
		uint index = 0;
		FOR_EACH(gn+1,gn+1,gn+1) {
			nodes[index++] = &nodeArray[i][j][k];
		} END_FOR
		
		// Perform fast march
		fastmarch3<FLOAT64>::fastMarch(nodes,extrapolate_dist,-extrapolate_dist,1);
		
		// Pickup extrapolated levelsets
		PARALLEL_FOR FOR_EACH(gn+1,gn+1,gn+1) {
			nodalLevelset[i][j][k] = nodeArray[i][j][k].levelset;
		} END_FOR
		
		dump("Done. Took %s.\n",stock("cellsurf_fastmarch"));
	}
}

static bool is_nan( vec3d p ) {
	return is_nan(p[0]) || is_nan(p[1]) || is_nan(p[2]);
}

void cellsurf3::removeNan( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, const std::vector<std::vector<uint> > &connections, const std::vector<FLOAT64> &areas ) {
	uint tri_cnt = 0;
	while( true ) {
		bool foundNan = false;
		for( uint n=0; n<vertices.size(); n++ ) {
			if( is_nan(vertices[n]) ) {
				foundNan = true;
				vec3d avg;
				FLOAT64 sum=0.0;
				for( uint m=0; m<connections[n].size(); m++ ) {
					uint idx = connections[n][m];
					if( ! is_nan(vertices[idx]) ) {
						FLOAT64 w = is_nan(areas[idx]) ? 0.0 : areas[idx]+1e-8;
						if( w ) {
							avg += w*vertices[idx];
							sum += w;
						}
					}
				}
				if( sum ) vertices[n] = avg / sum;
			}
		}
		for( uint n=0; n<vertices.size(); n++ ) {
			if( is_nan(normals[n]) ) {
				foundNan = true;
				vec3d avg;
				FLOAT64 sum=0.0;
				for( uint m=0; m<connections[n].size(); m++ ) {
					uint idx = connections[n][m];
					if( ! is_nan(normals[idx]) ) {
						FLOAT64 w = is_nan(areas[idx]) ? 0.0 : areas[idx]+1e-8;
						if( w ) {
							avg += w*normals[idx];
							sum += w;
						}
					}
				}
				if( sum ) normals[n] = avg.normal();
			}
		}
		if( tri_cnt++ > 20 ) {
			for( uint n=0; n<faces.size(); n++) {
				uint vidx0 = faces[n][0];
				uint vidx1 = faces[n][1];
				uint vidx2 = faces[n][2];
				vec3d fnormal = (vertices[vidx1]-vertices[vidx0]) ^ (vertices[vidx2]-vertices[vidx0]);
				if( is_nan(normals[vidx0]) ) normals[vidx0] = fnormal;
				if( is_nan(normals[vidx1]) ) normals[vidx1] = fnormal;
				if( is_nan(normals[vidx2]) ) normals[vidx2] = fnormal;
			}
		}
		if( ! foundNan ) break;
	}
}

void cellsurf3::smoothMesh( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, const levelset3 *solid, uint iterations ) {
	std::vector<vec3d> old_vertices = vertices;
	std::vector<vec3d> old_normals = normals;
	std::vector<std::vector<uint> > connections;
	std::vector<bool> at_inside;	
	at_inside.resize(vertices.size());
	
	std::vector<FLOAT64> areas;
	connections.resize(vertices.size());
	areas.resize(vertices.size());
	for( uint n=0; n<vertices.size(); n++ ) {
		areas[n] = 0.0;
	}
	for( uint n=0; n<faces.size(); n++ ) {
		for( uint i=0; i<3; i++ ) {
			connections[faces[n][i]].push_back(faces[n][(i+1)%3]);
		}
		FLOAT64 area = ((vertices[faces[n][1]]-vertices[faces[n][0]])^(vertices[faces[n][2]]-vertices[faces[n][0]])).len();
		for( uint j=0; j<3; j++ ) {
			areas[faces[n][j]] += area / 3.0;
		}
	}
	
	// Mark inside
	for( uint n=0; n<vertices.size(); n++ ) {
		bool inside = false;
		// Compute average connection length
		for( uint dim=0; dim<DIM; dim++) {
			if( vertices[n][dim] <= 1e-4 || vertices[n][dim] >= 1.0-1e-4 ) inside = true;
		}
		
		if( solid ) {
			FLOAT64 max_len = 0.0;
			for( uint m=0; m<connections[n].size(); m++ ) {
				uint idx = connections[n][m];
				max_len = fmax((vertices[idx]-vertices[n]).len(),max_len);
			}
			if( solid && solid->evalLevelset(vertices[n]) <= 0.1*max_len ) {
				inside = true;
			}
		}
		at_inside[n] = inside;
	}
	
	// Smooth out
	for( uint k=0; k<iterations; k++ ) {
		for( uint n=0; n<vertices.size(); n++ ) {
			if( at_inside[n] ) continue;
			
			FLOAT64 wsum = 0.0;
			vec3d pos;
			// Add self
			FLOAT64 w = is_nan(areas[n]) ? 0.0 : areas[n]+1e-8;
			if( ! is_nan(old_vertices[n]) && w ) {
				pos += w*old_vertices[n];
				wsum += w;
			}
			// Add neighbors
			for( uint m=0; m<connections[n].size(); m++ ) {
				uint idx = connections[n][m];
				if( ! is_nan(old_vertices[idx]) ) {
					FLOAT64 w = is_nan(areas[idx]) ? 0.0 : areas[idx]+1e-8;
					if( w ) {
						pos += w*old_vertices[idx];
						wsum += w;
					}
				}
			}
			if( wsum ) vertices[n] = pos / wsum;
		}
		old_vertices = vertices;
	}
	removeNan(vertices,normals,faces,connections,areas);
}

void cellsurf3::drawGL( const levelset3 *solid, uint kind ) const {
}
