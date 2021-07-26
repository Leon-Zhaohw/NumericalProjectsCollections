/*
 *	meshsurf3.cpp
 *
 *	Created by Ryoichi Ando on 2012/08/04
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "meshsurf3.h"
#include "util3.h"
#include "opengl.h"
#include "matutil.h"
#include "fastmarch3.h"
#include <vector>

meshsurf3::meshsurf3() {
	g = NULL;
	extrapolate_dist = 1.0;
}

meshsurf3::meshsurf3( const meshsurf3 &meshsurf ) {
	tick(); dump("Copying mesh surf...");
	*this = meshsurf;
	dump("Done. Took %s.\n",stock());
}

void meshsurf3::setReference(const mesher3 *g) {
	this->g = g;
}

void meshsurf3::setExtrapolationDist( FLOAT64 dist ) {
	extrapolate_dist = dist;
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

bool meshsurf3::isOpEdge( uint idx0, uint idx1 ) {
	if( idx0 == idx1 ) return false;
	const mesher3 &g = *this->g;
	for( uint i=0; i<2; i++ ) for( uint j=0; j<2; j++ ) if( g.edges[idx0][i] == g.edges[idx1][j] ) return false;
	return true;
}

void meshsurf3::fillHoles( std::vector<FLOAT64> &values, const levelset3 *solid, const std::vector<vec3d> &nodes, const std::vector<std::vector<uint> > &elements,
						   const std::vector<std::vector<uint> > &node2node, FLOAT64 dx ) {
	// Simple hack to fix check-mark artifacts
	tick(); dump("Filling holes on solid wall");
	uint count_filled = 0;
	for( uint k=0; k<3; k++ ) {
		std::vector<bool> adjacentFlag(nodes.size());
		std::vector<FLOAT64> old_levelset = values;
		PARALLEL_FOR for( uint n=0; n<values.size(); n++ ) {
			bool connected_wall = false;
			if( solid->evalLevelset(nodes[n]) > -dx ) {
				for( uint m=0; m<node2node[n].size(); m++ ) {
					if( solid->evalLevelset(nodes[node2node[n][m]]) < dx ) {
						connected_wall = true;
						break;
					}
				}
			}
			adjacentFlag[n] = connected_wall;
		}
		PARALLEL_FOR for( uint n=0; n<nodes.size(); n++ ) {
			if( adjacentFlag[n] && old_levelset[n] > 0.0 ) {
				uint fluid_count = 0;
				uint non_fluid_count = 0;
				FLOAT64 avg = 0.0;
				FLOAT64 sum = 0.0;
				for( uint m=0; m<node2node[n].size(); m++ ) {
					uint idx = node2node[n][m];
					if( adjacentFlag[idx] ) {
						if( old_levelset[idx] < 0.0 ) {
							fluid_count ++;
							avg += old_levelset[idx];
							sum ++;
						} else {
							non_fluid_count ++;
						}
					}
				}
				if( sum && fluid_count > 2*non_fluid_count ) {
					values[n] = avg / sum;
					count_filled ++;
				}
			}
		}
		dump(".");
	}
	dump("Done. Took %s. Filled %d holes.\n",stock("meshsurf_fill_hole"), count_filled );
	writeNumber("meshsurf_fill_hole_number", count_filled);
}

void meshsurf3::buildSurface( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, const levelset3 *fluid, const levelset3 *solid, bool enclose, FLOAT64 dpx, uint iteration, bool doFit ) {
	// Save a reference to the mesher
	const mesher3 &g = *this->g;
	
	tick(); dump("Computing levelset for %d nodes...", g.nodes.size());
	values.clear();
	values.resize(g.nodes.size());
	if( values.empty() ) {
		email::print("Mesher Uninitialized !\n");
		email::send();
		exit(0);
	}
	// Compute levelset for each node
	PARALLEL_FOR for( uint n=0; n<g.nodes.size(); n++ ) {
		FLOAT64 phi = fluid->evalLevelset(g.nodes[n]);
		if( enclose ) phi = fmax(phi,-dpx-solid->evalLevelset(g.nodes[n]));
		values[n] = phi;
	}
	dump("Done. Took %s.\n",stock("meshsurf_levelset_sample"));
	
	// Fill holes
	fillHoles(values,solid,g.nodes,g.elements,g.node2node,fluid->dx);
	
	// Calculate intersections on edges
	tick(); dump("Computing cut edges for each tet...");
	std::vector<vec3d> cutpoints(g.edges.size());
	std::vector<int> cutIndices(g.edges.size());
	uint index = 0;
	for( int n=0; n<g.edges.size(); n++ ) {
		FLOAT64 levelsets[2];
		for( uint m=0; m<2; m++ ) levelsets[m] = values[g.edges[n][m]];
		if( copysign(1.0,levelsets[0]) * copysign(1.0,levelsets[1]) <= 0 ) {
			FLOAT64 det = levelsets[0]-levelsets[1];
			if( fabs(det) < 1e-16 ) det = copysign(1e-16,det);
			FLOAT64 a = fmin(0.99,fmax(0.01,levelsets[0]/det));
			vec3d p = (1.0-a)*g.nodes[g.edges[n][0]]+a*g.nodes[g.edges[n][1]];
			cutpoints[n] = p;
			cutIndices[n] = index++;
		} else {
			cutIndices[n] = -1;
		}
	}
	
	// Copy them to packed array
	vertices.clear();
	vertices.resize(index);
	index = 0;
	for( uint n=0; n<cutIndices.size(); n++ ) {
		if( cutIndices[n] >= 0 ) vertices[index++] = cutpoints[n];
	}
	dump("Done. Took %s.\n",stock());
	
	// Stitch faces (marching tet)
	tick(); dump("Stitching edges...");
	std::vector<std::vector<uint> > node_faces;
	node_faces.resize(g.nodes.size());
	faces.clear();
	for( uint n=0; n<g.elements.size(); n++ ) {
		std::vector<uint> nv;
		for( uint m=0; m<g.element_edges[n].size(); m++ ) {
			uint idx = g.element_edges[n][m];
			if( cutIndices[idx] >= 0 ) {
				nv.push_back(cutIndices[idx]);
			}
		}
		if( nv.size() == 4 ) {
			// Triangulate the veritces
			int idx0 = -1;
			int idx1 = -1;
			for( uint m=0; m<g.element_edges[n].size(); m++ ) {
				if( cutIndices[g.element_edges[n][m]] >= 0 ) {
					idx0 = g.element_edges[n][m];
					break;
				}
			}
			// If found one intersection edge
			if( idx0 >= 0 ) {
				for( uint m=0; m<g.element_edges[n].size(); m++ ) {
					uint opID = g.element_edges[n][m];
					if( cutIndices[opID] >= 0 && isOpEdge(idx0,opID)) {
						idx1 = opID;
						break;
					}
				}
				// If found opposite intersection edge
				if( idx1 >= 0 ) {
					vec3i v;
					v[0] = cutIndices[idx0];
					v[1] = cutIndices[idx1];
					for( uint m=0; m<g.element_edges[n].size(); m++ ) {
						uint idx2 = g.element_edges[n][m];
						if( cutIndices[idx2] >= 0 && idx2 != idx0 && idx2 != idx1 ) {
							v[2] = cutIndices[idx2];
							faces.push_back(v);
							for( uint m=0; m<g.elements[n].size(); m++ ) {
								node_faces[g.elements[n][m]].push_back(faces.size()-1);
							}
						}
					}
				} else {
					dump( "Opposite edge was not found !\n");
					exit(0);
				}
			} else {
				dump( "Ground edge was not found !\n");
				exit(0);
			}
		} else if( nv.size() == 3 ) {
			faces.push_back(vec3i(nv[0],nv[1],nv[2]));
			for( uint m=0; m<g.elements[n].size(); m++ ) {
				node_faces[g.elements[n][m]].push_back(faces.size()-1);
			}
		} else if( nv.size() == 0 ) {
			// Empty surface. Do nothing...
		} else {
			dump( "\nBad stitch encounrtered: Element - edge count = %d.\n", g.element_edges[n].size() );
			dump( "Unknown set found N(%e,%e,%e,%e) = %d.\n",
				   values[g.elements[n][0]], values[g.elements[n][1]], values[g.elements[n][2]], values[g.elements[n][3]], nv.size());
			for( uint k=0; k<g.element_edges[n].size(); k++ ) {
				uint vidx[2] = { g.edges[g.element_edges[n][k]][0], g.edges[g.element_edges[n][k]][1] };
				dump( "Edge info[%d] = edge(%d,%d) = (%e,%e) = (%f,%f)\n", k, vidx[0], vidx[1], values[vidx[0]], values[vidx[1]], copysign(1.0,values[vidx[0]]), copysign(1.0,values[vidx[1]]) );
			}
		}
	}
	dump("Done. Took %s.\n",stock("meshsurf_stitch_mesh"));
	
	// Find closest position
	tick(); dump("Finding closest surface position at surface tets...");
	surfacePos.clear();
	surfacePos.resize(g.nodes.size());
	for( uint n=0; n<g.nodes.size(); n++ ) {
		FLOAT64 dist = 1e8;
		vec3d p = g.nodes[n];
		for( uint m=0; m<node_faces[n].size(); m++ ) {
			vec3d p0 = vertices[faces[node_faces[n][m]][0]];
			vec3d p1 = vertices[faces[node_faces[n][m]][1]];
			vec3d p2 = vertices[faces[node_faces[n][m]][2]];
			vec3d out = g.nodes[n];
			vec3d src = out;
			FLOAT64 d = fmax(1e-8,point_triangle_distance(src,p0,p1,p2,out));
			if( d < dist ) {
				p = out;
				dist = d;
			}
		}
		if( dist < 1.0 ) {
			values[n] = (values[n]>=0 ? 1.0 : -1.0) * dist;
			surfacePos[n].push_back(g.nodes[n]);
			surfacePos[n].push_back(p);
		} else {
			values[n] = (values[n]>=0 ? 1.0 : -1.0) * 1e8;
		}
	}
	dump("Done. Took %s.\n",stock("meshsurf_find_close_pos1"));
	
	if( extrapolate_dist ) {
		// Fast march
		tick();
		std::vector<fastmarch3<FLOAT64>::node3 *> fnodes(g.nodes.size());
		for( uint n=0; n<g.nodes.size(); n++ ) fnodes[n] = new fastmarch3<FLOAT64>::node3;
		for( uint n=0; n<g.nodes.size(); n++ ) {
			vec3d p = g.nodes[n];
			bool fixed = fabs(values[n]) < 1.0;
			fnodes[n]->p = p;
			fnodes[n]->fixed = fixed;
			fnodes[n]->levelset = values[n];
			fnodes[n]->value = 0.0;
			fnodes[n]->p2p.resize(g.node2node[n].size());
			for( uint m=0; m<g.node2node[n].size(); m++ ) {
				fnodes[n]->p2p[m] = fnodes[g.node2node[n][m]];
			}
		}
		fastmarch3<FLOAT64>::fastMarch(fnodes,extrapolate_dist,-extrapolate_dist,1);
		
		// Pick up values
		for( uint n=0; n<g.nodes.size(); n++ ) {
			values[n] = fnodes[n]->levelset;
			delete fnodes[n];
		}
		stock("meshsurf_fastmarch");
	}
	
	// Compute normal
	tick(); dump("Computing normals...");
	computeNormals(vertices,normals,faces,cutIndices);
	dump("Done. Took %s.\n",stock("meshsurf_normal"));

	// Flip facet rotation if necessary
	flipFacet(vertices,normals,faces,0.1*fluid->dx);
	
#if 1
	// Fit surface
	if( doFit ) {
		tick(); dump("Fitting %d surface vertices...", vertices.size() );
		for( uint k=0; k<1; k++ ) {
			smoothMesh(vertices,normals,faces,solid,dpx,1);
			PARALLEL_FOR for( uint n=0; n<vertices.size(); n++ ) {
				vec3d out = vertices[n];
				if( solid->evalLevelset(out) > dpx && fluid->getClosestSurfacePos(out) ) {
					vertices[n] = out;
				}				
			}
		}
		dump("Done. Took %s.\n",stock("meshsurf_fit"));
	}
#endif
	
	// Smooth mesh
#if 1
	tick(); dump("Smoothing faces...");
	smoothMesh(vertices,normals,faces,solid,dpx,iteration);
	dump("Done. Took %s.\n",stock("meshsurf_smooth"));
#endif
	
	// Find closest position again
	tick(); dump("Finding closest surface position at surface tets again...");
	surfacePos.clear();
	surfacePos.resize(g.nodes.size());
	for( uint n=0; n<g.nodes.size(); n++ ) {
		FLOAT64 dist = 1e8;
		vec3d p = g.nodes[n];
		for( uint m=0; m<node_faces[n].size(); m++ ) {
			vec3d p0 = vertices[faces[node_faces[n][m]][0]];
			vec3d p1 = vertices[faces[node_faces[n][m]][1]];
			vec3d p2 = vertices[faces[node_faces[n][m]][2]];
			vec3d out = g.nodes[n];
			vec3d src = out;
			FLOAT64 d = fmax(1e-8,point_triangle_distance(src,p0,p1,p2,out));
			if( d < dist ) {
				p = out;
				dist = d;
			}
		}
		if( dist < 1.0 ) {
			values[n] = (values[n]>=0 ? 1.0 : -1.0) * dist;
		}
	}
	dump("Done. Took %s.\n",stock("meshsurf_find_close_pos2"));
	
	// Flip again
	flipFacet(vertices,normals,faces,0.1*fluid->dx);
}

static bool is_nan( vec3d p ) {
	return is_nan(p[0]) || is_nan(p[1]) || is_nan(p[2]);
}

void meshsurf3::flipFacet( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, FLOAT64 min_dx ) {
	// Flip facet rotation if necessary
	for( uint n=0; n<faces.size(); n++) {
		uint vidx0 = faces[n][0];
		uint vidx1 = faces[n][1];
		uint vidx2 = faces[n][2];
		vec3d advnormal;
		FLOAT64 count = 0.0;
		if( ! is_nan(normals[vidx0]) ) { advnormal += normals[vidx0]; count++; }
		if( ! is_nan(normals[vidx1]) ) { advnormal += normals[vidx1]; count++; }
		if( ! is_nan(normals[vidx2]) ) { advnormal += normals[vidx2]; count++; }
		advnormal = advnormal / count;
		vec3d fnormal = (vertices[vidx1]-vertices[vidx0]) ^ (vertices[vidx2]-vertices[vidx0]);
		if( advnormal * fnormal > 0.0 ) {
			faces[n][0] = vidx2;
			faces[n][2] = vidx0;
		}
		if( ! fnormal.len2()) {
			vertices[vidx0] += 1e-8 * vec3d(nrand(),nrand(),nrand());
		}
	}
}

void meshsurf3::computeNormals( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, const std::vector<int> &cutIndices ) {
	// Save a reference to the mesher
	const mesher3 &g = *this->g;
	// Compute normals
	normals.resize(vertices.size());
	for( uint n=0; n<g.edges.size(); n++ ) {
		int vidx = cutIndices[n];
		if( vidx >= 0 ) {
			uint idx[2] = { g.edges[n][0], g.edges[n][1] };
			FLOAT64 dist[2] = { (vertices[vidx]-g.nodes[idx[0]]).len(), (vertices[vidx]-g.nodes[idx[1]]).len() };
			FLOAT64 dw[2] = { dist[1]/fmax(1e-8,dist[0]+dist[1]), dist[0]/fmax(1e-8,dist[0]+dist[1]) };
			vec3d normal;
			for( uint i=0; i<2; i++ ) {
				for( uint m=0; m<g.node_elements[idx[i]].size(); m++) {
					uint elm = g.node_elements[idx[i]][m];
					FLOAT64 w = dw[i]*g.volumes[elm];
					normal += -w*getElementGradient(elm);
				}
			}
			normals[vidx] = normal.normal();
		}
	}
}

void meshsurf3::smoothMesh( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, const levelset3 *solid, FLOAT64 dpx, uint iterations ) {
	std::vector<vec3d> old_vertices = vertices;
	std::vector<vec3d> old_normals = normals;
	std::vector<std::vector<uint> > connections;

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
			
	// Smooth out
	for( uint k=0; k<iterations; k++ ) {
		for( uint n=0; n<vertices.size(); n++ ) {
			if( solid->evalLevelset(vertices[n]) < dpx ) continue;
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
		
		for( uint n=0; n<normals.size(); n++ ) {
			if( solid->evalLevelset(vertices[n]) < dpx ) continue;
			FLOAT64 wsum = 0.0;
			vec3d normal;
			// Add self
			FLOAT64 w = is_nan(areas[n]) ? 0.0 : areas[n]+1e-8;
			if( ! is_nan(old_normals[n]) && w ) {
				normal += w*old_normals[n];
				wsum += w;
			}
			// Add neighbors
			for( uint m=0; m<connections[n].size(); m++ ) {
				uint idx = connections[n][m];
				if( ! is_nan(old_normals[idx]) ) {
					FLOAT64 w = is_nan(areas[idx]) ? 0.0 : areas[idx]+1e-8;
					if( w ) {
						normal += w*old_normals[idx];
						wsum += w;
					}
				}
			}
			if( wsum ) normals[n] = normal.normal();
		}
		old_normals = normals;
	}
	removeNan(vertices,normals,faces,connections,areas);
}

void meshsurf3::removeNan( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, const std::vector<std::vector<uint> > &connections, const std::vector<FLOAT64> &areas ) {
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

FLOAT64 meshsurf3::evalLevelset(vec3d p) const {
	const mesher3 &g = *this->g;	
	int n = g.hitElements(p);
	if( n >= 0 ) {
		FLOAT64 t[NUM_VERT];
		FLOAT64 b[] = { p[0], p[1], p[2], 1.0 };
		for( uint i=0; i<NUM_VERT; i++ ) {
			t[i] = 0.0;
			for( uint k=0; k<NUM_VERT; k++ ) t[i] += g.matrix[n].m[i][k]*b[k];
		}
		FLOAT64 res = 0.0;
		for( uint i=0; i<NUM_VERT; i++ ) {
			if( g.elements[n][i] >= g.nodes.size() ) {
				printf( "Mesh seems to be broken !!\n" );
				exit(0);
			}
			res += t[i]*values[g.elements[n][i]];
		}
		return res;
	}
	return 1.0;
}

vec3d meshsurf3::getElementGradient(uint elm) const {
	const mesher3 &g = *this->g;
	vec3d res;
	for( uint dim=0; dim<DIM; dim++ ) {
		res[dim] = 0.0;
		for( uint k=0; k<NUM_VERT; k++ ) res[dim] += g.matrix[elm].m[k][dim]*values[g.elements[elm][k]];
	}
	return res;
}

vec3d meshsurf3::evalGradient(vec3d p) const {
	const mesher3 &g = *this->g;
	int n = g.hitElements(p);
	if( n >= 0 ) {
		return getElementGradient(n);
	}
	return vec3d();
}

void meshsurf3::drawGL( const levelset3 *solid, uint kind ) const {
}
