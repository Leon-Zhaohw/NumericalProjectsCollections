/*
 *	meshsurf2.cpp
 *
 *	Created by Ryoichi Ando on 2012/08/04
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "meshsurf2.h"
#include "util2.h"
#include "opengl.h"
#include "matutil.h"
#include "fastmarch2.h"
#include "kernel.h"
#include <vector>

meshsurf2::meshsurf2() {
	g = NULL;
	extrapolate_dist = 1.0;
}

meshsurf2::meshsurf2( const meshsurf2 &meshsurf ) {
	*this = meshsurf;
}

void meshsurf2::setReference(const mesher2 *g) {
	this->g = g;
}

void meshsurf2::setExtrapolationDist( FLOAT64 dist ) {
	extrapolate_dist = dist;
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

void meshsurf2::fillHoles( std::vector<FLOAT64> &values, const levelset2 *solid, const std::vector<vec2d> &nodes, const std::vector<std::vector<uint> > &elements,
						   const std::vector<std::vector<uint> > &node2node, FLOAT64 dx ) {
	// Simple hack to fix check-mark artifacts
	for( uint k=0; k<3; k++ ) {
		std::vector<bool> adjacentFlag(nodes.size());
		std::vector<FLOAT64> old_levelset = values;
		PARALLEL_FOR for( uint n=0; n<values.size(); n++ ) {
			bool connected_wall = false;
			if( solid->evalLevelset(nodes[n]) > dx ) {
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
				}
			}
		}
	}
}

void meshsurf2::buildSurface( std::vector<vec2d> &vertices, std::vector<vec2d> &normals, std::vector<vec2i> &faces, const levelset2 *fluid, const levelset2 *solid, bool enclose, FLOAT64 dpx, uint iteration, bool doFit) {
	// Save a reference to the mesher
	const mesher2 &g = *this->g;
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
	
	// Fill holes
	fillHoles(values,solid,g.nodes,g.elements,g.node2node,fluid->dx);
	
	// Calculate intersections on edges
	std::vector<vec2d> cutpoints(g.edges.size());
	std::vector<int> cutIndices(g.edges.size());
	uint index = 0;
	for( int n=0; n<g.edges.size(); n++ ) {
		FLOAT64 levelsets[2];
		for( uint m=0; m<2; m++ ) levelsets[m] = values[g.edges[n][m]];
		if( copysign(1.0,levelsets[0]) * copysign(1.0,levelsets[1]) <= 0 ) {
			FLOAT64 det = levelsets[0]-levelsets[1];
			if( fabs(det) < 1e-16 ) det = copysign(1e-16,det);
			FLOAT64 a = fmin(0.99,fmax(0.01,levelsets[0]/det));
			vec2d p = (1.0-a)*g.nodes[g.edges[n][0]]+a*g.nodes[g.edges[n][1]];
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
	
	// Stitch faces
	std::vector<std::vector<uint> > node_faces;
	node_faces.resize(g.nodes.size());
	faces.clear();
	for( uint n=0; n<g.elements.size(); n++ ) {
		std::vector<uint> pair;
		for( uint m=0; m<g.element_edges[n].size(); m++ ) {
			uint idx = g.element_edges[n][m];
			if( cutIndices[idx] >= 0 ) {
				pair.push_back(cutIndices[idx]);
			}
		}
		// Append a surface patch
		if( pair.size() == 2 ) {
			faces.push_back(vec2i(pair[0],pair[1]));
			for( uint m=0; m<g.elements[n].size(); m++ ) {
				node_faces[g.elements[n][m]].push_back(faces.size()-1);
			}
		}
	}
	
	// Find closest position
	surfacePos.clear();
	surfacePos.resize(g.nodes.size());
	for( uint n=0; n<g.nodes.size(); n++ ) {
		FLOAT64 dist = 1e8;
		vec2d p = g.nodes[n];
		for( uint m=0; m<node_faces[n].size(); m++ ) {
			vec2d p0 = vertices[faces[node_faces[n][m]][0]];
			vec2d p1 = vertices[faces[node_faces[n][m]][1]];
			vec2d out = g.nodes[n];
			FLOAT64 d = fmax(1e-8,distance(p0,p1,out));
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
	
	if( extrapolate_dist ) {
		// Fast march
		std::vector<fastmarch2<FLOAT64>::node2 *> fnodes(g.nodes.size());
		for( uint n=0; n<g.nodes.size(); n++ ) fnodes[n] = new fastmarch2<FLOAT64>::node2;
		for( uint n=0; n<g.nodes.size(); n++ ) {
			vec2d p = g.nodes[n];
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
		fastmarch2<FLOAT64>::fastMarch(fnodes,extrapolate_dist,-extrapolate_dist,1);
		// Pick up values
		for( uint n=0; n<g.nodes.size(); n++ ) {
			values[n] = fnodes[n]->levelset;
			delete fnodes[n];
		}
	}
	
	// Compute normals
	normals.resize(vertices.size());
	for( uint n=0; n<g.edges.size(); n++ ) {
		int vidx = cutIndices[n];
		if( vidx >= 0 ) {
			uint idx[2] = { g.edges[n][0], g.edges[n][1] };
			FLOAT64 dist[2] = { (vertices[vidx]-g.nodes[idx[0]]).len(), (vertices[vidx]-g.nodes[idx[1]]).len() };
			FLOAT64 dw[2] = { dist[1]/(dist[0]+dist[1]), dist[0]/(dist[0]+dist[1]) };
			vec2d normal;
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
	
	// Fix facet direction
	for( uint n=0; n<faces.size(); n++ ) {
		vec2d avg_normal = 0.5*(normals[faces[n][0]]+normals[faces[n][1]]);
		vec2d dir_vec = vertices[faces[n][1]]-vertices[faces[n][0]];
		if( avg_normal * dir_vec.rotate() < 0.0 ) {
			uint save = faces[n][0];
			faces[n][0] = faces[n][1];
			faces[n][1] = save;
		}
	}
	
#if 1
	// Fit surface
	if( doFit ) {
		for( uint k=0; k<1; k++ ) {
			smoothMesh(vertices,normals,faces,solid,dpx,1);
			PARALLEL_FOR for( uint n=0; n<vertices.size(); n++ ) {
				vec2d out = vertices[n];
				if( solid->evalLevelset(out) > dpx && fluid->getClosestSurfacePos(out) ) {
					vertices[n] = out;
				}
			}
		}
	}
#endif
	
	// Smooth surface
	smoothMesh(vertices,normals,faces,solid,dpx,iteration);
	
	// Find closest position again !
	surfacePos.clear();
	surfacePos.resize(g.nodes.size());
	for( uint n=0; n<g.nodes.size(); n++ ) {
		FLOAT64 dist = 1e8;
		vec2d p = g.nodes[n];
		for( uint m=0; m<node_faces[n].size(); m++ ) {
			vec2d p0 = vertices[faces[node_faces[n][m]][0]];
			vec2d p1 = vertices[faces[node_faces[n][m]][1]];
			vec2d out = g.nodes[n];
			FLOAT64 d = fmax(1e-8,distance(p0,p1,out));
			if( d < dist ) {
				p = out;
				dist = d;
			}
		}
		if( dist < 1.0 ) {
			values[n] = (values[n]>=0 ? 1.0 : -1.0) * dist;
		}
	}
	
	// Fix facet direction
	for( uint n=0; n<faces.size(); n++ ) {
		vec2d avg_normal = 0.5*(normals[faces[n][0]]+normals[faces[n][1]]);
		vec2d dir_vec = vertices[faces[n][1]]-vertices[faces[n][0]];
		if( avg_normal * dir_vec.rotate() < 0.0 ) {
			uint save = faces[n][0];
			faces[n][0] = faces[n][1];
			faces[n][1] = save;
		}
	}
	
	// Try plot into matlab
#if 0
	std::vector<FLOAT64> values_export = values;
	for( uint n=0; n<values.size(); n++ ) {
		values_export[n] = values[n];
	}
	g.write_matlab("levelset.m",values_export);
	exit(0);
#endif
}

static bool is_nan( vec2d p ) {
	return is_nan(p[0]) || is_nan(p[1]) || is_nan(p[2]);
}

void meshsurf2::smoothMesh( std::vector<vec2d> &vertices, std::vector<vec2d> &normals, std::vector<vec2i> &faces, const levelset2 *solid, FLOAT64 dpx, uint iterations ) {
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
			if( solid->evalLevelset(vertices[n]) < dpx ) continue;
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
		
		for( uint n=0; n<normals.size(); n++ ) {
			if( solid->evalLevelset(vertices[n]) < dpx ) continue;
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
	}
}

FLOAT64 meshsurf2::evalLevelset(vec2d p) const {
	const mesher2 &g = *this->g;
	int n = g.hitElements(p);
	if( n >= 0 ) {
		FLOAT64 t[NUM_VERT];
		FLOAT64 b[] = { p[0], p[1], 1.0 };
		for( uint i=0; i<NUM_VERT; i++ ) {
			t[i] = 0.0;
			for( uint k=0; k<NUM_VERT; k++ ) t[i] += g.matrix[n].m[i][k]*b[k];
		}
		FLOAT64 res = 0.0;
		for( uint i=0; i<NUM_VERT; i++ ) res += t[i]*values[g.elements[n][i]];
		return res;
	}
	return 1.0;
}

vec2d meshsurf2::getElementGradient(uint elm) const {
	const mesher2 &g = *this->g;
	vec2d res;
	for( uint dim=0; dim<DIM; dim++ ) {
		res[dim] = 0.0;
		for( uint k=0; k<NUM_VERT; k++ ) res[dim] += g.matrix[elm].m[k][dim]*values[g.elements[elm][k]];
	}
	return res;
}

vec2d meshsurf2::evalGradient(vec2d p) const {
	const mesher2 &g = *this->g;
	int n = g.hitElements(p);
	if( n >= 0 ) {
		return getElementGradient(n);
	}
	return vec2d();
}

void meshsurf2::drawGL( const levelset2 *solid, uint kind ) const {
	const mesher2 &g = *this->g;
	// Fill this levelset
	if( kind == 0 ) {
		glColor4d(0.3,0.3,0.8,0.6);
		for( uint n=0; n<g.elements.size(); n++ ) {
#if 1
			std::vector<vec2d> nodes(3);
			std::vector<FLOAT64> levelsets(3);
			for( uint k=0; k<3; k++ ) {
				uint idx = g.elements[n][k];
				nodes[k] = g.nodes[idx];
				levelsets[k] = fmax(values[idx],solid?-solid->evalLevelset(nodes[k]):-1e9 );
			}
			std::vector<vec2d> points = util2::marchPoints(nodes,levelsets);
			glBegin(GL_TRIANGLE_FAN);
			for( int m=0; m<points.size(); m++ ) glVertex2dv(points[m].v);
			glEnd();
#else
			glBegin(GL_TRIANGLES);
			for( uint k=0; k<g.elements[n].size(); k++ ) {
				uint idx = g.elements[n][k];
				FLOAT64 levelset = values[idx];
				glColor4d(levelset>0,0.0,levelset<0.0,7.0*fabs(levelset));
				glVertex2dv(g.nodes[idx].v);
			}
			glEnd();
#endif
		}
	} else {
#if 0
		// Draw closest pos
		for( uint n=0; n<surfacePos.size(); n++ ) {
			if( surfacePos[n].size() ) {
				FLOAT64 levelset = values[n];
				glColor4d(levelset>0,0.5,levelset<=0.0,1.0);
				glBegin(GL_LINES);
				glVertex2dv(surfacePos[n][0].v);
				glVertex2dv(surfacePos[n][1].v);
				glEnd();
			}
		}
#endif
	}
}