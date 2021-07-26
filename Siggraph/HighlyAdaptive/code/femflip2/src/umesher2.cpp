/*
 *	umesher2.cpp
 *	
 *	Created by Ryoichi Ando on 1/16/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "opengl.h"
#include "umesher2.h"
#include "levelset2.h"
using namespace std;

void umesher2::generateElements( std::vector<vec2d> &nodes, std::vector<std::vector<uint> > &elements, const levelset2 *hint, uint gn ) {
	
	// Generate base nodes and elements
	int w = 1.17*gn; if(w%2) w++;
	int h = 1.02*(gn)*4+1;
	if( this->w != w || this->h != h ) {
		this->w = w;
		this->h = h;
		
		FLOAT64 dx = 0.5/gn;
		FLOAT64 sqrt3 = sqrt(3.0);
		
		// Clear all vector arrays
		nodes_bank.clear();
		elements_bank.clear();
		nodes.clear();
		elements.clear();
		subdivided.clear();
		nodes_bank.clear();
		local_indices.clear();
		local_neighbors.clear();
		local_elements.clear();
		directions.clear();
		enabled.clear();
		status.clear();
		numBig = 0;
		
		vec2d maxp(0.0,0.0);
		for( uint j=0; j<h; j++ ) for( uint i=0; i<w+1; i++ ) {
			if(j%2==0) {
				vec2d p(sqrt3*dx*i,dx*j/2);
				nodes_bank.push_back(p);
				maxp[0] = fmax(p[0],maxp[0]);
				maxp[1] = fmax(p[1],maxp[1]);
			} else {
				vec2d p(sqrt3*dx*(i+0.5),dx*j/2);
				if(i<w) {
					nodes_bank.push_back(p);
					maxp[0] = fmax(p[0],maxp[0]);
					maxp[1] = fmax(p[1],maxp[1]);
				}
			}
		}
		
		// Centerize node positions
		vec2d shift(maxp[0]/2-0.5,maxp[1]/2-0.5);
		for( uint n=0; n<nodes_bank.size(); n++ ) {
			nodes_bank[n] -= shift;
		}
		
		// Build local_indices
		uint a = 0;
		uint idx = 0;
		int type = 1;
		while( a+(4*w+2) < nodes_bank.size() ) {
			type = 1-type;
			vector<int> indices;
			if(type==0) {
				indices.push_back(a);
				indices.push_back(a+w+1);
				indices.push_back(a+2*w+1);
				indices.push_back(a+2*w+2);
				indices.push_back(a+3*w+2);
				indices.push_back(a+4*w+2);
			} else {
				indices.push_back(a);
				indices.push_back(a+w);
				indices.push_back(a+2*w+1);
				indices.push_back(a+2*w);
				indices.push_back(a+3*w+1);
				indices.push_back(a+4*w+2);
			}
			local_indices.push_back(indices);
			directions.push_back(type);
			if( idx%w==w-1 ) {
				a += w+2;
				type = 1-type;
			} else {
				a += (type==0)*2;
			}
			idx ++;
		}
		
		// Build neighbor indices
		numBig = local_indices.size();
		local_neighbors.resize(numBig);
		
		for( int n=0; n<numBig; n++ ) {
			vector<int> neighbors(3);
			for( uint k=0; k<3; k++ ) neighbors[k] = -1;
			if( n-w >= 0) neighbors[0] = n-w;	// Bottom
			if( n-1>=0 && directions[n]==0 && directions[n-1]==1 ) {
				neighbors[1] = n-1;	// Neighbor (left)
			}
			if( n+1<numBig && directions[n]==1 && directions[n+1]==0 ) {
				neighbors[1] = n+1;	// Neighbor (right)
			}
			if(n+w < numBig) neighbors[2] = n+w; // Top
			local_neighbors[n] = neighbors;
		}
		
		// Build local elements
		idx = 0;
		local_elements.resize(numBig);
		local_status.resize(numBig);
		for( int n=0; n<numBig; n++ ) {
			vector<int> &indices = local_indices[n];
			vector<uint> element(3);
			element[0] = indices[0];
			element[1] = indices[3];
			element[2] = indices[5];
			elements_bank.push_back(element); // Big element
			local_elements[n].push_back(idx++);
			local_status[n].push_back(BIG);
			
			element[0] = indices[0];
			element[1] = indices[1];
			element[2] = indices[2];
			elements_bank.push_back(element); // Down small element
			local_elements[n].push_back(idx++);
			local_status[n].push_back(SUBDIV);
			
			element[0] = indices[1];
			element[1] = indices[4];
			element[2] = indices[2];
			elements_bank.push_back(element); // Middle small element
			local_elements[n].push_back(idx++);
			local_status[n].push_back(SUBDIV);
			
			element[0] = indices[1];
			element[1] = indices[3];
			element[2] = indices[4];
			elements_bank.push_back(element); // Neighbor small element
			local_elements[n].push_back(idx++);
			local_status[n].push_back(SUBDIV);
			
			element[0] = indices[2];
			element[1] = indices[4];
			element[2] = indices[5];
			elements_bank.push_back(element); // Up small element
			local_elements[n].push_back(idx++);
			local_status[n].push_back(SUBDIV);
			
			element[0] = indices[0];
			element[1] = indices[1];
			element[2] = indices[5];
			elements_bank.push_back(element); // Down connection element (0)
			local_elements[n].push_back(idx++);
			local_status[n].push_back(DOWN);
			
			element[0] = indices[1];
			element[1] = indices[3];
			element[2] = indices[5];
			elements_bank.push_back(element); // Down connection element (1)
			local_elements[n].push_back(idx++);
			local_status[n].push_back(DOWN);
			
			element[0] = indices[0];
			element[1] = indices[3];
			element[2] = indices[2];
			elements_bank.push_back(element); // Neighbor connection element (0)
			local_elements[n].push_back(idx++);
			local_status[n].push_back(NEIGHBOR);
			
			element[0] = indices[3];
			element[1] = indices[5];
			element[2] = indices[2];
			elements_bank.push_back(element); // Neighbor connection element (1)
			local_elements[n].push_back(idx++);
			local_status[n].push_back(NEIGHBOR);
			
			element[0] = indices[0];
			element[1] = indices[3];
			element[2] = indices[4];
			elements_bank.push_back(element); // Upper connection element (0)
			local_elements[n].push_back(idx++);
			local_status[n].push_back(UP);
			
			element[0] = indices[4];
			element[1] = indices[5];
			element[2] = indices[0];
			elements_bank.push_back(element); // Upper connection element (1)
			local_elements[n].push_back(idx++);
			local_status[n].push_back(UP);
		}
	}
	
	// Set status
	subdivided.resize(numBig);
	for( int n=0; n<numBig; n++ ) {
		vec2d center;
		uint elm = local_elements[n][0];
		for( int m=0; m<3; m++ ) {
			center += nodes_bank[elements_bank[elm][m]]/3.0;
		}
		subdivided[n] = hint ? hint->evalLevelset(center)<0.0 : false;
	}
	
	enabled.resize(elements_bank.size());
	for( int n=0; n<elements_bank.size(); n++ ) {
		enabled[n] = 0;
	}
	
	// Update elements enabled variable
	std::vector<uint> changedElements;
	updateElements(changedElements);
	
	nodes_enabled.resize(nodes_bank.size());
	for( int n=0; n<nodes_bank.size(); n++ ) nodes_enabled[n] = false;
	for( int n=0; n<elements_bank.size(); n++ ) {
		for( int m=0; m<elements_bank[n].size(); m++ ) {
			nodes_enabled[elements_bank[n][m]] = nodes_enabled[elements_bank[n][m]] || enabled[n];
		}
	}

	// Put those elements and nodes as result
	int numNodeEnabled = 0;
	int numElementEnabled = 0;
	for( int n=0; n<elements_bank.size(); n++ ) {
		if( enabled[n] ) numElementEnabled++;
	}
	for( int n=0; n<nodes_bank.size(); n++ ) {
		if( nodes_enabled[n] ) numNodeEnabled++;
	}
	elements.resize(numElementEnabled);
	nodes.resize(numNodeEnabled);
	vector<int>  nodeMap(nodes_bank.size());
	
	uint eidx = 0;
	uint nidx = 0;
	for( int n=0; n<nodes_bank.size(); n++ ) {
		if( nodes_enabled[n] ) {
			nodes[nidx] = nodes_bank[n];
			nodeMap[n] = nidx;
			nidx++;
		} else {
			nodeMap[n] = -1;
		}
	}
	
	elementMap.resize(numElementEnabled);
	for( int n=0; n<elements_bank.size(); n++ ) {
		if( enabled[n] ) {
			elements[eidx].clear();
			elements[eidx].resize(elements_bank[n].size());
			elementMap[eidx] = n;
			for( uint m=0; m<elements_bank[n].size(); m++ ) elements[eidx][m] = nodeMap[elements_bank[n][m]];
			eidx++;
		}
	}
}

void umesher2::updateElements( std::vector<uint> &changedElements ) {
	status.resize(numBig);
	for( int n=0; n<numBig; n++ ) {
		status[n] = subdivided[n] ? SUBDIV : BIG;
	}
	
	// Change the state according to the neighborhood mesh status
	bool changed = true;
	while(changed) {
		changed = false;
		vector<int> save_status = status;
		for( uint n=0; n<numBig; n++ ) {
			if( status[n] != SUBDIV ) {
				uint cnt = 0;
				for( uint m=0; m<3; m++ ) {
					int nidx = local_neighbors[n][m];
					if( nidx >= 0 && save_status[nidx] == SUBDIV ) cnt++;
				}
				if( cnt > 1 ) {
					status[n] = SUBDIV;
					changed = true;
				}
			}
		}
		
		for( uint n=0; n<numBig; n++ ) {
			if( status[n] != SUBDIV ) {
				for( uint m=0; m<3; m++ ) {
					int nidx = local_neighbors[n][m];
					if( nidx >= 0 && status[nidx] == SUBDIV ) {
						int new_status;
						if( m == 0 ) new_status = DOWN;
						else if( m == 1 ) new_status = NEIGHBOR;
						else if( m == 2 ) new_status = UP;
						if( new_status != status[n] ) {
							status[n] = new_status;
							changed = true;
						}
					}
				}
			}
		}
	}
	
	// Gather changed elements
	changedElements.clear();
	for( uint n=0; n<numBig; n++ ) {
		for( uint m=0; m<local_elements[n].size(); m++ ) {
			uint idx = local_elements[n][m];
			bool new_enabled = local_status[n][m]==status[n];
			if( new_enabled != enabled[idx] ) {
				enabled[idx] = new_enabled;
				changedElements.push_back(idx);
			}
		}
	}
}

static void drawBitmapString( const char *string, void *font=NULL ) {
	if( ! font ) font = GLUT_BITMAP_HELVETICA_10;
	while (*string) glutBitmapCharacter(font, *string++);
}

void umesher2::drawBigMesh(uint n) const {
	uint elm = local_elements[n][0];
	int idx0 = elements_bank[elm][0];
	int idx1 = elements_bank[elm][1];
	int idx2 = elements_bank[elm][2];
	glBegin(GL_POLYGON);
	glVertex2dv(nodes_bank[idx0].v);
	glVertex2dv(nodes_bank[idx1].v);
	glVertex2dv(nodes_bank[idx2].v);
	glEnd();
}

void umesher2::drawMesh() const {
	mesher2::drawMesh();
}