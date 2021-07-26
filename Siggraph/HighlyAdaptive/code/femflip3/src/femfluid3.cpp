/*
 *	femfluid3.cpp
 *	
 *	Created by Ryoichi Ando on 2/3/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "femfluid3.h"
#include "pcgsolver/pcg_solver.h"
#include "pcgsolver/matutil.h"
#include "levelset3.h"
#include "opengl.h"
#include "util3.h"
#include "fastmarch3.h"
#include <float.h>
#include <algorithm>
using namespace std;

femfluid3::femfluid3() {
	gn = 0;
	dx = 0.0;
    dt = 0.1;
	tension = 0.0;
    solid = NULL;
    fluid = NULL;
	surf_order = 2;
	variation = true;
	extrapolate_dist = 1.0;
	g.setGenFacet(false);
	g.setCenterType(g.BARYCENTRIC);
}

void femfluid3::init( uint gn, const levelset3 *hint ) {
	gn ++;	// Translate into nodal grid size
	this->gn = gn;
	dx = 1.0/(gn-1);
    initMesh(hint);
}

void femfluid3::setParameter( int name, FLOAT64 value ) {
	switch(name) {
		case SURF_ORDER:
			surf_order = value;
			break;
		case EXTRAPOLATE_DIST:
			extrapolate_dist = value;
			break;
		case VARIATION:
			variation = value;
			break;
	}
}

void femfluid3::initMesh(const levelset3 *hint) {
	tick(); dump( ">>> Building tetrahedra started...\n" );
	
	// Generate nodal points
	tick();
	g.generateMesh(gn,hint);
	dump( "Generated %d tets with %d vertices. BCC generation took %s.\n", g.elements.size(), g.nodes.size(), stock("FEM_BCC"));
	
	tick(); dump( "Allocating mesh array..." );
	// Allocate elements
	elements.clear();
	elements.resize(g.elements.size());
	
	// Allocate nodes
	nodes.clear();
	nodes.resize(g.nodes.size());
	
	// For each element...
    PARALLEL_FOR for( uint n=0; n<g.elements.size(); n++ ) {
		for( uint i=0; i<NUM_VERT; i++ ) for( uint j=0; j<NUM_VERT; j++ ) {
			elements[n].VBtB[i][j] = 0.0;
		}
		elements[n].velocity = vec3d();
		elements[n].gradient = vec3d();
		elements[n].mapped = false;
		elements[n].frac = 0.0;
    }
	
	// For each node...
	for( uint n=0; n<g.nodes.size(); n++ ) {
		nodes[n].fluidLS = 1.0;
		nodes[n].solidLS = 1.0;
		nodes[n].pressure = 0.0;
		nodes[n].divergence = 0.0;
		nodes[n].strain = 0.0;
	}
	dump( "Done. Took %s.\n", stock());
	
	// Compute node volume
	tick(); dump( "Computing tet volumes..." );
	for( uint n=0; n<g.nodes.size(); n++ ) {
		nodes[n].volume = 0.0;
		for ( uint m=0; m<g.node_elements[n].size(); m++ ) {
			nodes[n].volume += g.volumes[g.node_elements[n][m]] / g.node_elements[n].size();
		}
	}
	dump( "Done. Took %s.\n", stock());
	dump( "<<< Building tetrahedra done. Took %s\n", stock("FEM_tet") );
}

void femfluid3::setTimestep( FLOAT64 dt ) {
    this->dt = dt;
}

void femfluid3::setupSolidLevelset( const levelset3 *solid ) {
	this->solid = solid;
	PARALLEL_FOR for( uint n=0; n<nodes.size(); n++ ) {
		nodes[n].solidLS = solid->evalLevelset(g.nodes[n]);
	}
	PARALLEL_FOR for( uint n=0; n<g.elements.size(); n++ ) {
		elements[n].solidLS = solid->evalLevelset(g.centers[n]);
	}
	
	// Compute solid volume
	PARALLEL_FOR for( uint n=0; n<g.elements.size(); n++ ) {
		FLOAT64 LS[4];
		for( uint i=0; i<4; i++ ) {
			LS[i] = nodes[g.elements[n][i]].solidLS;
		}
		FLOAT64 volume = g.getLevelsetVolume(n,LS);
		elements[n].frac = 1.0-fmax(0.0,fmin(1.0,volume/(g.volumes[n]+LDBL_MIN)));
		if( ! variation ) elements[n].frac = elements[n].frac < 1.0 ? 0.0 : 1.0;
	}
}

void femfluid3::computeElementLiquidFraction() {
	// Then compute fluid fraction (density)
	for( uint n=0; n<elements.size(); n++ ) {
		FLOAT64 LS[NUM_VERT];
		if( g.elements[n].size() != NUM_VERT ) {
			dump( "elements size did not match ! %d != %d\n", g.elements[n].size(), NUM_VERT);
		}
		for( uint i=0; i<NUM_VERT; i++ ) {
			uint idx = g.elements[n][i];
			LS[i] = nodes[idx].fluidLS;
		}
		FLOAT64 volume = g.getLevelsetVolume(n,LS);
		elements[n].rho = fmax(0.0,fmin(1.0,volume/(g.volumes[n]+LDBL_MIN)));
	}
}

void femfluid3::setupFluidLevelset( const levelset3 *fluid ) {
	this->fluid = fluid;
	PARALLEL_FOR for( uint n=0; n<g.nodes.size(); n++ ) {
		nodes[n].fluidLS = fluid->evalLevelset(g.nodes[n]);
	}
	PARALLEL_FOR for( uint n=0; n<g.elements.size(); n++ ) {
		elements[n].fluidLS = fluid->evalLevelset(g.centers[n]);
	}
	
	// Compute liquid fraction for each element
	tick(); dump( "Computing fluid volume fraction..." );
	computeElementLiquidFraction();
	dump( "Done. Took %s.\n", stock() );
}

void femfluid3::getSamplePoints( std::vector<vec3d> &pos ) {
	pos = g.centers;
}

void femfluid3::setupVelocity( const std::vector<vec3d> &vel, const std::vector<bool> &mapped ) {
	for( uint n=0; n<elements.size(); n++ ) {
		if( elements[n].frac ) {
			elements[n].velocity = vel[n];
			elements[n].mapped = mapped[n];
		} else {
			elements[n].velocity = vec3d();
			elements[n].mapped = false;
		}
	}
}

static FLOAT64 get_dx( vec3d pos, const octree3 &octree ) {
	FLOAT64 base_dx = 1.0;
	int hit = octree.hitTest(pos);
	if( hit >= 0 ) {
		base_dx = fmin(base_dx,octree.terminals[hit]->dx/(FLOAT64)octree.resolution);
	} else {
		dump( "femfluid2 octree failed to fetch at (%f,%f,%f)\n", pos[0], pos[1], pos[2] );
	}
	return base_dx;
}

bool femfluid3::computeGhostPressure( uint elm_idx, uint index, FLOAT64 coef[NUM_VERT] ) const {
	FLOAT64 theta[NUM_VERT];
	
	uint node_idx = g.elements[elm_idx][index];
	for( uint n=0; n<NUM_VERT; n++ ) {
		coef[n] = 0.0;
		theta[n] = 0.0;
	}
	uint num_found = 0;
	for( uint n=0; n<NUM_VERT; n++ ) {
		uint neigh_idx = g.elements[elm_idx][n];
		if( node_idx != neigh_idx && nodes[neigh_idx].fluidLS < 0.0 ) {
			theta[n] = elements[elm_idx].VBtB[n][index];
			num_found ++;
		}
	}
	if( num_found == 1 ) {
		for( uint n=0; n<NUM_VERT; n++ ) {
			uint neigh_idx = g.elements[elm_idx][n];
			if( node_idx != neigh_idx && nodes[neigh_idx].fluidLS < 0.0 ) {
				theta[n] = 1.0;
				break;
			}
		}
	}
	
	FLOAT64 sum = 0.0;
	for( uint n=0; n<NUM_VERT; n++ ) sum += theta[n];
	if( sum ) {
		for( uint n=0; n<NUM_VERT; n++ ) {
			theta[n] /= sum;
		}
	} else {
		return false;
	}
	
	FLOAT64 det=0.0;
	vec3d pos;
	for( uint n=0; n<NUM_VERT; n++ ) {
		uint neigh_idx = g.elements[elm_idx][n];
		if( theta[n] ) {
			det += theta[n]*nodes[neigh_idx].fluidLS;
			pos += theta[n]*g.nodes[neigh_idx];
		}
	}
	if( det < 0.0 ) {
		det = fmin(det,-1e-2*get_dx(pos,g.octree));
		for( uint n=0; n<NUM_VERT; n++ ) {
			uint neigh_idx = g.elements[elm_idx][n];
			if( node_idx != neigh_idx && nodes[neigh_idx].fluidLS < 0.0 ) {
				FLOAT64 t = nodes[node_idx].fluidLS / det;
				coef[n] = t * theta[n];
			}
		}
		FLOAT64 scale_factor = 1.0;
		uint num_ghost_node = 0;
		for( uint n=0; n<NUM_VERT; n++ ) {
			if( nodes[g.elements[elm_idx][n]].fluidLS < 0.0 ) num_ghost_node ++;
		}
		if( num_ghost_node ) {
			for( uint n=0; n<NUM_VERT; n++ ) {
				if( nodes[g.elements[elm_idx][n]].fluidLS < 0.0 ) {
					FLOAT64 diag = elements[elm_idx].VBtB[n][n];
					FLOAT64 tol = 0.25;
					FLOAT64 entry = elements[elm_idx].VBtB[n][index];
					if( diag+coef[n]*entry < diag*tol ) {
						FLOAT64 c = coef[n]*entry;
						if( c ) {
							scale_factor = fmin(scale_factor,(tol*num_ghost_node-1.0)*diag/c);
						}
					}
				}
			}
			for( uint n=0; n<NUM_VERT; n++ ) {
				coef[n] *= scale_factor;
				theta[n] *= scale_factor;
			}
		}
	}
	return true;
}

void femfluid3::buildMatrix() {
	tick(); dump( "Building matrix...");
	
	LhsMatrix.clear();
	LhsMatrix.resize(nodes.size());
	RhsMatrix.clear();
	RhsMatrix.resize(nodes.size());
	
	FLOAT64 min_volume = 1e9;
	for( uint n=0; n<elements.size(); n++ ) {
		min_volume = fmin(min_volume,g.volumes[n]);
	}
	
	for( uint n=0; n<g.elements.size(); n++ ) {
		FLOAT64 B[DIM][NUM_VERT];
		FLOAT64 Bt[NUM_VERT][DIM];
		FLOAT64 BtB[NUM_VERT][NUM_VERT];
		FLOAT64 volume = elements[n].frac*(g.volumes[n]/min_volume);
		FLOAT64 rho = elements[n].rho;
		if( ! rho || ! volume ) continue;
		
		// B = [ dN_1/dx dN_2/dx dN_3/dx dN_4/dx ]
		//	   [ dN_1/dy dN_2/dy dN_3/dy dN_4/dy ]
		//	   [ dN_1/dz dN_2/dz dN_3/dz dN_4/dz ]
		for( uint dim=0; dim<DIM; dim++ ) {
			for( uint i=0; i<NUM_VERT; i++ ) B[dim][i] = g.matrix[n].m[i][dim];
		}
		// BtB = B^t * B
		bool foundBadTerm = false;
		for( uint i=0; i<NUM_VERT; i++ ) for( uint j=0; j<NUM_VERT; j++ ) {
			BtB[i][j] = 0.0;
			for( uint k=0; k<DIM; k++ ) BtB[i][j] += B[k][i]*B[k][j];
			if( i < j && BtB[i][j] > 1e-8 ) foundBadTerm = true;
		}
		// Save this to the laplace matrix
		for( uint i=0; i<NUM_VERT; i++ ) for( uint j=0; j<NUM_VERT; j++ ) {
			elements[n].VBtB[i][j] = volume*BtB[i][j];
		}
		// Bt = B^t
		for( uint i=0; i<NUM_VERT; i++ ) for( uint j=0; j<DIM; j++ ) {
			Bt[i][j] = B[j][i];
		}
		// Fill the global matrix...
		for( uint m=0; m<NUM_VERT; m++ ) {
			for( uint k=0; k<NUM_VERT; k++ ) {
				LhsMatrix.add_to_element(g.elements[n][m],g.elements[n][k],volume*BtB[m][k]);
			}
			for( uint dim=0; dim<DIM; dim++ ) {
				RhsMatrix.add_to_element(g.elements[n][m],n+dim*g.elements.size(),volume*Bt[m][dim]);
			}
		}
	}
	
	dump( "Done. Took %s.\n", stock("FEM_global_matrix"));
	SparseMatrix<FLOAT64> LhsMatrixOld = LhsMatrix;
	// Count the number of surface node
	uint num_surface_node = 0;
	for( uint n=0; n<LhsMatrix.n; n++ ) {
		if( nodes[n].fluidLS < 0.0 ) {
			for( uint m=0; m<LhsMatrix.index[n].size(); m++ ) {
				uint idx = LhsMatrix.index[n][m];
				if( nodes[idx].fluidLS > 0.0 ) {
					// Plug pressure surface boundary condition
					num_surface_node ++;
					break;
				}
			}
		}
	}
	
	tick(); dump( "Embedding dirichlet boundary condition with %d nodes", num_surface_node );
	uint bad_terms = 0;
	// Embed boundary condition
	LhsMatrix = LhsMatrixOld;
	if( surf_order == 2 ) {
		// Now embed second order accuracy
		for( uint n=0; n<elements.size(); n++ ) {
			for( uint m=0; m<NUM_VERT; m++ ) {
				uint idx = g.elements[n][m];
				// For ghost nodes
				if( nodes[idx].fluidLS > 0.0 ) {
					// Compute weight coefficients
					FLOAT64 coef[NUM_VERT];
					computeGhostPressure(n,m,coef);
					
					// Embed this coefs to adjacent fluid nodes
					for( uint k=0; k<NUM_VERT; k++ ) {
						if( k != m ) {
							uint neigh_idx = g.elements[n][k];
							if( nodes[neigh_idx].fluidLS < 0.0 ) {
								// Embed there
								FLOAT64 theta = elements[n].VBtB[m][k];
								LhsMatrix.set_element(neigh_idx,idx,0.0);
								for( uint i=0; i<NUM_VERT; i++ ) {
									if( i != m ) {
										// Embed into matrix
										LhsMatrix.add_to_element(neigh_idx,g.elements[n][i],theta*coef[i]);
									}
								}
							}
						}
					}
				}
			}
		}
		dump(".");
	} else {
		PARALLEL_FOR for( uint n=0; n<nodes.size(); n++ ) {
			for( uint m=0; m<LhsMatrix.index[n].size(); m++ ) {
				uint idx = LhsMatrix.index[n][m];
				if( nodes[n].fluidLS < 0.0 && nodes[idx].fluidLS > 0.0 ) {
					// Plug pressure surface boundary condition
					LhsMatrix.set_element(n,idx,0.0);
				}
			}
		}
		dump(".");
	}
	// Ensure that the matrix is nicely conditioned (should be synmmetric and positive defnite)
	std::vector<bool> flags(LhsMatrix.n);
	for( uint n=0; n<LhsMatrix.n; n++ ) {
		flags[n] = false;
	}
	if( LhsMatrix.n > nodes.size() ) {
		dump( "Nodes and matrix size unmatch !\n");
		exit(0);
	}
	for( uint n=0; n<LhsMatrix.n; n++ ) {
		for( uint m=0; m<LhsMatrix.index[n].size(); m++ ) {
			uint idx0 = n;
			uint idx1 = LhsMatrix.index[n][m];
			if( idx1 < LhsMatrix.n && idx0 <= idx1 && nodes[idx0].fluidLS < 0.0 ) {
				if( idx0 == idx1 ) {
					// If found bad matrix term, dump them and exit
					if( LhsMatrix.value[n][m] < 0.0 ) {
						dump( "WARNING: Matrix diagonal term at row %d is negative at (%f,%f,%f)!\n", n, g.nodes[n][0], g.nodes[n][1], g.nodes[n][2] );
						dump( "----- Matrix row: (+) out-of-water term (*) diagonal term ---- \n" );
						for( uint j=0; j<LhsMatrix.index[n].size(); j++ ) {
							if( n == LhsMatrix.index[n][j] ) dump( "*" );
							else if( nodes[LhsMatrix.index[n][j]].fluidLS > 0.0 ) dump( "+" );
							else printf( " " );
							dump( "M[][%d] = %f\n", LhsMatrix.index[n][j], LhsMatrix.value[n][j] );
						}
						dump( "Matrix row before embedding boundary conditions: \n" );
						for( uint j=0; j<LhsMatrixOld.index[n].size(); j++ ) {
							if( n == LhsMatrixOld.index[n][j] ) dump( "*" );
							else if( nodes[LhsMatrixOld.index[n][j]].fluidLS > 0.0 ) dump( "+" );
							else dump( " " );
							dump( "M[][%d] = %f\n", LhsMatrixOld.index[n][j], LhsMatrixOld.value[n][j] );
						}
						dump( "----- List of neighboring small VBtB matrix of the node %d -----\n", n );
						for( uint k=0; k<g.node_elements[n].size(); k++ ) {
							uint elm = g.node_elements[n][k];
							dump( "BtB[%d] = \n", k );
							if( n == g.elements[elm][0] ) dump( "*" );
							else if( nodes[g.elements[elm][0]].fluidLS > 0.0 ) dump( "+" );
							else email::print( " " );
							dump( "%d [%f %f %f %f;\n", g.elements[elm][0], elements[elm].VBtB[0][0], elements[elm].VBtB[0][1], elements[elm].VBtB[0][2], elements[elm].VBtB[0][3] );
							if( n == g.elements[elm][1] ) dump( "*" );
							else if( nodes[g.elements[elm][1]].fluidLS > 0.0 ) dump( "+" );
							else dump( " " );
							dump( "%d  %f %f %f %f;\n", g.elements[elm][1], elements[elm].VBtB[1][0], elements[elm].VBtB[1][1], elements[elm].VBtB[1][2], elements[elm].VBtB[1][3] );
							if( n == g.elements[elm][2] ) dump( "*" );
							else if( nodes[g.elements[elm][2]].fluidLS > 0.0 ) dump( "+" );
							else dump( " " );
							dump( "%d  %f %f %f %f;\n", g.elements[elm][2], elements[elm].VBtB[2][0], elements[elm].VBtB[2][1], elements[elm].VBtB[2][2], elements[elm].VBtB[2][3] );
							if( n == g.elements[elm][3] ) dump( "*" );
							else if( nodes[g.elements[elm][3]].fluidLS > 0.0 ) dump( "+" );
							else dump( " " );
							dump( "%d  %f %f %f %f]\n", g.elements[elm][3], elements[elm].VBtB[3][0], elements[elm].VBtB[3][1], elements[elm].VBtB[3][2], elements[elm].VBtB[3][3] );
							bool bad = false;
							FLOAT64 badone;
							for( uint i=0; i<4; i++ ) for( uint j=0; j<4; j++ ) {
								if( i < j ) {
									if( elements[elm].VBtB[i][j] > 1e-8 ) {
										badone = elements[elm].VBtB[i][j];
										bad = true;
									}
								}
							}
							if( bad ) {
								dump( "Bad tet. %f\n", badone );
							} else {
								dump( "Good tet.\n" );
							}
							dump( "---\n" );
						}
						dump("Bad matrix term !\n");
						flags[idx0] = true;
						bad_terms ++;
						exit(0);
					}
				} else if( nodes[idx1].fluidLS < 0.0 ) {
					FLOAT64 diff = fabs(LhsMatrix(idx1,idx0)-LhsMatrix.value[n][m]);
					if( diff > 1e-5 ) {
						dump( "WARNING: Matrix might not be symmetric ! (Residual = %e)\n", diff);
					}
				}
			}
		}
	}
	dump( "...Done. Took %s.\n", stock() );
}

void femfluid3::embedGhostPressure( uint elm_idx, FLOAT64 press[NUM_VERT] ) const {
	for( uint n=0; n<NUM_VERT; n++ ) {
		uint node_idx = g.elements[elm_idx][n];
		if( nodes[node_idx].fluidLS > 0.0 ) {
			FLOAT64 coef[NUM_VERT];
			FLOAT64 sum = 0.0;
			computeGhostPressure(elm_idx,n,coef);
			for( uint m=0; m<NUM_VERT; m++ ) {
				sum += press[m]*coef[m];
			}
			press[n] = sum;
		}
	}
}

void femfluid3::resample( uint type, uint kind ) {
	if( type == 0 ) {
		// Downsample to vertices
		PARALLEL_FOR for( uint n=0; n<g.nodes.size(); n++ ) {
			if( ! nodes[n].mapped ) {
				vec3d velocity;
				FLOAT64 wsum=0;
				if( kind == 0 ) nodes[n].velocity = vec3d();
				else nodes[n].gradient = vec3d();
				nodes[n].mapped = false;
				for( uint m=0; m<g.node_elements[n].size(); m++ ) {
					uint idx = g.node_elements[n][m];
					if( elements[idx].mapped ) {
						FLOAT64 w = g.volumes[idx];
						if( kind == 0 ) velocity += w*elements[idx].velocity;
						else velocity += w*elements[idx].gradient;
						wsum += w;
					}
				}
				if( wsum ) {
					if( kind == 0 ) nodes[n].velocity = velocity / wsum;
					else nodes[n].gradient = velocity / wsum;
					nodes[n].mapped = true;
				}
			}
		}
	} else if( type == 1 ) {
		// Upsample to elements
		PARALLEL_FOR for( uint n=0; n<g.elements.size(); n++ ) {
			if( ! elements[n].mapped ) {
				vec3d velocity;
				FLOAT64 wsum=0;
				if( kind == 0 ) elements[n].velocity = vec3d();
				else elements[n].gradient = vec3d();
				elements[n].mapped = false;
				for( uint m=0; m<g.elements[n].size(); m++ ) {
					uint idx = g.elements[n][m];
					if( nodes[idx].mapped ) {
						FLOAT64 w = nodes[idx].volume;
						if( kind == 0 ) velocity += w*nodes[idx].velocity;
						else velocity += w*nodes[idx].gradient;
						wsum += w;
					}
				}
				if( wsum ) {
					if( kind == 0 ) elements[n].velocity = velocity / wsum;
					else elements[n].gradient = velocity / wsum;
					elements[n].mapped = true;
				}
			}
		}
	}
}

void femfluid3::project() {
	tick(); dump( ">>> Extrapolating velocity...\n" );
	// Cut out solid wall velocity
	for( uint n=0; n<elements.size(); n++ ) {
		if( ! elements[n].frac ) {
			elements[n].mapped = false;
			elements[n].velocity = vec3d();
		}
	}
	
	FLOAT64 max_vel0 = 0.0;
	for( uint n=0; n<g.elements.size(); n++ ) {
		if( elements[n].mapped ) max_vel0 = fmax(max_vel0,elements[n].velocity.len());
	}
	
    // Downsample to nodes
	for( uint n=0; n<nodes.size(); n++ ) nodes[n].mapped = false;
	resample(0,0);
	
	// Extrapolate element velocity
	extrapolate(g,nodes,0);
	
	// Upsample to elements
	resample(1,0);
	
	FLOAT64 max_vel1 = 0.0;
	for( uint n=0; n<g.elements.size(); n++ ) {
		if( elements[n].mapped )  max_vel1 = fmax(max_vel1,elements[n].velocity.len());
	}
	
	dump( "<<< Done. Took %s.\n", stock("FEM_vel_extrapolate") );
	
	// Build matrix
	tick(); dump( ">>> Building matrix started...\n" );
	buildMatrix();
	dump( "<<< Done. Took %s.\n", stock("FEM_matrix"));
	
	std::vector<FLOAT64> rhs(g.nodes.size());
	std::vector<FLOAT64> elmv(g.elements.size()*DIM);
	FLOAT64 max_vel = 0.0;
	for( uint n=0; n<g.elements.size(); n++ ) {
		elmv[n+0*g.elements.size()] = elements[n].velocity[0];
		elmv[n+1*g.elements.size()] = elements[n].velocity[1];
		elmv[n+2*g.elements.size()] = elements[n].velocity[2];
		max_vel = fmax(max_vel,elements[n].velocity.len());
	}
	
	tick(); dump( "Computing divergence (right hand side)..." );
	
	// Compute divergence
	multiply<FLOAT64>(RhsMatrix,elmv,rhs);
	
	// Copy divergence
	for( uint n=0; n<nodes.size(); n++ ) {
		nodes[n].divergence = rhs[n];
	}
	
	SparseMatrix<FLOAT64> LhsMatrix = this->LhsMatrix;

	// Check Nan value (for debug)
#if 1
	for( uint n=0; n<rhs.size(); n++ ) {
		if( is_nan(rhs[n]) ) {
			email::print( "Right hand side contains NaN value !\n" );
			email::send();
			exit(0);
		}
	}
#endif
	
	// Cut out air nodes
	for( uint n=0; n<nodes.size(); n++ ) {
		if( nodes[n].fluidLS > 0.0 ) {
			LhsMatrix.index[n].clear();
			LhsMatrix.value[n].clear();
		} 
	}
	dump( "Done. Took %s.\n", stock());
	tick(); dump( "Organizing the matrix (packing the matrix)..." );
	
	// Compact the matrix
	std::vector<int> idx_map;
	SparseMatrix<FLOAT64> rawLhsMatrix = LhsMatrix;
	compress( LhsMatrix, rhs, idx_map );
	
	// Solve by incomplete cholesky factorizaion preconditioned conjugate gradient method
    PCGSolver<FLOAT64> solver;
	solver.set_solver_parameters(1e-11, 150000, 0.97, 0.25);
    FLOAT64 residual_out;
	std::vector<FLOAT64> result(LhsMatrix.n);
	for( uint n=0; n<nodes.size(); n++ ) {
		if( idx_map[n] >= 0 ) result[idx_map[n]] = nodes[n].pressure;
	}
	
	dump( "Done. Took %s.\n", stock());
	dump( "Solving FEM pressure..." );
    int iterations;
	FLOAT64 msec;
	bool converged = solver.solve( LhsMatrix, rhs, result, residual_out, iterations, msec );
	if( ! converged ) {
		email::print("PCG did not converge. Residual = %e\n", residual_out );
		
		// Dump this matrix
		ofstream file;
		file.open(format_str("%s/Lhs_FEM_matrix.m",root_path));
		LhsMatrix.write_matlab(file,LhsMatrix.n,LhsMatrix.n,"Lhs_FEM");
		file.close();
		file.open(format_str("%s/Rhs_FEM_vector.m",root_path));
		write_matlab1d(file,rhs,"Rhs_FEM");
		file.close();
		
		email::send();
		exit(0);
	} else {
		dump("Done. Took %d iterations with %d unknowns. Residual = %e. Took %s.\n", iterations, LhsMatrix.n, residual_out, tstr(msec));
		writeNumber("FEM_PCG_iterations",iterations);
		writeNumber("FEM_PCG_time",msec);
		writeNumber("FEM_PCG_unknowns", LhsMatrix.n);
		writeNumber("FEM_PCG_residual",residual_out);
    }

	// Set pressure
	PARALLEL_FOR for( uint n=0; n<nodes.size(); n++ ) {
		nodes[n].pressure = idx_map[n] >= 0 ? result[idx_map[n]] : 0.0;
	}
	
	tick(); dump( "Computing pressure gradient..." );
	
	// Take gradient and subtract
	FLOAT64 max_grad = 0.0;
	PARALLEL_FOR for( uint n=0; n<elements.size(); n++ ) {
		elements[n].gradient = vec3d();
		elements[n].mapped = false;
		FLOAT64 rho = elements[n].rho;
		if( ! rho || ! elements[n].frac ) continue;
		
		uint idx[NUM_VERT] = { g.elements[n][0], g.elements[n][1], g.elements[n][2], g.elements[n][3] };
		FLOAT64 press[NUM_VERT] = { nodes[idx[0]].pressure, nodes[idx[1]].pressure, nodes[idx[2]].pressure, nodes[idx[3]].pressure };
		if( surf_order == 2 ) embedGhostPressure(n,press);
		
		vec3d value;
		for( uint i=0; i<NUM_VERT; i++ ) {
			for( uint dim=0;dim<DIM;dim++) {
				value[dim] += g.matrix[n].m[i][dim]*press[i];
			}
		}
		elements[n].gradient = value/dt;
		elements[n].mapped = true;
		max_grad = fmax(max_grad,elements[n].gradient.len());
	}
	
	dump( "Done. Took %s.\n", stock() );
	
	tick(); dump( ">>> Extrapolating pressure gradient...\n" );
	
	FLOAT64 max_grad0 = 0.0;
	for( uint n=0; n<g.elements.size(); n++ ) {
		if( elements[n].mapped ) max_grad0 = fmax(max_grad0,elements[n].gradient.len());
	}
	
	// Downsample to nodes
	for( uint n=0; n<nodes.size(); n++ ) nodes[n].mapped = false;
	resample(0,1);
	
	// Extrapolate element gradient
	extrapolate(g,nodes,1);
	
	// Upsample to elements
	resample(1,1);
	
	FLOAT64 max_grad1 = 0.0;
	for( uint n=0; n<g.elements.size(); n++ ) {
		if( elements[n].mapped ) max_grad1 = fmax(max_grad1,elements[n].gradient.len());
	}
	dump( "<<< Done. Took %s.\n", stock("FEM_gradient_extrapolate"));
	
	// Compute strain
	tick(); dump( "Computing strain..." );
	computeStrain();
	dump( "Done. Took %s.\n", stock("FEM_strain"));
}

void femfluid3::extrapolate( const mesher3 &g, std::vector<node3> &nodes, int kind ) {
	FLOAT64 max_vel0 = 0.0;
	for( uint n=0; n<g.nodes.size(); n++ ) {
		if( nodes[n].mapped ) max_vel0 = fmax(max_vel0, kind==0 ? nodes[n].velocity.len() : nodes[n].gradient.len());
	}
	
	std::vector<fastmarch3<vec3d>::node3 *> fnodes(g.nodes.size());
	for( uint n=0; n<g.nodes.size(); n++ ) fnodes[n] = new fastmarch3<vec3d>::node3;
	for( uint n=0; n<g.nodes.size(); n++ ) {
		vec3d p = g.nodes[n];
		bool fixed = nodes[n].mapped;
		fnodes[n]->p = p;
		fnodes[n]->fixed = fixed;
		fnodes[n]->levelset = fmax(nodes[n].fluidLS,-nodes[n].solidLS);
		fnodes[n]->value = (kind==0) ? nodes[n].velocity : nodes[n].gradient;
		fnodes[n]->p2p.resize(g.node2node[n].size());
		for( uint m=0; m<g.node2node[n].size(); m++ ) {
			fnodes[n]->p2p[m] = fnodes[g.node2node[n][m]];
		}
	}
	fastmarch3<vec3d>::fastMarch(fnodes,extrapolate_dist,-1.0,0);
	// Pick up extrapolated values
	for( uint n=0; n<g.nodes.size(); n++ ) {
		if( fnodes[n]->fixed ) {
			if (kind==0) nodes[n].velocity = fnodes[n]->value;
			else nodes[n].gradient = fnodes[n]->value;
			nodes[n].mapped = true;
		}
		delete fnodes[n];
	}
	
	FLOAT64 max_vel1 = 0.0;
	for( uint n=0; n<g.nodes.size(); n++ ) {
		if( nodes[n].mapped ) max_vel1 = fmax(max_vel1, kind==0 ? nodes[n].velocity.len() : nodes[n].gradient.len());
	}
	if( max_vel1 > 1.1*max_vel0 ) {
		dump( "Extrapolated quantify has larger value (%g > %g). Seems like you have bugs in extrapolation code. Now clamping...\n", max_vel1, max_vel0 );
		for( uint n=0; n<g.nodes.size(); n++ ) {
			if( nodes[n].mapped ) {
				if( kind == 0 ) {
					if( nodes[n].velocity.len() > max_vel1 ) nodes[n].velocity = max_vel0 * nodes[n].velocity.normal();
				} else {
					if( nodes[n].gradient.len() > max_vel1 ) nodes[n].gradient = max_vel0 * nodes[n].gradient.normal();
				}
			}
		}
	}
}

vec3d femfluid3::getVector( vec3d p, uint kind ) const {
	int hit = g.hitElements(p);
	if( hit >= 0 ) {
		for( uint n1=0; n1<NUM_VERT; n1++ ) {
			for( uint n2=n1+1; n2<NUM_VERT; n2++ ) {
				for( uint n3=n2+1; n3<NUM_VERT; n3++ ) {
					// Put this subdivision triangle
					std::vector<vec3d> tet(NUM_VERT);
					tet[0] = g.centers[hit];
					tet[1] = g.nodes[g.elements[hit][n1]];
					tet[2] = g.nodes[g.elements[hit][n2]];
					tet[3] = g.nodes[g.elements[hit][n3]];
					
					// Precompute a shape function here
					FLOAT64 M[NUM_VERT][NUM_VERT];
					FLOAT64 SM[NUM_VERT][NUM_VERT];
					for( uint k=0; k<NUM_VERT; k++ ) {
						for( uint dim=0; dim<DIM; dim++ ) {
							M[dim][k] = tet[k][dim];
						}
						M[DIM][k] = 1.0;
					}
					if( ! invert4x4( M, SM )) {
						dump( "Failed inverse 4x4! (subdividion part)\n" );
					}
					
					vec3d vec;
					bool skip = false;
					FLOAT64 x[NUM_VERT] = { p[0], p[1], p[2], 1.0 };
					vec3d v[NUM_VERT];
					if( kind == 0 ) {
						v[0] = elements[hit].velocity;
						v[1] = nodes[g.elements[hit][n1]].velocity;
						v[2] = nodes[g.elements[hit][n2]].velocity;
						v[3] = nodes[g.elements[hit][n3]].velocity;
					} else {
						v[0] = elements[hit].gradient;
						v[1] = nodes[g.elements[hit][n1]].gradient;
						v[2] = nodes[g.elements[hit][n2]].gradient;
						v[3] = nodes[g.elements[hit][n3]].gradient;
					}
					for( uint i=0; i<NUM_VERT; i++ ) {
						FLOAT64 t = 0.0;
						for( uint j=0;j<NUM_VERT;j++) {
							t += SM[i][j]*x[j];
						}
						if( t < -1e-4 || t > 1.0+1e-4 ) {
							skip = true;
							break;
						}
						vec += t*v[i];
					}
					if( ! skip ) return vec;
				}
			}
		}
		dump( "getVector failed kind = %d hit = %d (%f,%f,%f)\n", kind, hit, p[0], p[1], p[2]);
		return kind==0 ? elements[hit].velocity : elements[hit].gradient;
	}
	dump( "getVector failed kind = %d hit = %d (%f,%f,%f)\n", kind, hit, p[0], p[1], p[2]);
	return vec3d();
}

vec3d femfluid3::getPressureGradient( vec3d p ) const {
	return getVector(p,1);
}

vec3d femfluid3::getVelocity( vec3d p ) const {
	return getVector(p,0);
}

FLOAT64 femfluid3::getDivergence( vec3d p ) const {
	int n = g.hitElements(p);
	if( n >= 0 ) {
		FLOAT64 x[NUM_VERT] = { p[0], p[1], p[2], 1.0 };
		FLOAT64 div = 0.0;
		for( uint i=0; i<NUM_VERT; i++ ) {
			FLOAT64 t = 0.0;
			for( uint j=0; j<NUM_VERT; j++ ) {
				t += g.matrix[n].m[i][j]*x[j];
			}
			div += t*nodes[g.elements[n][i]].divergence;
		}
		return div;
	}
	return 0.0;
}

void femfluid3::computeStrain() {
	// Compute strain per tet
	std::vector<FLOAT64> strain_tet(elements.size());
	PARALLEL_FOR for( uint n=0; n<elements.size(); n++ ) {
		FLOAT64 strain_tensor[DIM][DIM];
		bool mapped = true;
		for( uint m=0; m<NUM_VERT; m++) {
			if( ! nodes[g.elements[n][m]].mapped ) {
				mapped = false;
				break;
			}
		}
		if( mapped ) {
			FOR_EACH2(DIM,DIM) {
				strain_tensor[i][j] = 0.0;
			} END_FOR
			for( uint dim_i=0; dim_i<DIM; dim_i++ ) {
				FLOAT64 u[DIM];
				for( uint m=0; m<NUM_VERT; m++ ) u[m] = 0.0;
				for( uint i=0; i<NUM_VERT; i++ ) {
					for( uint dim_j=0;dim_j<DIM;dim_j++) {
						u[dim_j] += g.matrix[n].m[i][dim_j]*(nodes[g.elements[n][i]].velocity[dim_i]-dt*nodes[g.elements[n][i]].gradient[dim_i]);
					}
				}
				for( uint dim_j=0;dim_j<DIM;dim_j++) {
					strain_tensor[dim_i][dim_j] += 0.5*u[dim_j];
					strain_tensor[dim_j][dim_i] += 0.5*u[dim_j];
				}
			}
			FLOAT64 strain2 = 0.0;
			FOR_EACH2(DIM,DIM) {
				strain2 += sqr(strain_tensor[i][j]);
			} END_FOR
			strain_tet[n] = sqrtf(strain2);
		} else {
			strain_tet[n] = 0.0;
		}
	}
	
	// Upsample to nodes and compute strain difference there
	PARALLEL_FOR for( uint n=0; n<nodes.size(); n++ ) {
		FLOAT64 avg_strain = 0.0;
		uint num = g.node_elements[n].size();
		for( uint m=0; m<num; m++ ) {
			avg_strain += strain_tet[g.node_elements[n][m]];
		}
		avg_strain /= num;
		nodes[n].strain = avg_strain;
	}
}

FLOAT64 femfluid3::getStrain( vec3d p ) const {
	if( nodes.empty()) return 0.0;
	int n = g.hitElements(p);
	if( n >= 0 ) {
		FLOAT64 x[NUM_VERT] = { p[0], p[1], 1.0 };
		FLOAT64 strain = 0.0;
		for( uint i=0; i<NUM_VERT; i++ ) {
			FLOAT64 t = 0.0;
			for( uint j=0; j<NUM_VERT; j++ ) {
				t += g.matrix[n].m[i][j]*x[j];
			}
			strain += t*nodes[g.elements[n][i]].strain;
		}
		return strain;
	}
	return 0.0;
}

static void drawBitmapString( const char *string, void *font=NULL ) {
	if( ! font ) font = GLUT_BITMAP_HELVETICA_10;
	while (*string) glutBitmapCharacter(font, *string++);
}

void femfluid3::render( int name ) const {
	switch(name) {
		case MESH:
			glColor4f(1.0,1.0,1.0,1.0);
			g.drawMesh();
			break;
		case NAME: {
			glColor4f(1.0,1.0,1.0,1.0);
			glRasterPos2d(0.07,0.92);
			drawBitmapString(getName(),GLUT_BITMAP_HELVETICA_18);
			break;
		}
	}
}
