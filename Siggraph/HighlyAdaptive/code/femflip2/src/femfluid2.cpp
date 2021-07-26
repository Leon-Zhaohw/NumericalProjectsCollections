/*
 *	femfluid2.cpp
 *	
 *	Created by Ryoichi Ando on 12/9/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "femfluid2.h"
#include "pcgsolver/pcg_solver.h"
#include "pcgsolver/matutil.h"
#include "levelset2.h"
#include "opengl.h"
#include "util2.h"
#include "fastmarch2.h"
#include <float.h>
#include <algorithm>
using namespace std;

femfluid2::femfluid2() {
	gn = 0;
	dx = 0.0;
    dt = 0.1;
    solid = NULL;
    fluid = NULL;
	surf_order = 2;
	variation = true;
	extrapolate_dist = 1.0;
	g.setGenFacet(false);
	g.setCenterType(g.BARYCENTRIC);
}

void femfluid2::init( uint gn, const levelset2 *hint ) {
	gn ++;	// Translate into nodal grid size
	this->gn = gn;
	dx = 1.0/(gn-1);
    initMesh(hint);
}

void femfluid2::setParameter( int name, FLOAT64 value ) {
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

void femfluid2::initMesh(const levelset2 *hint) {
	// Generate nodal points
	g.generateMesh(gn,hint);
	
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
		elements[n].velocity = vec2d();
		elements[n].gradient = vec2d();
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
	
	// Compute node volume
	for( uint n=0; n<g.nodes.size(); n++ ) {
		nodes[n].volume = 0.0;
		for ( uint m=0; m<g.node_elements[n].size(); m++ ) {
			nodes[n].volume += g.volumes[g.node_elements[n][m]] / g.node_elements[n].size();
		}
	}
}

void femfluid2::setTimestep( FLOAT64 dt ) {
    this->dt = dt;
}

void femfluid2::setupSolidLevelset( const levelset2 *solid ) {
	this->solid = solid;
	// Set solid levelset for each node
	PARALLEL_FOR for( uint n=0; n<nodes.size(); n++ ) {
		nodes[n].solidLS = solid->evalLevelset(g.nodes[n]);
	}
	PARALLEL_FOR for( uint n=0; n<g.elements.size(); n++ ) {
		elements[n].solidLS = solid->evalLevelset(g.centers[n]);
	}
	
	// Compute solid volume
	for( uint n=0; n<g.elements.size(); n++ ) {
		uint num = g.elements[n].size();
		vector<vec2d> points(num);
		vector<FLOAT64> LS(num);
		for( uint i=0; i<num; i++ ) {
			uint idx = g.elements[n][i];
			points[i] = g.nodes[idx];
			LS[i] = nodes[idx].solidLS;
		}
		FLOAT64 volume = util2::computeVolume(util2::marchPoints(points,LS));
		elements[n].frac = 1.0-fmax(0.0,fmin(1.0,volume/(g.volumes[n]+LDBL_MIN)));
		if( ! variation ) elements[n].frac = elements[n].frac < 1.0 ? 0.0 : 1.0;
	}
}

static FLOAT64 fraction( FLOAT64 phi0, FLOAT64 phi1 ) {
	return phi0*phi1 > 0 ? 1 - (phi0 > 0) : -fmin(phi0,phi1)/fmax(fabs(phi1-phi0),1.0e-18);
}

void femfluid2::computeElementLiquidFraction() {
	// Then compute fluid fraction (density)
	for( uint n=0; n<elements.size(); n++ ) {
		uint num = g.elements[n].size();
		vector<vec2d> points(num);
		vector<FLOAT64> LS(num);
		for( uint i=0; i<num; i++ ) {
			uint idx = g.elements[n][i];
			points[i] = g.nodes[idx];
			LS[i] = nodes[idx].fluidLS;
		}
		FLOAT64 volume = util2::computeVolume(util2::marchPoints(points,LS));
		elements[n].rho = fmax(0.0,fmin(1.0,volume/(g.volumes[n]+LDBL_MIN)));
	}
}

void femfluid2::setupFluidLevelset( const levelset2 *fluid ) {
	this->fluid = fluid;
	// Set fluid levelset for each node
	PARALLEL_FOR for( uint n=0; n<nodes.size(); n++ ) {
		nodes[n].fluidLS = fluid->evalLevelset(g.nodes[n]);
	}
	PARALLEL_FOR for( uint n=0; n<g.elements.size(); n++ ) {
		elements[n].fluidLS = fluid->evalLevelset(g.centers[n]);
	}
	// Compute fraction
	computeElementLiquidFraction();
}

void femfluid2::getSamplePoints( std::vector<vec2d> &pos ) {
	pos = g.centers;
}

void femfluid2::setupVelocity( const std::vector<vec2d> &vel, const std::vector<bool> &mapped ) {
	for( uint n=0; n<elements.size(); n++ ) {
		if( elements[n].frac ) {
			elements[n].velocity = vel[n];
			elements[n].mapped = mapped[n];
		} else {
			elements[n].velocity = vec2d();
			elements[n].mapped = false;
		}
	}
}

static FLOAT64 get_dx( vec2d pos, const octree2 &octree ) {
	FLOAT64 base_dx = 1.0;
	int hit = octree.hitTest(pos);
	if( hit >= 0 ) {
		base_dx = fmin(base_dx,octree.terminals[hit]->dx/(FLOAT64)octree.resolution);
	} else {
		dump( "femfluid2 octree failed to fetch at (%f,%f)\n", pos[0], pos[1] );
	}
	return base_dx;
}

bool femfluid2::computeGhostPressure( uint elm_idx, uint index, FLOAT64 coef[NUM_VERT] ) const {
	FLOAT64 theta[NUM_VERT];
	
	uint node_idx = g.elements[elm_idx][index];
	for( uint n=0; n<NUM_VERT; n++ ) {
		coef[n] = 0.0;
		theta[n] = 0.0;
	}
	uint num_found = 0;
	for( uint n=0; n<NUM_VERT; n++ ) {
		uint neigh_idx = g.elements[elm_idx][n];
		if( node_idx != neigh_idx && nodes[neigh_idx].fluidLS <= 0.0 ) {
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
	vec2d pos;
	for( uint n=0; n<NUM_VERT; n++ ) {
		uint neigh_idx = g.elements[elm_idx][n];
		if( theta[n] ) {
			det += theta[n]*nodes[neigh_idx].fluidLS;
			pos += theta[n]*g.nodes[neigh_idx];
		}
	}
	if( det < 0.0 ) {
		FLOAT64 base_dx;
		if( getOctree() ) base_dx = get_dx(pos,*getOctree());
		else base_dx = dx;
		det = fmin(det,-1e-2*base_dx);
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
		}		
		for( uint n=0; n<NUM_VERT; n++ ) {
			coef[n] *= scale_factor;
			theta[n] *= scale_factor;
		}
	}
	return true;
}

void femfluid2::buildMatrix() {
	LhsMatrix.clear();
	LhsMatrix.resize(nodes.size());
	RhsMatrix.clear();
	RhsMatrix.resize(nodes.size());
	
	FLOAT64 min_volume = 1e9;
	for( uint n=0; n<elements.size(); n++ ) {
		min_volume = fmin(min_volume,g.volumes[n]);
	}
	
	for( uint n=0; n<elements.size(); n++ ) {
		FLOAT64 B[DIM][NUM_VERT];
		FLOAT64 Bt[NUM_VERT][DIM];
		FLOAT64 BtB[NUM_VERT][NUM_VERT];
		FLOAT64 area = elements[n].frac*(g.volumes[n]/min_volume);
		FLOAT64 rho = elements[n].rho;
		if( ! rho || ! area ) continue;
		
		// B = [ dN_1/dx dN_2/dx dN_3/dx ]
		//     [ dN_1/dy dN_2/dy dN_3/dy ]
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
			elements[n].VBtB[i][j] = area*BtB[i][j];
		}
		// Bt = B^t
		for( uint i=0; i<NUM_VERT; i++ ) for( uint j=0; j<DIM; j++ ) {
			Bt[i][j] = B[j][i];
		}
		// Fill the global matrix...
		for( uint m=0; m<NUM_VERT; m++ ) {
			for( uint k=0; k<NUM_VERT; k++ ) {
				LhsMatrix.add_to_element(g.elements[n][m],g.elements[n][k],area*BtB[m][k]);
			}
			for( uint dim=0; dim<DIM; dim++ ) {
				RhsMatrix.add_to_element(g.elements[n][m],n+dim*elements.size(),area*Bt[m][dim]);
			}
		}
	}
	
	// Embed boundary condition
	SparseMatrix<FLOAT64> LhsMatrixOld = LhsMatrix;
	LhsMatrix = LhsMatrixOld;
	
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
	} else {
		for( uint n=0; n<LhsMatrix.n; n++ ) {
			for( uint m=0; m<LhsMatrix.index[n].size(); m++ ) {
				uint idx = LhsMatrix.index[n][m];
				if( nodes[n].fluidLS < 0.0 && nodes[idx].fluidLS > 0.0 ) {
					// Plug pressure surface boundary condition
					LhsMatrix.set_element(n,idx,0.0);
				}
			}
		}
	}
	
	// Ensure that the matrix is nicely conditioned (should be synmmetric and positive defnite)
	std::vector<bool> flags(LhsMatrix.n);
	for( uint n=0; n<LhsMatrix.n; n++ ) {
		flags[n] = false;
	}
	for( uint n=0; n<LhsMatrix.n; n++ ) {
		for( uint m=0; m<LhsMatrix.index[n].size(); m++ ) {
			uint idx0 = n;
			uint idx1 = LhsMatrix.index[n][m];
			if( idx1 < LhsMatrix.n && idx0 <= idx1 && nodes[idx0].fluidLS < 0.0 ) {
				if( idx0 == idx1 ) {
					if( LhsMatrix.value[n][m] < 0.0 ) {
						dump( "WARNING: Matrix diagonal term is negative ! n=%d m=%d\n", n, m );
						dump( "----- Matrix row: (+) out-of-water term (*) diagonal term ---- \n" );
						for( uint j=0; j<LhsMatrix.index[n].size(); j++ ) {
							if( n == LhsMatrix.index[n][j] ) dump( "*" );
							else if( nodes[LhsMatrix.index[n][j]].fluidLS > 0.0 ) dump( "+" );
							else printf( " " );
							dump( "M[][%d] = %f\n", LhsMatrix.index[n][j], LhsMatrix.value[n][j] );
						}
						dump( "Matrix row before embedding boundary conditions: \n" );
						for( uint j=0; j<LhsMatrixOld.index[n].size(); j++ ) {
							if( n == LhsMatrix.index[n][j] ) dump( "*" );
							else if( nodes[LhsMatrix.index[n][j]].fluidLS > 0.0 ) dump( "+" );
							else dump( " " );
							dump( "M[][%d] = %f\n", LhsMatrixOld.index[n][j], LhsMatrixOld.value[n][j] );
						}
						flags[idx0] = true;
						break;
					}
				} else if( nodes[idx1].fluidLS < 0.0 ) {
					FLOAT64 diff = fabs(LhsMatrix(idx1,idx0)-LhsMatrix.value[n][m]);
					if( diff > 1e-5 ) {
						dump( "Matrix is not symmetric ! (Residual = %e)\n", diff);
						exit(0);
					}
				}
			}
		}
		
	}
	for( uint n=0; n<LhsMatrix.n; n++ ) {
		if( flags[n] ) {
			nodes[n].fluidLS = 0.0;
		}
	}
}

void femfluid2::embedGhostPressure( uint elm_idx, FLOAT64 press[NUM_VERT] ) const {
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

void femfluid2::resample( uint type, uint kind ) {
	if( type == 0 ) {
		// Downsample to vertices
		PARALLEL_FOR for( uint n=0; n<g.nodes.size(); n++ ) {
			if( ! nodes[n].mapped ) {
				vec2d velocity;
				FLOAT64 wsum=0;
				if( kind == 0 ) nodes[n].velocity = vec2d();
				else nodes[n].gradient = vec2d();
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
				vec2d velocity;
				FLOAT64 wsum=0;
				if( kind == 0 ) elements[n].velocity = vec2d();
				else elements[n].gradient = vec2d();
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

void femfluid2::project() {
	// Cut out solid wall velocity
    for( uint n=0; n<elements.size(); n++ ) {
        if( ! elements[n].frac ) {
            elements[n].mapped = false;
			elements[n].velocity = vec2d();
        }
    }
	
	// Downsample to nodes
	for( uint n=0; n<nodes.size(); n++ ) nodes[n].mapped = false;
	resample(0,0);

	// Extrapolate element velocity
	extrapolate(g,nodes,0);
	
	// Upsample to elements
	resample(1,0);
	
	// Build matrix
	buildMatrix();
	
	// Compute divergence
	vector<FLOAT64> rhs(nodes.size());
	vector<FLOAT64> elmv(elements.size()*DIM);
	for( uint n=0; n<elements.size(); n++ ) {
		elmv[n] = elements[n].velocity[0];
		elmv[n+elements.size()] = elements[n].velocity[1];
	}
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

	// Compact the matrix
	vector<int> idx_map;
	SparseMatrix<FLOAT64> rawLhsMatrix = LhsMatrix;
	compress( LhsMatrix, rhs, idx_map );
	
	// Solve by incomplete cholesky factorizaion preconditioned conjugate gradient method
	PCGSolver<FLOAT64> solver;
	solver.set_solver_parameters(1e-11, 1500000, 0.97, 0.25);
	FLOAT64 residual_out;
	vector<FLOAT64> result(LhsMatrix.n);
	for( uint n=0; n<nodes.size(); n++ ) {
		if( idx_map[n] >= 0 ) result[idx_map[n]] = nodes[n].pressure;
	}
		
#if 0
	ofstream file;
	file.open("Lhs_FEM_matrix.m");
	LhsMatrix.write_matlab(file,LhsMatrix.n,LhsMatrix.n,"Lhs_FEM");
	file.close();
	file.open("Rhs_FEM_vector.m");
	write_matlab1d(file,rhs,"Rhs_FEM");
	file.close();
#endif
	
	int iterations;
	FLOAT64 msec;
	bool converged = solver.solve( LhsMatrix, rhs, result, residual_out, iterations, msec );
	if( converged && residual_out > 1.0 ) {
		printf( "Something is wrong with PCG\n");
		exit(0);
		converged = false;
	}
	if( ! converged ) {
		email::print("PCG did not converge.\n" );
		email::send();
		exit(0);
	} else {
		dump( "FEM: PCG Converged ! %d iterations and took %.3f msec with %d unknowns.\n", iterations, msec, LhsMatrix.n);
	}
	if( residual_out > 1.0 ) {
		email::print("PCG did not converge. (Residual too large)\n" );
		email::send();
		exit(0);
	}
	// Set pressure
	PARALLEL_FOR for( uint n=0; n<nodes.size(); n++ ) {
		nodes[n].pressure = idx_map[n] >= 0 ? result[idx_map[n]] : 0.0;
	}
	
	// Compute gradient
	PARALLEL_FOR for( uint n=0; n<elements.size(); n++ ) {
		elements[n].gradient = vec2d();
		elements[n].mapped = false;
		FLOAT64 rho = elements[n].rho;
		if( ! rho || ! elements[n].frac ) continue;
		uint idx[NUM_VERT] = { g.elements[n][0], g.elements[n][1], g.elements[n][2] };
		FLOAT64 press[NUM_VERT] = { nodes[idx[0]].pressure, nodes[idx[1]].pressure, nodes[idx[2]].pressure };
		if( surf_order == 2 ) embedGhostPressure(n,press);
	
		vec2d value;
		for( uint i=0; i<NUM_VERT; i++ ) {
			for( uint dim=0;dim<DIM;dim++) {
				value[dim] += g.matrix[n].m[i][dim]*press[i];
			}
		}
		elements[n].gradient = value/dt;
		elements[n].mapped = true;
	}
	
	// Downsample to nodes
	for( uint n=0; n<nodes.size(); n++ ) nodes[n].mapped = false;
	resample(0,1);
	
	// Extrapolate element gradient
	extrapolate(g,nodes,1);
	
	// Upsample to elements
	resample(1,1);
	
	// Compute strain
	computeStrain();
}

void femfluid2::extrapolate(const mesher2 &g, std::vector<node2> &nodes, int kind) {
	FLOAT64 max_vel0 = 0.0;
	for( uint n=0; n<g.nodes.size(); n++ ) {
		if( nodes[n].mapped ) max_vel0 = fmax(max_vel0, kind==0 ? nodes[n].velocity.len() : nodes[n].gradient.len());
	}
	
	std::vector<fastmarch2<vec2d>::node2 *> fnodes(g.nodes.size());
	for( uint n=0; n<g.nodes.size(); n++ ) fnodes[n] = new fastmarch2<vec2d>::node2;
	for( uint n=0; n<g.nodes.size(); n++ ) {
		vec2d p = g.nodes[n];
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
	fastmarch2<vec2d>::fastMarch(fnodes,extrapolate_dist,-1.0,0);
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

vec2d femfluid2::getVector( vec2d p, uint kind ) const {
	int hit = g.hitElements(p);
	if( hit >= 0 ) {
		for( uint n1=0; n1<NUM_VERT; n1++ ) {
			for( uint n2=n1+1; n2<NUM_VERT; n2++ ) {
				vec2d vec;
				bool skip = false;
				FLOAT64 x[NUM_VERT] = { p[0], p[1], 1.0 };
				vec2d v[NUM_VERT];
				if( kind == 0 ) {
					v[0] = elements[hit].velocity;
					v[1] = nodes[g.elements[hit][n1]].velocity;
					v[2] = nodes[g.elements[hit][n2]].velocity;
				} else {
					v[0] = elements[hit].gradient;
					v[1] = nodes[g.elements[hit][n1]].gradient;
					v[2] = nodes[g.elements[hit][n2]].gradient;
				}
				
				// Put this subdivision triangle
				std::vector<vec2d> triangle(NUM_VERT);
				triangle[0] = g.centers[hit];
				triangle[1] = g.nodes[g.elements[hit][n1]];
				triangle[2] = g.nodes[g.elements[hit][n2]];
				
				// Compute a shape function here
				FLOAT64 M[NUM_VERT][NUM_VERT];
				FLOAT64 SM[NUM_VERT][NUM_VERT];
				for( uint k=0; k<NUM_VERT; k++ ) {
					for( uint dim=0; dim<DIM; dim++ ) {
						M[dim][k] = triangle[k][dim];
					}
					M[DIM][k] = 1.0;
				}
				invert3x3( M, SM );
				
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
		dump( "getVector failed kind = %d hit = %d (%f,%f)\n", kind, hit, p[0], p[1]);
		return kind==0 ? elements[hit].velocity : elements[hit].gradient;
	}
	dump( "getVector failed kind = %d hit = %d (%f,%f)\n", kind, hit, p[0], p[1]);
	return vec2d();
}

vec2d femfluid2::getPressureGradient( vec2d p ) const {
	return getVector(p,1);
}

vec2d femfluid2::getVelocity( vec2d p ) const {
	return getVector(p,0);
}

FLOAT64 femfluid2::getDivergence( vec2d p ) const {
	int n = g.hitElements(p);
	if( n >= 0 ) {
		FLOAT64 x[NUM_VERT] = { p[0], p[1], 1.0 };
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

void femfluid2::computeStrain() {
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
			FOR_EACH(DIM,DIM) {
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
			FOR_EACH(DIM,DIM) {
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

FLOAT64 femfluid2::getStrain( vec2d p ) const {
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

void femfluid2::render( int name, vec2d mousePos ) const {
    switch(name) {
        case SOLID: {
			std::vector<FLOAT64> solidLS(nodes.size());
			for( uint n=0; n<nodes.size(); n++ ) solidLS[n] = nodes[n].solidLS;
			glColor4d(0.7,0.7,0.3,0.6);
			g.drawLevelset(solidLS);
			break;
		}
		case FLUID: {
			std::vector<FLOAT64> fluidLS(nodes.size());
			for( uint n=0; n<nodes.size(); n++ ) fluidLS[n] = nodes[n].fluidLS;
			glColor4d(0.3,0.3,0.8,0.6);
			g.drawLevelset(fluidLS);
			break;
        }
		case MESH: {
			glColor4f(1.0,1.0,1.0,0.2);
			g.drawMesh();
			break;
		}
		case MATRIX_CONNECTION: {	
			FLOAT64 mind = 999.0;
			uint index = 0;
			for( uint n=0; n<g.nodes.size(); n++ ) {
				FLOAT64 d = (g.nodes[n]-mousePos).len();
				if( d < dx && d < mind ) {
					mind = d;
					index = n;
				}
			}
			if( mind < 1.0 && index<LhsMatrix.n ) {
				glColor4d(1.0,1.0,1.0,0.75);
				glBegin(GL_LINES);
				for( uint n=0; n<LhsMatrix.index[index].size(); n++ ) {
					uint idx = LhsMatrix.index[index][n];
					if( idx != index ) {
						glVertex2f(g.nodes[index][0],g.nodes[index][1]);
						glVertex2f(g.nodes[idx][0],g.nodes[idx][1]);
					}
				}
				glEnd();
				FLOAT64 diag = LhsMatrix(index,index);
				if( fabs(diag) > 0 ) {
					for( uint n=0; n<LhsMatrix.index[index].size(); n++ ) {
						uint idx = LhsMatrix.index[index][n];
						if( idx != index ) {
							vec2d vdir = (g.nodes[idx]-g.nodes[index]).normal();
							FLOAT64 r = 1.0/16.0;
							glRasterPos2d(g.nodes[index][0]-0.25*r+r*vdir[0],g.nodes[index][1]+r*vdir[1]);
							drawBitmapString(format_str("%.2f",-LhsMatrix.value[index][n]/diag));
						}
					}
				}
				glColor4d(1.0,1.0,1.0,1.0);
				glPointSize(3.0);
				glBegin(GL_POINTS);
				glVertex2f(g.nodes[index][0],g.nodes[index][1]);
				glEnd();
				glPointSize(1.0);
			}
			break;
		}
		case VELOCITY: {
			std::vector<vec2d> vel(nodes.size());
			for( uint n=0; n<nodes.size(); n++ ) vel[n] = nodes[n].velocity-dt*nodes[n].gradient;
			glColor4d(1.0,1.0,0.0,0.75);
			g.drawVector(vel,g.NODAL,dx);
			break;
		}
		case PRESSURE: {
			std::vector<bool> mask(nodes.size());
			std::vector<FLOAT64> pressure(nodes.size());
			for( uint n=0; n<nodes.size(); n++ ) {
				pressure[n] = nodes[n].pressure;
				mask[n] = false;
				if( nodes[n].fluidLS<0.0 ) {
					for( uint m=0; m<g.node_elements[n].size(); m++ ) {
						if( elements[g.node_elements[n][m]].frac ) mask[n] = true;
					}
				}
			}
			g.drawScalar(pressure,g.NODAL,mask);
			// Plot zero pressure point
			glPointSize(3.0);
			glBegin(GL_POINTS);
			glColor4d(1.0,1.0,1.0,1.0);
			for( uint n=0; n<elements.size(); n++ ) {
				FLOAT64 rho = elements[n].rho;
				if( rho && elements[n].frac ) {
					for( uint m=0; m<NUM_VERT; m++ ) {
						uint node_idx = g.elements[n][m];
						if( nodes[node_idx].fluidLS > 0.0 ) {
							if( surf_order == 2 ) {
								FLOAT64 coef[NUM_VERT];
								computeGhostPressure(n,m,coef);
								vec2d op;
								FLOAT64 phi = 0.0;
								FLOAT64 sum = 0.0;
								for( uint k=0; k<NUM_VERT; k++ ) sum += coef[k];
								for( uint k=0; k<NUM_VERT; k++ ) coef[k] /= sum;
								for( uint k=0; k<NUM_VERT; k++ ) {
									op += g.nodes[g.elements[n][k]]*coef[k];
									phi += nodes[g.elements[n][k]].fluidLS*coef[k];
								}
								FLOAT64 theta = fraction(nodes[node_idx].fluidLS,phi);
								vec2d pos = theta*g.nodes[node_idx]+(1.0-theta)*op;
								glVertex2dv(pos.v);
							} else {
								glVertex2dv(g.nodes[node_idx].v);
							}
						}
					}
				}
			}
			glEnd();
			glPointSize(1.0);
			break;
		}
		case NAME: {
			glColor4f(1.0,1.0,1.0,1.0);
			glRasterPos2d(0.07,0.92);
			drawBitmapString("FEM",GLUT_BITMAP_HELVETICA_18);
			break;
		}
    }
}