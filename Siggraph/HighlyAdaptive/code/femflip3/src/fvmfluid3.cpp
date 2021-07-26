/*
 *	fvmfluid3.cpp
 *	
 *	Created by Ryoichi Ando on 2/3/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "fvmfluid3.h"
#include "pcgsolver/pcg_solver.h"
#include "levelset3.h"
#include "fastmarch3.h"
#include "opengl.h"
#include "pcgsolver/matutil.h"
#include "util3.h"
#include <stdlib.h>
#include <fstream>
using namespace std;

fvmfluid3::fvmfluid3() {
	dt = 0.1;
	solid=NULL;
	fluid=NULL;
	extrapolate_dist = 1.0;
	surf_order = 2;
	tension = 0.0;
	variation = true;
}

void fvmfluid3::init( uint gn, const levelset3 *hint ) {
	gn ++;	// Translate into nodal grid size
	this->gn = gn;
	dx = 1.0/(gn-1);
	initMesh(hint);
}

void fvmfluid3::initMesh(const levelset3 *hint) {
	// Clear all
	facets.clear();
	cells.clear();
	nodes.clear();
	
	tick(); dump( ">>> Building tetrahedra started...\n" );
	
	// Generate nodal points
	g.setCenterType(g.CIRCUMCENTRIC);
	g.generateMesh(gn,hint);
	
	dump( "<<< Done. Took %s.\n", stock("FVM_tet"));
	tick(); dump( "Computing FVM tet info..." );
	
	// Allocate array
	facets.resize(g.facets.size());
	for( uint n=0; n<g.facets.size(); n++ ) {
		facets[n].velocity = 0.0;
		facets[n].gradient = 0.0;
		facets[n].mapped = false;
	}
	cells.resize(g.elements.size());
	nodes.resize(g.nodes.size());
	
	// Compute facet centers and its normals
	for( uint n=0; n<g.facets.size(); n++ ) {
		vec3d midpos = (g.nodes[g.facets[n][0]]+g.nodes[g.facets[n][1]]+g.nodes[g.facets[n][2]])/3.0;
		vec3d normal = (g.nodes[g.facets[n][1]]-g.nodes[g.facets[n][0]])^(g.nodes[g.facets[n][2]]-g.nodes[g.facets[n][0]]);
		facets[n].center = midpos;
		facets[n].normal = normal.normal();
	}
	for( uint n=0; n<g.elements.size(); n++ ) {
		cells[n].pressure = 0.0;
		cells[n].divergence = 0.0;
	}
	
	// Compute node volume
	for( uint n=0; n<g.nodes.size(); n++ ) {
		nodes[n].volume = 0.0;
		for ( uint m=0; m<g.node_elements[n].size(); m++ ) {
			nodes[n].volume += g.volumes[g.node_elements[n][m]] / g.node_elements[n].size();
		}
	}
	dump( "Done. Took %s.\n", stock("FVM_tet_info"));
}

void fvmfluid3::setParameter( int name, FLOAT64 value ) {
	switch(name) {
		case EXTRAPOLATE_DIST:
			extrapolate_dist = value;
			break;
		case SURF_ORDER:
			surf_order = value;
			break;
		case VARIATION:
			variation = value;
			break;
	}
}

void fvmfluid3::setTimestep( FLOAT64 dt ) {
	this->dt = dt;
}

static FLOAT64 fraction( FLOAT64 phi0, FLOAT64 phi1 ) {
	return phi0*phi1 > 0 ? 1 - (phi0 > 0) : -fmin(phi0,phi1)/fmax(fabs(phi1-phi0),1.0e-6);
}

void fvmfluid3::setupSolidLevelset( const levelset3 *solid ) {
	this->solid = solid;
	
	// Compute nodal solid levelset
	for( uint n=0; n<g.nodes.size(); n++ ) {
		nodes[n].solidLS = solid->evalLevelset(g.nodes[n]);
	}
	
	// Compute center solid levelset
	for( uint n=0; n<g.centers.size(); n++ ) {
		cells[n].solidLS = solid->evalLevelset(g.centers[n]);
	}
	
	// Compute fractions
	for( uint n=0; n<g.facets.size(); n++ ) {
		std::vector<vec3d> pos(3);
		std::vector<FLOAT64> levelsets(3);
		for( uint m=0; m<3; m++ ) {
			pos[m] = g.nodes[g.facets[n][m]];
			levelsets[m] = nodes[g.facets[n][m]].solidLS;
		}
		std::vector<vec3d> marchPos = util3::march2DPoints(pos,levelsets);
		if( ! marchPos.empty() ) {
			vec3d vec0 = (g.nodes[g.facets[n][1]]-g.nodes[g.facets[n][0]]).normal();
			vec3d vec1 = (g.nodes[g.facets[n][1]]-g.nodes[g.facets[n][0]]).normal();
			vec3d n_vec = (vec1 ^ vec0).normal();
			vec3d x_vec = vec0;
			vec3d y_vec = (n_vec ^ x_vec).normal();
			for( uint k=0; k<marchPos.size(); k++ ) marchPos[k] = vec3d(marchPos[k]*x_vec,marchPos[k]*y_vec,0.0);
			FLOAT64 area = util3::compute2DArea(marchPos);
			facets[n].fraction = fmax(0.0,fmin(0.0,1.0-area/g.facetArea[n]));
			if( ! variation ) facets[n].fraction = facets[n].fraction < 1.0 ? 0.0 : 1.0;
		} else {
			facets[n].fraction = 1.0;
		}
	}
	
	// Compute solid volume
	tick(); dump( "Computing solid volume..." );
	for( uint n=0; n<g.elements.size(); n++ ) {
		FLOAT64 LS[4];
		for( uint i=0; i<4; i++ ) {
			LS[i] = nodes[g.elements[n][i]].solidLS;
		}
		cells[n].solidVolume = g.getLevelsetVolume(n,LS);
	}
	dump( "Done. Took %s.\n", stock());
}

void fvmfluid3::setupFluidLevelset( const levelset3 *fluid ) {
	this->fluid = fluid;
	
	// Compute nodal fluid levelset
	for( uint n=0; n<g.nodes.size(); n++ ) {
		nodes[n].fluidLS = fluid->evalLevelset(g.nodes[n]);
	}
	// Compute center fluid levelset
	for( uint n=0; n<g.centers.size(); n++ ) {
		cells[n].fluidLS = fluid->evalLevelset(g.centers[n]);
	}
	this->fluid = fluid;
}

void fvmfluid3::extrapolate(const mesher3 &g, std::vector<node3> &nodes, int kind) {
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
}

void fvmfluid3::getSamplePoints( std::vector<vec3d> &pos ) {
	pos.resize(facets.size());
	for( uint n=0; n<facets.size(); n++ ) {
		pos[n] = facets[n].center;
	}
}

void fvmfluid3::setupVelocity( const std::vector<vec3d> &vel, const std::vector<bool> &mapped ) {
	for( uint n=0; n<facets.size(); n++ ) {
		facets[n].velocity = vel[n] * facets[n].normal;
		facets[n].mapped = facets[n].fraction ? mapped[n] : false;
	}
	computeElementVelocity(g,facets,cells,0);
}

static FLOAT64 computeTriArea( vec3d p0, vec3d p1, vec3d p2 ) {
	return 0.5*((p1-p0)^(p2-p0)).len();
}

void fvmfluid3::buildMatrix() {
	tick(); dump( "Building left and right hand side matrix..." );
	// Build divergence matrix operator
	RhsMatrix.clear();
	RhsMatrix.resize(g.elements.size());
	for( uint n=0; n<g.element_facets.size(); n++ ) {
		for( uint m=0; m<g.element_facets[n].size(); m++ ) {
			uint fidx = g.element_facets[n][m];
			uint eidx = g.element_elements[n][m];
			vec3d fn = facets[fidx].normal;
			// Flip normal if necessary
			FLOAT64 sgn = 1.0;
			if( fn * (g.centers[eidx]-g.centers[n]) < 0 ) sgn = -1.0;
			FLOAT64 value = -sgn*g.facetArea[fidx]*facets[fidx].fraction;
			if( is_nan(value)) {
				dump( "Tried to add NaN value (%g) at right hand side matrix !\n", value);
				exit(0);
			} else {
				RhsMatrix.add_to_element(n,fidx,-sgn*g.facetArea[fidx]*facets[fidx].fraction);
			}
		}
	}
	
	// Build pressure matrix operator
	LhsMatrix.clear();
	LhsMatrix.resize(g.elements.size());
	
	for( uint n=0; n<g.elements.size(); n++ ) {
		for( uint m=0; m<g.element_facets[n].size(); m++ ) {
			uint fid = g.element_facets[n][m];
			FLOAT64 area = g.facetArea[fid] * facets[fid].fraction;
			uint n2 = g.element_elements[n][m];
			if( n==n2 ) {
				email::print( "Same ID. Should not happen !\n" );
				email::send();
				exit(0);
			}
			//FLOAT64 dist = (g.centers[n2]-facetCenters[fid]).len()+(g.centers[n]-facetCenters[fid]).len();
			FLOAT64 dist = fmax((g.centers[n2]-g.centers[n]).len(),1e-2*dx);
			FLOAT64 rho = cells[n].fluidLS < 0 || cells[n2].fluidLS < 0;
			if( surf_order == 2 ) rho = fraction(cells[n].fluidLS,cells[n2].fluidLS);
			if( rho ) {
				if( area && fabs(dist) < 1e-8 ) {
					printf( "Length zero ID: %d (%.2f,%.2f,%.2f) - ID: %d (%.2f,%.2f,%.2f)!\n",
						   n2, g.centers[n2][0], g.centers[n2][1], g.centers[n2][2],
						   n, g.centers[n][0], g.centers[n2][1], g.centers[n2][2] );
					printf( "Face Fraction %f Face Polygon: (%f,%f,%f,%f) - (%f,%f,%f,%f) - (%f,%f,%f,%f)\n",
						   facets[fid].fraction,
						   g.nodes[g.facets[fid][0]][0], g.nodes[g.facets[fid][0]][1], g.nodes[g.facets[fid][0]][2], nodes[g.facets[fid][0]].solidLS,
						   g.nodes[g.facets[fid][1]][0], g.nodes[g.facets[fid][1]][1], g.nodes[g.facets[fid][1]][2], nodes[g.facets[fid][1]].solidLS,
						   g.nodes[g.facets[fid][2]][0], g.nodes[g.facets[fid][2]][1], g.nodes[g.facets[fid][2]][2], nodes[g.facets[fid][2]].solidLS );
					printf( "Polygon: (%d,%d,%d,%d), (%d,%d,%d,%d).\n",
						   g.elements[n2][0], g.elements[n2][1], g.elements[n2][2], g.elements[n2][3],
						   g.elements[n][0], g.elements[n][1], g.elements[n][2], g.elements[n][3] );
					exit(0);
				} else {
					FLOAT64 scale = area ? area/dist : 0.0;
					if( scale ) {
						LhsMatrix.add_to_element(n,n2,-scale/rho);
						LhsMatrix.add_to_element(n,n,scale/rho);
						if( scale/rho < 0.0 ) {
							dump( "Tried to add negative value (%g/%g = %g) at diagonal matrix ! dist=%g, area=%g\n", area, dist, scale/rho, dist, area);
							dump( "FacetArea = %g FacetFraction = %g", g.facetArea[fid], facets[fid].fraction );
							exit(0);
						}
					}
				}
			}
		}
	}
	dump( "Done. Took %s.\n", stock("FVM_matrix"));
	
	tick(); dump( "Embedding dirichlet boundary condition..." );
	// Embed boundary conditions
	for( uint n=0; n<cells.size(); n++ ) {
		if( cells[n].fluidLS > 0.0 ) {
			LhsMatrix.index[n].clear();
			LhsMatrix.value[n].clear();
		} else {
			for( uint m=0; m<LhsMatrix.index[n].size(); m++ ) {
				uint idx = LhsMatrix.index[n][m];
				if( cells[idx].fluidLS > 0.0 ) {
					// Plug pressure surface boundary condition
					LhsMatrix.set_element(n,idx,0.0);
				}
			}
		}
	}
	dump( "Done. Took %s.\n", stock("FVM_embed"));
}

void fvmfluid3::resample( uint type, uint kind ) {
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
					if( cells[idx].mapped ) {
						FLOAT64 w = g.volumes[idx];
						if( kind == 0 ) velocity += w*cells[idx].velocity;
						else velocity += w*cells[idx].gradient;
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
			if( ! cells[n].mapped ) {
				vec3d velocity;
				FLOAT64 wsum=0;
				if( kind == 0 ) cells[n].velocity = vec3d();
				else cells[n].gradient = vec3d();
				cells[n].mapped = false;
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
					if( kind == 0 ) cells[n].velocity = velocity / wsum;
					else cells[n].gradient = velocity / wsum;
					cells[n].mapped = true;
				}
			}
		}
	}
}

void fvmfluid3::project() {
	tick(); dump( "Extrapolting velocity..." );
	for( uint n=0; n<nodes.size(); n++ ) nodes[n].mapped = false;
	resample(0,0);
	
	// Extrapolate element velocity
	extrapolate(g,nodes,0);
	
	// Upsample to elements
	resample(1,0);
	dump( "Done. Took %s.\n", stock("FVM_velocity_extrapolate"));
	
	// Compute facet velocity where still unmapped
	for( uint n=0; n<g.facets.size(); n++ ) {
		if( ! facets[n].mapped ) {
			vec3d velocity;
			FLOAT64 sum = 0.0;
			for( uint m=0; m<g.facet_elements[n].size(); m++ ) {
				uint idx = g.facet_elements[n][m];
				if( cells[idx].mapped ) {
					sum ++;
					velocity += cells[idx].velocity;
				}
			}
			if( sum ) {
				facets[n].mapped = true;
				facets[n].velocity = velocity*facets[n].normal/sum;
			}
		}
	}
	
	// Build matrix
	buildMatrix();
	
	tick(); dump( "Computing divergence..." );
	// Compute divergence
	std::vector<FLOAT64> divergence(g.elements.size());
	std::vector<FLOAT64> facetVelocity(g.facets.size());
	for( uint n=0; n<g.facets.size(); n++ ) facetVelocity[n] = facets[n].velocity;
	multiply<FLOAT64>(RhsMatrix,facetVelocity,divergence);
	dump( "Done. Took %s.\n", stock("divergence"));
	
	// Copy divergence
	for( uint n=0; n<cells.size(); n++ ) {
		cells[n].divergence = divergence[n];
	}

	// Check Nan value (for debug)
#if 1
	for( uint n=0; n<divergence.size(); n++ ) {
		if( is_nan(divergence[n]) ) {
			email::print( "Right hand side contains NaN value !\n" );
			email::send();
			exit(0);
		}
	}
#endif
	
	// Compact the matrix
	tick(); dump( "Compressing matrix..." );
	std::vector<int> idx_map;
	std::vector<FLOAT64> rhs = divergence;
	SparseMatrix<FLOAT64> rawLhsMatrix = LhsMatrix;
	compress( LhsMatrix, rhs, idx_map );
	dump( "Done. Took %s.\n", stock());
	
	// Solve by incomplete cholesky factorizaion preconditioned conjugate gradient method
    PCGSolver<FLOAT64> solver;
	solver.set_solver_parameters(1e-10, 50000, 0.97, 0.25);
    FLOAT64 residual_out;
	std::vector<FLOAT64> result(LhsMatrix.n);
	for( uint n=0; n<cells.size(); n++ ) {
		if( idx_map[n] >= 0 ) result[idx_map[n]] = cells[n].pressure;
	}
	
	int iterations;
	FLOAT64 msec;
	dump( "Solving pressure by FVM..." );
	bool converged = solver.solve( LhsMatrix, rhs, result, residual_out, iterations, msec );
	if( ! converged ) {
		// Dump this matrix
		ofstream file;
		file.open(format_str("%s/Lhs_FVM_matrix.m",root_path));
		LhsMatrix.write_matlab(file,LhsMatrix.n,LhsMatrix.n,"Lhs_FVM");
		file.close();
		file.open(format_str("%s/Rhs_FVM_vector.m",root_path));
		write_matlab1d(file,rhs,"Rhs_FVM");
		file.close();
		
		email::print( "PCG did not converge. Residual = %e\n", residual_out );
		email::send();
		exit(0);
	} else {
        dump("Done. Took %d iterations with %d unknowns. Residual = %e. Took %s.\n", iterations, LhsMatrix.n, residual_out, tstr(msec));
		writeNumber("FVM_PCG_iterations",iterations);
		writeNumber("FVM_PCG_time",msec);
		writeNumber("FVM_PCG_unknowns", LhsMatrix.n);
		writeNumber("FVM_PCG_residual",residual_out);
    }
	
	// Set pressure
	for( uint n=0; n<g.elements.size(); n++ ) {
		cells[n].pressure = idx_map[n] >= 0 ? result[idx_map[n]] : 0.0;
	}
	
	// Compute pressure gradient
	tick(); dump( "Computing facet gradient..." );
	for( uint n=0; n<g.facet_elements.size(); n++ ) {
		facets[n].mapped = false;
		facets[n].gradient = 0.0;
		vec3d fn = facets[n].normal;
		if( g.facet_elements[n].size() == 2 ) {
			uint idx0 = g.facet_elements[n][0];
			uint idx1 = g.facet_elements[n][1];
			
			// Flip normal if necessary
			FLOAT64 sgn = 1.0;
			if( fn * (g.centers[idx1]-g.centers[idx0]) < 0 ) sgn = -1.0;
			//FLOAT64 dist = (g.centers[idx0]-facetCenters[n]).len()+(g.centers[idx1]-facetCenters[n]).len();
			FLOAT64 dist = fmax((g.centers[idx0]-g.centers[idx1]).len(),1e-2*dx);
			FLOAT64 scale = sgn/dist;
			FLOAT64 rho = cells[idx0].fluidLS < 0 || cells[idx1].fluidLS < 0;
			if( surf_order == 2 ) rho = fraction(cells[idx0].fluidLS,cells[idx1].fluidLS);
			if( rho ) {
				// Compute pressure
				FLOAT64 press0 = cells[idx0].pressure;
				FLOAT64 press1 = cells[idx1].pressure;
				
				// Set pressure gradient
				facets[n].gradient = scale*(press1-press0)/(dt*rho);
				facets[n].mapped = true;
			}
		}
	}
	dump( "Done. Took %s.\n", stock("FVM_facet_gradient"));
	
	// Compute element center pressure gradient
	tick(); dump( "Computing element center gradient by the least square..." );
	computeElementVelocity(g,facets,cells,1);
	dump( "Done. Took %s.\n", stock("FVM_least_square_velocity"));
	
	tick(); dump( "Extrapolating element gradient..." );
	// Downsample to nodes
	for( uint n=0; n<nodes.size(); n++ ) nodes[n].mapped = false;
	resample(0,1);
	
	// Extrapolate element gradient
	extrapolate(g,nodes,1);
	
	// Upsample to elements
	resample(1,1);
	dump( "Done. Took %s.\n", stock("FVM_element_gradient_extrapolate"));
}

void fvmfluid3::computeElementVelocity( const mesher3 &g, const std::vector<facet3> &facets, std::vector<cell3> &cells, int kind ) {
	PARALLEL_FOR for( uint n=0; n<g.element_facets.size(); n++ ) {
		if( kind == 0 )	cells[n].velocity = vec3d();
		else cells[n].gradient = vec3d();
		cells[n].mapped = false;
		if( g.volumes[n] != cells[n].solidVolume ) {
			vector<vector<FLOAT64> > N;
			vector<FLOAT64> c;
			for( uint m=0; m<g.element_facets[n].size(); m++ ) {
				vector<FLOAT64> nv(3);
				uint fidx = g.element_facets[n][m];
				if( facets[fidx].mapped ) {
					nv[0] = facets[fidx].normal[0];
					nv[1] = facets[fidx].normal[1];
					nv[2] = facets[fidx].normal[2];
					N.push_back(nv);
					if( kind == 0 ) c.push_back(facets[fidx].velocity);
					else c.push_back(facets[fidx].gradient);
				}
			}
			if( N.size() >= DIM ) {
				FLOAT64 NtN[3][3];
				FLOAT64 Ntc[3];
				for( uint i=0; i<3; i++ ) for( uint j=0; j<3; j++ ) {
					NtN[i][j] = 0.0;
					for( uint k=0; k<N.size(); k++ ) NtN[i][j] += N[k][i]*N[k][j];
				}
				for( uint i=0; i<3; i++ ) {
					Ntc[i] = 0.0;
					for( uint k=0; k<c.size(); k++ ) Ntc[i] += N[k][i]*c[k];
				}
				FLOAT64 NtNinv[3][3];
				if( invert3x3(NtN,NtNinv) ) {
					FLOAT64 NtNinv_c[3];
					for( uint i=0; i<3; i++ ) {
						NtNinv_c[i] = 0.0;
						for( uint k=0; k<3; k++ ) NtNinv_c[i] += NtNinv[i][k]*Ntc[k];
					}
					if( kind == 0 ) cells[n].velocity = vec3d(NtNinv_c[0],NtNinv_c[1],NtNinv_c[2]);
					else cells[n].gradient = vec3d(NtNinv_c[0],NtNinv_c[1],NtNinv_c[2]);
					cells[n].mapped = true;
				}
			}
		}
	}
}

vec3d fvmfluid3::getVector( vec3d p, uint kind ) const {
	int hit = g.hitElements(p);
	if( hit >= 0 ) {
		uint idx = 0;
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
						// dump( "Failed inverse 4x4! (subdividion part)\n" );
					} else {					
						vec3d vec;
						bool skip = false;
						FLOAT64 x[NUM_VERT] = { p[0], p[1], p[2], 1.0 };
						vec3d v[NUM_VERT];
						if( kind == 0 ) {
							v[0] = cells[hit].velocity;
							v[1] = nodes[g.elements[hit][n1]].velocity;
							v[2] = nodes[g.elements[hit][n2]].velocity;
							v[3] = nodes[g.elements[hit][n3]].velocity;
						} else {
							v[0] = cells[hit].gradient;
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
						idx ++;
						if( ! skip ) return vec;
					}
				}
			}
		}
		dump( "getVector failed kind = %d hit = %d (%f,%f,%f)\n", kind, hit, p[0], p[1], p[2]);
		return kind==0 ? cells[hit].velocity : cells[hit].gradient;
	}
	dump( "getVector failed kind = %d hit = %d (%f,%f,%f)\n", kind, hit, p[0], p[1], p[2]);
	return vec3d();
}

vec3d fvmfluid3::getPressureGradient( vec3d p ) const {
	return getVector(p,1);
}

vec3d fvmfluid3::getVelocity( vec3d p ) const {
	return getVector(p,0);
}

FLOAT64 fvmfluid3::getDivergence( vec3d p ) const {
	int n = g.hitElements(p);
	if( n >= 0 ) {
		return cells[n].divergence/g.volumes[n];
	}
	return 0.0;
}

static void drawBitmapString( const char *string, void *font=NULL ) {
	if( ! font ) font = GLUT_BITMAP_HELVETICA_10;
	while (*string) glutBitmapCharacter(font, *string++);
}

void fvmfluid3::render( int name ) const {
	switch(name) {
		case MESH:
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
