/*
 *	fvmfluid2.cpp
 *	
 *	Created by Ryoichi Ando on 12/27/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "pcgsolver/pcg_solver.h"
#include "fvmfluid2.h"
#include "levelset2.h"
#include "fastmarch2.h"
#include "opengl.h"
#include "pcgsolver/matutil.h"
#include "util2.h"
#include <fstream>
using namespace std;

fvmfluid2::fvmfluid2() {
	dt = 0.1;
	solid=NULL;
	fluid = NULL;
	extrapolate_dist = 1.0;
	surf_order = 2;
	variation = true;
}

void fvmfluid2::init( uint gn, const levelset2 *hint ) {
	gn ++;	// Translate into nodal grid size
	this->gn = gn;
	dx = 1.0/(gn-1);
	initMesh(hint);
}

void fvmfluid2::initMesh(const levelset2 *hint) {
	// Clear all
	facets.clear();
	cells.clear();
	nodes.clear();
	
	// Generate nodal points
	g.setCenterType(g.CIRCUMCENTRIC);
	g.generateMesh(gn,hint);
	
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
		vec2d midpos = (g.nodes[g.facets[n][0]]+g.nodes[g.facets[n][1]])/2.0;
		vec2d normal = (g.nodes[g.facets[n][1]]-g.nodes[g.facets[n][0]]).rotate().normal();
		facets[n].center = midpos;
		facets[n].normal = normal;
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
}

void fvmfluid2::setParameter( int name, FLOAT64 value ) {
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

void fvmfluid2::setTimestep( FLOAT64 dt ) {
	this->dt = dt;
}

static FLOAT64 fraction( FLOAT64 phi0, FLOAT64 phi1 ) {
	return phi0*phi1 > 0 ? 1 - (phi0 > 0) : -fmin(phi0,phi1)/fmax(fabs(phi1-phi0),1.0e-18);
}

void fvmfluid2::setupSolidLevelset( const levelset2 *solid ) {
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
		uint idx0 = g.facets[n][0];
		uint idx1 = g.facets[n][1];
		facets[n].fraction = 1.0-fraction(nodes[idx0].solidLS,nodes[idx1].solidLS);
		if( ! variation ) facets[n].fraction = facets[n].fraction < 1.0 ? 0.0 : 1.0;
	}
	
	// Compute solid volume
	for( uint n=0; n<g.elements.size(); n++ ) {
		vector<vec2d> points(g.elements[n].size());
		vector<FLOAT64> LS(g.elements[n].size());
		for( uint i=0; i<g.elements[n].size(); i++ ) {
			points[i] = g.nodes[g.elements[n][i]];
			LS[i] = nodes[g.elements[n][i]].solidLS;
		}
		cells[n].solidVolume = util2::computeVolume(util2::marchPoints(points,LS));
	}
}

void fvmfluid2::setupFluidLevelset( const levelset2 *fluid ) {
	this->fluid = fluid;
	
	// Compute nodal fluid levelset
	for( uint n=0; n<g.nodes.size(); n++ ) {
		nodes[n].fluidLS = fluid->evalLevelset(g.nodes[n]);
	}
	// Compute center fluid levelset
	for( uint n=0; n<g.centers.size(); n++ ) {
		cells[n].fluidLS = fluid->evalLevelset(g.centers[n]);
	}
}

void fvmfluid2::extrapolate(const mesher2 &g, std::vector<node2> &nodes, int kind) {
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
}

void fvmfluid2::getSamplePoints( std::vector<vec2d> &pos ) {
	pos.resize(facets.size());
	for( uint n=0; n<facets.size(); n++ ) {
		pos[n] = facets[n].center;
	}
}

void fvmfluid2::setupVelocity( const std::vector<vec2d> &vel, const std::vector<bool> &mapped ) {
	for( uint n=0; n<facets.size(); n++ ) {
		facets[n].velocity = vel[n] * facets[n].normal;
		facets[n].mapped = facets[n].fraction ? mapped[n] : false;
	}
	computeElementVelocity(g,facets,cells,0);
}

void fvmfluid2::buildMatrix() {
	// Build divergence matrix operator
	RhsMatrix.clear();
	RhsMatrix.resize(g.elements.size());
	for( uint n=0; n<g.element_facets.size(); n++ ) {
		for( uint m=0; m<g.element_facets[n].size(); m++ ) {
			uint fidx = g.element_facets[n][m];
			uint eidx = g.element_elements[n][m];
			vec2d fn = facets[fidx].normal;
			// Flip normal if necessary
			FLOAT64 sgn = 1.0;
			if( fn * (g.centers[eidx]-g.centers[n]) < 0 ) sgn = -1.0;
			RhsMatrix.add_to_element(n,fidx,-sgn*g.facetArea[fidx]*facets[fidx].fraction);
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
				email::print("Same ID. Should not happen !\n" );
				email::send();
				exit(0);
			}
			//FLOAT64 dist = (g.centers[n2]-facetCenters[fid]).len()+(g.centers[n]-facetCenters[fid]).len();
			FLOAT64 dist = fmax((g.centers[n2]-g.centers[n]).len(),dx*1e-2);
			FLOAT64 rho = cells[n].fluidLS < 0 || cells[n2].fluidLS < 0;
			if( surf_order == 2 ) rho = fraction(cells[n].fluidLS,cells[n2].fluidLS);
			if( rho && area && fabs(dist) < 1e-8 ) {
				printf( "Length zero ID: %d (%.2f,%.2f) - ID: %d (%.2f,%.2f)!\n",
					   n2, g.centers[n2][0], g.centers[n2][1], 
					   n, g.centers[n][0], g.centers[n2][1] );
				printf( "Face Fraction %f Face Polygon: (%f,%f,%f) - (%f,%f,%f) - (%f,%f,%f)\n",
					   facets[fid].fraction,
					   g.nodes[g.facets[fid][0]][0], g.nodes[g.facets[fid][0]][1], nodes[g.facets[fid][0]].solidLS,
					   g.nodes[g.facets[fid][1]][0], g.nodes[g.facets[fid][1]][1], nodes[g.facets[fid][1]].solidLS,
					   g.nodes[g.facets[fid][2]][0], g.nodes[g.facets[fid][2]][1], nodes[g.facets[fid][2]].solidLS );
				printf( "Polygon: (%d,%d,%d), (%d,%d,%d).\n",
					   g.elements[n2][0], g.elements[n2][1], g.elements[n2][2],
					   g.elements[n][0], g.elements[n][1], g.elements[n][2]);
				exit(0);
			} else {
				FLOAT64 scale = area ? area/dist : 0.0;
				if( scale ) {
					LhsMatrix.add_to_element(n,n2,-scale/rho);
					LhsMatrix.add_to_element(n,n,scale/rho);
				}
			}
		}
	}
	
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
}

void fvmfluid2::resample( uint type, uint kind ) {
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
				vec2d velocity;
				FLOAT64 wsum=0;
				if( kind == 0 ) cells[n].velocity = vec2d();
				else cells[n].gradient = vec2d();
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

void fvmfluid2::project() {
	// Downsample to nodes
	for( uint n=0; n<nodes.size(); n++ ) nodes[n].mapped = false;
	resample(0,0);
	
	// Extrapolate element velocity
	extrapolate(g,nodes,0);
	
	// Upsample to elements
	resample(1,0);
	
	// Compute facet velocity where still unmapped
	for( uint n=0; n<g.facets.size(); n++ ) {
		if( ! facets[n].mapped ) {
			vec2d velocity;
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

	// Build Matrix
	buildMatrix();

	// Compute divergence
	std::vector<FLOAT64> divergence(g.elements.size());
	std::vector<FLOAT64> facetVelocity(g.facets.size());
	for( uint n=0; n<g.facets.size(); n++ ) facetVelocity[n] = facets[n].velocity;
	multiply<FLOAT64>(RhsMatrix,facetVelocity,divergence);
	
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
	SparseMatrix<FLOAT64> LhsMatrix = this->LhsMatrix;
	std::vector<int> idx_map;
	std::vector<FLOAT64> rhs = divergence;
	SparseMatrix<FLOAT64> rawLhsMatrix = LhsMatrix;
	compress( LhsMatrix, rhs, idx_map );
	
	// Solve by incomplete cholesky factorizaion preconditioned conjugate gradient method
    PCGSolver<FLOAT64> solver;	
    FLOAT64 residual_out;
	std::vector<FLOAT64> result(LhsMatrix.n);
	for( uint n=0; n<cells.size(); n++ ) {
		if( idx_map[n] >= 0 ) result[idx_map[n]] = cells[n].pressure;
	}
    int iterations;
	FLOAT64 msec;
	bool converged = solver.solve( LhsMatrix, rhs, result, residual_out, iterations, msec );
	if( residual_out > 1.0 ) converged = false;
	if( ! converged ) {
		email::print( "PCG did not converge.\n" );
		email::send();
#if 0
		ofstream file;
		file.open("Lhs_FVM_matrix.m");
		LhsMatrix.write_matlab(file,LhsMatrix.n,LhsMatrix.n,"Lhs_FVM");
		file.close();
		file.open("Rhs_FVM_vector.m");
		write_matlab1d(file,rhs,"Rhs_FVM");
		file.close();
#endif
		exit(0);
	} else {
        dump( "FVM: PCG Converged ! %d iterations and took %.3f msec with %d unknowns.\n", iterations, msec, LhsMatrix.n);
    }
	
	// Set pressure
	for( uint n=0; n<g.elements.size(); n++ ) {
		cells[n].pressure = idx_map[n] >= 0 ? result[idx_map[n]] : 0.0;
	}
	
	// Compute pressure gradient	
	for( uint n=0; n<g.facet_elements.size(); n++ ) {
		facets[n].mapped = false;
		facets[n].gradient = 0.0;
		vec2d fn = facets[n].normal;
		if( facets[n].fraction && g.facet_elements[n].size() == 2 ) {
			uint idx0 = g.facet_elements[n][0];
			uint idx1 = g.facet_elements[n][1];
			
			// Flip normal if necessary
			FLOAT64 sgn = 1.0;
			if( fn * (g.centers[idx1]-g.centers[idx0]) < 0 ) sgn = -1.0;
			//FLOAT64 dist = (g.centers[idx0]-facetCenters[n]).len()+(g.centers[idx1]-facetCenters[n]).len();
			FLOAT64 dist = fmax((g.centers[idx0]-g.centers[idx1]).len(),dx*1e-2);
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
	
	// Compute element center pressure gradient
	computeElementVelocity(g,facets,cells,1);
	
	// Downsample to nodes
	for( uint n=0; n<nodes.size(); n++ ) nodes[n].mapped = false;
	resample(0,1);
	
	// Extrapolate element gradient
	extrapolate(g,nodes,1);
	
	// Upsample to elements
	resample(1,1);
}

void fvmfluid2::computeElementVelocity( const mesher2 &g, const std::vector<facet2> &facets, std::vector<cell2> &cells, int kind ) {
	PARALLEL_FOR for( uint n=0; n<g.element_facets.size(); n++ ) {
		if( kind == 0 )	cells[n].velocity = vec2d();
		else cells[n].gradient = vec2d();
		cells[n].mapped = false;
		if( g.volumes[n] != cells[n].solidVolume ) {
			std::vector<vector<FLOAT64> > N;
			std::vector<FLOAT64> c;
			for( uint m=0; m<g.element_facets[n].size(); m++ ) {
				std::vector<FLOAT64> nv(2);
				uint fidx = g.element_facets[n][m];
				if( facets[fidx].mapped ) {
					nv[0] = facets[fidx].normal[0];
					nv[1] = facets[fidx].normal[1];
					N.push_back(nv);
					if( kind == 0 )	c.push_back(facets[fidx].velocity);
					else c.push_back(facets[fidx].gradient);
				}
			}
			if( N.size() > DIM ) {
				FLOAT64 NtN[2][2];
				FLOAT64 Ntc[2];
				for( uint i=0; i<2; i++ ) for( uint j=0; j<2; j++ ) {
					NtN[i][j] = 0.0;
					for( uint k=0; k<N.size(); k++ ) NtN[i][j] += N[k][i]*N[k][j];
				}
				for( uint i=0; i<2; i++ ) {
					Ntc[i] = 0.0;
					for( uint k=0; k<c.size(); k++ ) Ntc[i] += N[k][i]*c[k];
				}
				FLOAT64 NtNinv[2][2];
				if( invert2x2(NtN,NtNinv) ) {
					FLOAT64 NtNinv_c[2];
					for( uint i=0; i<2; i++ ) {
						NtNinv_c[i] = 0.0;
						for( uint k=0; k<2; k++ ) NtNinv_c[i] += NtNinv[i][k]*Ntc[k];
					}
					if( kind == 0 ) cells[n].velocity = vec2d(NtNinv_c[0],NtNinv_c[1]);
					else cells[n].gradient = vec2d(NtNinv_c[0],NtNinv_c[1]);
					cells[n].mapped = true;
				}
			}
		}
	}
}

vec2d fvmfluid2::getVector( vec2d p, uint kind ) const {
	int hit = g.hitElements(p);
	if( hit >= 0 ) {
		for( uint n1=0; n1<NUM_VERT; n1++ ) {
			for( uint n2=n1+1; n2<NUM_VERT; n2++ ) {
				
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
				vec2d vec;
				bool skip = false;
				FLOAT64 x[NUM_VERT] = { p[0], p[1], 1.0 };
				vec2d v[NUM_VERT];
				if( kind == 0 ) {
					v[0] = cells[hit].velocity;
					v[1] = nodes[g.elements[hit][n1]].velocity;
					v[2] = nodes[g.elements[hit][n2]].velocity;
				} else {
					v[0] = cells[hit].gradient;
					v[1] = nodes[g.elements[hit][n1]].gradient;
					v[2] = nodes[g.elements[hit][n2]].gradient;
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
		dump( "getVector failed kind = %d hit = %d (%f,%f)\n", kind, hit, p[0], p[1]);
		return kind==0 ? cells[hit].velocity : cells[hit].gradient;
	}
	dump( "getVector failed kind = %d hit = %d (%f,%f)\n", kind, hit, p[0], p[1]);
	return vec2d();
}

vec2d fvmfluid2::getPressureGradient( vec2d p ) const {
	return getVector(p,1);
}

vec2d fvmfluid2::getVelocity( vec2d p ) const {
	return getVector(p,0);
}

FLOAT64 fvmfluid2::getDivergence( vec2d p ) const {
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

void fvmfluid2::render( int name, vec2d mousePos ) const {
	switch(name) {
		case MESH: {
			glColor4f(1.0,1.0,1.0,0.2);
			g.drawMesh();
			glPointSize(2.0);
			glColor4f(1.0,1.0,1.0,0.4);
			g.drawCenters();
			glPointSize(1.0);	
			break;
		}
		case SOLID: {
			glColor4d(0.7,0.7,0.3,0.6);
			vector<FLOAT64> nodalSolidLS(g.nodes.size());
			for( uint n=0; n<g.nodes.size(); n++ ) nodalSolidLS[n] = nodes[n].solidLS;
			g.drawLevelset(nodalSolidLS);
			break;
		}
		case FLUID: {
			glColor4d(0.3,0.3,0.8,0.6);
			vector<FLOAT64> nodalFluidLS(g.nodes.size());
			for( uint n=0; n<g.nodes.size(); n++ ) nodalFluidLS[n] = nodes[n].fluidLS;
			g.drawLevelset(nodalFluidLS);
			break;
		}
		case PRESSURE: {
			vector<bool> mask(g.elements.size());
			vector<FLOAT64> pressure(g.elements.size());
			for( uint n=0; n<g.elements.size(); n++ ) {
				mask[n] = cells[n].fluidLS<0.0 && g.volumes[n]!=cells[n].solidVolume;
				pressure[n] = cells[n].pressure;
			}
			g.drawScalar(pressure,g.ELEMENT,mask);
			
			// Draw estimated zero pressure positions
			glColor4d(1.0,1.0,1.0,1.0);
			glPointSize(3.0);
			glBegin(GL_POINTS);
			for( uint n=0; n<cells.size(); n++ ) {
				uint idx0 = n;
				for( uint m=0; m<g.element_elements[n].size(); m++ ) {
					uint idx1 = g.element_elements[n][m];
					if( cells[idx0].fluidLS * cells[idx1].fluidLS < 0.0 ) {
						FLOAT64 t = fraction(cells[idx0].fluidLS,cells[idx1].fluidLS);
						vec2d pos;
						if( surf_order == 2 ) {
							vec2d p0, p1;
							if( cells[idx0].fluidLS < 0.0 ) {
								p0 = g.centers[idx0];
								p1 = g.centers[idx1];
							} else {
								p0 = g.centers[idx1];
								p1 = g.centers[idx0];
							}
							pos = (1.0-t)*p0+t*p1;
						} else {
							if( cells[idx0].fluidLS < 0.0 ) pos = g.centers[idx1];
							else pos = g.centers[idx0];
						}						
						if( solid->evalLevelset(pos) > 0.0 ) {
							glVertex2dv(pos.v);
						}
					}
				}
			}
			glEnd();
			glPointSize(1.0);
			break;
		}
		case VELOCITY: {
			std::vector<vec2d> vel(cells.size());
			for( uint n=0; n<cells.size(); n++ ) vel[n] = cells[n].velocity-dt*cells[n].gradient;
			glColor4d(1.0,1.0,0.0,0.75);
			g.drawVector(vel,g.ELEMENT,dx);
			break;
		}
		case MATRIX_CONNECTION: {
			FLOAT64 mind = 999.0;
			uint index = 0;
			for( uint n=0; n<g.elements.size(); n++ ) {
				FLOAT64 d = (g.centers[n]-mousePos).len();
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
						glVertex2f(g.centers[index][0],g.centers[index][1]);
						glVertex2f(g.centers[idx][0],g.centers[idx][1]);
					}
				}
				glEnd();
				FLOAT64 diag = LhsMatrix(index,index);
				if( fabs(diag) > 0 ) {
					for( uint n=0; n<LhsMatrix.index[index].size(); n++ ) {
						uint idx = LhsMatrix.index[index][n];
						if( idx != index ) {
							vec2d vdir = (g.centers[idx]-g.centers[index]).normal();
							FLOAT64 r = 1.0/16.0;
							glRasterPos2d(g.centers[index][0]-0.25*r+r*vdir[0],g.centers[index][1]+r*vdir[1]);
							drawBitmapString(format_str("%.2f", -LhsMatrix.value[index][n]/diag));
						}
					}
				}
				glColor4d(1.0,1.0,1.0,1.0);
				glPointSize(3.0);
				glBegin(GL_POINTS);
				glVertex2f(g.centers[index][0],g.centers[index][1]);
				glEnd();
				glPointSize(1.0);
			}
			break;
		}
		case NAME: {
			glColor4f(1.0,1.0,1.0,1.0);
			glRasterPos2d(0.07,0.92);
			drawBitmapString("FVM",GLUT_BITMAP_HELVETICA_18);
			break;
		}
	}
}