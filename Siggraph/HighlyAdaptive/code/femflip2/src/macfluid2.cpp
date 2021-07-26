/*
 *	macfluid2.cpp
 *	
 *	Created by Ryoichi Ando on 2012/04/30
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"
#include "macfluid2.h"
#include "levelset2.h"
#include "fastmarch2.h"
#include "opengl.h"
#include "util2.h"
#include "pcgsolver/matutil.h"
#include <stdio.h>
#include <fstream>

using namespace std;

macfluid2::macfluid2() {
	init(0,NULL);
	extrapolate_dist = 1.0;
	surf_order = 2;
	variation = 1;
	dt=0.1;
	tension = 0.0;
	frac_eps = 1e-8;
}

void macfluid2::init( uint gn, const levelset2 *hint ) {
	// Init variables
	solid = NULL;
	fluid = NULL;
	this->gn = gn;
	dx = 1.0/gn;
	
	for(uint dim=0; dim<DIM; dim++ ) {
		u[dim].resize(gn+(dim==0),gn+(dim==1));
		u[dim].clear();
		area[dim].resize(gn+(dim==0),gn+(dim==1));
		area[dim].clear();
		flag[dim].resize(gn+(dim==0),gn+(dim==1));
		flag[dim].clear();
		gradp[dim].resize(gn+(dim==0),gn+(dim==1));
		gradp[dim].clear();
		facetfluidLS[dim].resize(gn+(dim==0),gn+(dim==1));
		facetfluidLS[dim].clear();
		facetsolidLS[dim].resize(gn+(dim==0),gn+(dim==1));
		facetsolidLS[dim].clear();
	}
	
	volume.resize(gn,gn);
	volume.clear();
    p.resize(gn,gn);
	p.clear();
	div.resize(gn,gn);
	div.clear();
	cellfluidLS.resize(gn,gn);
	cellfluidLS.clear();
	cellCurvature.resize(gn,gn);
	cellCurvature.clear();
}

void macfluid2::setParameter( int name, FLOAT64 value ) {
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

void macfluid2::setTimestep( FLOAT64 dt ) {
	this->dt = dt;
}

void macfluid2::setupSolidLevelset( const levelset2 *solid ) {
	this->solid = solid;
	// Compute solid fractions and volumes
	computeFractions( area, volume, solid );
	
	// Compute facet centered solid levelset
	for( uint dim=0; dim<DIM; dim++ ) {
		FOR_EACH(gn+(dim==0),gn+(dim==1)) {
			facetsolidLS[dim][i][j] = solid->evalLevelset(vec2d(i*dx+0.5*(dim!=0)*dx,j*dx+0.5*(dim!=1)*dx));
		} END_FOR
	}
}

void macfluid2::setupFluidLevelset( const levelset2 *fluid ) {
	this->fluid = fluid;
	// Build cell levelset
	FOR_EACH(gn,gn) {
		cellfluidLS[i][j] = fluid->evalLevelset(vec2d(dx*(i+0.5),dx*(j+0.5)));
	} END_FOR
	
	// Build facet centered fluid levelset
	for( uint dim=0; dim<DIM; dim++ ) {
		FOR_EACH(gn+(dim==0),gn+(dim==1)) {
			facetfluidLS[dim][i][j] = fluid->evalLevelset(vec2d(i*dx+0.5*(dim!=0)*dx,j*dx+0.5*(dim!=1)*dx));
		} END_FOR
	}
}

void macfluid2::getSamplePoints( std::vector<vec2d> &pos ) {
	// Set up cell-face positions
	uint usize = (gn+1)*gn;
	pos.resize(DIM*usize);
	for(uint dim=0; dim<DIM; dim++ ) {
		FOR_EACH(gn+(dim==0),gn+(dim==1)) {
			uint idx = i+j*(gn+(dim==0))+dim*usize;
			pos[idx] = vec2d(i*dx+0.5*(dim!=0)*dx,j*dx+0.5*(dim!=1)*dx);
		} END_FOR
	}
}

void macfluid2::setupVelocity( const std::vector<vec2d> &vel, const std::vector<bool> &mapped ) {
	uint usize = (gn+1)*gn;
	for(uint dim=0; dim<DIM; dim++ ) {
		FOR_EACH(gn+(dim==0),gn+(dim==1)) {
			uint idx = i+j*(gn+(dim==0))+dim*usize;
			u[dim][i][j] = vel[idx][dim];
			flag[dim][i][j] = area[dim][i][j] ? mapped[idx] : false;
		} END_FOR
	}
}

void macfluid2::project() {	
	extrapolate(u,flag);
	project(u,cellfluidLS,solid,variation,surf_order,dt);
	extrapolate(gradp,flag);
}

static FLOAT64 fraction( FLOAT64 phi0, FLOAT64 phi1 ) {
	return phi0*phi1 > 0 ? 1 - (phi0 > 0) : -fmin(phi0,phi1)/fmax(fabs(phi1-phi0),1.0e-6);
}

void macfluid2::computeFractions( array2<FLOAT64> area[], array2<FLOAT64> &volume, const levelset2 *solid ) {
    // Nodal Levelset
    array2<FLOAT64> nodalSolidLS(gn+1,gn+1);
    
    // Setup solid levelset
	PARALLEL_FOR FOR_EACH(gn+1,gn+1) {	
		vec2d p( dx*i, dx*j );
		nodalSolidLS[i][j] = solid->evalLevelset(p);
	} END_FOR
    
	// For each dimension
	for( int dim=0; dim<DIM; dim++ ) {
		PARALLEL_FOR FOR_EACH(gn+(dim==0),gn+(dim==1)) {
			FLOAT64 frac = 1.0 - fraction(nodalSolidLS[i][j],nodalSolidLS[i+(dim!=0)][j+(dim!=1)]);
			if( frac ) frac = fmax(frac,1e-6);
			area[dim][i][j] = frac;
		} END_FOR
	}
	
	// Compute volume fraction using marching cubed mesh
	PARALLEL_FOR FOR_EACH(gn,gn) {
		vec2d p[8];
		FLOAT64 v = 0.0;
		int pnum;
		levelset2::marchPoints( i, j, nodalSolidLS, p, pnum );
		for( int m=0; m<pnum; m++ ) {
			v += 0.5*(p[m][0]*p[(m+1)%pnum][1]-p[m][1]*p[(m+1)%pnum][0]);
		}
		volume[i][j] = 1.0-v/util2::square(dx);
	} END_FOR
}

vec2d macfluid2::getVelocity( vec2d p ) const {
	return getVelocity( p, u );
}

FLOAT64 macfluid2::getDivergence( vec2d p ) const {
	return util2::interp<FLOAT64>(vec2d(gn*p[0]-0.5,gn*p[1]-0.5),div);
}

vec2d macfluid2::getVelocity( vec2d p, const array2<FLOAT64> u[] ) const {
	return vec2d(util2::interp<FLOAT64>(vec2d(gn*p[0],gn*p[1]-0.5),u[0]),
				 util2::interp<FLOAT64>(vec2d(gn*p[0]-0.5,gn*p[1]),u[1]));
}

vec2d macfluid2::getPressureGradient( vec2d p ) const {
	return vec2d(util2::interp<FLOAT64>(vec2d(gn*p[0],gn*p[1]-0.5),gradp[0]),
				 util2::interp<FLOAT64>(vec2d(gn*p[0]-0.5,gn*p[1]),gradp[1]));
}

void macfluid2::copyGridVelocity( array2<FLOAT64> dest[], const array2<FLOAT64> src[] ) {
	for( int dim=0; dim<DIM; dim++ ) {
		dest[dim] = src[dim];
	}
}

void macfluid2::extrapolate( array2<FLOAT64> u[], const array2<bool> flag[] ) {
	// For Each Dimension
	for( int dim=0; dim<DIM; dim++ ) {
		// Setup a network
		FLOAT64 dx = 1.0/gn;
		array2<fastmarch2<FLOAT64>::node2> nodeArray;
		nodeArray.resize(gn+(dim==0),gn+(dim==1));
		PARALLEL_FOR FOR_EACH(gn+(dim==0),gn+(dim==1)) {
			vec2d p = vec2d(i*dx,j*dx);
			nodeArray[i][j].p = p;
			nodeArray[i][j].fixed = flag[dim][i][j];
			nodeArray[i][j].levelset = fmax(facetfluidLS[dim][i][j],-facetsolidLS[dim][i][j]);
			nodeArray[i][j].value = u[dim][i][j];
			int q[][2] = {{i-1,j},{i+1,j},{i,j-1},{i,j+1}};
			for( uint n=0; n<4; n++ ) {
				int ni = q[n][0];
				int nj = q[n][1];
				if( ni>=0 && ni<gn+(dim==0) && nj>=0 && nj<gn+(dim==1) ) {
					nodeArray[i][j].p2p.push_back(&nodeArray[ni][nj]);
				}
			}
		} END_FOR
		std::vector<fastmarch2<FLOAT64>::node2 *> nodes(gn*(gn+1));
		uint index = 0;
		FOR_EACH(gn+(dim==0),gn+(dim==1)) {
			nodes[index++] = &nodeArray[i][j];
		} END_FOR
		
		// Perform fast march
		fastmarch2<FLOAT64>::fastMarch(nodes,extrapolate_dist,-1.0,0);
		
		// Pickup extrapolated velocity
		PARALLEL_FOR FOR_EACH(gn+(dim==0),gn+(dim==1)) {
			u[dim][i][j] = nodeArray[i][j].fixed ? nodeArray[i][j].value : 0.0;
		} END_FOR
	}
}

void macfluid2::project( array2<FLOAT64> u[], const array2<FLOAT64> &cellfluidLS, const levelset2* solid,
					bool variation, uint surf_order, FLOAT64 dt ) {
	// Allocate matrix and right hand side
    SparseMatrix<FLOAT64> matrix;
	std::vector<FLOAT64> rhs(gn*gn);            // Right hand side vector
	std::vector<FLOAT64> pressure(gn*gn);       // Pressure in 1D
	
    // Copy to one dimension vector
    FOR_EACH(gn,gn) {
        rhs[i+j*gn] = 0.0;
    } END_FOR

	// Clear all except for pressure
    matrix.clear();
	matrix.resize(gn*gn);
	
	// Resample curvature on a coarse grid
	FOR_EACH(gn,gn) {
		cellCurvature[i][j] = tension*fluid->evalCurvature(vec2d(dx*(i+0.5),dx*(j+0.5)));
	} END_FOR
	
	// Compute the right hand side and the matrix for pressure solve
	PARALLEL_FOR FOR_EACH(gn,gn) {
		FLOAT64 phi0 = cellfluidLS[i][j];	// Fluid levelset at this cell
		FLOAT64 dx2 = util2::square(dx);
		if( phi0 <= 0 ) {					// If in the liquid...
			// Cell indices: Left, Bottom, Right, Top...
			int q[][2] = { {i-1,j}, {i,j-1}, {i+1,j}, {i,j+1} };
			// Set face indices (Direcion,i-index,j-index) followed by the order above
			int f[][3] = { {0,i,j}, {1,i,j}, {0,i+1,j}, {1,i,j+1} };
			for( int n=0; n<4; n++ ) {
                if( q[n][0] < 0 || q[n][0] > gn-1 || q[n][1] < 0 || q[n][1] > gn-1 ) continue;
				FLOAT64 phi1 = cellfluidLS.clampFetch(q[n][0],q[n][1]);		// Fluid levelset outside
				FLOAT64 rho = fmax(frac_eps,fraction(phi0,phi1));			// Fluid density
				if( surf_order==1 ) rho = 1.0;
				FLOAT64 frac = area[f[n][0]][f[n][1]][f[n][2]];				// Cell face fraction
				FLOAT64 scale = frac*dt/dx2;
				FLOAT64 solidLS[2] = {	solid->evalLevelset(vec2d(dx*(i+0.5),dx*(j+0.5))),
                                        solid->evalLevelset(vec2d(dx*(q[n][0]+0.5),dx*(q[n][1]+0.5))) };
				// If not use variational method, set fraction 0 or 1 on cell faces
				if( ! variation ) frac = (solidLS[0] < 0 || solidLS[1] < 0 ) ? 0.0 : 1.0;
				if( phi1 <= 0 ) {											// If adjacent to fluid cell
                    matrix.add_to_element(i+j*gn,i+j*gn,scale);
                    matrix.add_to_element(i+j*gn,q[n][0]+q[n][1]*gn,-scale);
				} else {													// If adjacent to empty cell
                    matrix.add_to_element(i+j*gn,i+j*gn,scale/rho);
				}
				FLOAT64 sgn = ( i < q[n][0] || j < q[n][1] ) ? -1.0 : 1.0;
                rhs[i+j*gn] += sgn*frac*u[f[n][0]][f[n][1]][f[n][2]] / dx;  // Compute negative divergence
				
				// Give surface tension force
				if( tension ) {
					if( phi0 * phi1 < 0 ) {
						vec2d spos = (1.0-rho)*vec2d(dx*(i+0.5),dx*(j+0.5))+rho*vec2d(dx*(q[n][0]+0.5),dx*(q[n][1]+0.5));
						FLOAT64 tensF = tension*util2::interp<FLOAT64>(vec2d(gn*spos[0]-0.5,gn*spos[1]-0.5),cellCurvature);
						rhs[i+j*gn] -= scale*tensF/rho;
					}
				}
            }
		}
	} END_FOR
	
	// Copy divergence in 2D array for visualization
	div.resize(gn,gn);
    FOR_EACH(gn,gn) {
        div[i][j] = rhs[i+j*gn];
    } END_FOR
	
	// Compact the matrix
	vector<int> idx_map;
	compress( matrix, rhs, idx_map );
	
#if 0
	ofstream file;
	file.open("Lhs_MAC_matrix.m");
	matrix.write_matlab(file,matrix.index.size(),matrix.index.size(),"Lhs_MAC");
	file.close();
	file.open("Rhs_MAC_vector.m");
	write_matlab1d(file,rhs,"Rhs_MAC");
	file.close();
#endif
	
	// Set initial guess
    FOR_EACH(gn,gn) {
		uint idx = i+gn*j;
        if( idx_map[idx] >= 0 ) pressure[idx_map[idx]] = p[i][j];
    } END_FOR
	
	// Solve by incomplete cholesky factorizaion preconditioned conjugate gradient method
    PCGSolver<FLOAT64> solver;
    FLOAT64 residual_out, msec;
    int iterations;
	bool converged = solver.solve( matrix, rhs, pressure, residual_out, iterations, msec );
	if( residual_out > 1.0 ) converged = false;
	if( ! converged ) {
		email::print( "PCG did not converge.\n" );
		email::send();
		exit(0);
	} else {
        dump( "MAC: PCG Converged ! %d iterations and took %.3f msec with %d unknowns.\n", iterations, msec, matrix.index.size());
    }
    
    // Copy the result
    FOR_EACH(gn,gn) {
		uint idx = i+gn*j;
        p[i][j] = idx_map[idx] >= 0 ? pressure[idx_map[idx]] : 0.0;
    } END_FOR

	// Subtract pressure gradient ( We must take into account ghost pressures here... )
	for( int dim=0; dim<DIM; dim++ ) {
		PARALLEL_FOR FOR_EACH(gn+(dim==0),gn+(dim==1)) {
			int ni = i-(dim==0);
			int nj = j-(dim==1);

			// Gather solid levelset info
			FLOAT64 solidLS[2] = {	solid->evalLevelset(vec2d(dx*(i+0.5),dx*(j+0.5))),
									solid->evalLevelset(vec2d(dx*(ni+0.5),dx*(nj+0.5))) };
			// If not use variational method, do not subtract pressure gradient on solid faces, but instead set zero
			if( ! variation ) {
				if( solidLS[0] * solidLS[1] < 0 ) {
					u[dim][i][j] = 0.0;
					if( solidLS[0] < 0 ) p.safeSet(i,j,p[ni][nj]);
					else p.safeSet(ni,nj,p.clampFetch(i,j));
					// Nothing to do for this face, so skip to next
					continue;
				}
			}
			// Fetch forward and backward pressures
			FLOAT64 press_f = p.clampFetch(i,j);
			FLOAT64 press_b = p.clampFetch(ni,nj);
			
			// Fetch liquid levelset here
			FLOAT64 fluidLS[2] = { cellfluidLS.clampFetch(i,j), cellfluidLS.clampFetch(ni,nj) };
			
			// Compute density
			FLOAT64 rho = fmax(frac_eps,fraction(fluidLS[0],fluidLS[1]));
			if( surf_order==1 ) rho = 1.0;
			
			// Now compute the pressure gradient
			gradp[dim][i][j] = area[dim][i][j] ? (press_f-press_b)/(dx*rho) : 0.0;
			
			// Give surface tension force
			if( tension ) {
				if( fluidLS[0] * fluidLS[1] < 0 ) {
					FLOAT64 sgn = fluidLS[0] < 0.0 ? -1.0 : 1.0;
					vec2d spos;
					if( sgn < 0.0 ) spos = (1.0-rho)*vec2d(dx*(i+0.5),dx*(j+0.5))+rho*vec2d(dx*(ni+0.5),dx*(nj+0.5));
					else spos = rho*vec2d(dx*(i+0.5),dx*(j+0.5))+(1.0-rho)*vec2d(dx*(ni+0.5),dx*(nj+0.5));
					FLOAT64 tensF = tension*util2::interp<FLOAT64>(vec2d(gn*spos[0]-0.5,gn*spos[1]-0.5),cellCurvature);
					gradp[dim][i][j] -= sgn*tensF/(rho*dx);
				}
			}
		} END_FOR
	}
}

static void drawMarchingCube( const array2<FLOAT64> &L ) {	
	uint w = L.size().w;
	uint h = L.size().h;
	FOR_EACH(w-1,h-1) {
		vec2d p[8];
		int pnum;
		glBegin(GL_TRIANGLE_FAN);
		levelset2::marchPoints( i, j, L, p, pnum );
		for( int m=0; m<pnum; m++ ) {
			glVertex2f(p[m][0],p[m][1]);
		}
		glEnd();
	} END_FOR;
}

static void drawBitmapString( const char *string, void *font=NULL ) {
	if( ! font ) font = GLUT_BITMAP_HELVETICA_10;
	while (*string) glutBitmapCharacter(font, *string++);
}

void macfluid2::render( int name, vec2d mousePos ) const {
	switch( name ) {
		case DIVERGENCE:
		case PRESSURE: {
			const array2<FLOAT64> &q = (name == DIVERGENCE) ? div : p;
			FLOAT64 maxv = util2::arrayMax<FLOAT64>(q);
			FLOAT64 minv = util2::arrayMin<FLOAT64>(q);
			
			FOR_EACH(gn,gn) {
				if( fabs(fluid->evalLevelset(vec2d(dx*(i+0.5),dx*(j+0.5)))) < dx ) {
					maxv = fmax(maxv,q[i][j]);
					minv = fmin(minv,q[i][j]);
				}
			} END_FOR
			
			FLOAT64 det = maxv - minv;
			if( det > 0 ) {
				FOR_EACH(gn,gn) {
					if( cellfluidLS[i][j] <= 0.0 && volume[i][j] > 0.0 ) {
						FLOAT64 normp = 2.0*(q[i][j]-minv)/det-1.0;
						glColor4d(normp>0,0.3,normp<=0,fabs(normp));
						glBegin(GL_QUADS);
						glVertex2d(i*dx,j*dx);
						glVertex2d((i+1)*dx,j*dx);
						glVertex2d((i+1)*dx,(j+1)*dx);
						glVertex2d(i*dx,(j+1)*dx);
						glEnd();
					}
				} END_FOR
			}
			break;
		}
		case SOLID: {
			// Draw solid levelset
			glColor4d(0.7,0.7,0.3,0.6);
			FOR_EACH(gn,gn) {
				vector<vec2d> nodes(4);
				nodes[0] = vec2d(dx*i,dx*j);
				nodes[1] = vec2d(dx*(i+1),dx*j);
				nodes[2] = vec2d(dx*(i+1),dx*(j+1));
				nodes[3] = vec2d(dx*i,dx*(j+1));
				
				vector<FLOAT64> levelsets(4);
				for( uint k=0; k<4; k++ ) levelsets[k] = solid->evalLevelset(nodes[k]);

				vector<vec2d> points = util2::marchPoints(nodes,levelsets);
				glBegin(GL_POLYGON);
					for( int m=0; m<points.size(); m++ ) glVertex2f(points[m][0],points[m][1]);
				glEnd();
			} END_FOR
			break;
		}
		case FLUID: {
			// Draw fluid levelset
			glColor4d(0.3,0.3,0.8,0.6);
			glPushMatrix();
			glScaled(gn/(FLOAT64)(gn+1),gn/(FLOAT64)(gn+1),1.0);
			glTranslated(0.5*dx,0.5*dx,0.0);
			drawMarchingCube(cellfluidLS);
			glPopMatrix();
			break;
		}
		case VELOCITY: {
			// Draw combined flow
			glColor4d(1.0,1.0,0.0,0.75);
			FOR_EACH(gn,gn) {
				FLOAT64 p[2] = {i*dx+dx/2.0,j*dx+dx/2.0};
				FLOAT64 vel[2] = {0.5*u[0][i][j]+0.5*u[0][i+1][j],0.5*u[1][i][j]+0.5*u[1][i][j+1]};
				glBegin(GL_LINES);
				glVertex2d(p[0],p[1]);
				glVertex2d(p[0]+dx*vel[0],p[1]+dx*vel[1]);
				glEnd();
			} END_FOR
			break;
		}
		case MESH: {
			// Draw grid
			glColor4f(1.0,1.0,1.0,0.2);
			glBegin(GL_LINES);
			for( int i=0; i<=gn; i++ ) {
				glVertex2f(i*dx,0.0);
				glVertex2f(i*dx,1.0);
			}
			for( int j=0; j<=gn; j++ ) {
				glVertex2f(0.0,j*dx);
				glVertex2f(1.0,j*dx);
			}
			glEnd();
			break;
		}
		
		case MATRIX_CONNECTION: {
			FLOAT64 mind = 999.0;
			uint index_i = 0;
			uint index_j = 0;
			FOR_EACH(gn,gn) {
				vec2d p = vec2d(i*dx+dx/2.0,j*dx+dx/2.0);
				FLOAT64 d = (p-mousePos).len();
				if( d < dx && d < mind ) {
					mind = d;
					index_i = i;
					index_j = j;
				}
			} END_FOR
			
			// Cell indices: Left, Bottom, Right, Top... (Should match to the order defined in macros.h !)
			int q[][2] = { {index_i-1,index_j}, {index_i,index_j-1}, {index_i+1,index_j}, {index_i,index_j+1} };
			FLOAT64 diag = 0.0;
			FLOAT64 values[4] = { 0.0, 0.0, 0.0, 0.0 };
			// Set face indices (Direcion,i-index,j-index) followed by the order above
			int f[][3] = { {0,index_i,index_j}, {1,index_i,index_j}, {0,index_i+1,index_j}, {1,index_i,index_j+1} };
			for( int n=0; n<4; n++ ) {
                if( q[n][0] < 0 || q[n][1] > gn-1 ) continue;
				FLOAT64 scale = dt/(dx*dx);
				FLOAT64 frac = area[f[n][0]][f[n][1]][f[n][2]];				// Cell face fraction
				FLOAT64 solidLS[2] = {	solid->evalLevelset(vec2d(dx*(index_i+0.5),dx*(index_j+0.5))),
										solid->evalLevelset(vec2d(dx*(q[n][0]+0.5),dx*(q[n][1]+0.5))) };
				// If not use variational method, set fraction 0 or 1 on cell faces
				if( ! variation ) frac = (solidLS[0] < 0 || solidLS[1] < 0 ) ? 0.0 : 1.0;
				FLOAT64 term = frac*scale;
				diag += term;
				values[n] -= term;
            }
			glColor4d(1.0,1.0,1.0,0.75);
			glBegin(GL_LINES);
			for( int n=0; n<4; n++ ) {
				if( values[n] && diag ) {
					glVertex2f((index_i+0.5)*dx,(index_j+0.5)*dx);
					glVertex2f((q[n][0]+0.5)*dx,(q[n][1]+0.5)*dx);
				}
			}
			glEnd();
			for( int n=0; n<4; n++ ) {
				if( values[n] && diag ) {
					vec2d p0((index_i+0.5)*dx,(index_j+0.5)*dx);
					vec2d p1((q[n][0]+0.5)*dx,(q[n][1]+0.5)*dx);
					vec2d vdir = (p1-p0).normal();
					FLOAT64 r = 1.0/16.0;
					glRasterPos2d(p0[0]-0.25*r+r*vdir[0],p0[1]+r*vdir[1]);
					drawBitmapString(format_str("%.2f", -values[n]/diag));
				}
			}
			glColor4d(1.0,1.0,1.0,1.0);
			glPointSize(3.0);
			glBegin(GL_POINTS);
			glVertex2f((index_i+0.5)*dx,(index_j+0.5)*dx);
			glEnd();
			glPointSize(1.0);
			break;
		}
		case NAME: {
			glColor4f(1.0,1.0,1.0,1.0);
			glRasterPos2d(0.07,0.92);
			drawBitmapString("MAC",GLUT_BITMAP_HELVETICA_18);
			break;
		}
	}
}
