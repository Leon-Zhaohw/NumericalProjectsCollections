/*
 *	macfluid3.cpp
 *	
 *	Created by Ryoichi Ando on 2012/05/03
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"
#include "macfluid3.h"
#include "levelset3.h"
#include "fastmarch3.h"
#include "opengl.h"
#include "util3.h"
#include "pcgsolver/matutil.h"
#include <stdio.h>
#include <fstream>

using namespace std;

macfluid3::macfluid3() {
	extrapolate_dist = 1.0;
	wall_separate = true;
	surf_order = 2;
	variation = 1;
	dt=0.1;
}

void macfluid3::init( uint gn, const levelset3 *hint ) {
	// Init variables
	solid = NULL;
	fluid = NULL;
	this->gn = gn;
	dx = 1.0/gn;
	tick(); dump( "Building MAC mesh..." );
	for(uint dim=0; dim<DIM; dim++ ) {
		u[dim].resize(gn+(dim==0),gn+(dim==1),gn+(dim==2));
		u[dim].clear();
		area[dim].resize(gn+(dim==0),gn+(dim==1),gn+(dim==2));
		area[dim].clear();
		flag[dim].resize(gn+(dim==0),gn+(dim==1),gn+(dim==2));
		flag[dim].clear();
		gradp[dim].resize(gn+(dim==0),gn+(dim==1),gn+(dim==2));
		gradp[dim].clear();
		facetfluidLS[dim].resize(gn+(dim==0),gn+(dim==1),gn+(dim==2));
		facetfluidLS[dim].clear();
		facetsolidLS[dim].resize(gn+(dim==0),gn+(dim==1),gn+(dim==2));
		facetsolidLS[dim].clear();
	}
	
	volume.resize(gn,gn,gn);
	volume.clear();
    p.resize(gn,gn,gn);
	p.clear();
	div.resize(gn,gn,gn);
	div.clear();
	
	cellfluidLS.resize(gn,gn,gn);
	cellfluidLS.clear();
	dump( "Done. Took %s.\n", stock("MAC_mesh"));
}

void macfluid3::setParameter( int name, FLOAT64 value ) {
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

void macfluid3::setTimestep( FLOAT64 dt ) {
	this->dt = dt;
}

void macfluid3::setupSolidLevelset( const levelset3 *solid ) {
	this->solid = solid;
	// Compute solid fractions and volumes
	computeFractions( area, volume, solid );
	
	// Compute facet centered solid levelset
	for( uint dim=0; dim<DIM; dim++ ) {
		FOR_EACH(gn+(dim==0),gn+(dim==1),gn+(dim==2)) {
			facetsolidLS[dim][i][j][k] = solid->evalLevelset(vec3d(i*dx+0.5*(dim!=0)*dx,
																   j*dx+0.5*(dim!=1)*dx,
																   k*dx+0.5*(dim!=2)*dx));
		} END_FOR
	}
}

void macfluid3::setupFluidLevelset( const levelset3 *fluid ) {
	this->fluid = fluid;
	// Build cell levelset
	FOR_EACH(gn,gn,gn) {
		cellfluidLS[i][j][k] = fluid->evalLevelset(vec3d(dx*(i+0.5),dx*(j+0.5),dx*(k+0.5)));
	} END_FOR
	
	// Build facet centered fluid levelset
	for( uint dim=0; dim<DIM; dim++ ) {
		FOR_EACH(gn+(dim==0),gn+(dim==1),gn+(dim==2)) {
			facetfluidLS[dim][i][j][k] = fluid->evalLevelset(vec3d(i*dx+0.5*(dim!=0)*dx,
																   j*dx+0.5*(dim!=1)*dx,
																   k*dx+0.5*(dim!=2)*dx));
		} END_FOR
	}
}

void macfluid3::getSamplePoints( std::vector<vec3d> &pos ) {
	// Set up cell-face positions
	pos.clear();
	uint usize = (gn+1)*gn*gn;
	pos.resize(DIM*usize);
	for(uint dim=0; dim<DIM; dim++ ) {
		FOR_EACH(gn+(dim==0),gn+(dim==1),gn+(dim==2)) {
			uint idx = i+j*(gn+(dim==0))+k*(gn+(dim==0))*(gn+(dim==1))+dim*usize;
			pos[idx] = vec3d(i*dx+0.5*(dim!=0)*dx,j*dx+0.5*(dim!=1)*dx,k*dx+0.5*(dim!=2)*dx);
		} END_FOR
	}
}

void macfluid3::setupVelocity( const std::vector<vec3d> &vel, const std::vector<bool> &mapped ) {
	uint usize = (gn+1)*gn*gn;
	for(uint dim=0; dim<DIM; dim++ ) {
		FOR_EACH(gn+(dim==0),gn+(dim==1),gn+(dim==2)) {
			uint idx = i+j*(gn+(dim==0))+k*(gn+(dim==0))*(gn+(dim==1))+dim*usize;
			u[dim][i][j][k] = vel[idx][dim];
			flag[dim][i][j][k] = area[dim][i][j][k] ? mapped[idx] : false;
		} END_FOR
	}
}

void macfluid3::project() {	
	tick(); dump( "Extrapolating velocity..." );
	extrapolate(u,flag);
	dump( "Done. Took %s.\n", stock("MAC_velocity_extrapolate"));
	
	project(u,cellfluidLS,solid,variation,surf_order,wall_separate,dt);
	
	tick(); dump( "Extrapolating gradient..." );
	extrapolate(gradp,flag);
	dump( "Done. Took %s.\n", stock("MAC_gradient_extrapolate"));
}

static FLOAT64 fraction( FLOAT64 phi0, FLOAT64 phi1 ) {
	return phi0*phi1 > 0 ? 1 - (phi0 > 0) : -fmin(phi0,phi1)/fmax(fabs(phi1-phi0),1.0e-6);
}

void macfluid3::computeFractions( array3<FLOAT64> area[], array3<FLOAT64> &volume, const levelset3 *solid ) {
    // Nodal Levelset
    array3<FLOAT64> nodalSolidLS(gn+1,gn+1,gn+1);
    
    // Setup solid levelset
	PARALLEL_FOR FOR_EACH(gn+1,gn+1,gn+1) {	
		vec3d p( dx*i, dx*j, dx*k );
		nodalSolidLS[i][j][k] = solid->evalLevelset(p);
	} END_FOR
	
	// Compute cell face fractions
	for(uint dim=0; dim<DIM; dim++ ) {
		PARALLEL_FOR FOR_EACH(gn+(dim==0),gn+(dim==1),gn+(dim==2)) {
			// Collect face 4 levelset value
			std::vector<FLOAT64> lv(4);
			std::vector<vec3d> p(4);
			int q[4][2] = { {0,0}, {1,0}, {1,1}, {0,1} };
			for( uint n=0; n<4; n++ ) {
				uint xi = q[n][0];
				uint yi = q[n][1];
				if(dim==0) lv[n] = nodalSolidLS[i][j+xi][k+yi];
				if(dim==1) lv[n] = nodalSolidLS[i+yi][j][k+xi];
				if(dim==2) lv[n] = nodalSolidLS[i+xi][j+yi][k];
				p[n][0] = xi;
				p[n][1] = yi;
				p[n][2] = 0.0;
			}
			
			// Compute Area
			std::vector<vec3d> mesh = util3::march2DPoints(p,lv);
			area[dim][i][j][k] = 1.0-util3::compute2DArea(mesh);
		} END_FOR
	}
	
	// Compute Approximated volume fraction by sub-voxel sum
	PARALLEL_FOR FOR_EACH(gn,gn,gn) {
		uint nRes = 5;
		uint sum = 0;
		FLOAT64 w = 1.0/(nRes-1);
		for( uint si=0; si<nRes; si++ ) for( uint sj=0; sj<nRes; sj++ ) for( uint sk=0; sk<nRes; sk++ ) {
			vec3d shift = vec3d(si*w,sj*w,sk*w)+0.5*vec3d(w,w,w);
			sum += (solid->evalLevelset(vec3d(i*dx,j*dx,k*dx)+dx*shift)>0.0);
		}
		volume[i][j][k] = sum*(w*w*w);
	} END_FOR
}

vec3d macfluid3::getVelocity( vec3d p ) const {
	return getVelocity( p, u );
}

vec3d macfluid3::getVelocity( vec3d p, const array3<FLOAT64> u[] ) const {
	return vec3d(util3::interp<FLOAT64>(vec3d(gn*p[0],gn*p[1]-0.5,gn*p[2]-0.5),u[0]),
				 util3::interp<FLOAT64>(vec3d(gn*p[0]-0.5,gn*p[1],gn*p[2]-0.5),u[1]),
				 util3::interp<FLOAT64>(vec3d(gn*p[0]-0.5,gn*p[1]-0.5,gn*p[2]),u[2]));
}

FLOAT64 macfluid3::getDivergence( vec3d p ) const {
	return util3::interp<FLOAT64>(vec3d(gn*p[0]-0.5,gn*p[1]-0.5,gn*p[2]-0.5),div);
}

vec3d macfluid3::getPressureGradient( vec3d p ) const {
	return vec3d(util3::interp<FLOAT64>(vec3d(gn*p[0],gn*p[1]-0.5,gn*p[2]-0.5),gradp[0]),
				 util3::interp<FLOAT64>(vec3d(gn*p[0]-0.5,gn*p[1],gn*p[2]-0.5),gradp[1]),
				 util3::interp<FLOAT64>(vec3d(gn*p[0]-0.5,gn*p[1]-0.5,gn*p[2]),gradp[2]));
}

void macfluid3::copyGridVelocity( array3<FLOAT64> dest[], const array3<FLOAT64> src[] ) {
	for( int dim=0; dim<DIM; dim++ ) {
		dest[dim] = src[dim];
	}
}

void macfluid3::extrapolate( array3<FLOAT64> u[], const array3<bool> flag[] ) {
	// For Each Dimension
	for( int dim=0; dim<DIM; dim++ ) {
		// Setup a network
		FLOAT64 dx = 1.0/gn;
		array3<fastmarch3<FLOAT64>::node3> nodeArray;
		nodeArray.resize(gn+(dim==0),gn+(dim==1),gn+(dim==2));
		PARALLEL_FOR FOR_EACH(gn+(dim==0),gn+(dim==1),gn+(dim==2)) {
			vec3d p = vec3d(i*dx,j*dx,k*dx);
			nodeArray[i][j][k].p = p;
			nodeArray[i][j][k].fixed = flag[dim][i][j][k];
			nodeArray[i][j][k].levelset = fmax(facetfluidLS[dim][i][j][k],-facetsolidLS[dim][i][j][k]);
			nodeArray[i][j][k].value = u[dim][i][j][k];
			int q[][3] = {{i-1,j,k},{i+1,j,k},{i,j-1,k},{i,j+1,k},{i,j,k-1},{i,j,k+1}};
			for( uint n=0; n<6; n++ ) {
				int ni = q[n][0];
				int nj = q[n][1];
				int nk = q[n][2];
				if( ni>=0 && ni<gn+(dim==0) && nj>=0 && nj<gn+(dim==1) && nk>=0 && nk<gn+(dim==2)) {
					nodeArray[i][j][k].p2p.push_back(&nodeArray[ni][nj][nk]);
				}
			}
		} END_FOR
		std::vector<fastmarch3<FLOAT64>::node3 *> nodes(gn*gn*(gn+1));
		uint index = 0;
		FOR_EACH(gn+(dim==0),gn+(dim==1),gn+(dim==2)) {
			nodes[index++] = &nodeArray[i][j][k];
		} END_FOR
		
		// Perform fast march
		fastmarch3<FLOAT64>::fastMarch(nodes,extrapolate_dist,-1.0,0);
		
		// Pickup extrapolated levelsets
		PARALLEL_FOR FOR_EACH(gn+(dim==0),gn+(dim==1),gn+(dim==2)) {
			u[dim][i][j][k] = nodeArray[i][j][k].fixed ? nodeArray[i][j][k].value : 0.0;
		} END_FOR
	}
}

#define X(i,j,k) ((i)+(j)*(gn)+(k)*(gn*gn))
void macfluid3::project( array3<FLOAT64> u[], const array3<FLOAT64> &cellfluidLS, const levelset3* solid,
						 bool variation, uint surf_order, bool wallSepration, FLOAT64 dt ) {
	// Allocate matrix and right hand side
    SparseMatrix<FLOAT64> matrix;
	vector<FLOAT64> rhs(gn*gn*gn);            // Right hand side vector
    vector<FLOAT64> pressure(gn*gn*gn);       // Pressure in 1D
	FLOAT64 eps = 1e-3;
	
    // Copy to one dimension vector
    FOR_EACH(gn,gn,gn) {
        rhs[X(i,j,k)] = 0.0;
    } END_FOR
	
	// Clear all except for pressure
    matrix.clear();
	matrix.resize(gn*gn*gn);
	
	tick(); dump( "Building matrix and computing right hand side..." );
	// Compute the right hand side and the matrix for pressure solve
	PARALLEL_FOR FOR_EACH(gn,gn,gn) {
		FLOAT64 phi0 = cellfluidLS[i][j][k];	// Fluid levelset at this cell
		FLOAT64 dx2 = util3::square(dx);
		if( phi0 <= 0 ) {					// If in the liquid...
			// Cell indices: Left, Bottom, Right, Top...
			int q[][3] = { {i-1,j,k}, {i,j-1,k}, {i,j,k-1}, {i+1,j,k}, {i,j+1,k}, {i,j,k+1} };
			// Set face indices (Direcion,i-index,j-index) followed by the order above
			int f[][4] = { {0,i,j,k}, {1,i,j,k}, {2,i,j,k}, {0,i+1,j,k}, {1,i,j+1,k}, {2,i,j,k+1} };
			for( int n=0; n<6; n++ ) {
                if( q[n][0] < 0 || q[n][0] > gn-1 || q[n][1] < 0 || q[n][1] > gn-1 || q[n][2] < 0 || q[n][2] > gn-1 ) continue;
				FLOAT64 phi1 = cellfluidLS.clampFetch(q[n][0],q[n][1],q[n][2]);		// Fluid levelset outside
				FLOAT64 rho = fmax(eps,fraction(phi0,phi1));				// Fluid density
				if( surf_order==1 ) rho = 1.0;
				FLOAT64 frac = area[f[n][0]][f[n][1]][f[n][2]][f[n][3]];			// Cell face fraction
				FLOAT64 scale = frac*dt/dx2;
				FLOAT64 solidLS[2] = {	solid->evalLevelset(vec3d(dx*(i+0.5),dx*(j+0.5),dx*(k+0.5))),
										solid->evalLevelset(vec3d(dx*(q[n][0]+0.5),dx*(q[n][1]+0.5),dx*(q[n][2]+0.5))) };
				// If not use variational method, set fraction 0 or 1 on cell faces
				if( ! variation ) frac = (solidLS[0] < 0 || solidLS[1] < 0 ) ? 0.0 : 1.0;
				if( phi1 <= 0 ) {											// If adjacent to fluid cell
                    matrix.add_to_element(X(i,j,k),X(i,j,k),scale);
                    matrix.add_to_element(X(i,j,k),X(q[n][0],q[n][1],q[n][2]),-scale);
				} else {													// If adjacent to empty cell
                    matrix.add_to_element(X(i,j,k),X(i,j,k),scale/rho);
				}
				FLOAT64 sgn = ( i < q[n][0] || j < q[n][1] || k < q[n][2] ) ? -1.0 : 1.0;
                rhs[X(i,j,k)] += sgn*frac*u[f[n][0]][f[n][1]][f[n][2]][f[n][3]] / dx;  // Compute negative divergence
            }
		}
	} END_FOR
	
	// Copy divergence in 3D array for visualization
	div.resize(gn,gn,gn);
    FOR_EACH(gn,gn,gn) {
        div[i][j][k] = rhs[X(i,j,k)];
    } END_FOR
	
	dump( "Done. Took %s.\n", stock("MAC_matrix"));
	tick(); dump( "Organizing matrix (compressing matrix)..." );
	
	// Compact the matrix
	vector<int> idx_map;
	compress( matrix, rhs, idx_map );
	dump( "Done. Took %s.\n", stock());
	
#if 0
	ofstream file;
	file.open("Lhs_matrix.m");
	matrix.write_matlab(file,matrix.index.size(),matrix.index.size(),"Lhs");
	file.close();
#endif
	
	// Set initial guess
    FOR_EACH(gn,gn,gn) {
		uint idx = X(i,j,k);
        if( idx_map[idx] >= 0 ) pressure[idx_map[idx]] = p[i][j][k];
    } END_FOR
	
	// Solve by incomplete cholesky factorizaion preconditioned conjugate gradient method
    PCGSolver<FLOAT64> solver;
    FLOAT64 residual_out, msec;
    int iterations;
	dump( "Solving pressure by MAC..." );
	bool converged = solver.solve( matrix, rhs, pressure, residual_out, iterations, msec );
	if( residual_out > 1.0 ) {
		dump( "Too large Residual... ! Residual = %e", residual_out );
		converged = false;
	}
	if( ! converged ) {
		email::print( "PCG did not converge. Residual = %e\n", residual_out );
		email::send();
	} else {
        dump("Done. Took %d iterations with %d unknowns. Took %s.\n", iterations, matrix.index.size(), tstr(msec));
		writeNumber("MAC_PCG_iterations",iterations);
		writeNumber("MAC_PCG_time",msec);
		writeNumber("MAC_PCG_unknowns", matrix.index.size());
		writeNumber("MAC_PCG_residual",residual_out);
    }
    
    // Copy the result
    FOR_EACH(gn,gn,gn) {
		uint idx = X(i,j,k);
        p[i][j][k] = idx_map[idx] >= 0 ? pressure[idx_map[idx]] : 0.0;
    } END_FOR
	
	tick(); dump( "Computing gradient..." );
	
	// Subtract pressure gradient ( We must take into account ghost pressures here... )
	for( int dim=0; dim<DIM; dim++ ) {
		PARALLEL_FOR FOR_EACH(gn+(dim==0),gn+(dim==1),gn+(dim==2)) {
			int ni = i-(dim==0);
			int nj = j-(dim==1);
			int nk = k-(dim==2);
			// Gather solid levelset info
			FLOAT64 solidLS[2] = {	solid->evalLevelset(vec3d(dx*(i+0.5),dx*(j+0.5),dx*(k+0.5))),
									solid->evalLevelset(vec3d(dx*(ni+0.5),dx*(nj+0.5),dx*(nk+0.5))) };
			// If not use variational method, do not subtract pressure gradient on solid faces, but instead set zero
			if( ! variation ) {
				if( solidLS[0] * solidLS[1] < 0 ) {
					if( wallSepration )
						u[dim][i][j][k] = (solid->evalGradient(vec3d(dx*(i+0.5),dx*(j+0.5),dx*(k+0.5)))[dim] * u[dim][i][j][k] > 0.0 ) * u[dim][i][j][k];
					else 
						u[dim][i][j][k] = 0.0;
					if( solidLS[0] < 0 ) p.safeSet(i,j,k,p[ni][nj][nk]);
					else p.safeSet(ni,nj,nk,p.clampFetch(i,j,k));
					// Nothing to do for this face, so skip to next
					continue;
				}
			}
			// Fetch forward and backward pressures
			FLOAT64 press_f = p.clampFetch(i,j,k);
			FLOAT64 press_b = p.clampFetch(ni,nj,nk);
			
			// Fetch liquid levelset here
			FLOAT64 fluidLS[2] = { cellfluidLS.clampFetch(i,j,k), cellfluidLS.clampFetch(ni,nj,nk) };
			
			// Compute density
			FLOAT64 rho = fmax(eps,fraction(fluidLS[0],fluidLS[1]));
			if( surf_order==1 ) rho = 1.0;
			
			// Instant wall separating condition
			if( wall_separate ) {
				if( solidLS[0] < 0.5*dx ) press_f = fmax(0.0,press_f);
				if( solidLS[1] < 0.5*dx ) press_b = fmax(0.0,press_b);
			}
			
			// Now compute the pressure gradient
			gradp[dim][i][j][k] = area[dim][i][j][k] ? (press_f-press_b)/(dx*rho) : 0.0;
		} END_FOR
	}
	dump( "Done. Took %s.\n", stock("MAC_project"));
}

static void drawBitmapString( const char *string, void *font=NULL ) {
	if( ! font ) font = GLUT_BITMAP_HELVETICA_10;
	while (*string) glutBitmapCharacter(font, *string++);
}

void macfluid3::render( int name ) const {
	switch(name) {
		case MESH:
			break;
		case NAME: {
			glColor4f(1.0,1.0,1.0,1.0);
			glRasterPos2d(0.07,0.92);
			drawBitmapString(getName(),GLUT_BITMAP_HELVETICA_18);
			break;
		}
	}
}
