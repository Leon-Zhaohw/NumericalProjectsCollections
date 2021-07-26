/*
 *	octfluid3.cpp
 *
 *	Created by Ryoichi Ando on 6/29/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "octfluid3.h"
#include "levelset3.h"
#include "fastmarch3.h"
#include "opengl.h"
#include "util3.h"
#include "matutil.h"
#include "pcgsolver/pcg_solver.h"
#include <math.h>
#include <queue>
#include <fstream>

octfluid3::octfluid3() {
	maxdepth = 2;
	dt = 0.1;
	extrapolate_dist = 1.0;
	surf_order = 2;
	variation = true;
	tension = 0.0;
	solid = NULL;
	fluid = NULL;
	free_e = 1e-18;
}

uint64 octfluid3::computeIndex( vec3i p ) {
	for( uint dim=0; dim<DIM; dim++ ) {
		if( p[dim] < 0 || p[dim] > octree.resolution ) return 0;
	}
	uint64 R = octree.resolution;
	return p[0]+p[1]*R+p[2]*R*R+1;
}

static void add_to_element( std::vector<uint> &indices, std::vector<FLOAT64> &values, uint index, FLOAT64 value ) {
	for( uint n=0; n<indices.size(); n++ ) {
		if( indices[n] == index ) {
			values[n] += value;
			return;
		}
	}
	indices.push_back(index);
	values.push_back(value);
}

void octfluid3::init( uint gn, const levelset3 *hint ) {
	// Build octree
	uint maxdepth = log2(gn)+1;
	octree.buildOctree(hint,maxdepth);

	// Create terminal cells
	for( uint n=0; n<cells.size(); n++ ) delete cells[n];
	cells.clear();
	cells.resize(octree.terminals.size());
	for( uint n=0; n<octree.terminals.size(); n++ ) {
		const octree3::leaf3 &leaf = *octree.terminals[n];
		cells[n] = new cell3;
		cells[n]->dx = leaf.dx;
		cells[n]->p = leaf.center;
		cells[n]->pressure = 0.0;
		cells[n]->divergence = 0.0;
		cells[n]->fluidLS = 0.0;
		cells[n]->index = n;
	}
	
	// Create possible facets
	std::map<uint64,facet3 *> facetDictionary;
	for( uint n=0; n<octree.terminals.size(); n++ ) {
		for( uint dim=0; dim<DIM; dim++ ) {
			for( int dir=-1; dir<=1; dir+=2 ) {
				const octree3::leaf3 &leaf = *octree.terminals[n];
				int dx = leaf.dx;
				vec3i fcenter = leaf.center+0.5*dx*vec3i(dir*(dim==0),dir*(dim==1),dir*(dim==2));
				uint64 idx = computeIndex(fcenter);
				// If this turns to be a new facet, just append to the dictionary
				facet3 *facet;
				if( idx && facetDictionary.find(idx) == facetDictionary.end() ) {
					facet = new facet3;
					facet->p = fcenter;
					facet->dx = leaf.dx;
					facet->dim = dim;
					facet->velocity = 0.0;
					facet->gradient = 0.0;
					facet->fraction = 1.0;
					facet->mapped = false;
					facet->Lcell = NULL;
					facet->ds = 0.0;
					facet->index = 0;
					bool edge = false;
					if( facet->p[dim] == 0 || facet->p[dim] == octree.resolution ) edge = true;
					facet->edge = edge;
					facetDictionary[idx] = facet;
				}
			}
		}
	}
	
	// Connect facets and cells for facets
	for( uint n=0; n<cells.size(); n++ ) {
		cell3 &cell = *cells[n];
		// For all center directions...
		FOR_EACH(3,3,3) {
			int dx = cell.dx;
			vec3i fcenter = cell.p-0.5*dx*vec3i(1,1,1)+0.5*dx*vec3i(i,j,k);
			uint64 idx = computeIndex(fcenter);
			if( facetDictionary.find(idx) != facetDictionary.end() ) {
				facet3 &facet = *facetDictionary[idx];
				facet.cells.push_back(&cell);
			}
		} END_FOR
	}
			
	// Remove some facets
	std::map<uint64,facet3 *> newFacetDictionary;
	for( std::map<uint64,facet3 *>::iterator it=facetDictionary.begin(); it!=facetDictionary.end(); it++ ) {
		facet3 &facet = *(*it).second;
		if( facet.cells.size() <= 1 && facet.edge == false ) {
			delete (*it).second;
		} else {
			newFacetDictionary[(*it).first] = (*it).second;
		}
	}
	facetDictionary = newFacetDictionary;
	
	// Connect facet and cells for cells
	for( uint n=0; n<cells.size(); n++ ) {
		cell3 &cell = *cells[n];
		// For all center directions...
		FOR_EACH(3,3,3) {
			int dx = cell.dx;
			vec3i fcenter = cell.p-0.5*dx*vec3i(1,1,1)+0.5*dx*vec3i(i,j,k);
			uint64 idx = computeIndex(fcenter);
			if( facetDictionary.find(idx) != facetDictionary.end() ) {
				facet3 *facet = facetDictionary[idx];
				cell.facets.push_back(facet);
			}
		} END_FOR
	}
	
	// Pack facets into the vector
	for( uint n=0; n<facets.size(); n++ ) {
		delete facets[n];
	}
	facets.resize(facetDictionary.size());
	uint index = 0;
	for( std::map<uint64,facet3 *>::iterator it=facetDictionary.begin(); it!=facetDictionary.end(); it++ ) {
		facet3 *facet = (*it).second;
		facet->index = index;
		facets[index++] = facet;
	}
	
	// Build nodes
	for( std::map<uint64,node3 *>::iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		node3 *node = (*it).second;
		delete node;
	}
	nodes.clear();
	for( uint n=0; n<cells.size(); n++ ) {
		cell3 *cell = cells[n];
		int dx = cell->dx;
		vec3i corner = cell->p-0.5*dx*vec3i(1,1,1);
		FOR_EACH(2,2,2) {
			vec3i p = corner+dx*vec3i(i,j,k);
			uint idx = computeIndex(p);
			if( idx ) {
				node3 *node = NULL;
				if( nodes.find(idx) == nodes.end() ) {
					node = new node3();
					node->p = p;
					node->dx = cell->dx;
					node->fluidLS = 0.0;
					node->solidLS = 0.0;
					node->mapped = false;
					for( uint dim=0; dim<DIM; dim++ ) {
						node->weights[dim].clear();
						node->samples[dim].clear();
						if( node->p[dim] == 0 || node->p[dim] == octree.resolution ) node->edge[dim] = true;
						else node->edge[dim] = false;
					}
					nodes[idx] = node;
				} else {
					node = nodes[idx];
					node->dx = min(node->dx,cell->dx);
				}
				cell->nodes[i][j][k] = node;
			}
		} END_FOR
	}
	
	// Scatter samples beforehand
	for( uint dim=0; dim<DIM; dim++ ) {
		for( std::map<uint64,sample3 *>::iterator it=samples[dim].begin(); it!=samples[dim].end(); it++ ) {
			sample3 *sample = (*it).second;
			delete sample;
		}
		samples[dim].clear();
	}
	for( uint n=0; n<facets.size(); n++ ) {
		facet3 *facet = facets[n];
		sample3 *sample = new sample3();
		sample->p = facet->p;
		sample->facets[0] = facet;
		sample->weights[0] = 1.0;
		sample->num = 1;
		uint idx = computeIndex(facet->p);
		samples[facet->dim][idx] = sample;
	}
	for( uint dim=0; dim<DIM; dim++ ) {
		for( uint n=0; n<cells.size(); n++ ) {
			cell3 *cell = cells[n];
			int dx = cell->dx;
			bool findJunction = false;
			for( int dir=-1; dir<=1; dir+=2 ) {
				vec3i p = cell->p+0.5*dx*dir*vec3i(dim==0,dim==1,dim==2);
				uint idx = computeIndex(p);
				if( facetDictionary.find(idx) != facetDictionary.end() ) {
					findJunction = true;
					break;
				}
			}
			if( findJunction ) {
				for( uint dim2=0; dim2<DIM; dim2++ ) {
					if( dim != dim2 ) {
						uint idx = computeIndex(cell->p);
						if( samples[dim2].find(idx) == samples[dim2].end() ) {
							sample3 *sample = new sample3();
							sample->p = cell->p;
							int count = 0;
							for( int dir=-1; dir<=1; dir+=2 ) {
								vec3i p = cell->p+0.5*dx*dir*vec3i(dim2==0,dim2==1,dim2==2);
								uint idx = computeIndex(p);
								if( facetDictionary.find(idx) != facetDictionary.end() && facetDictionary[idx]->dim == dim2 ) {
									facet3 *facet = facetDictionary[idx];
									sample->facets[count] = facet;
									sample->weights[count] = 1.0;
									count ++;
								}
							}
							sample->num = count;
							if( count ) {
								for( uint m=0; m<count; m++ ) {
									sample->weights[m] /= (FLOAT64)count;
								}
							}
							samples[dim2][idx] = sample;
						}
					}
				}
			}
		}
	}
	
	// Connect nodes and samples
	for( std::map<uint64,node3 *>::iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		node3 *node = (*it).second;
		int dx = node->dx;
		for( uint dim=0; dim<DIM; dim++ ) {
			// If happens to be on the T-junctoin, facet may be found just there
			vec3i p = node->p;
			uint idx = computeIndex(p);
			if( samples[dim].find(idx) != samples[dim].end() ) {
				node->samples[dim].push_back(samples[dim][idx]);
			} else {
				// Otherwise probe for other directions
				vec3i vec0, vec1;
				getTangentVectors(dim,vec0,vec1);
				vec3i dirv[] = { vec0+vec1, -1*vec0+vec1, -1*vec0-1*vec1,vec0-1*vec1, vec0, vec1, -1*vec0, -1*vec1 };
				FLOAT64 arms[] = { 0.5, 1.0, 1.5 };
				// Try find corner ones
				for( uint i=0; i<3; i++ ) {
					if( i == 2 && node->samples[dim].size() ) break;
					for( int dir=0; dir<8; dir++ ) {
						vec3i p = node->p+dx*arms[i]*dirv[dir];	
						uint idx = computeIndex(p);
						if( samples[dim].find(idx) != samples[dim].end() ) {
							sample3 *sample = samples[dim][idx];
							node->samples[dim].push_back(sample);
						}
					}
				}
			}
			// Compute weights for each samples
			if( ! computeWeights( node->samples[dim], node->weights[dim], dim, node->p )) {
				writeConnections("octfluid.obj",dim,node);
				exit(0);
			}
		}
	}
	
	// Connect facets and nodes
	for( uint n=0; n<facets.size(); n++ ) {
		facet3 *facet = facets[n];
		uint dim = facet->dim;
		int dx = facet->dx;
		// On T-junction, node can exist just under the facet
		uint idx = computeIndex(facet->p);
		if( nodes.find(idx) != nodes.end() ) {
			facet->nodes.push_back(nodes[idx]);
		} else {
			// If not found, just search for both directions
			vec3i vec0, vec1;
			getTangentVectors(dim,vec0,vec1);
			vec3i dirv[] = { vec0+vec1, -1*vec0+vec1, -1*vec0-1*vec1,vec0-1*vec1 };
			for( int dir=0; dir<4; dir++ ) {
				vec3i p = facet->p+dx*0.5*dirv[dir];
				uint idx = computeIndex(p);
				if( nodes.find(idx) != nodes.end() ) {
					facet->nodes.push_back(nodes[idx]);
				}
			}
		}
	}
	
	// Connect node and node
	for( std::map<uint64,node3 *>::iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		node3 *node = (*it).second;
		int dx = 1.1*node->dx;
		for( uint dim=0; dim<DIM; dim++ ) {
			for( int dir=-1; dir<=1; dir+=2 ) {
				for( int arm=1; arm<=2; arm++ ) {
					vec3i p = node->p+dir*arm*dx*vec3i(dim==0,dim==1,dim==2);
					uint idx = computeIndex(p);
					if( nodes.find(idx) != nodes.end() ) {
						node->nodes.push_back(nodes[idx]);
						break;
					}
				}
			}
		}
	}
	
	// Clear matrix (against viewer functions)
	LhsMatrix.clear();
	RhsMatrix.clear();
}

void octfluid3::setParameter( int name, FLOAT64 value ) {
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

void octfluid3::setTimestep( FLOAT64 dt ) {
	this->dt = dt;
}

static FLOAT64 fraction( FLOAT64 phi0, FLOAT64 phi1 ) {
	return phi0*phi1 > 0 ? 1 - (phi0 > 0) : -fmin(phi0,phi1)/fmax(fabs(phi1-phi0),1.0e-6);
}

void octfluid3::setupSolidLevelset( const levelset3 *solid ) {
	this->solid = solid;
	for( std::map<uint64,node3 *>::iterator it=this->nodes.begin(); it!=this->nodes.end(); it++ ) {
		node3 *node = (*it).second;
		vec3d p = vec3d(node->p)/octree.resolution;
		node->solidLS = solid->evalLevelset(p);
	}
	// Compute fractions of facets
	for( uint n=0; n<facets.size(); n++ ) {
		facet3 *facet = facets[n];
		vec3d p = vec3d(facet->p)/(FLOAT64)octree.resolution;
		FLOAT64 dx = facet->dx/(FLOAT64)octree.resolution;
		uint dim = facet->dim;
		vec3i vec0, vec1;
		getTangentVectors(dim,vec0,vec1);
		std::vector<vec3d> q(4);
		q[0] = p-0.5*dx*vec3d(vec0)-0.5*dx*vec3d(vec1);
		q[1] = p+0.5*dx*vec3d(vec0)-0.5*dx*vec3d(vec1);
		q[2] = p+0.5*dx*vec3d(vec0)+0.5*dx*vec3d(vec1);
		q[3] = p-0.5*dx*vec3d(vec0)+0.5*dx*vec3d(vec1);
		std::vector<vec3d> lq(4);
		lq[0] = vec3d(-0.5,-0.5,0.0);
		lq[1] = vec3d(+0.5,-0.5,0.0);
		lq[2] = vec3d(+0.5,+0.5,0.0);
		lq[3] = vec3d(-0.5,+0.5,0.0);
		
		std::vector<FLOAT64> ls(4);
		for( uint i=0; i<4; i++ ) {
			ls[i] = solid->evalLevelset(q[i]);
		}
		FLOAT64 area = util3::compute2DArea(util3::march2DPoints(lq,ls));
		facet->fraction = 1.0-area;
		if( ! variation ) facet->fraction = ceil(facet->fraction);
	}
}

void octfluid3::setupFluidLevelset( const levelset3 *fluid ) {
	for( std::map<uint64,node3 *>::iterator it=this->nodes.begin(); it!=this->nodes.end(); it++ ) {
		node3 *node = (*it).second;
		vec3d p = vec3d(node->p)/octree.resolution;
		node->fluidLS = fluid->evalLevelset(p);
	}
	for( uint n=0; n<cells.size(); n++ ) {
		cells[n]->fluidLS = fluid->evalLevelset(vec3d(cells[n]->p)/octree.resolution);
	}
}

void octfluid3::getSamplePoints( std::vector<vec3d> &pos ) {
	// Gather facet position
	pos.resize(nodes.size());
	uint n=0;
	for( std::map<uint64,node3 *>::const_iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		const node3 *node = (*it).second;
		pos[n++] = vec3d(node->p)/octree.resolution;
	}
}

void octfluid3::setupVelocity( const std::vector<vec3d> &vel, const std::vector<bool> &mapped ) {
	// Set whether node is mapped
	uint index=0;
	for( std::map<uint64,node3 *>::iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		node3 *node = (*it).second;
		node->velocity = vel[index];
		for( uint dim=0; dim<DIM; dim++ ) {
			if( node->edge[dim] ) node->velocity[dim] = 0.0;
		}
		node->mapped = mapped[index];
		index ++;
	}
}

void octfluid3::resample( int kind, int type ) { // type 0: node -> facet, 1: facet -> node
	if( type == 0 ) {
		// Resample node -> faccet
		for( uint n=0; n<facets.size(); n++ ) {
			facet3 *facet = facets[n];
			facet->mapped = false;
			uint dim = facet->dim;
			FLOAT64 value=0.0;
			uint sum=0;
			for( std::list<node3 *>::iterator it=facet->nodes.begin(); it!=facet->nodes.end(); it++ ) {
				node3 *node = *it;
				if( node->mapped ) {
					if( kind == 0 ) value += node->velocity[dim];
					else value += node->gradient[dim];
					sum++;
				}
			}
			if( sum ) {
				if( kind == 0 ) facet->velocity = value / sum;
				else facet->gradient = value / sum;
				facet->mapped = true;
			}
		}
	} else {
		// Resample facet -> node
		for( std::map<uint64,node3 *>::iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
			node3 *node = (*it).second;
			if( kind == 0 ) node->velocity = vec3d();
			else node->gradient = vec3d();
			bool mapped[DIM] = { false, false };
			for( uint dim=0; dim<DIM; dim++ ) {
				FLOAT64 wsum = 0.0;
				for( uint n=0; n<node->samples[dim].size(); n++ ) {
					sample3 *sample = node->samples[dim][n];
					FLOAT64 w = node->weights[dim][n];
					for( uint m=0; m<sample->num; m++ ) {
						facet3 *facet = sample->facets[m];
						if( facet->mapped ) {
							FLOAT64 w2 = w*sample->weights[m];
							if( kind == 0 ) node->velocity[dim] += w2*facet->velocity;
							else node->gradient[dim] += w2*facet->gradient;
							wsum += w2;
						}
					}
				}
				if( wsum ) {
					if( kind == 0 ) node->velocity[dim] /= wsum;
					else node->gradient[dim]  /= wsum;
					if( node->solidLS > 0.0 ) mapped[dim] = true;
				}
			}
			if( DIM == 2 ) node->mapped = mapped[0] && mapped[1];
			else if( DIM == 3 ) node->mapped = mapped[0] && mapped[1] && mapped[2];
		}
	}
}

void octfluid3::buildMatrix() {
	// Compute a facet distance "ds" and find an adjacent biggest cell
	for( uint n=0; n<facets.size(); n++ ) {
		facet3 *facet = facets[n];
		uint dim = facet->dim;
		// Find a cell with maximal cell area
		FLOAT64 max_area = 0.0;
		cell3 *Lcell = NULL;
		for( std::list<cell3 *>::iterator it=facet->cells.begin(); it!=facet->cells.end(); it++ ) {
			cell3 *cell = *it;
			FLOAT64 area = cell->dx;
			if( max_area < area ) {
				max_area = area;
				Lcell = cell;
			}
		}
		facet->Lcell = Lcell;
		
		// Compute average (signed) distance between small and large cells
		FLOAT64 ds=0.0;
		for( std::list<cell3 *>::iterator it=facet->cells.begin(); it!=facet->cells.end(); it++ ) {
			cell3 *Scell = *it;
			if( Scell != Lcell ) {
				FLOAT64 dist = fabs(Scell->p[dim]-Lcell->p[dim]);
				FLOAT64 area = Scell->dx;
				ds += area/max_area * dist;
			}
		}
		facet->ds = ds;
	}
	
	// Fill global matrix
	LhsMatrix.clear();
	LhsMatrix.resize(cells.size());
	RhsMatrix.clear();
	RhsMatrix.resize(cells.size());
	
	for( uint n=0; n<cells.size(); n++ ) {
		cell3 *cell = cells[n];
		if( cell->fluidLS > 0.0 ) continue;
		for( std::list<facet3 *>::iterator it=cell->facets.begin(); it!=cell->facets.end(); it++ ) {
			facet3 *facet = *it;	
			uint dim = facet->dim;
			int sgn = cell->p[dim] < facet->p[dim] ? 1 : -1;
			FLOAT64 ds = facet->ds;
			cell3 *Lcell = facet->Lcell;
			if( cell == Lcell ) {
				for( std::list<cell3 *>::iterator it2=facet->cells.begin(); it2!=facet->cells.end(); it2++ ) {
					cell3 *Scell = *it2;
					if( Scell == Lcell ) continue;		
					FLOAT64 frac = facet->fraction;
					if( ! variation ) frac = ceil(frac);
					if( frac ) {
						FLOAT64 area = Scell->dx;
						FLOAT64 rho = fraction(Scell->fluidLS,Lcell->fluidLS);
						if( surf_order == 1 ) rho = ceil(rho);
						if( rho ) {
							FLOAT64 value = frac*area/(ds*fmax(free_e,rho));
							if( Scell->fluidLS < 0.0 ) LhsMatrix.add_to_element(n,Scell->index,-value);	// p_small
							if( Lcell->fluidLS < 0.0 ) LhsMatrix.add_to_element(n,Lcell->index,value);	// p_large (diagonal term)
							RhsMatrix.add_to_element(n,facet->index,-sgn*frac*area);
						}
					}
				}
			} else {
				FLOAT64 frac = facet->fraction;
				if( ! variation ) frac = ceil(frac);
				if( frac ) {
					cell3 *Scell = cell;
					FLOAT64 area = Scell->dx;
					FLOAT64 rho = fraction(Scell->fluidLS,Lcell->fluidLS);
					if( surf_order == 1 ) rho = ceil(rho);
					if( rho ) {						
						FLOAT64 value = frac*area/(ds*fmax(free_e,rho));
						if( Scell->fluidLS < 0.0 ) LhsMatrix.add_to_element(n,Scell->index,value);	// p_small (diagonal term)
						if( Lcell->fluidLS < 0.0 ) LhsMatrix.add_to_element(n,Lcell->index,-value);	// p_large
						RhsMatrix.add_to_element(n,facet->index,-sgn*frac*area);
					}
				}
			}
		}
	}
}

void octfluid3::project() {	
	
	// Extrapolate
	tick(); dump( "Velocity extrapolation started..." );
	extrapolate(0);
	dump( "Done. Took %s.\n", stock("OCT_velocity_extrapolate"));
	
	// Resample to facets
	tick(); dump( "Resampling velocity node -> facet..." );
	resample(0,0);
	dump( "Done. Took %s.\n", stock("OCT_resample_velocity_facet"));
	
	// Build matrix
	tick(); dump( "Building matrix started..." );
	buildMatrix();
	dump( "Done. Took %s.\n", stock("OCT_matrix"));
	
	// Compute divergence
	tick(); dump( "Computing divergence started..." );
	std::vector<FLOAT64> rhs(cells.size());
	std::vector<FLOAT64> fvel(facets.size());
	for( uint n=0; n<facets.size(); n++ ) {
		fvel[n] = facets[n]->velocity;
	}
	multiply<FLOAT64>(RhsMatrix,fvel,rhs);
	dump( "Done. Took %s.\n", stock("OCT_divergence"));
	
	// Copy divergence
	for( uint n=0; n<cells.size(); n++ ) {
		cells[n]->divergence = rhs[n];
	}
	
	tick(); dump( "Organizing matrix started..." );
	SparseMatrix<FLOAT64> LhsMatrix = this->LhsMatrix;
	
	// Compact the matrix
	std::vector<int> idx_map;
	compress( LhsMatrix, rhs, idx_map );
	dump( "Done. Took %s.\n", stock("OCT_matrix_organize"));
	
#if 0
	std::ofstream file;
	file.open("Lhs_OCT_matrix.m");
	LhsMatrix.write_matlab(file,LhsMatrix.n,LhsMatrix.n,"Lhs_OCT");
	file.close();
	file.open("Rhs_OCT_vector.m");
	write_matlab1d(file,rhs,"Rhs_OCT");
	file.close();
	exit(0);
#endif
	
	// Solve by incomplete cholesky factorizaion preconditioned conjugate gradient method
	PCGSolver<FLOAT64> solver;
	FLOAT64 residual_out;
	std::vector<FLOAT64> result(LhsMatrix.n);
	for( uint n=0; n<cells.size(); n++ ) {
		if( idx_map[n] >= 0 ) result[idx_map[n]] = cells[n]->pressure;
	}
	
	int iterations;
	FLOAT64 msec;
	bool converged = solver.solve( LhsMatrix, rhs, result, residual_out, iterations, msec );
	if( ! converged ) {
		email::print("PCG did not converge. Residual = %e\n", residual_out );
		email::send();
		
		// Dump this matrix
		std::ofstream file;
		file.open(format_str("%s/Lhs_OCT_matrix.m",root_path));
		LhsMatrix.write_matlab(file,LhsMatrix.n,LhsMatrix.n,"Lhs_FEM");
		file.close();
		file.open(format_str("%s/Rhs_OCT_vector.m",root_path));
		write_matlab1d(file,rhs,"Rhs_FEM");
		file.close();
		exit(0);
	} else {
		dump( "OCT: PCG Converged ! %d iterations and took %.3f msec with %d unknowns.\n", iterations, msec, LhsMatrix.n);
		writeNumber("OCT_iterations",iterations);
		writeNumber("OCT_PCG_iterations",iterations);
		writeNumber("OCT_PCG_time",msec);
		writeNumber("OCT_PCG_unknowns", LhsMatrix.n);
		writeNumber("OCT_PCG_residual",residual_out);
	}
	
	// Set pressure
	tick(); dump( "Computing pressure gradient started..." );
	for( uint n=0; n<cells.size(); n++ ) {
		cells[n]->pressure = idx_map[n] >= 0 ? result[idx_map[n]] : 0.0;
	}
	
	// Compute gradient
	for( uint n=0; n<facets.size(); n++ ) {
		facet3 *facet = facets[n];
		uint dim = facet->dim;
		FLOAT64 ds = facet->ds;
		cell3 *Lcell = facet->Lcell;
		int sgn = Lcell->p[dim] < facet->p[dim] ? 1 : -1;
		FLOAT64 gradient = 0.0;
		facet->mapped = false;
		if( facet->fraction ) {
			for( std::list<cell3 *>::iterator it2=facet->cells.begin(); it2!=facet->cells.end(); it2++ ) {
				cell3 *Scell = *it2;
				if( Scell == Lcell ) continue;
				FLOAT64 rho = fraction(Scell->fluidLS,Lcell->fluidLS);
				if( surf_order == 1 ) rho = ceil(rho);
				if( rho ) {
					gradient += -sgn * (Scell->dx/(FLOAT64)Lcell->dx) * (Lcell->pressure - Scell->pressure) / (ds*fmax(free_e,rho));
					facet->mapped = true;
				}
			}
		} else {
			facet->mapped = false;
		}
		facet->gradient = gradient;
	}
	dump( "Done. Took %s.\n", stock("OCT_pressure_gradient"));
	
	// Resample to nodes
	tick(); dump( "Resampling gradients from facet -> node..." );
	resample(1,1);
	dump( "Done. Took %s.\n", stock("OCT_resample_facet_node"));
	
	// Extrapolate
	tick(); dump( "Gradient extrapolation started..." );
	extrapolate(1);
	dump( "Done. Took %s.\n", stock("OCT_gradient_extrapolate"));
}


vec3d octfluid3::getPressureGradient( vec3d p ) const {
	return sampleVector(p,1);
}

vec3d octfluid3::getVelocity( vec3d p ) const {
	return sampleVector(p,0);
}

FLOAT64 octfluid3::getDivergence( vec3d p ) const {
	int n = octree.hitTest(p);
	if( n >= 0 ) {
		return cells[n]->divergence;
	}
	return 0.0;
}

vec3d octfluid3::sampleVector( vec3d p, char kind ) const {
	int hit = octree.hitTest(p);
	if( hit >= 0 ) {
		// Do sample here...
		FLOAT64 dx = cells[hit]->dx/(FLOAT64)octree.resolution;
		FLOAT64 x = (p[0]-cells[hit]->nodes[0][0][0]->p[0]/(FLOAT64)octree.resolution)/dx;
		FLOAT64 y = (p[1]-cells[hit]->nodes[0][0][0]->p[1]/(FLOAT64)octree.resolution)/dx;
		FLOAT64 z = (p[2]-cells[hit]->nodes[0][0][0]->p[2]/(FLOAT64)octree.resolution)/dx;
		vec3d v[8];
		if( kind == 1 ) {
			FLOAT64 L[8];
			L[0] = cells[hit]->nodes[0][0][0]->fluidLS;
			L[1] = cells[hit]->nodes[1][0][0]->fluidLS;
			L[2] = cells[hit]->nodes[1][1][0]->fluidLS;
			L[3] = cells[hit]->nodes[0][1][0]->fluidLS;
			L[4] = cells[hit]->nodes[0][0][1]->fluidLS;
			L[5] = cells[hit]->nodes[1][0][1]->fluidLS;
			L[6] = cells[hit]->nodes[1][1][1]->fluidLS;
			L[7] = cells[hit]->nodes[0][1][1]->fluidLS;
			FLOAT64 fluidLS = (1.-z)*(L[0]*(1.-x)*(1.-y)+L[1]*x*(1.-y)+L[2]*x*y+L[3]*(1.-x)*y)+
								z*(L[4]*(1.-x)*(1.-y)+L[5]*x*(1.-y)+L[7]*x*y+L[7]*(1.-x)*y);
			if( fluidLS > 0.0 ) kind = 0;
		}
		if( kind == 0 ) {
			v[0] = cells[hit]->nodes[0][0][0]->velocity;
			v[1] = cells[hit]->nodes[1][0][0]->velocity;
			v[2] = cells[hit]->nodes[1][1][0]->velocity;
			v[3] = cells[hit]->nodes[0][1][0]->velocity;
			v[4] = cells[hit]->nodes[0][0][1]->velocity;
			v[5] = cells[hit]->nodes[1][0][1]->velocity;
			v[6] = cells[hit]->nodes[1][1][1]->velocity;
			v[7] = cells[hit]->nodes[0][1][1]->velocity;
		} else {
			v[0] = cells[hit]->nodes[0][0][0]->gradient/dt;
			v[1] = cells[hit]->nodes[1][0][0]->gradient/dt;
			v[2] = cells[hit]->nodes[1][1][0]->gradient/dt;
			v[3] = cells[hit]->nodes[0][1][0]->gradient/dt;
			v[4] = cells[hit]->nodes[0][0][1]->gradient/dt;
			v[5] = cells[hit]->nodes[1][0][1]->gradient/dt;
			v[6] = cells[hit]->nodes[1][1][1]->gradient/dt;
			v[7] = cells[hit]->nodes[0][1][1]->gradient/dt;
		}
		return (1.-z)*(v[0]*(1.-x)*(1.-y)+v[1]*x*(1.-y)+v[2]*x*y+v[3]*(1.-x)*y)+z*(v[4]*(1.-x)*(1.-y)+v[5]*x*(1.-y)+v[6]*x*y+v[7]*(1.-x)*y);
	}
	return vec3d();
}

void octfluid3::extrapolate( char kind ) {
	// Assign index
	uint n=0;
	for( std::map<uint64,node3 *>::iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		node3 *node = (*it).second;
		node->index = n++;
	}
	std::vector<fastmarch3<vec3d>::node3 *> fnodes(nodes.size());
	for( uint n=0; n<nodes.size(); n++ ) fnodes[n] = new fastmarch3<vec3d>::node3;
	for( std::map<uint64,node3 *>::iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		node3 *node = (*it).second;
		uint n = node->index;
		vec3d p = node->p;
		bool fixed = node->mapped;
		fnodes[n]->p = p;
		fnodes[n]->fixed = fixed;
		fnodes[n]->levelset = fmax(node->fluidLS,-node->solidLS);
		fnodes[n]->value = (kind==0) ? node->velocity : node->gradient;
		fnodes[n]->p2p.resize(node->nodes.size());
		uint m = 0;
		for( std::list<node3 *>::iterator it=node->nodes.begin(); it!=node->nodes.end(); it++ ) {
			fnodes[n]->p2p[m++] = fnodes[(*it)->index];
		}
	}
	fastmarch3<vec3d>::fastMarch(fnodes,extrapolate_dist,-1.0,0);
	// Pick up extrapolated values
	for( std::map<uint64,node3 *>::iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		node3 *node = (*it).second;
		uint n = node->index;
		if( fnodes[n]->fixed ) {
			if (kind==0) node->velocity = fnodes[n]->value;
			else node->gradient = fnodes[n]->value;
		}
		delete fnodes[n];
	}
}

bool octfluid3::computeWeights( const std::vector<sample3 *> &samples, std::vector<FLOAT64> &weights, int dim, vec3i p ) {
	uint num = samples.size();
	weights.resize(num);
	for( uint n=0; n<num; n++ ) weights[n] = 0.0;
	if( num == 0 ) {
		vec3d pos = vec3d(p)/octree.resolution;
		printf( "Unexpected size ! dim = %d, num = %d, p=(%f,%f,%f)\n", dim, num, pos[0], pos[1], pos[2] );
		return false;
	} else if( num == 1 ) {
		weights[0] = 1.0;
	} else if( num == 2 ) {
		FLOAT64 dist0 = (p-samples[0]->p).len();
		FLOAT64 dist1 = (p-samples[1]->p).len();
		weights[0] = dist1/(dist0+dist1);
		weights[1] = dist0/(dist0+dist1);
	} else {
		// Collect low dimension points
		std::vector<vec3i> pos;
		vec3i center;
		pos.resize(num);
		uint slot=0;
		for( uint d=0; d<DIM; d++ ) {
			if( d != dim ) {
				for( uint n=0; n<num; n++ ) {
					pos[n][slot] = samples[n]->p[d];
				}
				center[slot] = p[d];
				slot++;
			}
		}
		if( num == 3 ) {
			FLOAT64 A[3][3];
			for( uint i=0; i<3; i++ ) for( uint j=0; j<3; j++ ) {
				if( i < 2 ) A[i][j] = pos[j][i];
				else A[i][j] = 1.0;
			}
			FLOAT64 C[3][3];
			if( invert3x3(A,C) ) {
				FLOAT64 x[3] = { center[0], center[1], 1.0 };
				for( uint i=0; i<3; i++ ) {
					weights[i] = 0.0;
					for( uint j=0; j<3; j++ ) weights[i] += C[i][j]*x[j];
				}
			} else {
				printf( "Failed to invert3x3 when computing weights.\n" );
				exit(0);
			}
		} else if( num == 4 ) {
			FLOAT64 A[4][4];
			for( uint i=0; i<4; i++ ) for( uint j=0; j<4; j++ ) {
				if(i<2) A[i][j] = pos[j][i];
				else if( i==2 ) A[i][j] = pos[j][0]*pos[j][1];
				else A[i][j] = 1.0;
			}
			FLOAT64 C[4][4];
			if( invert4x4(A,C) ) {
				FLOAT64 x[4] = { center[0], center[1], center[0]*center[1], 1.0 };
				for( uint i=0; i<4; i++ ) {
					weights[i] = 0.0;
					for( uint j=0; j<4; j++ ) weights[i] += C[i][j]*x[j];
				}
			} else {
				printf( "Failed to invert4x4 when computing weights.\n" );
				return false;
			}
		} else {
			for( uint n=0; n<num; n++ ) weights[n] = 1.0/num;
		}
	}
	return true;
}

void octfluid3::getTangentVectors( int dim, vec3i &vec0, vec3i &vec1 ) {
	switch( dim ) {
		case 0:
			vec0 = vec3i(0,1,0);
			vec1 = vec3i(0,0,1);
			break;
		case 1:
			vec0 = vec3i(1,0,0);
			vec1 = vec3i(0,0,1);
			break;
		case 2:
			vec0 = vec3i(1,0,0);
			vec1 = vec3i(0,1,0);
			break;
	}
}

static void drawBitmapString( const char *string, void *font=NULL ) {
	if( ! font ) font = GLUT_BITMAP_HELVETICA_10;
	while (*string) glutBitmapCharacter(font, *string++);
}

void octfluid3::render( int name ) const {
	switch(name) {
		case MESH:
			//glColor4d(1.0,1.0,1.0,0.05);
			//octree.drawOctree();
			break;
		case NAME: {
			glColor4f(1.0,1.0,1.0,1.0);
			glRasterPos2d(0.07,0.92);
			drawBitmapString(getName(),GLUT_BITMAP_HELVETICA_18);
			break;
		}
	}
}

void octfluid3::writeConnections( const char *path, uint dim, node3 *target_node ) {
	FILE *fp = fopen(path,"w");
	if( fp ) {
		FLOAT64 scale = 1.0/octree.resolution;
	
		// Write facets
		uint counter = 0;
		for( uint n=0; n<facets.size(); n++ ) {
			facet3 *facet = facets[n];
			if( (facet->p-target_node->p).len() < target_node->dx*5 && facet->p[dim] == target_node->p[dim] ) {
				vec3i vec0, vec1;
				getTangentVectors(facet->dim, vec0, vec1);
				int dx = facet->dx;
				int q[][2] = { {-1,-1}, {1,-1}, {1,1}, {-1,1} };
				for( uint m=0; m<4; m++ ) {
					vec3i corner = facet->p+0.5*dx*q[m][0]*vec0+0.5*dx*q[m][1]*vec1;
					uint64 idx = computeIndex(corner);
					if( nodes.find(idx) == nodes.end() ) {
						printf( "Unexpected facet encountered.\n" );
					} else {
						fprintf(fp, "v %f %f %f\n", scale*corner[0], scale*corner[1], scale*corner[2] );
						counter++;
					}
				}
				fprintf( fp, "f %d %d %d %d\n", counter-3, counter-2, counter-1, counter );
			}
		}
		
		// Write sample points
		for( std::map<uint64,sample3 *>::iterator it=samples[dim].begin(); it!=samples[dim].end(); it++ ) {
			sample3 *sample = (*it).second;
			if( (sample->p-target_node->p).len() < target_node->dx*5 ) {
				vec3i vec0, vec1;
				getTangentVectors(dim,vec0,vec1);
				int q[][2] = { {-1,-1}, {1,-1}, {1,1}, {-1,1} };
				for( uint m=0; m<4; m++ ) {
					vec3i corner = sample->p+1*q[m][0]*vec0+1*q[m][1]*vec1;
					fprintf(fp, "v %f %f %f\n", scale*corner[0], scale*corner[1], scale*corner[2] );
					counter++;
				}
				fprintf( fp, "f %d %d %d %d\n", counter-3, counter-2, counter-1, counter );
			}
		}
		
		// Write node points
		if( target_node ) {
			vec3i vec0, vec1;
			getTangentVectors(dim,vec0,vec1);
			int q[][2] = { {-1,-1}, {1,-1}, {1,1}, {-1,1} };
			for( uint m=0; m<4; m++ ) {
				vec3d corner = target_node->p+q[m][0]*vec0+q[m][1]*vec1;
				fprintf(fp, "v %f %f %f\n", scale*corner[0], scale*corner[1], scale*corner[2] );
			}
			
			vec3i dirv[] = { vec0+vec1, -1*vec0+vec1, -1*vec0-1*vec1,vec0-1*vec1, vec0, vec1, -1*vec0, -1*vec1 };
			for( int dir=0; dir<8; dir++ ) {
				FLOAT64 arms[] = { 0.5, 1.0, 1.5 };
				// Try find corner ones
				for( uint i=0; i<3; i++ ) {
					int dx = target_node->dx;
					vec3i p = target_node->p+dx*arms[i]*dirv[dir];	
					fprintf(fp, "v %f %f %f\n", scale*p[0], scale*p[1], scale*p[2] );
					uint idx = computeIndex(p);
					if( samples[dim].find(idx) != samples[dim].end() ) {
						printf( "Found out 1 !!\n" );
						break;
					}
				}
			}
		}
		fclose(fp);
	}
}
