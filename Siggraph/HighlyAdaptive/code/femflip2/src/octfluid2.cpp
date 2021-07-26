/*
 *	octfluid2.cpp
 *
 *	Created by Ryoichi Ando on 6/9/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "octfluid2.h"
#include "levelset2.h"
#include "fastmarch2.h"
#include "opengl.h"
#include "util2.h"
#include "matutil.h"
#include "pcgsolver/pcg_solver.h"
#include <math.h>
#include <queue>

octfluid2::octfluid2() {
	maxdepth = 2;
	dt = 0.1;
	extrapolate_dist = 1.0;
	surf_order = 2;
	variation = true;
	solid = NULL;
	fluid = NULL;
	free_e = 1e-18;
}

uint64 octfluid2::computeIndex( vec2i p ) {
	for( uint dim=0; dim<DIM; dim++ ) {
		if( p[dim] < 0 || p[dim] > octree.resolution ) return 0;
	}
	uint64 R = octree.resolution;
	return p[0]+p[1]*R+1;
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

void octfluid2::init( uint gn, const levelset2 *hint ) {
	// Build octree
	uint maxdepth = log2(gn)+1;
	octree.buildOctree(hint,maxdepth);
	
	// Create terminal cells
	for( uint n=0; n<cells.size(); n++ ) delete cells[n];
	cells.clear();
	cells.resize(octree.terminals.size());
	for( uint n=0; n<octree.terminals.size(); n++ ) {
		const octree2::leaf2 &leaf = *octree.terminals[n];
		cells[n] = new cell2;
		cells[n]->dx = leaf.dx;
		cells[n]->p = leaf.center;
		cells[n]->pressure = 0.0;
		cells[n]->divergence = 0.0;
		cells[n]->fluidLS = 0.0;
		cells[n]->index = n;
	}
	
	// Create possible facets
	std::map<uint64,facet2 *> facetDictionary;
	for( uint n=0; n<octree.terminals.size(); n++ ) {
		for( uint dim=0; dim<DIM; dim++ ) {
			for( char dir=-1; dir<=1; dir+=2 ) {
				const octree2::leaf2 &leaf = *octree.terminals[n];
				int dx = leaf.dx;
				vec2i fcenter = leaf.center+0.5*dx*vec2i(dir*(dim==0),dir*(dim==1));
				uint64 idx = computeIndex(fcenter);
				// If this turns to be a new facet, just append to the dictionary
				facet2 *facet;
				if( idx && facetDictionary.find(idx) == facetDictionary.end() ) {
					facet = new facet2;
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
		cell2 &cell = *cells[n];
		// For all center directions...
		FOR_EACH(3,3) {
			int dx = cell.dx;
			vec2i fcenter = cell.p-0.5*dx*vec2i(1,1)+0.5*dx*vec2i(i,j);
			uint64 idx = computeIndex(fcenter);
			if( facetDictionary.find(idx) != facetDictionary.end() ) {
				facet2 &facet = *facetDictionary[idx];
				facet.cells.push_back(&cell);
			}
		} END_FOR
	}
	
	// Remove some facets
	std::map<uint64,facet2 *> newFacetDictionary;
	for( std::map<uint64,facet2 *>::iterator it=facetDictionary.begin(); it!=facetDictionary.end(); it++ ) {
		facet2 &facet = *(*it).second;
		if( facet.cells.size() <= 1 && facet.edge == false ) {
			delete (*it).second;
		} else {
			newFacetDictionary[(*it).first] = (*it).second;
		}
	}
	facetDictionary = newFacetDictionary;
	
	// Connect facet and cells for cells
	for( uint n=0; n<cells.size(); n++ ) {
		cell2 &cell = *cells[n];
		// For all center directions...
		FOR_EACH(3,3) {
			int dx = cell.dx;
			vec2i fcenter = cell.p-0.5*dx*vec2i(1,1)+0.5*dx*vec2i(i,j);
			uint64 idx = computeIndex(fcenter);
			if( facetDictionary.find(idx) != facetDictionary.end() ) {
				facet2 *facet = facetDictionary[idx];
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
	for( std::map<uint64,facet2 *>::iterator it=facetDictionary.begin(); it!=facetDictionary.end(); it++ ) {
		facet2 *facet = (*it).second;
		facet->index = index;
		facets[index++] = facet;
	}
	
	// Build nodes
	for( std::map<uint64,node2 *>::iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		node2 *node = (*it).second;
		delete node;
	}
	nodes.clear();
	for( uint n=0; n<cells.size(); n++ ) {
		cell2 *cell = cells[n];
		int dx = cell->dx;
		vec2i corner = cell->p-0.5*dx*vec2i(1,1);
		FOR_EACH(2,2) {
			vec2i p = corner+dx*vec2i(i,j);
			uint idx = computeIndex(p);
			if( idx ) {
				node2 *node = NULL;
				if( nodes.find(idx) == nodes.end() ) {
					node = new node2();
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
				cell->nodes[i][j] = node;
			}
		} END_FOR
	}
	
	// Scatter samples beforehand
	for( uint dim=0; dim<DIM; dim++ ) {
		for( std::map<uint64,sample2 *>::iterator it=samples[dim].begin(); it!=samples[dim].end(); it++ ) {
			sample2 *sample = (*it).second;
			delete sample;
		}
		samples[dim].clear();
	}
	for( uint n=0; n<facets.size(); n++ ) {
		facet2 *facet = facets[n];
		sample2 *sample = new sample2();
		sample->p = facet->p;
		sample->facets[0] = facet;
		sample->weights[0] = 1.0;
		sample->num = 1;
		uint idx = computeIndex(facet->p);
		samples[facet->dim][idx] = sample;
	}
	for( uint dim=0; dim<DIM; dim++ ) {
		for( uint n=0; n<cells.size(); n++ ) {
			cell2 *cell = cells[n];
			int dx = cell->dx;
			bool findJunction = false;
			for( int dir=-1; dir<=1; dir+=2 ) {
				vec2i p = cell->p+0.5*dx*dir*vec2i(dim==0,dim==1);
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
							sample2 *sample = new sample2();
							sample->p = cell->p;
							int count = 0;
							for( int dir=-1; dir<=1; dir+=2 ) {
								vec2i p = cell->p+0.5*dx*dir*vec2i(dim2==0,dim2==1);
								uint idx = computeIndex(p);
								if( facetDictionary.find(idx) != facetDictionary.end() && facetDictionary[idx]->dim == dim2 ) {
									facet2 *facet = facetDictionary[idx];
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
	for( std::map<uint64,node2 *>::iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		node2 *node = (*it).second;
		int dx = node->dx;
		for( uint dim=0; dim<DIM; dim++ ) {
			// If on the T-junctoin, facet may be found just there
			vec2i p = node->p;
			uint idx = computeIndex(p);
			std::vector<FLOAT64> dists;	
			if( samples[dim].find(idx) != samples[dim].end() ) {
				node->samples[dim].push_back(samples[dim][idx]);
				dists.push_back(1.0);
			} else {
				// Otherwise probe for other directions
				FLOAT64 arms[] = { 0.5, 1.0, 1.5 };
				for( uint i=0; i<3; i++ ) {
					if( i == 2 && dists.size() ) break;
					for( int dir=-1; dir<=1; dir+=2 ) {
						// Try find corner ones
						vec2i p = node->p+dx*dir*arms[i]*vec2i(dim!=0,dim!=1);	
						uint idx = computeIndex(p);
						if( samples[dim].find(idx) != samples[dim].end() ) {
							sample2 *sample = samples[dim][idx];
							node->samples[dim].push_back(sample);
							dists.push_back(arms[i]);
						}
					}
				}
			}
			node->weights[dim].resize(dists.size());
			if( dists.size() == 2 ) {
				node->weights[dim][0] = dists[1]/(dists[0]+dists[1]);
				node->weights[dim][1] = dists[0]/(dists[0]+dists[1]);
			} else if( dists.size() == 1 ) {
				node->weights[dim][0] = 1.0;
			} else {
				if( node->samples[dim].size() != dists.size() ) {
					email::print( "Unexpected dists size.\n" );
					email::send();
					exit(0);
				}
			}
		}
	}
	
	// Connect facets and nodes
	for( uint n=0; n<facets.size(); n++ ) {
		facet2 *facet = facets[n];
		uint dim = facet->dim;
		int dx = facet->dx;
		// On T-junction, node can exist just under the facet
		uint idx = computeIndex(facet->p);
		if( nodes.find(idx) != nodes.end() ) {
			facet->nodes.push_back(nodes[idx]);
		} else {
			// If not found, just search for both directions
			for( int dir=-1; dir<=1; dir+=2 ) {
				vec2i p = facet->p+0.5*dir*dx*vec2i(dim!=0,dim!=1);
				uint idx = computeIndex(p);
				if( nodes.find(idx) != nodes.end() ) {
					facet->nodes.push_back(nodes[idx]);
				}
			}
		}
	}
	
	// Connect node and node
	for( std::map<uint64,node2 *>::iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		node2 *node = (*it).second;
		int dx = node->dx;
		for( uint dim=0; dim<DIM; dim++ ) {
			for( int dir=-1; dir<=1; dir+=2 ) {
				for( int arm=1; arm<=2; arm++ ) {
					vec2i p = node->p+dir*arm*dx*vec2i(dim==0,dim==1);
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

void octfluid2::setParameter( int name, FLOAT64 value ) {
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

void octfluid2::setTimestep( FLOAT64 dt ) {
	this->dt = dt;
}

static FLOAT64 fraction( FLOAT64 phi0, FLOAT64 phi1 ) {
	return phi0*phi1 > 0 ? 1 - (phi0 > 0) : -fmin(phi0,phi1)/fmax(fabs(phi1-phi0),1.0e-6);
}

void octfluid2::setupSolidLevelset( const levelset2 *solid ) {
	this->solid = solid;
	for( std::map<uint64,node2 *>::iterator it=this->nodes.begin(); it!=this->nodes.end(); it++ ) {
		node2 *node = (*it).second;
		vec2d p = vec2d(node->p)/octree.resolution;
		node->solidLS = solid->evalLevelset(p);
	}
	
	// Compute fractions of facets
	for( uint n=0; n<facets.size(); n++ ) {
		facet2 *facet = facets[n];
		FLOAT64 dx = facet->dx/(FLOAT64)octree.resolution;
		uint dim = facet->dim;
		vec2d p = vec2d(facet->p)/octree.resolution;
		vec2d up = p+0.5*dx*vec2d(dim!=0,dim!=1);
		vec2d down = p-0.5*dx*vec2d(dim!=0,dim!=1);
		facet->fraction = 1.0-fraction(solid->evalLevelset(up),solid->evalLevelset(down));
		if( ! variation ) facet->fraction = ceil(facet->fraction);
	}
}

void octfluid2::setupFluidLevelset( const levelset2 *fluid ) {
	this->fluid = fluid;
	for( std::map<uint64,node2 *>::iterator it=this->nodes.begin(); it!=this->nodes.end(); it++ ) {
		node2 *node = (*it).second;
		vec2d p = vec2d(node->p)/octree.resolution;
		node->fluidLS = fluid->evalLevelset(p);
	}
	for( uint n=0; n<cells.size(); n++ ) {
		cells[n]->fluidLS = fluid->evalLevelset(vec2d(cells[n]->p)/octree.resolution);
	}
}

void octfluid2::getSamplePoints( std::vector<vec2d> &pos ) {
	// Gather facet position
	pos.resize(nodes.size());
	uint n=0;
	for( std::map<uint64,node2 *>::const_iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		const node2 *node = (*it).second;
		pos[n++] = vec2d(node->p)/octree.resolution;
	}
}

void octfluid2::setupVelocity( const std::vector<vec2d> &vel, const std::vector<bool> &mapped ) {
	// Set whether node is mapped
	uint index=0;
	for( std::map<uint64,node2 *>::iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		node2 *node = (*it).second;
		node->velocity = vel[index];
		for( uint dim=0; dim<DIM; dim++ ) {
			if( node->edge[dim] ) node->velocity[dim] = 0.0;
		}
		node->mapped = mapped[index];
		index ++;
	}
}

void octfluid2::resample( int kind, int type ) { // type 0: node -> facet, 1: facet -> node
	if( type == 0 ) {
		// Resample node -> faccet
		for( uint n=0; n<facets.size(); n++ ) {
			facet2 *facet = facets[n];
			facet->mapped = false;
			uint dim = facet->dim;
			FLOAT64 value=0.0;
			uint sum=0;
			for( std::list<node2 *>::iterator it=facet->nodes.begin(); it!=facet->nodes.end(); it++ ) {
				node2 *node = *it;
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
		for( std::map<uint64,node2 *>::iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
			node2 *node = (*it).second;
			if( kind == 0 ) node->velocity = vec2d();
			else node->gradient = vec2d();
			bool mapped[DIM] = { false, false };
			for( uint dim=0; dim<DIM; dim++ ) {
				FLOAT64 wsum = 0.0;
				for( uint n=0; n<node->samples[dim].size(); n++ ) {
					sample2 *sample = node->samples[dim][n];
					FLOAT64 w = node->weights[dim][n];
					for( uint m=0; m<sample->num; m++ ) {
						facet2 *facet = sample->facets[m];
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

void octfluid2::buildMatrix() {
	// Compute a facet distance "ds" and find an adjacent biggest cell
	for( uint n=0; n<facets.size(); n++ ) {
		facet2 *facet = facets[n];
		uint dim = facet->dim;
		// Find a cell with maximal cell area
		FLOAT64 max_area = 0.0;
		cell2 *Lcell = NULL;
		for( std::list<cell2 *>::iterator it=facet->cells.begin(); it!=facet->cells.end(); it++ ) {
			cell2 *cell = *it;
			FLOAT64 area = cell->dx;
			if( max_area < area ) {
				max_area = area;
				Lcell = cell;
			}
		}
		facet->Lcell = Lcell;
		
		// Compute average (signed) distance between small and large cells
		FLOAT64 ds=0.0;
		for( std::list<cell2 *>::iterator it=facet->cells.begin(); it!=facet->cells.end(); it++ ) {
			cell2 *Scell = *it;
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
		cell2 *cell = cells[n];
		if( cell->fluidLS > 0.0 ) continue;
		for( std::list<facet2 *>::iterator it=cell->facets.begin(); it!=cell->facets.end(); it++ ) {
			facet2 *facet = *it;
			uint dim = facet->dim;
			int sgn = cell->p[dim] < facet->p[dim] ? 1 : -1;
			FLOAT64 ds = facet->ds;
			cell2 *Lcell = facet->Lcell;
			if( cell == Lcell ) {
				for( std::list<cell2 *>::iterator it2=facet->cells.begin(); it2!=facet->cells.end(); it2++ ) {
					cell2 *Scell = *it2;
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
					cell2 *Scell = cell;
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

void octfluid2::project() {	
	
	// Extrapolate
	extrapolate(0);

	// Resample to facets
	resample(0,0);
	
	// Build matrix
	buildMatrix();
	
	// Compute divergence
	std::vector<FLOAT64> rhs(cells.size());
	std::vector<FLOAT64> fvel(facets.size());
	for( uint n=0; n<facets.size(); n++ ) {
		fvel[n] = facets[n]->velocity;
	}
	multiply<FLOAT64>(RhsMatrix,fvel,rhs);
	
	// Copy divergence
	for( uint n=0; n<cells.size(); n++ ) {
		cells[n]->divergence = rhs[n];
	}
	SparseMatrix<FLOAT64> LhsMatrix = this->LhsMatrix;
	
	// Compact the matrix
	std::vector<int> idx_map;
	compress( LhsMatrix, rhs, idx_map );
	
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
	if( residual_out > 1.0 ) converged = false;
	if( ! converged ) {
		email::print("PCG did not converge.\n" );
		email::send();
		exit(0);
	} else {
		dump( "OCT: PCG Converged ! %d iterations and took %.3f msec with %d unknowns.\n", iterations, msec, LhsMatrix.n);
	}
	
	// Set pressure
	for( uint n=0; n<cells.size(); n++ ) {
		cells[n]->pressure = idx_map[n] >= 0 ? result[idx_map[n]] : 0.0;
	}
	
	// Compute gradient
	for( uint n=0; n<facets.size(); n++ ) {
		facet2 *facet = facets[n];
		uint dim = facet->dim;
		FLOAT64 ds = facet->ds;
		cell2 *Lcell = facet->Lcell;
		int sgn = Lcell->p[dim] < facet->p[dim] ? 1 : -1;
		FLOAT64 gradient = 0.0;
		facet->mapped = false;
		if( facet->fraction ) {
			for( std::list<cell2 *>::iterator it2=facet->cells.begin(); it2!=facet->cells.end(); it2++ ) {
				cell2 *Scell = *it2;
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
	
	// Resample to nodes
	resample(1,1);

	// Extrapolate
	extrapolate(1);
}

vec2d octfluid2::getPressureGradient( vec2d p ) const {
	return sampleVector(p,1);
}

vec2d octfluid2::getVelocity( vec2d p ) const {
	return sampleVector(p,0);
}

FLOAT64 octfluid2::getDivergence( vec2d p ) const {
	int n = octree.hitTest(p);
	if( n >= 0 ) {
		return cells[n]->divergence;
	}
	return 0.0;
}

vec2d octfluid2::sampleVector( vec2d p, char kind ) const {
	int hit = octree.hitTest(p);
	if( hit >= 0 ) {
		// Do sample here...
		FLOAT64 dx = cells[hit]->dx/(FLOAT64)octree.resolution;
		FLOAT64 x = (p[0]-cells[hit]->nodes[0][0]->p[0]/(FLOAT64)octree.resolution)/dx;
		FLOAT64 y = (p[1]-cells[hit]->nodes[0][0]->p[1]/(FLOAT64)octree.resolution)/dx;
		vec2d v[4];
		if( kind == 1 ) {
			FLOAT64 L[4];
			L[0] = cells[hit]->nodes[0][0]->fluidLS;
			L[1] = cells[hit]->nodes[1][0]->fluidLS;
			L[2] = cells[hit]->nodes[1][1]->fluidLS;
			L[3] = cells[hit]->nodes[0][1]->fluidLS;
			FLOAT64 fluidLS = L[0]*(1.-x)*(1.-y)+L[1]*x*(1.-y)+L[2]*x*y+L[3]*(1.-x)*y;
			if( fluidLS > 0.0 ) kind = 0;
		}
		if( kind == 0 ) {
			v[0] = cells[hit]->nodes[0][0]->velocity;
			v[1] = cells[hit]->nodes[1][0]->velocity;
			v[2] = cells[hit]->nodes[1][1]->velocity;
			v[3] = cells[hit]->nodes[0][1]->velocity;
		} else {
			v[0] = cells[hit]->nodes[0][0]->gradient/dt;
			v[1] = cells[hit]->nodes[1][0]->gradient/dt;
			v[2] = cells[hit]->nodes[1][1]->gradient/dt;
			v[3] = cells[hit]->nodes[0][1]->gradient/dt;
		}
		return v[0]*(1.-x)*(1.-y)+v[1]*x*(1.-y)+v[2]*x*y+v[3]*(1.-x)*y;
	}
	return vec2d();
}

void octfluid2::extrapolate( char kind ) {
	// Assign index
	uint n=0;
	for( std::map<uint64,node2 *>::iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		node2 *node = (*it).second;
		node->index = n++;
	}
	std::vector<fastmarch2<vec2d>::node2 *> fnodes(nodes.size());
	for( uint n=0; n<nodes.size(); n++ ) fnodes[n] = new fastmarch2<vec2d>::node2;
	for( std::map<uint64,node2 *>::iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		node2 *node = (*it).second;
		uint n = node->index;
		vec2d p = node->p;
		bool fixed = node->mapped;
		fnodes[n]->p = p;
		fnodes[n]->fixed = fixed;
		fnodes[n]->levelset = fmax(node->fluidLS,-node->solidLS);
		fnodes[n]->value = (kind==0) ? node->velocity : node->gradient;
		fnodes[n]->p2p.resize(node->nodes.size());
		uint m = 0;
		for( std::list<node2 *>::iterator it=node->nodes.begin(); it!=node->nodes.end(); it++ ) {
			fnodes[n]->p2p[m++] = fnodes[(*it)->index];
		}
	}
	fastmarch2<vec2d>::fastMarch(fnodes,extrapolate_dist,-1.0,0);
	// Pick up extrapolated values
	for( std::map<uint64,node2 *>::iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		node2 *node = (*it).second;
		uint n = node->index;
		if( fnodes[n]->fixed ) {
			if (kind==0) node->velocity = fnodes[n]->value;
			else node->gradient = fnodes[n]->value;
		}
		delete fnodes[n];
	}
}

void octfluid2::drawLevelset( int kind ) const {
	// Draw wall levelset
	for( uint n=0; n<cells.size(); n++ ) {
		cell2 *cell = cells[n];
		FLOAT64 dx = cell->dx/(FLOAT64)octree.resolution;
		vec2d origin = vec2d(cell->p)/(FLOAT64)octree.resolution-0.5*dx*vec2d(1.0,1.0);
		FLOAT64 solidLS[2][2];
		vec2d p[8];
		int pnum;
		FOR_EACH(2,2) {
			solidLS[i][j] = kind==0 ? cell->nodes[i][j]->fluidLS : cell->nodes[i][j]->solidLS;
		} END_FOR
		levelset2::marchPoints(solidLS,p,pnum);
		for( uint m=0; m<pnum; m++ ) {
			p[m] = p[m]*dx+origin;
		}
		glBegin(GL_TRIANGLE_FAN);
		for( uint m=0; m<pnum; m++ ) {
			glVertex2d(p[m][0],p[m][1]);
		}
		glEnd();
	}
}

static void drawBitmapString( const char *string, void *font=GLUT_BITMAP_HELVETICA_12 ) {
	if( ! font ) font = GLUT_BITMAP_HELVETICA_12;
	while (*string) glutBitmapCharacter(font, *string++);
}

void octfluid2::render( int name, vec2d mousePos ) const {
	// Find closest facet
	FLOAT64 dist = 1e9;
	facet2 *facet = NULL;
	for( uint n=0; n<facets.size(); n++ ) {		
		vec2d pos = vec2d(facets[n]->p)/octree.resolution;
		FLOAT64 d = (mousePos-pos).len();
		if( d < dist ) {
			facet = facets[n];
			dist = d;
		}
	}
	
	// Find closest node
	dist = 1e9;
	const node2 *node = NULL;
	for( std::map<uint64,node2 *>::const_iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
		vec2d pos = vec2d((*it).second->p)/octree.resolution;
		FLOAT64 d = (mousePos-pos).len();
		if( d < dist ) {
			node = (*it).second;
			dist = d;
		}
	}
	
	// Find closest cell
	dist = 1e9;
	const cell2 *cell = NULL;
	for( uint n=0; n<cells.size(); n++ ) {
		vec2d pos = vec2d(cells[n]->p)/octree.resolution;
		FLOAT64 d = (mousePos-pos).len();
		if( d < dist ) {
			cell = cells[n];
			dist = d;
		}
	}
	FLOAT64 dx = 1e9;
	for( uint n=0; n<cells.size(); n++ ) {
		dx = fmin(dx,cells[n]->dx/(FLOAT64)octree.resolution);
	}
	
	switch( name ) {
		case MESH: {
			// Draw octree
			glColor4d(1.0,1.0,1.0,0.05);
			octree.drawOctree();
			break;
		}
		case FLUID: {
			// Draw fluid levelset
			glColor4d(0.3,0.3,0.8,0.6);
			drawLevelset(0);
			break;
		}
		case SOLID: {
			// Draw wall levelset
			glColor4d(0.7,0.7,0.3,0.6);
			drawLevelset(1);
			break;
		}
		case MATRIX_CONNECTION: {
			// Render neighboring cells
			if( cell ) {
				uint i = cell->index;
				if( LhsMatrix.n ) {
					FLOAT64 diag = LhsMatrix(cell->index,cell->index);
					for( uint j=0; j<LhsMatrix.index[i].size(); j++ ) {
						uint index = LhsMatrix.index[i][j];
						const cell2 &cell = *cells[index];
						FLOAT64 dx = cell.dx/(FLOAT64)octree.resolution;
						vec2d p = vec2d(cell.p[0],cell.p[1])/(FLOAT64)octree.resolution;
						glColor4d(1.0,0.5,0.5,0.5);
						glRectf(p[0]-0.5*dx,p[1]-0.5*dx,p[0]+0.5*dx,p[1]+0.5*dx);
						
						if( index != i ) {
							glColor4d(1.0,1.0,1.0,1.0);
							glRasterPos2d(p[0],p[1]);
							drawBitmapString(format_str("%.3f", -LhsMatrix.value[i][j]/diag));
						}
					}
				}
			}
			break;
		}
		case VELOCITY: {
			// Draw corner velocity
			for( std::map<uint64,node2 *>::const_iterator it=nodes.begin(); it!=nodes.end(); it++ ) {
				const node2 *node = (*it).second;
				vec2d p = vec2d(node->p)/octree.resolution;
				vec2d vel = node->velocity;
				glColor4d(1.0,1.0,0.0,0.8);
				glLineWidth(2.0);
				glBegin(GL_LINES);
				glVertex2dv(p.v);
				glVertex2dv((p+dx*vel).v);
				glEnd();
				glLineWidth(1.0);
			}
			break;
		}
		case PRESSURE: {
			// Visualize pressure
			FLOAT64 minv = 1e9;
			FLOAT64 maxv = -1e9;
			for( uint n=0; n<cells.size(); n++ ) {
				minv = fmin(minv,cells[n]->pressure);
				maxv = fmax(maxv,cells[n]->pressure);
			}
			FLOAT64 det = maxv - minv;
			if( det ) {
				for( uint n=0; n<cells.size(); n++ ) {
					bool hasAir = false;
					for( std::list<facet2 *>::iterator it=cells[n]->facets.begin(); it!=cells[n]->facets.end(); it++ ) {
						facet2 *facet = *it;
						if( facet->fraction ) hasAir = true;
					}
					if( cells[n]->fluidLS < 0.0 && hasAir ) {
						vec2d p[4] = {	
							vec2d(cells[n]->nodes[0][0]->p)/octree.resolution,
							vec2d(cells[n]->nodes[1][0]->p)/octree.resolution,
							vec2d(cells[n]->nodes[1][1]->p)/octree.resolution,
							vec2d(cells[n]->nodes[0][1]->p)/octree.resolution };
						FLOAT64 normp = 2.0*(cells[n]->pressure-minv)/det-1.0;
						glColor4d(normp>0,0.3,normp<=0,fabs(normp));
						glBegin(GL_QUADS);
						for( uint m=0; m<4; m++ ) glVertex2dv(p[m].v);
						glEnd();
					}
				}
			}
			break;
		}
		case NAME: {
			glColor4f(1.0,1.0,1.0,1.0);
			glRasterPos2d(0.07,0.92);
			drawBitmapString("OCT",GLUT_BITMAP_HELVETICA_18);
			break;
		}
	}
}