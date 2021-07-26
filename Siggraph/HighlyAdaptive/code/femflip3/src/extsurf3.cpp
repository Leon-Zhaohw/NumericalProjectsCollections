/*
 *	extsurf3.cpp
 *
 *	Created by Ryoichi Ando on Dec 15 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <math.h>
#include <zlib.h>
#include "octlevelset3.h"
#include "meshsurf3.h"
#include "extsurf3.h"
#include "octree3.h"
#include "particle3.h"
#include "util.h"
#include "kernel.h"
#include "opengl.h"
#include "kernel.h"

extsurf3::extsurf3() {
	numSample = 20;
	method = 0;
	kernel_size = 4.0;
}

void extsurf3::writeFile( const char *path, const std::vector<particle3 *> &particles, FLOAT64 dpx, const octree3 *octree ) {
	gzFile fp = gzopen(path,"w");
	if( ! fp ) {
		printf( "\nCould not open path %s for writing particles.\n", path );
		exit(0);
	}
	// Load particles
	loadParticles(particles,dpx,octree);
	
	// First write the number of particles
	uint number = spheres.size();
	gzwrite(fp,&number,sizeof(uint));
	
	// And write the dpx
	gzwrite(fp,&dpx,sizeof(FLOAT64));
	
	// Write each sphere info
	for( uint n=0; n<spheres.size(); n++ ) {
		sphere3 &p = spheres[n];
		gzwrite(fp,&p.p,sizeof(vec3d));
		gzwrite(fp,&p.r,sizeof(FLOAT64));
		gzwrite(fp,&p.remesh_r,sizeof(FLOAT64));
		gzwrite(fp,&p.n,sizeof(uint));
		gzwrite(fp,&p.levelset,sizeof(FLOAT64));
		gzwrite(fp,&p.gradient,sizeof(vec3d));
		gzwrite(fp,&p.isolated,sizeof(bool));
		gzwrite(fp,&p.depth,sizeof(int));
	}
	gzclose(fp);
}

void extsurf3::loadFile( const char *path ) {
	gzFile fp = gzopen(path,"r");
	if( ! fp ) {
		printf( "Could not open path %s for reading particles.\n", path );
		exit(0);
	}

	// First read the number of particles
	uint number;
	gzread(fp,&number,sizeof(uint));
	
	tick(); dump( "Loading %d particle data...", number);
	
	// And read the dpx
	gzread(fp,&dpx,sizeof(FLOAT64));
	
	// Read each sphere info
	spheres.resize(number);
	for( uint n=0; n<spheres.size(); n++ ) {
		sphere3 &p = spheres[n];
		gzread(fp,&p.p,sizeof(vec3d));
		gzread(fp,&p.r,sizeof(FLOAT64));
		gzread(fp,&p.remesh_r,sizeof(FLOAT64));
		gzread(fp,&p.n,sizeof(uint));
		gzread(fp,&p.levelset,sizeof(FLOAT64));
		gzread(fp,&p.gradient,sizeof(vec3d));
		gzread(fp,&p.isolated,sizeof(bool));
		gzread(fp,&p.depth,sizeof(int));
	}
	gzclose(fp);
	dump( "Done. Took %s.\n", stock("extmesher_load"));

	for( uint n=0; n<spheres.size(); n++ ) {
		spheres[n].p[1] += 4*dpx;
	}
}

void extsurf3::setMethod( int method ) {
	this->method = method;
}

void extsurf3::loadParticles( const std::vector<particle3 *> &particles, FLOAT64 dpx, const octree3 *octree ) {
	spheres.resize(particles.size());
	FOR_EACH_PARTICLE(particles) {
		sphere3 sphere;
		sphere.p = p.p;
		sphere.r = 0.5*dpx*p.r;
		sphere.remesh_r = p.remesh_r;
		sphere.n = p.n;
		sphere.levelset = p.levelset;
		sphere.gradient = p.gradient;
		sphere.isolated = p.isolated;
		if( octree ) {
			sphere.depth = octree->terminals[octree->hitTest(sphere.p)]->depth;
		} else {
			sphere.depth = -1;
		}
		spheres[pcount] = sphere;
	} END_FOR
}

static FLOAT64 get_dx( vec3d pos, const octree3 &octree, FLOAT64 dpx ) {
	FLOAT64 base_dx = 1.0;
	std::vector<uint> array;
	if( octree.hitTest( pos, 0.5*dpx, array ) ) {
		for( uint n=0; n<array.size(); n++ ) {
			int hit = array[n];
			base_dx = fmin(base_dx,octree.terminals[hit]->dx/(FLOAT64)octree.resolution);
		}
	} else {
		dump( "extsurf3::evalLevelset base grid size failed to fetch at (%f,%f)\n", pos[0], pos[1] );
	}
	return base_dx;
}

class bcclevelset3 : public levelset3 {
public:
	bcclevelset3() {
		bcc = NULL;
		q = NULL;
	}
	virtual void setBCC( const bcc3 &bcc ) {
		this->bcc = &bcc;
	}
	virtual void setNodalLevelset( const std::vector<FLOAT64> &q ) {
		this->q = &q;
	}
	virtual FLOAT64 evalLevelset(vec3d p) const {
		int n = bcc->hitTest(p);
		if( n >= 0 ) {
			FLOAT64 x[NUM_VERT] = { p[0], p[1], p[2], 1.0 };
			FLOAT64 div = 0.0;
			for( uint i=0; i<NUM_VERT; i++ ) {
				FLOAT64 t = 0.0;
				for( uint j=0; j<NUM_VERT; j++ ) {
					t += bcc->matrix[n].m[i][j]*x[j];
				}
				div += t*q->at(bcc->cleanElements[n][i]);
			}
			return div;
		}
		return 1e9;
	}
	const bcc3 *bcc;
	const std::vector<FLOAT64> *q;
};

class solid_hint3 : public levelset3 {
public:
	solid_hint3( FLOAT64 dx, levelset3 *solid, const ann3 &sorter, const std::vector<vec3d> &positions, const std::vector<FLOAT64> &positions_r,
				const ann3 &inside_sorter, const std::vector<vec3d> &inside_positions, const std::vector<FLOAT64> &inside_r ) :
	sorter(sorter), positions(positions), positions_r(positions_r), inside_sorter(inside_sorter), inside_positions(inside_positions), inside_r(inside_r)
	{
		this->solid = solid;
		this->dx = dx;
	}
	virtual FLOAT64 evalLevelset(vec3d p) const {
		FLOAT64 min_dist = 1e9;
		FLOAT64 phi = solid->evalLevelset(p);
		if( fabs(phi) < dx ) {
			std::vector<ANNidx> closeOne = sorter.getkNeighbors(p,1);
			if(closeOne.size()) {
				vec3d cpos = positions[closeOne[0]];
				FLOAT64 r = positions_r[closeOne[0]];
				min_dist = fmin(min_dist,2*dx+(cpos-p).len()-r);
			}
			if( fabs(solid->evalCurvature(p)) > 10.0 ) {
				std::vector<ANNidx> closeOne = inside_sorter.getkNeighbors(p,1);
				if(closeOne.size()) {
					vec3d cpos = inside_positions[closeOne[0]];
					FLOAT64 r = inside_r[closeOne[0]];
					if( phi < 4.0*r && phi > -0.5*r ) {
						min_dist = fmin(min_dist,fabs(4*dx+phi));
					}
				}
			}
		}
		return fmax(dx,min_dist);
	}
	levelset3 *solid;
	FLOAT64 dx;
	const ann3 &sorter;
	const std::vector<vec3d> &positions;
	const std::vector<FLOAT64> &positions_r;
	
	const ann3 &inside_sorter;
	const std::vector<vec3d> &inside_positions;
	const std::vector<FLOAT64> &inside_r;
};

static void computeAnisotropy( const ann3& sorter, const std::vector<extsurf3::sphere3> &spheres, std::vector<svd3> &anisotropy, FLOAT64 r ) {
	anisotropy.resize(spheres.size());
	for( uint n=0; n<spheres.size(); n++ ) {
		anisotropy[n] = svd3();
	}	
    for( uint n=0; n<spheres.size(); n++ ) {
		// Compute SVD info
		svd3 svd;
		if( spheres[n].levelset > -3.0*spheres[n].r ) {
			// Let's make it small as default
			for( uint i=0; i<DIM; i++ ) svd.eig[i] *= 0.5;
			
			// ks is heuristically given
			FLOAT64 ks = 1.0/sqr(spheres[n].r);
			FLOAT64 C[DIM][DIM];
			
			// Gather neighbor particles
			std::vector<ANNidx> neighbors = sorter.getNeighbors(spheres[n].p,r*spheres[n].r);
			// If the number of neighbors is too small, skip SVD
			if( neighbors.size() >= 2 ) {
				// Compute center position first
				vec3d center;
				FLOAT64 wsum = 0.0;
				for( uint m=0; m<neighbors.size(); m++ ) {
					FLOAT64 w = kernel::smooth_kernel((spheres[n].p-spheres[neighbors[m]].p).len2(),r*spheres[n].r);
					center += w*spheres[neighbors[m]].p;
					wsum += w;
				}
				if( wsum ) center = center / wsum;
				else center = spheres[n].p;
				
				// Then compute covariance
				for( uint i=0; i<DIM; i++ ) for( uint j=0; j<DIM; j++ ) {
					FLOAT64 wsum = 0.0;
					C[i][j] = 0.0;
					if(i<=j) {
						for( uint m=0; m<neighbors.size(); m++ ) {
							FLOAT64 w = kernel::smooth_kernel((center-spheres[neighbors[m]].p).len2(),r*spheres[n].r);
							wsum += w;
							C[i][j] += w*(center[i]-spheres[neighbors[m]].p[i])*(center[j]-spheres[neighbors[m]].p[j]);
						}
						if( wsum ) C[i][j] /= wsum;
					}
				}
				for( uint i=0; i<DIM; i++ ) for( uint j=0; j<DIM; j++ ) {
					if(i>j) C[i][j]=C[j][i];
				}
				// Now run SVD !
				if( svd.run(C) ) {
					// If succeeded, multiply by some constant to scale eigen values.
					for( uint i=0; i<DIM; i++ )	svd.eig[i] *= ks;
				} else {
					printf( "Failed for some reason...\n" );
					exit(0);
				}
			}
		}
		anisotropy[n] = svd;
	}
}

void extsurf3::writeMesh( uint frame, levelset3 *solid ) {
	if( spheres.empty() ) {
		printf( "Empty data. Load first.\n");
		exit(0);
	}
	// Estimate the maximal depth
	FLOAT64 min_r = 1e9;
	for( uint n=0; n<spheres.size(); n++ ) {
		min_r = fmin(min_r,spheres[n].r);
	}
	solid->resize(1.0/(2.0*dpx));
	
	this->solid = solid;
	uint maxdepth = log2(1.0/min_r);
	
	// Sort spheres
	tick(); dump( "Removing floating particles..." );
	std::vector<vec3d> positions(spheres.size());
	for( uint n=0; n<spheres.size(); n++ ) {
		positions[n] = spheres[n].p;
	}
	sorter.sort(positions);
	
	// Remove floating spheres
	std::vector<bool> is_floating(spheres.size());
	PARALLEL_FOR for( uint n=0; n<spheres.size(); n++ ) {
		is_floating[n] = false;
		std::vector<ANNidx> neighbors = sorter.getkNeighbors(spheres[n].p,numSample);
		uint larger_count = 0;
		uint smaller_count = 0;
		for( uint m=0; m<neighbors.size(); m++ ) {
			uint idx = neighbors[m];
			if( idx != n ) {
				if( 1.9*spheres[n].r < spheres[idx].r ) {
					larger_count ++;
				} else {
					smaller_count ++;
				}
			}
			if( larger_count > smaller_count ) {
				is_floating[n] = true;
			}
		}
	}
	uint valid_count = 0;
	uint removed_count = 0;
	for( uint n=0; n<spheres.size(); n++ ) {
		if( ! is_floating[n] ) valid_count ++;
		else removed_count ++;
	}
	std::vector<sphere3> new_spheres(valid_count);
	uint index = 0;
	for( uint n=0; n<spheres.size(); n++ ) {
		if( ! is_floating[n] ) new_spheres[index++] = spheres[n];
	}
	spheres = new_spheres;
	
	// Sort spheres again
	positions.resize(spheres.size());
	for( uint n=0; n<spheres.size(); n++ ) {
		positions[n] = spheres[n].p;
	}
	sorter.sort(positions);
	dump("Done. Took %s. Removed %d floating particles.\n", stock(), removed_count );
	
	// Build a sphere array for octree
	std::vector<octree3::sphere3 *> oct_spheres(spheres.size());
	for( uint n=0; n<spheres.size(); n++ ) {
		bool aggresive = spheres[n].levelset > 0.0 || spheres[n].isolated;
		octree3::sphere3 *oct_sphere = new octree3::sphere3;
		oct_sphere->p = spheres[n].p;
		if( aggresive ) {
			oct_sphere->r = fmax(spheres[n].r,-spheres[n].levelset);
		} else {
			if( spheres[n].depth >= 0 ) oct_sphere->r = 1.0/powf(2,spheres[n].depth);
			else oct_sphere->r = 2.0 * fmax(spheres[n].r,-spheres[n].levelset);
		}
		// oct_sphere->r *= 0.5; for comparison
		oct_spheres[n] = oct_sphere;
	}
	
	// Sort surface particles
	std::vector<vec3d> surface_positions;
	surface_spheres.clear();
	for( uint n=0; n<spheres.size(); n++ ) {
		if( spheres[n].levelset > -4.0*spheres[n].r ) {
			surface_spheres.push_back(n);
			surface_positions.push_back(spheres[n].p);
		}
	}
	surf_sorter.sort(surface_positions);
	
	// Find surface particles that are close enough to surface
#if 1
	ann3 surf_solid_sorter;
	ann3 inside_sorter;
	std::vector<vec3d> surf_solid_positions;
	std::vector<FLOAT64> surf_solid_r;
	std::vector<vec3d> inside_positions;
	std::vector<FLOAT64> inside_r;
	
	for( uint n=0; n<spheres.size(); n++ ) {
		if( spheres[n].levelset > -2.0*spheres[n].r && spheres[n].levelset < -0.25*spheres[n].r &&
			solid->evalLevelset(spheres[n].p) < 2.0*spheres[n].r && ! spheres[n].isolated ) {
			surf_solid_positions.push_back(spheres[n].p);
			surf_solid_r.push_back(spheres[n].r);
		} else {
			inside_positions.push_back(spheres[n].p);
			inside_r.push_back(spheres[n].r);
		}
	}
	surf_solid_sorter.sort(surf_solid_positions);
	inside_sorter.sort(inside_positions);
	
	// Build octree for the final one
	solid_hint3 hint(2*min_r,solid,surf_solid_sorter,surf_solid_positions,surf_solid_r,inside_sorter,inside_positions,inside_r);
#endif
	octree.buildOctree(NULL,oct_spheres,maxdepth);

	// Generate BCC mesh
	{
		bcc3 bcc;
		bcc.buildBCC(octree);
		
		// Assign nodes
		mesher.nodes = bcc.getNodes();
		
		// Assign elements
		mesher.elements = bcc.getElements();
		
		// Free tets
		bcc.freeTets();
	}
	
	// Update connection
	mesher.setGenHash(false);
	mesher.setGenFacet(false);
	mesher.updateConnection();
#if 0
	// Precompute every levelset here on corase octree
	octree3 unfold_octree;
	unfold_octree.buildOctree(NULL,oct_spheres,maxdepth);
	
	// Delete the sphere array for octree generation
	for( uint n=0; n<spheres.size(); n++ ) delete oct_spheres[n];
	
	// Generate BCC mesh for unfoled mesh
	bcc3 unfold_bcc;
	unfold_bcc.buildBCC(unfold_octree);
	
	const std::vector<vec3d> &nodes = unfold_bcc.cleanNodes;
	std::vector<FLOAT64> levelsets(nodes.size());
	PARALLEL_FOR for( uint n=0; n<nodes.size(); n++ ) {
		levelsets[n] = evalLevelset(nodes[n]);
	}
	bcclevelset3 bccLevelset;
	bccLevelset.setBCC(unfold_bcc);
	bccLevelset.setNodalLevelset(levelsets);
#endif
	// Yu and turk method is specified. Should compute particle nisotropy.
	if( method == 3 ) {
		// Compute density of particles
		density.resize(spheres.size());
		PARALLEL_FOR for( uint n=0; n<spheres.size(); n++ ) {
			std::vector<ANNidx> neighbors = sorter.getkNeighbors(spheres[n].p,numSample);
			FLOAT64 dens = 0.0;
			for( uint m=0; m<neighbors.size(); m++ ) {
				FLOAT64 r = spheres[neighbors[m]].r;
				dens += kernel::smooth_kernel((spheres[n].p-spheres[neighbors[m]].p).len2(),kernel_size*r);
			}
			density[n] = dens / powf(dpx,DIM);
		}
		// Compute anisotropy
		computeAnisotropy(sorter,spheres,anisotropy,kernel_size);
	}
	
	// Generate a surface mesh
	meshsurf3 meshsurf;
	meshsurf.setReference(&mesher);
	meshsurf.setExtrapolationDist(0.0);
	meshsurf.buildSurface(vertices,normals,faces,this,solid,true,0.0,2,false);
	for( uint n=0; n<vertices.size(); n++ ) {
		vertices[n][1] -= 4*dpx;
	}
}

particle3 extsurf3::genParticle( const sphere3 &sphere ) const {
	particle3 particle;
	particle.p = sphere.p;
	particle.r = 2.0 * sphere.r / dpx;
	particle.levelset = sphere.levelset;
	particle.gradient = sphere.gradient;
	particle.isolated = sphere.isolated;
	return particle;
}

bool extsurf3::getClosestSurfacePos(vec3d &pos) const {
	FLOAT64 min_phi = 1e9;
	if( solid->evalLevelset(pos) < get_dx(pos,octree,dpx) ) {
		return false;
	}
	std::vector<ANNidx> surf_neighbors = surf_sorter.getkNeighbors(pos,numSample);
	std::vector<particle3 *> p_neighbors(surf_neighbors.size());
	for( uint n=0; n<surf_neighbors.size(); n++ ) {
		particle3 *particle = new particle3;
		*particle = genParticle(spheres[surface_spheres[surf_neighbors[n]]]);
		p_neighbors[n] = particle;
	}
	FLOAT64 r = 0.0;
	FLOAT64 rsum = 0.0;
	for( uint n=0; n<p_neighbors.size(); n++ ) {
		r += p_neighbors[n]->r;
		rsum += 1.0;
	}
	if( rsum ) r = r / rsum;
	vec3d p = pos;
	FLOAT64 shrink = 0.85;
	min_phi = fmin(min_phi,levelset3::getParticleLevelset(p_neighbors,dpx,p,pos,solid,shrink,false));
	for( uint n=0; n<p_neighbors.size(); n++ ) {
		delete p_neighbors[n];
	}
	return fabs(min_phi) < dpx*r;
}

FLOAT64 extsurf3::evalLevelset(vec3d pos) const {
	// Collect neighbors
	std::vector<ANNidx> neighbors = sorter.getkNeighbors(pos,numSample);
	std::vector<sphere3> p_neighbors(neighbors.size());
	for( uint n=0; n<neighbors.size(); n++ ) {
		p_neighbors[n] = spheres[neighbors[n]];
	}
	
	// Slide position
	int max_atempt = 5;
	FLOAT64 solid_phi = solid->evalLevelset(pos);
	while( solid_phi < 0.0 && max_atempt >= 0 ) {
		pos += -solid_phi*solid->evalGradient(pos);
		solid_phi = solid->evalLevelset(pos);
		max_atempt --;
	}
	
	// Slide position aggresively...
#if 1
	{
		int max_atempt = 5;
		FLOAT64 max_r = 0.0;
		for( uint n=0; n<p_neighbors.size(); n++ ) {
			max_r = fmax(p_neighbors[n].r,max_r);
		}
		FLOAT64 solid_phi = solid->evalLevelset(pos);
		while( solid_phi < 2*max_r && max_atempt >= 0 ) {
			pos += (2*max_r-solid_phi)*solid->evalGradient(pos);
			solid_phi = solid->evalLevelset(pos);
			max_atempt --;
		}
	}
#endif
	
	// Ours
	if( method == 0 ) {
		FLOAT64 max_dx = 0.0;
		FLOAT64 min_phi = 1e9;
		FLOAT64 shrink = 0.85;
		for( uint n=0; n<neighbors.size(); n++ ) {
			const sphere3 &p = spheres[neighbors[n]];
			FLOAT64 r = p.r;
			r = fmax(r,-shrink*p.levelset);
			min_phi = fmin(min_phi,(p.p-pos).len()-r);
			max_dx = fmax(1.0/powf(2,p.depth),max_dx);
		}
		if( fabs(min_phi) < 2*max_dx ) {
			std::vector<ANNidx> surf_neighbors = surf_sorter.getkNeighbors(pos,numSample);
			std::vector<particle3 *> p_neighbors(surf_neighbors.size());
			for( uint n=0; n<surf_neighbors.size(); n++ ) {
				particle3 *particle = new particle3;
				*particle = genParticle(spheres[surface_spheres[surf_neighbors[n]]]);
				p_neighbors[n] = particle;
			}
			vec3d outpos;
			min_phi = fmin(min_phi,levelset3::getParticleLevelset(p_neighbors,dpx,pos,outpos,solid,shrink,false));
			for( uint n=0; n<p_neighbors.size(); n++ ) {
				delete p_neighbors[n];
			}
		}
		return min_phi;
		// Adams et al. 07
	} else if( method == 1 ) {
		FLOAT64 d = 0.0;
		vec3d a;
		FLOAT64 wsum = 0.0;
		for( uint n=0; n<p_neighbors.size(); n++ ) {
			vec3d x = p_neighbors[n].p;
			FLOAT64 radius = p_neighbors[n].r;
			FLOAT64 w = kernel::smooth_kernel((pos-x).len2(),4*radius)+1e-30;
			d += w*fabs(p_neighbors[n].levelset);
			a += w*x;
			wsum += w;
		}
		if( wsum ) {
			d = d / wsum;
			a = a / wsum;
			return fabs((a-pos).len())-d;
		} else {
			return 1.0;
		}
		// Zhu and Brison 2005
	} else if( method == 2 ) {
		FLOAT64 r = 0.0;
		vec3d a;
		FLOAT64 wsum = 0.0;
		for( uint n=0; n<p_neighbors.size(); n++ ) {
			vec3d x = p_neighbors[n].p;
			FLOAT64 radius = p_neighbors[n].r;
			FLOAT64 w = kernel::smooth_kernel((pos-x).len2(),4*radius)+1e-18;
			r += w*radius;
			a += w*x;
			wsum += w;
		}
		if( wsum ) {
			r = r / wsum;
			a = a / wsum;
			return fabs((a-pos).len())-r;
		} else {
			return 1.0;
		}
		// Yu and Turk 2010
	} else if( method == 3 ) {
		if( anisotropy.empty()) {
			printf( "Compute anisotropy of particles to enable this feature.\n");
			exit(0);
		}		
		FLOAT64 sum = 0.0;
		for( uint n=0; n<p_neighbors.size(); n++ ) {
			vec3d r = p_neighbors[n].p-pos;
			FLOAT64 dens = density[neighbors[n]];
			if( dens ) {
				svd3 svd = anisotropy[neighbors[n]];
				r = util3::stretchPosition(svd,r,0.00001,10000.0);
				sum += (1.0/dens) * kernel::smooth_kernel(r.len2(), kernel_size*p_neighbors[n].r) / powf(dpx,DIM);
			}
		}
		return 0.2*kernel_size-sum;
	}
	return 1e9;
}

void extsurf3::draw() {
	glColor4d(1.0,1.0,1.0,0.1);
	//bcc.drawElements();
	for( int i=0; i<faces.size(); i++ ) {
		glBegin(GL_LINES);
		for( int idx=0; idx<3; idx++ ) {
			int face = faces[i][idx];
			glVertex3dv(vertices[face].v);
		}
		glEnd();
	}
}
