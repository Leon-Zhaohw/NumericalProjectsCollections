/*
 *	bccsurf3.cpp
 *
 *	Created by Ryoichi Ando on 2012/08/04
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */
 
#include "bccsurf3.h"
#include "util3.h"
#include "matutil.h"
#include "fastmarch3.h"
#include <vector>

bccsurf3::bccsurf3() {
	meshsurf.setReference(&this->g);
}

bccsurf3::bccsurf3( const bccsurf3 &bccsurf ) {
	*this = bccsurf;
}

void bccsurf3::operator=( const bccsurf3 &bccsurf ) {
	tick(); dump(">>> Copying BCC surface...\n");
	this->g = bccsurf.g;
	tick(); dump("Copying mesh surface...");
	this->meshsurf = bccsurf.meshsurf;
	meshsurf.setReference(&this->g);
	dump("Done. Took %s.\n",stock());
	dump("<<< Done. Took %s.\n",stock());
}

void bccsurf3::setBCC( const bccmesher3 &g ) {
	this->g = g;
}

void bccsurf3::setExtrapolationDist( FLOAT64 dist ) {
	meshsurf.setExtrapolationDist(dist);
}

void bccsurf3::buildSurface( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, const levelset3 *fluid, const levelset3 *solid, bool enclose ) {
	meshsurf.buildSurface(vertices,normals,faces,fluid,solid,enclose);
}

void bccsurf3::drawGL( const levelset3 *solid, uint kind ) const {
	meshsurf.drawGL(solid,kind);
}

FLOAT64 bccsurf3::evalLevelset(vec3d p) const {
	if( meshsurf.g != &this->g ) {
		printf( "Mesh reference is not correct.\n" );
		exit(0);
	}
	return meshsurf.evalLevelset(p);
}

vec3d bccsurf3::evalGradient(vec3d p) const {
	return meshsurf.evalGradient(p);
}