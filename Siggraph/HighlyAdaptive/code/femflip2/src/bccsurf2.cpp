/*
 *	bccsurf2.cpp
 *
 *	Created by Ryoichi Ando on 7/17/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "bccsurf2.h"
#include "util2.h"
#include "matutil.h"
#include "fastmarch2.h"
#include <vector>

bccsurf2::bccsurf2() {
	meshsurf.setReference(&this->g);
}

bccsurf2::bccsurf2( const bccsurf2 &bccsurf ) {
	*this = bccsurf;
}

void bccsurf2::operator=( const bccsurf2 &bccsurf ) {
	this->g = bccsurf.g;
	this->meshsurf = bccsurf.meshsurf;
	meshsurf.setReference(&this->g);
}

void bccsurf2::setBCC( const bccmesher2 &g ) {
	this->g = g;
}

void bccsurf2::setExtrapolationDist( FLOAT64 dist ) {
	meshsurf.setExtrapolationDist(dist);
}

void bccsurf2::buildSurface( std::vector<vec2d> &vertices, std::vector<vec2d> &normals, std::vector<vec2i> &faces, const levelset2 *fluid, const levelset2 *solid, bool enclose) {
	meshsurf.buildSurface(vertices,normals,faces,fluid,solid,enclose);
}

void bccsurf2::drawGL( const levelset2 *solid, uint kind ) const {
	meshsurf.drawGL(solid,kind);
}

FLOAT64 bccsurf2::evalLevelset(vec2d p) const {
	return meshsurf.evalLevelset(p);
}

vec2d bccsurf2::evalGradient(vec2d p) const {
	return meshsurf.evalGradient(p);
}

FLOAT64 bccsurf2::evalCurvature(vec2d p) const {
	return meshsurf.evalCurvature(p);
}