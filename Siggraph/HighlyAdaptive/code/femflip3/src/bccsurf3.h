/*
 *	bccsurf3.h
 *
 *	Created by Ryoichi Ando on 8/21/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "meshsurf3.h"
#include "bccmesher3.h"

#ifndef _BCCSURF3_H
#define _BCCSURF3_H

class bccsurf3 : public levelset3 {
public:
	bccsurf3();
	bccsurf3( const bccsurf3 &bccsurf );
	void operator=( const bccsurf3 &bccsurf );
	void setBCC( const bccmesher3 &g );
	void setExtrapolationDist( FLOAT64 dist );
	virtual void buildSurface( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, const levelset3 *fluid, const levelset3 *solid, bool enclose);
	virtual void drawGL( const levelset3 *solid, uint kind ) const;
	virtual FLOAT64 evalLevelset(vec3d p) const;
	virtual vec3d evalGradient(vec3d p) const;

	meshsurf3 meshsurf;
	bccmesher3 g;
};

#endif