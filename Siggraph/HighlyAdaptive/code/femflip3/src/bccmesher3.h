/*
 *	bccmesher3.h
 *
 *	Created by Ryoichi Ando on 5/29/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "mesher3.h"
#include "bcc3.h"

#ifndef _BCCMESHER3_H
#define _BCCMESHER3_H

class bccmesher3 : public mesher3 {
public:
	virtual void generateElements( std::vector<vec3d> &nodes, std::vector<std::vector<uint> > &elements, const levelset3 *hint, uint gn );
	virtual int	 hitElements(const vec3d p) const;
	octree3 octree;
	bcc3 bcc;
};

#endif