/*
 *	octlevelset3.h
 *
 *	Created by Ryoichi Ando on 9/7/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "octgrid3.h"
#include "levelset3.h"

#ifndef _OCTLEVELSET3_H
#define _OCTLEVELSET3_H

class octlevelset3 : public levelset3 {
public:
	octlevelset3() {}
	virtual void setOctree( const octree3 &octree ) {
		this->octree = octree;
		octgrid.setOctree(this->octree);
	}
	virtual void setNodalLevelset( const std::vector<FLOAT64> &q ) {
		octgrid.setValue(q);
	}
	virtual FLOAT64 evalLevelset(vec3d p) const {
		return octgrid.getValue(p);
	}
	octree3 octree;
	octgrid3<FLOAT64> octgrid;
};


#endif