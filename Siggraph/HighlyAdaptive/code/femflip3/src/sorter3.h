/*
 *	sorter3.h
 *
 *	Created by Ryoichi Ando on 8/21/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <vector>
#include "vec3.h"

#ifndef _SORTER3_H
#define _SORTER3_H

class particle3;
class sorter3 {
public:
	virtual void sortParticles( std::vector<particle3 *> &particle ) = 0;
	virtual std::vector<particle3 *> getNeighbors( vec3d p, FLOAT64 r ) const = 0;
	virtual std::vector<particle3 *> getkNeighbors( vec3d p, uint n ) const = 0;
	virtual void setDirty( bool isDirty=true ) = 0;
};

#endif