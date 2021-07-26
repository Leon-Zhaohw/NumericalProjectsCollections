/*
 *	sorter2.h
 *
 *	Created by Ryoichi Ando on 7/22/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <vector>
#include "vec2.h"

#ifndef _SORTER2_H
#define _SORTER2_H

class particle2;
class sorter2 {
public:
	virtual void sortParticles( std::vector<particle2 *> &particle ) = 0;
	virtual std::vector<particle2 *> getNeighbors( vec2d p, FLOAT64 r ) const = 0;
	virtual std::vector<particle2 *> getkNeighbors( vec2d p, uint n ) const = 0;
	virtual void setDirty( bool isDirty=true ) = 0;
};

#endif