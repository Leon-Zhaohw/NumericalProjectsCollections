/*
 *	adaptive2.h
 *	
 *	Created by Ryoichi Ando on 11/9/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "sorter2.h"
#include <vector>

#ifndef _ADAPTIVE2_H
#define _ADAPTIVE2_H

class particle2;
class levelset2;
class adaptive2 {
public:
	void merge( sorter2& sorter, const levelset2 *level, const levelset2 *surface, std::vector<particle2 *> &particles, FLOAT64 dpx, FLOAT64 dx );
    void split( sorter2& sorter, const levelset2 *level, const levelset2 *surface, std::vector<particle2 *> &particles, FLOAT64 dpx, FLOAT64 dx, bool splitAll=false );
private:
    particle2 *findClosestParticle( const sorter2& sorter, const particle2 &p, FLOAT64 r );
};

#endif