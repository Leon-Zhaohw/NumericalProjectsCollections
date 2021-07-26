/*
 *	annsorter2.h
 *
 *	Created by Ryoichi Ando on 9 Nov 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "sorter2.h"
#include "ann2.h"

#ifndef _ANNSORTER2_H
#define _ANNSORTER2_H

class particle2;
class annsorter2 : public sorter2 {
public:
	annsorter2();
	virtual void sortParticles( std::vector<particle2 *> &particle );
	virtual std::vector<particle2 *> getNeighbors( vec2d p, FLOAT64 r ) const;
	virtual std::vector<particle2 *> getkNeighbors( vec2d p, uint n ) const;
	virtual void setDirty( bool isDirty=true );
protected:
	ann2 ann;
	std::vector<particle2 *> particles;
	bool dirty;
};

#endif