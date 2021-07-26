/*
 *	annsorter2.h
 *
 *	Created by Ryoichi Ando on 9 Nov 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "sorter3.h"
#include "ANN.h"

#ifndef _ANNSORTER3_H
#define _ANNSORTER3_H

class particle3;
class annsorter3 : public sorter3 {
public:
	annsorter3();
	virtual void sortParticles( std::vector<particle3 *> &particle );
	virtual std::vector<particle3 *> getNeighbors( vec3d p, FLOAT64 r ) const;
	virtual std::vector<particle3 *> getkNeighbors( vec3d p, uint n ) const;
	virtual void setDirty( bool isDirty=true );
protected:
	uint numbers;
	ANNpointArray dataPts;
	ANNkd_tree * kdtree;
	std::vector<particle3 *> particles;
	bool dirty;
};

#endif