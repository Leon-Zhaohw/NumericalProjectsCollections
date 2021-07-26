/*
 *	fastmarch3.h
 *
 *	Created by Ryoichi Ando on 8/15/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "vec3.h"
#include <vector>

#ifndef _FASTMARCH3_H
#define _FASTMARCH3_H

template <class T> class fastmarch3 {
public:
	typedef struct _node3 {
		vec3d p;					// Position
		FLOAT64 levelset;			// Levelset (give sign info as input for levelset extrapolation)
		T value;					// Value (value to be extrapolated with given levelset value)
		bool fixed;					// Fixed flag
		std::vector<_node3 *> p2p;	// Connection
		
		///////// Just ignore below ///////////
		FLOAT64 new_levelset;
		T new_value;
		bool operator()(const _node3* lhs, const _node3* rhs) const {
			return fabs(lhs->levelset) < fabs(rhs->levelset);
		}
	} node3;
	
	// type = 0 -> quantity extrapolation
	// type = 1 -> levelset extrapolation
	static bool fastMarch( std::vector<node3 *> &nodes, FLOAT64 maxdist, FLOAT64 mindist, char type );
};

#endif