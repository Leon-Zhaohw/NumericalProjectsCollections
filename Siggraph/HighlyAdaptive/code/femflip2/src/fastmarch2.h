/*
 *	fastmarch2.h
 *
 *	Created by Ryoichi Ando on 7/15/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "vec2.h"
#include <vector>

#ifndef _FASTMARCH2_H
#define _FASTMARCH2_H

template <class T> class fastmarch2 {
public:
	typedef struct _node2 {
		vec2d p;					// Position
		FLOAT64 levelset;			// Levelset (give sign info as input for levelset extrapolation)
		T value;					// Value (value to be extrapolated with given levelset value)
		bool fixed;					// Fixed flag
		std::vector<_node2 *> p2p;	// Connection
		
		///////// Just ignore below ///////////
		FLOAT64 new_levelset;
		T new_value;
		bool operator()(const _node2* lhs, const _node2* rhs) const {
			return fabs(lhs->levelset) < fabs(rhs->levelset);
		}
	} node2;
	
	// type = 0 -> quantity extrapolation
	// type = 1 -> levelset extrapolation
	static bool fastMarch( std::vector<node2 *> &nodes, FLOAT64 maxdist, FLOAT64 mindist, char type );
};

#endif