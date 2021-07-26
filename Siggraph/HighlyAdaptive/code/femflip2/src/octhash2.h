/*
 *	octhash2.h
 *
 *	Created by Ryoichi Ando on 5/31/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "octree2.h"
#include <stdlib.h>
#include <vector>
#include <list>
#ifndef _OCTHASH2_H
#define _OCTHASH2_H

template <class T>
class octhash2 {
public:
	// Point structure
	typedef struct {
		T data;
		vec2d p;
	} point2;
	
	// Hash table
	std::vector<std::list<point2> > hash;
	
	// Octree
	const octree2 *octree;
	
	// Initialize the octree structure (can be called multiple times)
	void init( const octree2 *octree ) {	
		clear();
		hash.resize(octree->terminals.size());
		this->octree = octree;
	}
	
	// Clear hash table
	void clear() {
		for( uint n=0; n<hash.size(); n++ ) hash[n].clear();
	}
	
	// Add a point
	void push( vec2d p, T data ) {
		int index = octree->hitTest(p);
		if( index >= 0 ) {
			point2 point = { data, p };
			hash[index].push_back(point);
		}
	}
	
	// Remove a point
	void pop( vec2d p, T data ) {
		std::list<T> neighbors;
		int index = octree->hitTest(p);
		class std::list<point2>::iterator it;
		for( it=hash[index].begin(); it!=hash[index].end(); it++ ) {
			if( (*it).data == data ) {
				hash[index].erase(it);
				return;
			}
		}
	}
	
	// Collect neighbor points
	std::vector<T> getNeighbors( vec2d p, FLOAT64 r ) const {
		std::vector<T> neighbors;
		std::vector<uint> array;
		octree->hitTest(p,r,array);
		std::vector<uint>::iterator it;
		for( it=array.begin(); it!=array.end(); it++ ) {
			uint index = *it;
			typename std::list<point2>::const_iterator it2;
			for( it2=hash[index].begin(); it2!=hash[index].end(); it2++ ) {
				if( ((*it2).p-p).len2() < r*r ) {
					neighbors.push_back((*it2).data);
				}
			}
		}
		return neighbors;
	}
};

#endif