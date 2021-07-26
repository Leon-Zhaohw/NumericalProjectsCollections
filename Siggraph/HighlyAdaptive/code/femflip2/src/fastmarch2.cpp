/*
 *	fastmarch2.cpp
 *
 *	Created by Ryoichi Ando on 7/15/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "util2.h"
#include "fastmarch2.h"
#include "matutil.h"
#include <algorithm>
#include <list>
#include "pcgsolver/util.h"

template<typename T> inline FLOAT64 length (T value) {
	return fabs(value);
}

template <> inline FLOAT64 length (vec2d value) {
	return value.len();
}

template<typename T> inline T normalize (T value, FLOAT64 len) {
	return copysign(1.0,value) * len;
}

template <> inline vec2d normalize (vec2d value, FLOAT64 len) {
	return value.normal() * len;
}

template <class T> bool fastmarch2<T>::fastMarch( std::vector<node2 *> &nodes, FLOAT64 maxdist, FLOAT64 mindist, char type ) {
	// Check the sign of coefficients
	if( mindist > 0.0 ) {
		email::print("mindist must be negative.\n" );
		email::send();
		exit(0);
	}
	if( maxdist < 0.0 ) {
		email::print("maxdist must be positive.\n" );
		email::send();
		exit(0);
	}
	
	// Initialize
	for( uint n=0; n<nodes.size(); n++ ) {
		node2 *node = nodes[n];
		node->new_value = T();
		node->new_levelset = 0.0;
	}
	
	// Gather unfixed nodes
	std::list<node2 *> unfixed;
	for( uint n=0; n<nodes.size(); n++ ) {
		node2 *node = nodes[n];
		if( node->p2p.empty() ) continue;
		if( ! node->fixed ) {
			if( type == 1 ) node->levelset = (node->levelset > 0.0 ? maxdist : mindist);
			unfixed.push_back(node);
		}
	}
		
	// Now repeat the propagation...
	while( true ) {
		// Gather narrow band nodes
		std::list<node2 *> narrowList;
		for( typename std::list<node2 *>::iterator it=unfixed.begin(); it!=unfixed.end(); it++ ) {
			node2 *node = *it;
			if( node->p2p.empty() ) continue;
			if( ! node->fixed ) {
				// Sort order of connections
				std::sort(node->p2p.begin(),node->p2p.end(),node2());
				// Collect narrow bands
				node2 *neigh = node->p2p[0];
				if( neigh->fixed ) {
					if( neigh->levelset < maxdist && neigh->levelset > mindist ) {
						narrowList.push_back(node);
					}
				}
			}
		}
		
		// If not found, just leave the loop
		if( narrowList.empty() ) break;
		
		// Find the minimum edge length and min distance
		FLOAT64 ds = 1.0;
		FLOAT64 dist = 1.0;
		for( typename std::list<node2 *>::iterator it=narrowList.begin(); it!=narrowList.end(); it++ ) {
			node2 *node = *it;
			node2 *neigh = node->p2p[0];
			if( neigh->fixed ) {
				ds = fmin(ds,(node->p-neigh->p).len());
				dist = fmin(dist,fabs(neigh->levelset));
			}
		}
		
		// Cut out the narrow bands to regularize the propagation speed
		for( typename std::list<node2 *>::iterator it=narrowList.begin(); it!=narrowList.end(); ) {
			node2 *node = *it;
			if( node->p2p.empty() ) continue;
			node2 *neigh = node->p2p[0];
			if( fabs(neigh->levelset) > dist+ds ) {
				it = narrowList.erase(it);
			} else {
				it ++;
			}
		}
		
		// Tranfer to vector container
		std::vector<node2 *> narrowNodes;
		narrowNodes.insert(narrowNodes.end(),narrowList.begin(),narrowList.end());

		// Propagate once
		PARALLEL_FOR for( uint n=0; n<narrowNodes.size(); n++ ) {
			node2 *node = narrowNodes[n];
			if( node->p2p.empty() ) continue;
			
			// Pick neighboring nodes
			std::vector<node2 *> tri(node->p2p.size()+1); tri[0] = node;
			for( uint i=0; i<node->p2p.size(); i++ ) tri[1+i] = node->p2p[i];
			
			// Find the number of valid connections
			uint numValid=0;
			if( node->p2p.size() > 2 && tri[1]->fixed && tri[2]->fixed ) numValid = 2;
			else if( node->p2p.size() > 1 && tri[1]->fixed ) numValid = 1;
			
			// Compute shape function if necessary
			FLOAT64 M[NUM_VERT][NUM_VERT];
			if( numValid == 2 ) {
				FLOAT64 A[NUM_VERT][NUM_VERT];
				for( uint i=0; i<numValid+1; i++ ) for( uint j=0; j<numValid+1; j++ ) {
					if( i < numValid ) A[i][j] = tri[j]->p[i];
					else A[i][j] = 1.0;
				}
				if( ! invert3x3(A,M)) numValid = 1;
			}
			
			if( type == 0 ) {
				// Quantity extrapolation				
				// Compute neighbors absolute value
				FLOAT64 max_len = 0.0;
				FLOAT64 min_len = 1e18;
				for( uint m=0; m<node->p2p.size(); m++ ) {
					if( node->p2p[m]->fixed ) {
						max_len = fmax(max_len,length(node->p2p[m]->value));
						min_len = fmin(min_len,length(node->p2p[m]->value));
					}
				}
				
				if( numValid == 2 ) {
					// Compute gradient of levelset
					vec2d gradient;
					for( uint dim=0; dim<numValid; dim++ ) {
						for( uint k=0; k<numValid+1; k++ ) gradient[dim] += tri[k]->levelset*M[k][dim];
					}
					// Compute determinant
					FLOAT64 det = 0.0;
					for( uint dim=0; dim<numValid; dim++ ) {
						det += gradient[dim]*M[0][dim];
					}
					// Check that gradient is correct
					FLOAT64 deviation = fabs(gradient.len2()-1.0);
					if( deviation < 1e-3 ) { // If below the tolerance...
						// Compute right hand side
						T rhs = T();
						for( uint i=0; i<numValid; i++ ) for( uint j=1; j<numValid+1; j++ ) {
							rhs += gradient[i]*(-M[j][i])*tri[j]->value;
						}
						// Compute the value
						if( det ) {
							node->new_value = rhs / det;
						} else {
							node->new_value = 0.5*(tri[1]->value+tri[2]->value);
							printf( "Extrapolation failed. Determination was zero.\n" );
						}
					} else {
						// Just copy averaged one
						node->new_value = 0.5*(tri[1]->value+tri[2]->value);
					}
				} else if( numValid == 1 ) {
					// If only one neighbor is fixed, just copy that one
					node->new_value = tri[1]->value;
				} else {
					node->value = T();
				}
				
				// Normalize if necessary
				FLOAT64 len = length(node->new_value);
				if( len > max_len ) node->new_value = normalize(node->new_value, max_len);
				if( len < min_len ) node->new_value = normalize(node->new_value, min_len);
			} else if( type == 1 ) {
				// Levelset extrapolation
				int sgn = tri[0]->levelset > 0.0 ? 1 : -1;
				if( numValid == 2 ) {
					// Build quadric equation
					vec2d det;
					vec2d coef;
					for( uint dim=0; dim<numValid; dim++ ) {
						det[dim] = M[0][dim];
						coef[dim] = 0.0;
						for( uint k=1; k<numValid+1; k++ ) {
							coef[dim] += M[k][dim]*tri[k]->levelset;
						}
					}
					// Compute quadric coefficients
					FLOAT64 A = det.len2();
					FLOAT64 B = 2.*det*coef;
					FLOAT64 C = coef.len2()-1.0;
					if( A ) {
						FLOAT64 D = B/A;
						node->new_levelset = sgn*0.5*sqrtf(fmax(1e-8,D*D-4.0*C/A))-0.5*D;
					} else {
						email::print( "determinant was zero !\n" );
						email::send();
						exit(0);
					}
				} else if( numValid == 1 ) {
					// If only one neighbor is fixed, just copy that one
					node->new_levelset = tri[1]->levelset+sgn*(tri[1]->p-tri[0]->p).len();
				} else {
					node->new_levelset = 0.0;
				}
			}
		}
		// Fix the narrow bands
		PARALLEL_FOR for( uint n=0; n<narrowNodes.size(); n++ ) {
			node2 *node = narrowNodes[n];
			if( node->p2p.empty() ) continue;
			if( type == 0 ) node->value = node->new_value;
			else if( type == 1 ) node->levelset = fmax(mindist,fmin(maxdist,node->new_levelset));
			node->fixed = true;
		}
	}
	return true;
}

template class fastmarch2<FLOAT64>;
template class fastmarch2<vec2d>;