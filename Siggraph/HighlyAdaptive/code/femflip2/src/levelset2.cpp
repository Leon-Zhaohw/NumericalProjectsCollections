/*
 *	levelset2.cpp
 *	
 *	Created by Ryoichi Ando on 11/3/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <math.h>
#include "macros.h"
#include "util2.h"
#include "levelset2.h"
#include "particle2.h"
#include "pcgsolver/util.h"

void levelset2::resize( uint gsize ) {
	this->gsize = gsize;
	dx = 1.0/(gsize-1);
}

vec2d levelset2::evalGradient(vec2d p) const {
	FLOAT64 h = dx;
	return vec2d(0.5*(evalLevelset(vec2d(p[0]+h,p[1]))-evalLevelset(vec2d(p[0]-h,p[1])))/h,
				 0.5*(evalLevelset(vec2d(p[0],p[1]+h))-evalLevelset(vec2d(p[0],p[1]-h)))/h );
}

// Box levelset
FLOAT64 levelset2::box( vec2d p, vec2d p0, vec2d p1 ) {
	FLOAT64 sd = -9999.0;
	sd = fmax(sd,p0[0]-p[0]);
	sd = fmax(sd,p0[1]-p[1]);
	sd = fmax(sd,p[0]-p1[0]);
	sd = fmax(sd,p[1]-p1[1]);
// For outside of the box, following may be turned on for accurate
// outside levelset, which may result in weired animations in our applications
// due to the non-grid-aligned gradients in these area.	
#if 0
	if( p[0]<p0[0] && p[1]<p0[1] ) {
		sd = (p-vec2d(p0[0],p0[1])).len();
	} else if( p[0]<p0[0] && p[1]>p1[1] ) {
		sd = (p-vec2d(p0[0],p1[1])).len();
	} else if( p[0]>p1[0] && p[1]>p1[1] ) {
		sd = (p-vec2d(p1[0],p1[1])).len();
	} else if( p[0]>p1[0] && p[1]<p0[1] ) {
		sd = (p-vec2d(p1[0],p0[1])).len();
	}
#endif
	return sd;
}

// Sphere levelset
FLOAT64 levelset2::circle( vec2d p, vec2d c, FLOAT64 r ) {
	return (p-c).len() - r;
}

void levelset2::marchPoints( FLOAT64 L[2][2], vec2d p[8], int &pnum ) {
    pnum = 0;
	int quads[][2] = { {0, 0}, {1, 0}, {1, 1}, {0, 1} };
	for( int n=0; n<4; n++ ) {
		// Inside domain
		if( L[quads[n][0]][quads[n][1]] <= 0.0 ) {
			p[pnum][0] = quads[n][0];
			p[pnum][1] = quads[n][1];
			pnum ++;
		}
		// If encountered an intersection
		if( L[quads[n][0]][quads[n][1]] * L[quads[(n+1)%4][0]][quads[(n+1)%4][1]] < 0.0 ) {
			// Calculate cross position
			FLOAT64 y0 = L[quads[n][0]][quads[n][1]];
			FLOAT64 y1 = L[quads[(n+1)%4][0]][quads[(n+1)%4][1]];
			FLOAT64 a = y0/(y0-y1);
			FLOAT64 p0[2] = { quads[n][0], quads[n][1] };
			FLOAT64 p1[2] = { quads[(n+1)%4][0], quads[(n+1)%4][1] };
			p[pnum][0] = (1.0-a)*p0[0]+a*p1[0];
			p[pnum][1] = (1.0-a)*p0[1]+a*p1[1];
			pnum ++;
		}
	}	
}

// March levelset with set of points
void levelset2::marchPoints( int i, int j, const array2<FLOAT64> &L, vec2d p[8], int &pnum ) {
	pnum = 0;
	FLOAT64 dx = 1.0/(L.size().w-1);
	FLOAT64 dy = 1.0/(L.size().h-1);
    FLOAT64 LL[2][2];
    int quads_LL[][2] = { {0, 0}, {1, 0}, {1, 1}, {0, 1} };
	int quads_L[][2] = { {i, j}, {i+1, j}, {i+1, j+1}, {i, j+1} };
	for( int n=0; n<4; n++ ) {
		LL[quads_LL[n][0]][quads_LL[n][1]] = L[quads_L[n][0]][quads_L[n][1]];
	}
    marchPoints(LL,p,pnum);
    for( int n=0; n<pnum; n++ ) {
        p[n][0] *= dx;
        p[n][1] *= dy;
        p[n]+= vec2d(dx*i,dy*j);
    }
}

FLOAT64 levelset2::getParticleConvextHullLevelset( vec2d &pos, const particle2 &p0, FLOAT64 dpx, FLOAT64 shrink, bool consistent ) const {
	FLOAT64 r = 0.5*dpx*p0.r;
	r = fmax(r,-shrink*p0.levelset);
	FLOAT64 phi = (p0.p-pos).len()-r;
	vec2d dir = (pos-p0.p).normal();
	if( consistent && dir * p0.gradient < 0.0 ) {
		dir = -1.0 * dir;
		phi += r;
	}
	pos = r*dir+p0.p;
	return phi;
}

FLOAT64 levelset2::getParticleConvextHullLevelset( vec2d &pos, const particle2 &p0, const particle2 &p1, FLOAT64 dpx, FLOAT64 shrink, bool consistent ) const {
	// Line
	FLOAT64 max_phi = -1e9;
	bool out_of_bound = false;
	vec2d offset = (p0.p-p1.p).rotate();
	FLOAT64 x1 = (p0.p+offset)[0];
	FLOAT64 y1 = (p0.p+offset)[1];
	FLOAT64	x2 = (p1.p+offset)[0];
	FLOAT64 y2 = (p1.p+offset)[1];
	pos += offset;
	FLOAT64 r1 = 0.5*dpx*p0.r;
	FLOAT64 r2 = 0.5*dpx*p1.r;
	r1 = fmax(r1,-shrink*p0.levelset);
	r2 = fmax(r2,-shrink*p1.levelset);
	FLOAT64 DET = x2*y1-x1*y2;
	vec2d cpos;
	if( ! (p0.p-p1.p).len2() ) return 1e9;
	if( ! DET ) {
		return 1e9;
	} else {
		// Build a quadric equation
		FLOAT64 A1 = (y2-y1)/DET;
		FLOAT64 B1 = (r2*y1-r1*y2)/DET;
		FLOAT64 A2 = (x1-x2)/DET;
		FLOAT64 B2 = (r1*x2-r2*x1)/DET;
		FLOAT64 det = sqr(A1)*(1.0-sqr(B2))+sqr(A2)*(1.0-sqr(B1))+2.0*A1*A2*B1*B2;
		if( det > 0.0 ) {
			// Solve it !
			for( uint k=0; k<2; k++ ) {
				FLOAT64 d = ((k==0?1:-1)*sqrtf(det)-A1*B1-A2*B2)/(sqr(A1)+sqr(A2));
				vec2d normal(A1*d+B1,A2*d+B2);
				if( consistent && (normal * p0.gradient > 0.0 && normal * p1.gradient > 0.0) ) continue;
				vec2d head0 = (p0.p+offset)-normal*r1;
				vec2d head1 = (p1.p+offset)-normal*r2;
				FLOAT64 dist = -(normal*pos+d);
				vec2d out = pos+normal*dist;
				// Check if out is inside the line
				if( (out-head0)*(out-head1) < 0.0 ) {
					if( dist > max_phi ) {
						max_phi = dist;
						cpos = out;
					}
				} else {
					max_phi = -1e9;
					break;
				}
			}
		}
	}
	pos -= offset;
	cpos -= offset;
	if( max_phi < -1.0 ) out_of_bound = true;
	
	// Sphere
	if( out_of_bound ) {
		out_of_bound = false;
		max_phi = 1e9;
		const particle2 *p[] = { &p0, &p1 };
		for( uint n=0; n<2; n++ ) {
			vec2d fpos = pos;
			FLOAT64 phi = getParticleConvextHullLevelset(fpos,*p[n],dpx,shrink,consistent);
			if( phi < max_phi ) {
				max_phi = phi;
				cpos = fpos;
				out_of_bound = false;
			}
		}
	}
	pos = cpos;
	return out_of_bound ? 1e9 : max_phi;
}

FLOAT64 levelset2::getParticleLevelset( const std::vector<particle2 *> &neighbors, FLOAT64 dpx, vec2d pos, vec2d &outpos, const levelset2 *solid, FLOAT64 shrink, bool consistent ) const {
	// Final neighbors
	std::vector<particle2 *> fneighbors = neighbors;

	// Mirror particles
	std::vector<particle2> mirrors;
	
	// Put mirror particles
	FOR_EACH_PARTICLE(neighbors) {
		FLOAT64 phi_s = solid->evalLevelset(p.p);
		if( p.levelset < 0.0 && phi_s > 0.0 && ! p.isolated && phi_s < dpx*p.r ) {
			particle2 mirror_p = p;
			mirror_p.p = p.p-(0.5*dpx+phi_s)*solid->evalGradient(p.p);
			mirrors.push_back(mirror_p);
		}
	} END_FOR
	for( uint n=0; n<mirrors.size(); n++) fneighbors.push_back(&mirrors[n]);

	// Compute maximum radius
	FLOAT64 r = 0.0;
	for( uint n=0; n<fneighbors.size(); n++ ) {
		r = fmax(r,fneighbors[n]->r);
	}
	
	// Compute solid levelset there
	FLOAT64 solid_phi = solid->evalLevelset(pos);
	
	// Slide position
	int attempt=0;
	FLOAT64 min_dist = 0.0;
	while( solid_phi < min_dist ) {
		vec2d grad = solid->evalGradient(pos);
		if( ! grad.len2()) {
			pos += 1e-4*dpx*vec2d(nrand(),nrand());
		}
		pos += (1e-4*dpx+min_dist-solid_phi)*grad;
		solid_phi = solid->evalLevelset(pos);
		attempt ++;
		if( attempt > 5 ) {
			break;
		}
	}
	
	// Normalize
	vec2d avg_normal;
	FOR_EACH_PARTICLE(fneighbors) {
		avg_normal += p.gradient;
	} END_FOR
	avg_normal.normalize();
	
	FOR_EACH_PARTICLE(fneighbors) {
		if( p.gradient * avg_normal < 0 ) {
			consistent = false;
			break;
		}
	} END_FOR
	
	// Now compute convex hull levelset
	FLOAT64 min_phi = 1e9;
	for( uint n0=0; n0<fneighbors.size(); n0++ ) for( uint n1=n0+1; n1<fneighbors.size(); n1++ ) {
		particle2 &p0 = *fneighbors[n0];
		particle2 &p1 = *fneighbors[n1];
		FLOAT64 r0 = fmax(dpx*p0.r,-shrink*p0.levelset);
		FLOAT64 r1 = fmax(dpx*p1.r,-shrink*p1.levelset);
		
		if( ! p0.isolated && ! p1.isolated &&
			(p0.p-pos).len2() <= sqr(2*r0) &&
			(p1.p-pos).len2() <= sqr(2*r1) &&
		    (p1.p-p0.p).len2() <= sqr(1.24*(r0+r1))) {
			vec2d out = pos;
			FLOAT64 phi = levelset2::getParticleConvextHullLevelset(out,p0,p1,dpx,shrink,consistent);
			if( phi < min_phi ) {
				min_phi = phi;
				outpos = out;
			}
		}
	}
	return min_phi;
}
