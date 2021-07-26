/*
 *	levelset3.cpp
 *	
 *	Created by Ryoichi Ando on 11/26/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "levelset3.h"
#include "matutil.h"
#include "util3.h"
#include "pcgsolver/util.h"
#include "particle3.h"

void levelset3::resize( uint gsize ) {
	this->gsize = gsize;
	dx = 1.0/gsize;
}

vec3d levelset3::evalGradient(vec3d p) const {
	FLOAT64 h = 0.1*dx;
	return vec3d(0.5*(evalLevelset(vec3d(p[0]+h,p[1],p[2]))-evalLevelset(vec3d(p[0]-h,p[1],p[2])))/h,
				 0.5*(evalLevelset(vec3d(p[0],p[1]+h,p[2]))-evalLevelset(vec3d(p[0],p[1]-h,p[2])))/h,
				 0.5*(evalLevelset(vec3d(p[0],p[1],p[2]+h))-evalLevelset(vec3d(p[0],p[1],p[2]-h)))/h );
}

FLOAT64 levelset3::box( vec3d p, vec3d p0, vec3d p1 ) {
	FLOAT64 sd = -9999.0;
	sd = fmax(sd,p0[0]-p[0]);
	sd = fmax(sd,p0[1]-p[1]);
	sd = fmax(sd,p0[2]-p[2]);
	sd = fmax(sd,p[0]-p1[0]);
	sd = fmax(sd,p[1]-p1[1]);
	sd = fmax(sd,p[2]-p1[2]);
	// For outside of the box, following may be turned on for accurate
	// outside levelset, which may result in weired animations in our applications
	// due to the non-grid-aligned gradients in these area.	
#if 0
	if( p[2]<p0[2] ) {
		if( p[0]<p0[0] && p[1]<p0[1] ) {
			sd = (p-vec3d(p0[0],p0[1],p0[2])).len();
		} else if( p[0]<p0[0] && p[1]>p1[1] ) {
			sd = (p-vec3d(p0[0],p1[1],p0[2])).len();
		} else if( p[0]>p1[0] && p[1]>p1[1] ) {
			sd = (p-vec3d(p1[0],p1[1],p0[2])).len();
		} else if( p[0]>p1[0] && p[1]<p0[1] ) {
			sd = (p-vec3d(p1[0],p0[1],p0[2])).len();
		}
	} if( p[2]>p1[2] ) {
		if( p[0]<p0[0] && p[1]<p0[1] ) {
			sd = (p-vec3d(p0[0],p0[1],p1[2])).len();
		} else if( p[0]<p0[0] && p[1]>p1[1] ) {
			sd = (p-vec3d(p0[0],p1[1],p1[2])).len();
		} else if( p[0]>p1[0] && p[1]>p1[1] ) {
			sd = (p-vec3d(p1[0],p1[1],p1[2])).len();
		} else if( p[0]>p1[0] && p[1]<p0[1] ) {
			sd = (p-vec3d(p1[0],p0[1],p1[2])).len();
		}
	}
#endif
	return sd;
}

FLOAT64 levelset3::sphere( vec3d p, vec3d c, FLOAT64 r ) {
	return (p-c).len() - r;
}

FLOAT64 levelset3::getParticleConvextHullLevelset( vec3d &pos, const particle3 &p0, FLOAT64 dpx, FLOAT64 shrink, bool consistent ) const {
	FLOAT64 r = 0.5*dpx*p0.r;
	r = fmax(r,-shrink*p0.levelset);
	FLOAT64 phi = (p0.p-pos).len()-r;
	vec3d dir = (pos-p0.p).normal();
	if( consistent && dir * p0.gradient < 0.0 ) {
		dir = -1.0 * dir;
		phi += r;
	}
	pos = r*dir+p0.p;
	return phi;
}

FLOAT64 levelset3::getParticleConvextHullLevelset( vec3d &pos, const particle3 &p0, const particle3 &p1, FLOAT64 dpx, FLOAT64 shrink, bool consistent ) const {
	vec3d original_pos = pos;
	vec3d positions[5] = { p0.p, p1.p, pos, p0.gradient, p1.gradient };
	
	vec3d nm = ((positions[2]-positions[0])^(positions[1]-positions[0])).normal();
	vec3d e[2];
	e[0] = (positions[1]-positions[0]).normal();
	e[1] = nm ^ e[0];
	
	// Convert to 2D
	vec3d old_points[5];
	for( uint i=0; i<5; i++ ) old_points[i] = positions[i];
	for( uint i=0; i<5; i++ ) {
		positions[i][0] = e[0] * (old_points[i]-old_points[0]);
		positions[i][1] = e[1] * (old_points[i]-old_points[0]);
		positions[i][2] = 0.0;
	}
	
	FLOAT64 max_phi = -1e9;
	bool out_of_bound = false;
	vec3d offset = vec3d((positions[0]-positions[1])[1],-(positions[0]-positions[1])[0],0.0);
	FLOAT64 x1 = (positions[0]+offset)[0];
	FLOAT64 y1 = (positions[0]+offset)[1];
	FLOAT64	x2 = (positions[1]+offset)[0];
	FLOAT64 y2 = (positions[1]+offset)[1];
	pos += offset;
	FLOAT64 r1 = 0.5*dpx*p0.r;
	FLOAT64 r2 = 0.5*dpx*p1.r;
	r1 = fmax(r1,-shrink*p0.levelset);
	r2 = fmax(r2,-shrink*p1.levelset);
	FLOAT64 DET = x2*y1-x1*y2;
	vec3d cpos;
	if( ! (p0.p-p1.p).len2() ) return 1e9;
	if( ! DET ) {
		out_of_bound = true;
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
				FLOAT64 detn = sqr(A1)+sqr(A2);
				if( detn ) {
					FLOAT64 d = ((k==0?1:-1)*sqrtf(det)-A1*B1-A2*B2)/detn;
					vec3d normal(A1*d+B1,A2*d+B2,0.0);
					if( consistent && (normal * positions[3] >= 0.0 && normal * positions[4] >= 0.0) ) continue;
					
					vec3d head0 = (positions[0]+offset)-normal*r1;
					vec3d head1 = (positions[1]+offset)-normal*r2;
					FLOAT64 dist = normal*positions[2]+d;
					vec3d out = positions[2]+normal*dist;
					// Check if out is inside the line
					if( (out-head0)*(out-head1) <= 0.0 ) {
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
	}
	pos = original_pos;
	cpos -= offset;
	
	if( max_phi < -1.0 ) out_of_bound = true;
	else {
		// Convert 2D position to 3D
		cpos = p0.p + e[0]*cpos[0] + e[1]*cpos[1];
	}
	
	if( out_of_bound ) {
		// Sphere
		max_phi = 1e9;
		const particle3 *p[] = { &p0, &p1 };
		for( uint n=0; n<2; n++ ) {
			vec3d fpos = original_pos;
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

FLOAT64 levelset3::getParticleConvextHullLevelset( vec3d &pos, const particle3 &p0, const particle3 &p1, const particle3 &p2, FLOAT64 dpx, FLOAT64 shrink, bool consistent ) const {
	// Triangle volume
	FLOAT64 max_phi = -1e9;
	bool out_of_bound = false;
	vec3d original_pos = pos;
	vec3d offset = (p1.p-p0.p) ^ (p2.p-p0.p);
	FLOAT64 x1 = (p0.p+offset)[0];
	FLOAT64 y1 = (p0.p+offset)[1];
	FLOAT64 z1 = (p0.p+offset)[2];
	FLOAT64	x2 = (p1.p+offset)[0];
	FLOAT64 y2 = (p1.p+offset)[1];
	FLOAT64 z2 = (p1.p+offset)[2];
	FLOAT64	x3 = (p2.p+offset)[0];
	FLOAT64 y3 = (p2.p+offset)[1];
	FLOAT64 z3 = (p2.p+offset)[2];

	pos += offset;
	FLOAT64 r1 = 0.5*dpx*p0.r;
	FLOAT64 r2 = 0.5*dpx*p1.r;
	FLOAT64 r3 = 0.5*dpx*p2.r;
	r1 = fmax(r1,-shrink*p0.levelset);
	r2 = fmax(r2,-shrink*p1.levelset);
	r3 = fmax(r3,-shrink*p2.levelset);
	FLOAT64 DET = -x3*y2*z1+x2*y3*z1+x3*y1*z2-x1*y3*z2-x2*y1*z3+x1*y2*z3;
	vec3d cpos;
	if( ! (p0.p-p1.p).len2() || ! (p1.p-p2.p).len2() || ! (p2.p-p0.p).len2() ) return 1e9;
	if( ! DET ) {
#if 0
		dump( "getParticleConvextHullLevelset3: Denominator was zero !\n" );
		x1 -= offset[0]; x2 -= offset[0]; x3 -= offset[0];
		y1 -= offset[1]; y2 -= offset[1]; y3 -= offset[1];
		z1 -= offset[2]; z2 -= offset[2]; z3 -= offset[2];
		dump( "\nx1=%f y1=%f z1=%f\nx2=%f y2=%f z2=%f\nx3=%f y3=%f z3=%f\n", x1,y1,z1,x2,y2,z2,x3,y3,z3 );
		exit(0);
#endif
	} else {
		// Build a quadric equation
		FLOAT64 A1 = (-y3*z1-y1*z2+y3*z2+y2*(z1-z3)+y1*z3)/DET;
		FLOAT64 B1 = (-r3*y2*z1+r2*y3*z1+r3*y1*z2-r1*y3*z2-r2*y1*z3+r1*y2*z3)/DET;
		FLOAT64 A2 = (x3*z1+x1*z2-x3*z2-x1*z3+x2*(-z1+z3))/DET;
		FLOAT64 B2 = (r3*x2*z1-r2*x3*z1-r3*x1*z2+r1*x3*z2+r2*x1*z3-r1*x2*z3)/DET;
		FLOAT64 A3 = (-x3*y1-x1*y2+x3*y2+x2*(y1-y3)+x1*y3)/DET;
		FLOAT64 B3 = (-r3*x2*y1+r2*x3*y1+r3*x1*y2-r1*x3*y2-r2*x1*y3+r1*x2*y3)/DET;
		FLOAT64 detn = sqr(A1)+sqr(A2)+sqr(A3);
		FLOAT64 det = 4.0*sqr(A1*B1+A2*B2+A3*B3)-4.0*detn*(sqr(B1)+sqr(B2)+sqr(B3)-1.0);
		if( det > 0.0 && detn ) {
			// Solve it !
			for( uint k=0; k<2; k++ ) {
				FLOAT64 d = ((k==0?1:-1)*0.5*sqrtf(det)-A1*B1-A2*B2-A3*B3)/detn;
				vec3d normal(A1*d+B1,A2*d+B2,A3*d+B3);
				if( consistent && (normal * p0.gradient > 0.0 && normal * p1.gradient > 0.0 && normal * p2.gradient > 0.0) ) continue;
				vec3d head0 = (p0.p+offset)-normal*r1;
				vec3d head1 = (p1.p+offset)-normal*r2;
				vec3d head2 = (p2.p+offset)-normal*r3;
				FLOAT64 dist = -(normal*pos+d);
				vec3d out = pos+normal*dist;
				// Check if the out is inside the triangle
				FLOAT64 cross1 = ((head1-head0)^(out-head0)) * normal;
				FLOAT64 cross2 = ((head2-head1)^(out-head1)) * normal;
				FLOAT64 cross3 = ((head0-head2)^(out-head2)) * normal;
				if( (cross1 >= 0.0 && cross2 >= 0.0 && cross3 >= 0.0) || (cross1 <= 0.0 && cross2 <= 0.0 && cross3 <= 0.0) ) {
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
	pos = original_pos;
	cpos -= offset;
	if( max_phi < -1.0 ) out_of_bound = true;
	if( out_of_bound ) {
		// Cylineder
		max_phi = 1e9;
		const particle3 *p[] = { &p0, &p1, &p2 };
		for( uint n0=0; n0<3; n0++ ) for( uint n1=n0+1; n1<3; n1++ ) {
			vec3d fpos = pos;
			FLOAT64 phi = getParticleConvextHullLevelset(fpos,*p[n0],*p[n1],dpx,shrink,consistent);
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

static FLOAT64 adamKernel( FLOAT64 len2, FLOAT64 radius ) {
	return powf(fmax(0.0,1.0-(len2/(radius*radius))),3.0);
}

FLOAT64 levelset3::getAdams07ParticleLevelset( const std::vector<particle3 *> &neighbors, FLOAT64 dpx, vec3d pos, vec3d &outpos, const levelset3 *solid, FLOAT64 shrink ) const {
	FLOAT64 d = 0.0;
	vec3d gradient;
	vec3d a;
	FLOAT64 wsum = 0.0;
	for( uint n=0; n<neighbors.size(); n++ ) {
		vec3d x = neighbors[n]->p;
		FLOAT64 radius = 0.5*dpx*neighbors[n]->r;
		FLOAT64 w = adamKernel((pos-x).len2(),radius);
		d += fabs(neighbors[n]->levelset);
		a += w*x;
		gradient += w*neighbors[n]->gradient;
		wsum += w;
	}
	if( wsum ) {
		d = d / wsum;
		a = a / wsum;
		gradient.normalize();
		
		outpos = a + gradient*d;
		return d - (a-pos).len();
	} else {
		return 1e9;
	}
}

FLOAT64 levelset3::getParticleLevelset( const std::vector<particle3 *> &neighbors, FLOAT64 dpx, vec3d pos, vec3d &outpos, const levelset3 *solid, FLOAT64 shrink, bool consistent ) const {
	// Final neighbors
	std::vector<particle3 *> fneighbors = neighbors;
	
	// Mirror particles
	std::vector<particle3> mirrors;
	
	// Put mirror particles
	FOR_EACH_PARTICLE(neighbors) {
		FLOAT64 phi_s = solid->evalLevelset(p.p);
		if( p.levelset < dpx*p.r && phi_s > 0.0 && ! p.isolated && phi_s < dpx*p.r ) {
			particle3 mirror_p = p;
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
	FLOAT64 min_dist = 0.65*dpx*r;
	while( solid_phi < min_dist ) {
		vec3d grad = solid->evalGradient(pos);
		if( ! grad.len2()) {
			pos += 1e-4*dpx*vec3d(nrand(),nrand(),nrand());
		}
		pos += (1e-4*dpx+min_dist-solid_phi)*grad;
		solid_phi = solid->evalLevelset(pos);
		attempt ++;
		if( attempt > 5 ) {
			break;
		}
	}
	
	// Normalize
	vec3d avg_normal;
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
	FLOAT64 min_phi = 1e9;
	
	for( uint n0=0; n0<fneighbors.size(); n0++ ) for( uint n1=n0+1; n1<fneighbors.size(); n1++ ) for( uint n2=n1+1; n2<fneighbors.size(); n2++ ) {
		particle3 &p0 = *fneighbors[n0];
		particle3 &p1 = *fneighbors[n1];
		particle3 &p2 = *fneighbors[n2];
		FLOAT64 r0 = fmax(dpx*p0.r,-shrink*p0.levelset);
		FLOAT64 r1 = fmax(dpx*p1.r,-shrink*p1.levelset);
		FLOAT64 r2 = fmax(dpx*p2.r,-shrink*p2.levelset);
		particle3 *p[DIM] = { &p0, &p1, &p2 };
		FLOAT64 r[DIM] = { r0, r1, r2 };
		
		if( ! p0.isolated && ! p1.isolated && ! p2.isolated &&
		   (p0.p-pos).len2() <= sqr(2*r0) &&
		   (p1.p-pos).len2() <= sqr(2*r1) &&
		   (p2.p-pos).len2() <= sqr(2*r2) &&
		   (p0.p-p1.p).len2() <= sqr(r0+r1) &&
		   (p1.p-p2.p).len2() <= sqr(r1+r2) &&
		   (p2.p-p0.p).len2() <= sqr(r2+r0) ) {
			vec3d out = pos;
			FLOAT64 phi = levelset3::getParticleConvextHullLevelset(out,p0,p1,p2,dpx,shrink,consistent);
			if( phi < min_phi ) {
				min_phi = phi;
				outpos = out;
			}			
		} else {
			for( uint n=0; n<DIM; n++ ) {
				uint m = (n+1) % DIM;
				if( ! p[n]->isolated && ! p[m]->isolated &&
				    (p[n]->p-pos).len2() <= sqr(2*r[n]) &&
				    (p[m]->p-pos).len2() <= sqr(2*r[m]) &&
				    (p[m]->p-p[n]->p).len2() <= sqr(r[n]+r[m]) ) {
					vec3d out = pos;
					FLOAT64 phi = levelset3::getParticleConvextHullLevelset(out,*p[n],*p[m],dpx,shrink,consistent);
					if( phi < min_phi ) {
						min_phi = phi;
						outpos = out;
					}
				} else {
					particle3 *p_2[2] = { p[n], p[m] };
					FLOAT64 r_2[2] = { r[n], r[m] };
					for( uint k=0; k<2; k++ ) {
						if( ! p_2[k]->isolated && (p_2[k]->p-pos).len2() <= sqr(2*r_2[k]) ) {
							vec3d out = pos;
							FLOAT64 phi = levelset3::getParticleConvextHullLevelset(out,*p_2[k],dpx,shrink,consistent);
							if( phi < min_phi ) {
								min_phi = phi;
								outpos = out;
							}
						}
					}
				}
			}
		}
	}
	return min_phi;
}