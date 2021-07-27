/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Basic quaternion class
 *
 *****************************************************************************/

#ifndef UTIL_QUATERIONS_H
#define UTIL_QUATERIONS_H

#include "matrixbase.h"
#include "vectorbase.h"


#define smax(a,b) ((a>b)?(a):(b))
#define scpysign(a,b) ((a)*(b)>=0?(a):-(a))

namespace DDF {

template<class Scalar>
class Quaternion 
{
public:
	Scalar x,y,z,w;

	Quaternion () : x(0),y(0),z(0),w(0) {}
	
	Quaternion (Scalar nx, Scalar ny, Scalar nz, Scalar nw) : 
		x(nx),y(ny),z(nz),w(nw) {}
	
	Quaternion (const Quaternion &q)
	{
		x=q.x; y=q.y; z=q.z; w=q.w;
	}

	Quaternion (const ntlMatrix4x4<Scalar> &M)
	{
		w = sqrt( smax( 0, 1 + M.value[0][0] + M.value[1][1] + M.value[2][2] ) ) / 2;
		x = sqrt( smax( 0, 1 + M.value[0][0] - M.value[1][1] - M.value[2][2] ) ) / 2;
		y = sqrt( smax( 0, 1 - M.value[0][0] + M.value[1][1] - M.value[2][2] ) ) / 2;
		z = sqrt( smax( 0, 1 - M.value[0][0] - M.value[1][1] + M.value[2][2] ) ) / 2;
		x = scpysign( x, M.value[2][1] - M.value[1][2] );
		y = scpysign( y, M.value[0][2] - M.value[2][0] );
		z = scpysign( z, M.value[1][0] - M.value[0][1] );
	}

	Quaternion (const ntlVector3Dim<Scalar>& axis, Scalar angle)
	{
		Scalar mult = sin(angle*.5);
		x = axis.x * mult;
		y = axis.y * mult;
		z = axis.z * mult;
		w = cos (angle*.5);
	}

	Quaternion (Scalar rx, Scalar ry, Scalar rz)
	{
		Quaternion qx(ntlVector3Dim<Scalar> (1.,0.,0.),-rx);
		Quaternion qy(ntlVector3Dim<Scalar> (0.,1.,0.),-ry);
		Quaternion qz(ntlVector3Dim<Scalar> (0.,0.,1.),-rz);
		Quaternion q = qz*qy*qx;
		x=q.x; y=q.y; z=q.z; w=q.w;
	}

	ntlMatrix4x4<Scalar> getRotMat() const
	{
		const Scalar x2 = x*x;
		const Scalar y2 = y*y;
		const Scalar z2 = z*z;
		// const float w2 = w*w;
		
		ntlMatrix4x4<Scalar> M;
		M.initId();
		M.value[0][0] = 1 - 2*(y2+z2);
		M.value[0][1] = 2*(x*y - w*z);
		M.value[0][2] = 2*(x*z+w*y);
		M.value[1][0] = 2*(x*y+w*z);
		M.value[1][1] = 1-2*(x2+z2);
		M.value[1][2] = 2*(y*z-w*x);
		M.value[2][0] = 2*(x*z-w*y);
		M.value[2][1] = 2*(y*z+w*x);
		M.value[2][2] = 1-2*(x2+y2);
		return M;
	}

	ntlVector3Dim<Scalar> getAxis() const
	{
		Scalar phi2 = acos(w);
		ntlVector3Dim<Scalar> axis = ntlVector3Dim<Scalar> (x,y,z) * (1./sin(phi2));
		normalize(axis);
		return axis * 2.* phi2;
	}
	
	inline const Quaternion operator+(const Quaternion &q) const { return Quaternion(x+q.x,y+q.y,z+q.z,w+q.w);		};
	inline const Quaternion operator-(const Quaternion &q) const { return Quaternion(x-q.x,y-q.y,z-q.z,w-q.w);		};
	inline const Quaternion operator*(const Scalar m) const {	return Quaternion(x*m,y*m,z*m,w*m);		};
	inline const Quaternion operator-() const { return Quaternion(-x,-y,-z,-w);	 };
	inline Scalar dot(const Quaternion &q) const { return x*q.x+y*q.y+z*q.z+w*q.w; }
	inline Scalar normSq() const { return x*x+y*y+z*z+w*w; }
	inline Scalar norm() const { return sqrt(normSq()); }
	inline const Quaternion unit() const { Scalar d=1./norm(); return Quaternion(x*d,y*d,z*d,w*d); }
		
	inline const Quaternion operator*(const Quaternion &q) const
	{
		ntlVector3Dim<Scalar> v1(x,y,z), v2(q.x,q.y,q.z);
		ntlVector3Dim<Scalar> nv = v1*q.w + v2*w + cross(v2,v1);
		Scalar nw = w*q.w - (v1.x*v2.x+v1.y*v2.y+v1.z*v2.z);

		return Quaternion(nv.x,nv.y,nv.z,nw);
	}

	static const Quaternion slerp(Quaternion q, Quaternion p, Scalar t)
	{
		Scalar cosphi = q.dot(p);

		if(cosphi < 0.0f)
		{
			cosphi *= -1.0f;
			q = -q;
		}

		const Scalar DOT_THRESHOLD = 0.9995f;
		if (cosphi > DOT_THRESHOLD) {
			// interpolate linearly
			return (q + (p - q) * t).unit();
		}

		Scalar sinphi = sqrt(1. - cosphi * cosphi);
		Scalar phi = acos(cosphi);

		Quaternion res = q * (sin( phi * (1.-t) ) / sinphi) + p * (sin( phi * t) / sinphi);

		return res;
	}
};

typedef Quaternion<Real> Quat;

}


#endif // __Quaternion_h__


