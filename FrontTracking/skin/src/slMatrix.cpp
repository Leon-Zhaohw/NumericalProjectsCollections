// Copyright (c) 2011, Regents of the University of Utah
// Copyright (c) 2003-2005, Regents of the University of California.  
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the <organization> nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#include "slMatrix.H"
#include <cmath>
#include "slUtil.H"

//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//  Implementation for  SlMatrix3x3

std::istream &operator>>(std::istream &strm, SlMatrix3x3 &m) {
	std::ios::fmtflags orgFlags = strm.setf(std::ios::skipws);
	eatChar('[',strm);

	int r,c;

	for(r=0;r<3;r++) {
		eatChar('[',strm);
		for(c=0;c<3;c++) {
			strm >> m(r,c);
			if(c!=2) eatChar(',',strm);
		}
		eatChar(']',strm);
	}

	eatChar(']',strm);
	strm.flags(orgFlags);
	return strm;
}

std::ostream &operator<<(std::ostream &strm,const SlMatrix3x3 &m) {
	strm << "[";
	int r,c;

	for(r=0;r<3;r++) {
		strm << "[ ";
		for(c=0;c<3;c++) {
			strm << m(r,c);
			if(c!=2) strm << ",";
		}
		strm << " ]";
		if(r!=2) strm << "\n ";
	}

	strm << " ]\n";

	return strm;
}


//-------------------------------------------------------------------

SlMatrix3x3 &SlMatrix3x3::inplaceMultPre (const SlMatrix3x3 &that) {
	SlVector3 tmp;

	int i;

	for(i=0;i<3;i++) {

		tmp[0] = ( that.data[0][0] * data[0][i] +
				  that.data[0][1] * data[1][i] +
				  that.data[0][2] * data[2][i] );

		tmp[1] = ( that.data[1][0] * data[0][i] +
				  that.data[1][1] * data[1][i] +
				  that.data[1][2] * data[2][i] );

		tmp[2] = ( that.data[2][0] * data[0][i] +
				  that.data[2][1] * data[1][i] +
				  that.data[2][2] * data[2][i] );

		data[0][i] = tmp[0];
		data[1][i] = tmp[1];
		data[2][i] = tmp[2];
	}

	return (*this);
}


//-------------------------------------------------------------------

SlMatrix3x3 &SlMatrix3x3::inplaceMultPost(const SlMatrix3x3 &that) {
	SlVector3 tmp;

	int i;

	for(i=0;i<3;i++) {

		tmp[0] = ( data[i][0] * that.data[0][0] +
				  data[i][1] * that.data[1][0] +
				  data[i][2] * that.data[2][0] );

		tmp[1] = ( data[i][0] * that.data[0][1] +
				  data[i][1] * that.data[1][1] +
				  data[i][2] * that.data[2][1] );

		tmp[2] = ( data[i][0] * that.data[0][2] +
				  data[i][1] * that.data[1][2] +
				  data[i][2] * that.data[2][2] );

		data[i][0] = tmp[0];
		data[i][1] = tmp[1];
		data[i][2] = tmp[2];
	}

	return (*this);
}

//-------------------------------------------------------------------

#define SWAP(a,b) { tmp = (a) ; (a) = (b) ; (b) = tmp; }

static inline double smFabs(double x) {
	return (x>0)?(x):(-x);
}

static void rowswap(SlMatrix3x3 &m, int i, int j) {
	double tmp;
	SWAP(m(i,0),m(j,0));
	SWAP(m(i,1),m(j,1));
	SWAP(m(i,2),m(j,2));
}


SlMatrix3x3 inverse(const SlMatrix3x3 &a) {
	SlMatrix3x3 _a(a);
	SlMatrix3x3 _b   ; _b.setIdentity();
	int i,j,i1;
	for(j=0;j<3;j++) {
		i1 = j;
		for(i=j+1;i<3;i++)
			if (smFabs(_a(i,j)) > smFabs(_a(i1,j)))
				i1 = i;
		rowswap(_a,i1,j);
		rowswap(_b,i1,j);
		if (((_a(j,j) >= 0.0) && (_a(j,j) <  1e-6)) ||
			((_a(j,j) <= 0.0) && (_a(j,j) > -1e-6)) ){
			std::cerr << "Inverse of singular matrix\n" << std::flush;
			abort();
		}
		double div = 1.0/_a(j,j);
		_b(j,0) *= div;
		_b(j,1) *= div;
		_b(j,2) *= div;
		_a(j,0) *= div;
		_a(j,1) *= div;
		_a(j,2) *= div;
		for(i=0;i<3;i++) {
			if (i != j) {
				double tmp = _a(i,j);
				_b(i,0) -= tmp*_b(j,0);
				_b(i,1) -= tmp*_b(j,1);
				_b(i,2) -= tmp*_b(j,2);
				_a(i,0) -= tmp*_a(j,0);
				_a(i,1) -= tmp*_a(j,1);
				_a(i,2) -= tmp*_a(j,2);
			}
		}
	}
#ifdef DEBUG
	SlMatrix3x3 test = a*_b;
	for(i=0;i<3;i++) {
		for(j=0;j<3;j++) {
			if (smFabs(test(i,j)-((i==j)?(1):(0))) > 1e-6) {
				std::cerr << "There is a bug in the invese routine for SlMatrix3x3.\n";
				abort();
			}
		}
	}
#endif
	return _b;
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------

// Lots of support code for SlSymetricEigenDecomp (SlMatrix3x3)
// This code is also used by the SVD routines.

static const double EigenEps = 1e-10;

#define EIGEN_ROTATE(a,i,j,k,l) {		\
g=a(i,j);					\
h=a(k,l);					\
a(i,j)=g-s*(h+g*tau); 			\
a(k,l)=h+s*(g-h*tau);				\
}

static inline double eigenFabs(double x) {
	return (x>0)?(x):(-x);
}

static inline double eigenMax(double a, double b) {
	return (a>b)?(a):(b);
}

static inline double eigenEpsEq(double a, double b) {
	return (eigenFabs((a)-(b)) < EigenEps);
}

static inline double eigenEpsEqZero(double a) {
	return (eigenFabs(a) < EigenEps);
}

static inline double eigenSqr(double x) {
	return x * x;
}

static inline double eigenPythag(double a, double b) {
	double absa,absb;
	absa=eigenFabs(a);
	absb=eigenFabs(b);
	if (absa > absb) return absa*sqrt(1.0+eigenSqr(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+eigenSqr(absa/absb)));
}

static inline double eigenSign(double a, double b) {
	return ((b) >= 0.0 ? eigenFabs(a) : -eigenFabs(a));
}


static void eigenTred2( SlMatrix3x3 &a, SlVector3 &d, SlVector3 &e) {
	// Adapted from the Numerical Recipies in C routine "tred2"

	double scale,hh,h,g,f;

	scale = eigenFabs(a(2,0)) + eigenFabs(a(2,1));

	if (eigenEpsEqZero(scale)) {
		e[2]=a(2,1);
		h = 0.0;
	} else {
		a(2,0) /= scale;
		a(2,1) /= scale;
		h = a(2,0)*a(2,0) + a(2,1)*a(2,1);
		f = a(2,1);
		g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
		e[2] = scale*g;
		h -= f*g;
		a(2,1) = f-g;

		a(0,2) = a(2,0)/h;
		g = a(0,0)*a(2,0) + a(1,0)*a(2,1);
		e[0]=g/h;
		f = e[0]*a(2,0);

		a(1,2) = a(2,1)/h;
		g = a(1,0)*a(2,0) + a(1,1)*a(2,1);
		e[1] = g/h;
		f += e[1]*a(2,1);

		hh=f/(h+h);

		f = a(2,0);
		e[0]=g=e[0]-hh*f;
		a(0,0) -= (f*e[0]+g*a(2,0));

		f = a(2,1);
		e[1] = g = e[1]-hh*f;
		a(1,0) -= (f*e[0]+g*a(2,0));
		a(1,1) -= (f*e[1]+g*a(2,1));
	}
	d[2]=h;

	e[1]=a(1,0);
	e[0]=0.0;

	d[0]=a(0,0);
	a(0,0)=1.0;

	d[1]=a(1,1);
	a(1,1)=1.0;
	a(0,1)=a(1,0)=0.0;

	if (d[2]) {
		g = a(2,0)*a(0,0);
		a(0,0) -= g*a(0,2);
		a(1,0) -= g*a(1,2);

		g = a(2,1);
		a(0,1) -= g*a(0,2);
		a(1,1) -= g*a(1,2);
	}
	d[2]=a(2,2);
	a(2,2)=1.0;
	a(0,2)=a(2,0)=0.0;
	a(1,2)=a(2,1)=0.0;
}



static void eigenTqli(SlVector3 &d, SlVector3 &e, SlMatrix3x3 &z) {
	// Adapted from the Numerical Recipies in C routine "tqli"

	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;

	for (i=1;i<3;i++) e[i-1]=e[i];
	e[2]=0.0;

	for (l=0;l<3;l++) {

		iter=0;
		do {
			for (m=l;m<2;m++) {
				dd= eigenFabs(d[m])+eigenFabs(d[m+1]);
				if (eigenEpsEq((eigenFabs(e[m])+dd),dd)) break;
			}
			if (m != l) {
				if (iter++ >= 30) {
					std::cerr << "SlSymetricEigenDecomp:: Too many iterations in eigenTqli.\n";
					return;
				}
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=eigenPythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+eigenSign(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=eigenPythag(f,g));
					if (eigenEpsEqZero(r)) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					for (k=0;k<3;k++) {
						f=z(k,i+1);
						z(k,i+1)=s*z(k,i)+c*f;
						z(k,i)=c*z(k,i)-s*f;
					}
				}
				if (eigenEpsEqZero(r) && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}

void SlSymetricEigenDecomp(const SlMatrix3x3 &a,
							  /* */ SlVector3   &vals,
							  /* */ SlMatrix3x3 &vecs) {
	SlVector3 e;
	SlMatrix3x3 q(a);

#ifdef DEBUG
	if (!(eigenEpsEq( q(0,1),q(1,0) ) &&
		  eigenEpsEq( q(2,1),q(1,2) ) &&
		  eigenEpsEq( q(2,0),q(0,2) ) )) {
		std::cerr << "SlSymetricEigenDecomp called with non-symetric matrix.\n";
		abort();
	}
#endif

	vals = 0.0;
	eigenTred2(q,vals,e);
	vecs.setIdentity();
	eigenTqli(vals,e,vecs);
	vecs.inplaceMultPre(transpose(q));
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------

static
void SlSVDecompInternal(SlMatrix3x3 &a, SlVector3 &w, SlMatrix3x3 &v) {

	// Adapted from the Numerical Recipies in C routine "dsvdcmp"

	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z;

	SlVector3 rv1;
	g=scale=anorm=0.0;
	for (i=0;i<3;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < 3) {
			for (k=i;k<3;k++) scale += eigenFabs(a(k,i));
			if (scale) {
				for (k=i;k<3;k++) {
					a(k,i) /= scale;
					s += a(k,i)*a(k,i);
				}
				f=a(i,i);
				g = -eigenSign(sqrt(s),f);
				h=f*g-s;
				a(i,i)=f-g;
				for (j=l;j<3;j++) {
					for (s=0.0,k=i;k<3;k++) s += a(k,i)*a(k,j);
					f=s/h;
					for (k=i;k<3;k++) a(k,j) += f*a(k,i);
				}
				for (k=i;k<3;k++) a(k,i) *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i < 3 && i != 2 ) {
			for (k=l;k<3;k++) scale += eigenFabs(a(i,k));
			if (scale) {
				for (k=l;k<3;k++) {
					a(i,k) /= scale;
					s += a(i,k)*a(i,k);
				}
				f=a(i,l);
				g = -eigenSign(sqrt(s),f);
				h=f*g-s;
				a(i,l)=f-g;
				for (k=l;k<3;k++) rv1[k]=a(i,k)/h;
				for (j=l;j<3;j++) {
					for (s=0.0,k=l;k<3;k++) s += a(j,k)*a(i,k);
					for (k=l;k<3;k++) a(j,k) += s*rv1[k];
				}
				for (k=l;k<3;k++) a(i,k) *= scale;
			}
		}
		anorm=eigenMax(anorm,(eigenFabs(w[i])+eigenFabs(rv1[i])));
	}


	for (i=2;i>=0;i--) {
		if (i < 2) {
			if (g) {
				for (j=l;j<3;j++) v(j,i)=(a(i,j)/a(i,l))/g;
				for (j=l;j<3;j++) {
					for (s=0.0,k=l;k<3;k++) s += a(i,k)*v(k,j);
					for (k=l;k<3;k++) v(k,j) += s*v(k,i);
				}
			}
			for (j=l;j<3;j++) v(i,j)=v(j,i)=0.0;
		}
		v(i,i)=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=2;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<3;j++) a(i,j)=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<3;j++) {
				for (s=0.0,k=l;k<3;k++) s += a(k,i)*a(k,j);
				f=(s/a(i,i))*g;
				for (k=i;k<3;k++) a(k,j) += f*a(k,i);
			}
			for (j=i;j<3;j++) a(j,i) *= g;
		} else for (j=i;j<3;j++) a(j,i)=0.0;
		++a(i,i);
	}

	for (k=2;k>=0;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if ((double)(eigenFabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if(nm < 0) break; //avoid a memory bug when there are NaNs
				if ((double)(eigenFabs(w[nm])+anorm) == anorm) break;
			}

			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(eigenFabs(f)+anorm) == anorm) break;
					g=w[i];
					h=eigenPythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<3;j++) {
						y=a(j,nm);
						z=a(j,i);
						a(j,nm)=y*c+z*s;
						a(j,i)=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<3;j++) v(j,k) = -v(j,k);
				}
				break;
			}
			if (its == 30) {
				std::cerr << "SlSVDecomp:: Too many iterations in SlSVDecompInternal.\n";
				return;
			}
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=eigenPythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+eigenSign(g,f)))-h))/x;
			c=s=1.0;

			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=eigenPythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<3;jj++) {
					x=v(jj,j);
					z=v(jj,i);
					v(jj,j)=x*c+z*s;
					v(jj,i)=z*c-x*s;
				}
				z=eigenPythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<3;jj++) {
					y=a(jj,j);
					z=a(jj,i);
					a(jj,j)=y*c+z*s;
					a(jj,i)=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
}




void SlSVDecomp(const SlMatrix3x3 &a,
				   /* */ SlMatrix3x3 &u,
				   /* */ SlVector3   &w,
				   /* */ SlMatrix3x3 &v) {
	u = a;
	SlSVDecompInternal(u,w,v);
}


#define SMEIG_EPS  1.e-16
#define SMEIG_EPSQ 1.e-8   /* sqrt(EPS) */

#define SMEIG_DET3(m) ( m(0,0)*m(1,1)*m(2,2)-m(0,0)*m(2,1)*m(1,2)-m(0,1)*m(1,0)*m(2,2) \
+m(0,1)*m(2,0)*m(1,2)+m(0,2)*m(1,0)*m(2,1)-m(0,2)*m(2,0)*m(1,1) )

#define SMEIG_PI 3.14159265358979323846

#define SMEIG_MIN(a,b) (((a)>(b)) ? (b) : (a))

#define SMEIG_SWAP(x,y) (th=(x),(x)=(y),(y)=th)

#define SMEIG_CSWAP(i,j) (SMEIG_SWAP(vecs(0,i),vecs(0,j)),SMEIG_SWAP(vecs(1,i),vecs(1,j)),SMEIG_SWAP(vecs(2,i),vecs(2,j)))

void symeig_3(const SlMatrix3x3 &a, SlVector3 &vals, SlMatrix3x3 &vecs)
//, double *e , int dovec )
{
	double aa,bb,cc,dd,ee,ff ;
	double a1,a2,a3 , qq,rr, qs,th , lam1,lam2,lam3 ;
	double aba,abb,abc,abd,abe,abf , ann ;
	double d12,d13,d23 ;
	double u1,u2,u3 , v1,v2,v3 , w1,w2,w3 , t1,t2,t3 , tn ;

	aa = a(0,0) ; bb = a(0,1) ; cc = a(0,2) ;  // matrix is [ aa bb cc ]
	dd = a(1,1) ; ee = a(1,2) ; ff = a(2,2) ;  //           [ bb dd ee ]
	//           [ cc ee ff ]

	aba = fabs(aa) ; abb = fabs(bb) ; abc = fabs(cc) ;
	abd = fabs(dd) ; abe = fabs(ee) ; abf = fabs(ff) ;
	ann = aba+abb+abc+abd+abe+abf   ;                 /* matrix 'norm' */

	if( ann == 0.0 ){             /* matrix is all zero! */
		vals[0] = vals[1] = vals[2] = 0.0 ;
		vecs(0,0) = vecs(1,1) = vecs(2,2) = 1.0 ;
		vecs(0,1) = vecs(0,2) = vecs(1,0) = vecs(1,2) = vecs(2,0) = vecs(2,1) = 0.0;

		return ;
	}

	/*----- check for matrix that is essentially diagonal -----*/

	if( SMEIG_EPS*aba > (abb+abc) && SMEIG_EPS*abd > (abb+abe) &&
       SMEIG_EPS*abf > (abc+abe) ){

		lam1 = aa ; lam2 = dd ; lam3 = ff ;

		vecs(0,0) = vecs(1,1) = vecs(2,2) = 1.0 ;
		vecs(0,1) = vecs(0,2) = vecs(1,0) = vecs(1,2) = vecs(2,0) = vecs(2,1) = 0.0 ;

		if( lam1 > lam2 ){ SMEIG_SWAP(lam1,lam2) ; SMEIG_CSWAP(0,1) ; }
		if( lam1 > lam3 ){ SMEIG_SWAP(lam1,lam3) ; SMEIG_CSWAP(0,2) ; }
		if( lam2 > lam3 ){ SMEIG_SWAP(lam2,lam3) ; SMEIG_CSWAP(1,2) ; }
		if( SMEIG_DET3(a) < 0.0 ){
			vecs(0,2) = -vecs(0,2); vecs(1,2) = -vecs(1,2); vecs(2,2) = -vecs(2,2);
		}
		vals[0] = lam1 ; vals[1] = lam2 ; vals[2] = lam3 ;
		return ;
	}

	/*----- not diagonal ==> must solve cubic polynomial for eigenvalues -----*/
	/*      the cubic polynomial is x**3 + a1*x**2 + a2*x + a3 = 0            */

	a1 = -(aa+dd+ff) ;
	a2 =  (aa*ff+aa*dd+dd*ff - bb*bb-cc*cc-ee*ee) ;
	a3 =  ( aa*(ee*ee-dd*ff) + bb*(bb*ff-cc*ee) + cc*(cc*dd-bb*ee) ) ;

	qq = (a1*a1 - 3.0*a2) / 9.0 ;
	rr = (2.0*a1*a1*a1 - 9.0*a1*a2 + 27.0*a3) / 54.0 ;

	qs = sqrt(qq) ; rr = rr / (qs*qq) ;
	if( rr < -1.0 ) rr = -1.0 ; else if( rr > 1.0 ) rr = 1.0 ;
	th = acos(rr) ;

	lam1 = -2.0 * qs * cos(  th        /3.0 ) - a1 / 3.0 ;
	lam2 = -2.0 * qs * cos( (th+2.0*SMEIG_PI)/3.0 ) - a1 / 3.0 ;
	lam3 = -2.0 * qs * cos( (th+4.0*SMEIG_PI)/3.0 ) - a1 / 3.0 ;

	/*-- are doing eigenvectors; must do double root as a special case --*/

#undef  CROSS
#define CROSS(x1,x2,x3,y1,y2,y3,z1,z2,z3) \
( (z1)=(x2)*(y3)-(x3)*(y2), (z2)=(x3)*(y1)-(x1)*(y3), (z3)=(x1)*(y2)-(x2)*(y1) )

	d12 = fabs(lam1-lam2) ; d13 = fabs(lam1-lam3) ; d23 = fabs(lam2-lam3) ;
	rr  = SMEIG_MIN(d12,d13)    ; rr  = SMEIG_MIN(rr,d23)     ;

	if( rr > SMEIG_EPS*ann ){  /*---- not a double root ----*/

		if( lam1 > lam2 )  SMEIG_SWAP(lam1,lam2) ;  /* start by sorting eigenvalues */
		if( lam1 > lam3 )  SMEIG_SWAP(lam1,lam3) ;
		if( lam2 > lam3 )  SMEIG_SWAP(lam2,lam3) ;
		vals[0] = lam1 ; vals[1] = lam2 ; vals[2] = lam3 ;

		/* find eigenvector for lam1 by computing Ay-lam1*y for
		 vectors y1=[1,0,0], y2=[0,1,0], and y3=[0,0,1]; the eigenvector
		 is orthogonal to all of these, so use the cross product to get it */

		u1 = aa-lam1 ; u2 = bb      ; u3 = cc ;   /* A*y1 - lam1*y1 */
		v1 = bb      ; v2 = dd-lam1 ; v3 = ee ;   /* A*y2 - lam1*y2 */
		CROSS(u1,u2,u3 , v1,v2,v3 , t1,t2,t3 ) ;     tn = sqrt(t1*t1+t2*t2+t3*t3) ;
		if( tn < SMEIG_EPSQ*ann ){                      /* u and v were parallel? */
			w1 = cc ; w2 = ee ; w3 = ff-lam1 ;      /* A*y3 - lam1*y3 */
			CROSS(u1,u2,u3 , w1,w2,w3 , t1,t2,t3 ) ;   tn = sqrt(t1*t1+t2*t2+t3*t3) ;
			if( tn < SMEIG_EPSQ*ann ){                    /* u and w were parallel? */
				CROSS(v1,v2,v3 , w1,w2,w3 , t1,t2,t3 ) ; tn = sqrt(t1*t1+t2*t2+t3*t3) ;
			}
		}
		vecs(0,0) = t1/tn ; vecs(1,0) = t2/tn ; vecs(2,0) = t3/tn ;  /* normalize */

		/* do same for lam2 */

		u1 = aa-lam2 ; u2 = bb      ; u3 = cc ;
		v1 = bb      ; v2 = dd-lam2 ; v3 = ee ;
		CROSS(u1,u2,u3 , v1,v2,v3 , t1,t2,t3 ) ;     tn = sqrt(t1*t1+t2*t2+t3*t3) ;
		if( tn < SMEIG_EPSQ*ann ){
			w1 = cc ; w2 = ee ; w3 = ff-lam2 ;
			CROSS(u1,u2,u3 , w1,w2,w3 , t1,t2,t3 ) ;   tn = sqrt(t1*t1+t2*t2+t3*t3) ;
			if( tn < SMEIG_EPSQ*ann ){
				CROSS(v1,v2,v3 , w1,w2,w3 , t1,t2,t3 ) ; tn = sqrt(t1*t1+t2*t2+t3*t3) ;
			}
		}
		vecs(0,1) = t1/tn ; vecs(1,1) = t2/tn ; vecs(2,1) = t3/tn ;

		/* orthgonality of eigenvectors ==> can get last one by cross product */

#if 1
		CROSS( vecs(0,1),vecs(1,0),vecs(2,0) , vecs(0,1),vecs(1,1),vecs(2,1) ,
			  vecs(0,2),vecs(1,2),vecs(2,2) ) ;
#else
		u1 = aa-lam3 ; u2 = bb      ; u3 = cc ;
		v1 = bb      ; v2 = dd-lam3 ; v3 = ee ;
		CROSS(u1,u2,u3 , v1,v2,v3 , t1,t2,t3 ) ;     tn = sqrt(t1*t1+t2*t2+t3*t3) ;
		if( tn < SMEIG_EPSQ*ann ){
			w1 = cc ; w2 = ee ; w3 = ff-lam3 ;
			CROSS(u1,u2,u3 , w1,w2,w3 , t1,t2,t3 ) ;   tn = sqrt(t1*t1+t2*t2+t3*t3) ;
			if( tn < SMEIG_EPSQ*ann ){
				CROSS(v1,v2,v3 , w1,w2,w3 , t1,t2,t3 ) ; tn = sqrt(t1*t1+t2*t2+t3*t3) ;
			}
		}
		vecs(0,2) = t1/tn ; vecs(1,2) = t2/tn ; vecs(2,2) = t3/tn ;
#endif

		return ;

	} else { /*---- if here, we have a double root ----*/

		/* make sure that we have lam1=lam2 and lam3 is the outlier */

		if( d13 < d12 && d13 < d23 ) SMEIG_SWAP(lam2,lam3) ;
		else if( d23 < d12 && d23 < d13 ) SMEIG_SWAP(lam1,lam3) ;
		lam1 = lam2 = 0.5*(lam1+lam2) ;

		/* compute eigenvector for lam3 using method as above */

		u1 = aa-lam3 ; u2 = bb      ; u3 = cc ;
		v1 = bb      ; v2 = dd-lam3 ; v3 = ee ;
		CROSS(u1,u2,u3 , v1,v2,v3 , t1,t2,t3 ) ;     tn = sqrt(t1*t1+t2*t2+t3*t3) ;
		if( tn < SMEIG_EPSQ*ann ){
			w1 = cc ; w2 = ee ; w3 = ff-lam3 ;
			CROSS(u1,u2,u3 , w1,w2,w3 , t1,t2,t3 ) ;   tn = sqrt(t1*t1+t2*t2+t3*t3) ;
			if( tn < SMEIG_EPSQ*ann ){
				CROSS(v1,v2,v3 , w1,w2,w3 , t1,t2,t3 ) ; tn = sqrt(t1*t1+t2*t2+t3*t3) ;
			}
		}
		w1 = vecs(0,2) = t1/tn ; w2 = vecs(1,2) = t2/tn ; w3 = vecs(2,2) = t3/tn ;

		/* find a vector orthogonal to it */

		CROSS(w1,w2,w3 , 1,0,0 , t1,t2,t3 ) ; tn = sqrt(t1*t1+t2*t2+t3*t3) ;
		if( tn < SMEIG_EPSQ ){
			CROSS(w1,w2,w3 , 0,1,0 , t1,t2,t3 ) ; tn = sqrt(t1*t1+t2*t2+t3*t3) ;
			if( tn < SMEIG_EPSQ ){
				CROSS(w1,w2,w3 , 0,0,1 , t1,t2,t3 ) ; tn = sqrt(t1*t1+t2*t2+t3*t3) ;
			}
		}
		vecs(0,0) = t1/tn ; vecs(1,0) = t2/tn ; vecs(2,0) = t3/tn ;

		/* and the final vector is the cross product of these two */

		CROSS( w1,w2,w3 , vecs(0,0),vecs(0,1),vecs(0,2) ,
			  vecs(0,1),vecs(1,1),vecs(2,1) ) ;

		/* sort results (we know lam1==lam2) */

		if( lam1 > lam3 ){
			SMEIG_SWAP(lam1,lam3) ; SMEIG_CSWAP(0,2) ;
			if( SMEIG_DET3(a) < 0.0 ){
				vecs(0,2) = -vecs(0,2); vecs(1,2) = -vecs(1,2); vecs(2,2) = -vecs(2,2);
			}
		}

		vals[0] = lam1 ; vals[1] = lam2 ; vals[2] = lam3 ;
		return ;
	}
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//  Implementation for  SlMatrix2x2

std::istream &operator>>(std::istream &strm, SlMatrix2x2 &m) {
	std::ios::fmtflags orgFlags = strm.setf(std::ios::skipws);
	eatChar('[',strm);

	int r,c;

	for(r=0;r<2;r++) {
		eatChar('[',strm);
		for(c=0;c<2;c++) {
			strm >> m(r,c);
			if(c!=1) eatChar(',',strm);
		}
		eatChar(']',strm);
	}

	eatChar(']',strm);
	strm.flags(orgFlags);
	return strm;
}

std::ostream &operator<<(std::ostream &strm,const SlMatrix2x2 &m) {
	strm << "[";
	int r,c;

	for(r=0;r<2;r++) {
		strm << "[ ";
		for(c=0;c<2;c++) {
			strm << m(r,c);
			if(c!=1) strm << ",";
		}
		strm << " ]";
		if(r!=1) strm << "\n ";
	}
	strm << " ]\n";
	return strm;
}

//-------------------------------------------------------------------

SlMatrix2x2 &SlMatrix2x2::inplaceMultPre (const SlMatrix2x2 &that) {
	SlVector2 tmp;
	int i;

	for(i=0;i<2;i++) {
		tmp[0] = ( that.data[0][0] * data[0][i] +
				  that.data[0][1] * data[1][i] );
		tmp[1] = ( that.data[1][0] * data[0][i] +
				  that.data[1][1] * data[1][i] );
		data[0][i] = tmp[0];
		data[1][i] = tmp[1];
	}
	return (*this);
}


//-------------------------------------------------------------------

SlMatrix2x2 &SlMatrix2x2::inplaceMultPost(const SlMatrix2x2 &that) {
	SlVector2 tmp;
	int i;
	for(i=0;i<2;i++) {
		tmp[0] = ( data[i][0] * that.data[0][0] +
				  data[i][1] * that.data[1][0] );
		tmp[1] = ( data[i][0] * that.data[0][1] +
				  data[i][1] * that.data[1][1] );
		data[i][0] = tmp[0];
		data[i][1] = tmp[1];
	}
	return (*this);
}


//-------------------------------------------------------------------

SlMatrix2x2 inverse(const SlMatrix2x2 &a) {
	SlMatrix2x2 result;

	double invdet = determinant(a);
	if(invdet == 0.0) {
		std::cerr << "Inverse of singular matrix\n" << std::flush;
		abort();
	}
	invdet = 1.0 / invdet;

	result(1, 1) = invdet * a(0, 0);
	result(0, 1) = -invdet * a(0, 1);
	result(1, 0) = -invdet * a(1, 0);
	result(0, 0) = invdet * a(1, 1);

#ifdef DEBUG
	SlMatrix2x2 test = a * result;
	for(int i=0;i<2;i++) {
		for(int j=0;j<2;j++) {
			if (smFabs(test(i,j)-((i==j)?(1):(0))) > 1e-10) {
				std::cerr << "There is a bug in the inverse routine for SlMatrix2x2.\n";
				abort();
			}
		}
	}
#endif
	return result;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------

void SlSymetricEigenDecomp(const SlMatrix2x2 &a,
							  /* */ SlVector2   &vals,
							  /* */ SlMatrix2x2 &vecs) {
#ifdef DEBUG
    if (a(1,0) != a(0,1)) {
			std::cerr << "SlSymetricEigenDecomp :: Routine called with unsymetric matrix.\n" << flush;
			abort();
    }
    if (isnan(a(1,0))) {
			std::cerr << "SlSymetricEigenDecomp :: Routine called with NaN matrix.\n" << flush;
			abort();
    }
#endif

	double ct = sqrt(4*eigenSqr(a(0,1)) + eigenSqr(a(0,0) - a(1,1)));

	double val0 = 0.5*(a(0,0) + a(1,1) - ct);
	double val1 = 0.5*(a(0,0) + a(1,1) + ct);

	double d = (a(0,0) >= a(1,1)) ?
    sqrt(eigenSqr(a(0,1)) + 0.25*eigenSqr( a(0,0) - a(1,1) + ct )) :
    sqrt(eigenSqr(a(0,1)) + 0.25*eigenSqr(-a(0,0) + a(1,1) + ct )) ;

	double p;
	double q;

	if ((d+1.0) == 1.0) {
		p = 0.0;
		q = 0.0;
	}else{
		d = 1.0 / d;
		if (a(0,0) >= a(1,1)) {
			p = d * a(0,1);
			q = d * 0.5 * (-a(0,0) + a(1,1) - ct );
		}else{
			q = d * a(0,1);
			p = d * 0.5 * ( a(0,0) - a(1,1) - ct );
		}
	}

	vals[0] = val0;
	vals[1] = val1;

	vecs(0,0) = p;
	vecs(0,1) = q;

	vecs(1,0) = -q;
	vecs(1,1) =  p;
}





//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------

