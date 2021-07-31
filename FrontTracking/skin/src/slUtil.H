// Copyright (c) 2011, Regents of the University of Utah
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


#ifndef slUtil_H
#define slUtil_H

#include <iostream>
#include <fstream>
#include "slVector.H"
#include <vector>

///////////////////////////////////////////////////////////////////

inline std::istream &eatChar(char c,std::istream &buf) {
	char r;
	buf >> r;
	if (r!=c) {
		buf.clear(buf.rdstate() | std::ios::failbit);
	}
	return buf;
}

inline std::istream &eatStr(const char *s,std::istream &buf) {
  while (*s != '\0') {
    eatChar(*s,buf);
    s++;
  }
  return buf;
}

///////////////////////////////////////////////////////////////////

class SlTri {
public:
  int indices[3];
  inline int &operator[](const unsigned int &i) { return indices[i];};
  inline int operator[](const unsigned int &i) const { return indices[i];};
  inline void set(int x, int y, int z) {indices[0] = x; indices[1] = y; indices[2] = z;};
  inline SlTri(int x, int y, int z) {indices[0] = x; indices[1] = y; indices[2] = z;};
  inline SlTri() {};
  inline SlTri &operator=(const SlTri &that);
};
  
/////////////////////////////////////////////////////////////////

class SlInt3 {
public:
	inline SlInt3() { data[0] = data[1] = data[2] = 0; };
	inline SlInt3(int i) { data[0] = data[1] = data[2] = i; };
	inline SlInt3(int i1, int i2, int i3) { data[0] = i1; data[1] = i2; data[2] = i3; };
	inline SlInt3(const SlInt3 &da) { data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; };
	inline SlInt3(const int *da) { data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; };

	inline int &operator[](unsigned int i) { return data[i]; };
	inline int  operator[](unsigned int i) const { return data[i]; };

	inline int &operator()(unsigned int i) { return data[i]; };
	inline int  operator()(unsigned int i) const { return data[i]; };

	inline SlInt3 &set(int i) { data[0] = data[1] = data[2] = i; return *this; };
	inline SlInt3 &set(int i1, int i2, int i3) { data[0] = i1; data[1] = i2; data[2] = i3; return *this; };
	
	inline SlInt3 &set(const SlInt3 &da) { data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; return *this; };
	inline SlInt3 &set(const int *da) { data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; return *this; };
	
	inline SlInt3 &operator=(int i) { data[0] = data[1] = data[2] = i; return *this; };
	inline SlInt3 &operator=(const SlInt3 &da) { data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; return *this; };
	inline SlInt3 &operator=(const int *da) { data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; return *this; };

private:
	int data[3];
};

class SlInt2 {
public:
	inline SlInt2() { data[0] = data[1] = 0; };
	inline SlInt2(int i) { data[0] = data[1] = i; };
	inline SlInt2(int i1, int i2) { data[0] = i1; data[1] = i2; };
	inline SlInt2(const SlInt2 &da) { data[0] = da[0]; data[1] = da[1]; };
	inline SlInt2(const int *da) { data[0] = da[0]; data[1] = da[1]; };

	inline int &operator[](unsigned int i) { return data[i]; };
	inline int  operator[](unsigned int i) const { return data[i]; };

	inline int &operator()(unsigned int i) { return data[i]; };
	inline int  operator()(unsigned int i) const { return data[i]; };

	inline SlInt2 &set(int i) { data[0] = data[1] = i; return *this; };
	inline SlInt2 &set(int i1, int i2) { data[0] = i1; data[1] = i2; return *this; };
	
	inline SlInt2 &set(const SlInt3 &da) { data[0] = da[0]; data[1] = da[1]; return *this; };
	inline SlInt2 &set(const int *da) { data[0] = da[0]; data[1] = da[1]; return *this; };
	
	inline SlInt2 &operator=(int i) { data[0] = data[1] = i; return *this; };
	inline SlInt2 &operator=(const SlInt2 &da) { data[0] = da[0]; data[1] = da[1]; return *this; };
	inline SlInt2 &operator=(const int *da) { data[0] = da[0]; data[1] = da[1]; return *this; };

private:
	int data[2];
};

class SlInt6 {
public:
	inline SlInt6() { data[0] = data[1] = data[2] = data[3] = data[4] = data[5] = 0; };
	inline SlInt6(int i) { data[0] = data[1] = data[2] = data[3] = data[4] = data[5] = i; };
	inline SlInt6(int i0, int i1, int i2, int i3, int i4, int i5) { data[0] = i0; data[1] = i1; data[2] = i2; data[3] = i3; data[4] = i4; data[5] = i5;};
	inline SlInt6(const SlInt6 &da) { data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; data[3] = da[3]; data[4] = da[4]; data[5] = da[5];};
	inline SlInt6(const int *da) { data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; data[3] = da[3]; data[4] = da[4]; data[5] = da[5];};

	inline int &operator[](unsigned int i) { return data[i]; };
	inline int  operator[](unsigned int i) const { return data[i]; };

	inline int &operator()(unsigned int i) { return data[i]; };
	inline int  operator()(unsigned int i) const { return data[i]; };

	inline SlInt6 &set(int i) { data[0] = data[1] = data[2] = data[3] = data[4] = data[5] = i; return *this; };
	inline SlInt6 &set(int i0, int i1, int i2, int i3, int i4, int i5) { data[0] = i0; data[1] = i1; data[2] = i2; data[3] = i3; data[4] = i4; data[5] = i5; return *this; };
	
	inline SlInt6 &set(const SlInt6 &da) { data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; data[3] = da[3]; data[4] = da[4]; data[5] = da[5]; return *this;}; 
	inline SlInt6 &set(const int *da) { data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; data[3] = da[3]; data[4] = da[4]; data[5] = da[5]; return *this;};
	
	inline SlInt6 &operator=(int i) { data[0] = data[1] = data[2] = data[3] = data[4] = data[5] = i; return *this; };
	inline SlInt6 &operator=(const SlInt6 &da) { data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; data[3] = da[3]; data[4] = da[4]; data[5] = da[5]; return *this;}; 
	inline SlInt6 &operator=(const int *da) { data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; data[3] = da[3]; data[4] = da[4]; data[5] = da[5]; return *this;}; 
	void read(std::ifstream &in);
	void write(std::ofstream &out) const;

private:
	int data[12];
};
std::istream &operator>>(std::istream &strm, SlInt6 &v);
std::ostream &operator<<(std::ostream &strm,const SlInt6 &v);
std::istream &operator>>(std::istream &strm, std::vector<SlInt6> &l);
std::ostream &operator<<(std::ostream &strm,const std::vector<SlInt6> &l);

std::istream &operator>>(std::istream &strm, SlTri &t);
std::ostream &operator<<(std::ostream &strm, const SlTri &t);

std::istream &operator>>(std::istream &strm, std::vector<SlTri> &p);
std::ostream &operator<<(std::ostream &strm, const std::vector<SlTri> &p);
  
std::istream &operator>>(std::istream &strm, std::vector<SlVector3> &p);
std::ostream &operator<<(std::ostream &strm, const std::vector<SlVector3> &p);

std::istream &operator>>(std::istream &strm, SlInt3 &p);
std::ostream &operator<<(std::ostream &strm, const SlInt3 &p);

std::istream &operator>>(std::istream &strm, std::vector<SlInt3> &p);
std::ostream &operator<<(std::ostream &strm, const std::vector<SlInt3> &p);
  
std::istream &operator>>(std::istream &strm, std::vector<int> &p);
std::ostream &operator<<(std::ostream &strm, const std::vector<int> &p);

std::istream &operator>>(std::istream &strm, std::vector<double> &p);
std::ostream &operator<<(std::ostream &strm, const std::vector<double> &p);

namespace sl {
	void read(std::ifstream &in, std::vector<int> &l);
	void write(std::ofstream &out, const std::vector<int> &l);
	void read(std::ifstream &in, std::vector<double> &l);
	void write(std::ofstream &out, const std::vector<double> &l);
	void read(std::ifstream &in, SlTri &t);
	void write(std::ofstream &out, const SlTri &t);
	void read(std::ifstream &in, std::vector<SlTri> &l);
	void write(std::ofstream &out, const std::vector<SlTri> &l);
	void read(std::ifstream &in, std::vector<SlVector3> &l);
	void write(std::ofstream &out, const std::vector<SlVector3> &l);
	void read(std::ifstream &in, std::vector<SlInt3> &l);
	void write(std::ofstream &out, const std::vector<SlInt3> &l);
	void read(std::ifstream &in, std::vector<SlInt6> &l);
	void write(std::ofstream &out, const std::vector<SlInt6> &l);
};


void randomize(std::vector<int> &list);
void computeNormals(const std::vector<SlVector3> &meshPts, const std::vector<SlTri> &triangles, 
										std::vector<SlVector3> &vertexNormals, std::vector<SlVector3> &faceNormals);
void computeFaceNormals(const std::vector<SlVector3> &meshPts, const std::vector<SlTri> &triangles, 
												std::vector<SlVector3> &faceNormals);
void bvtRead(char *fname, std::vector<SlVector3> &meshPts, std::vector<SlTri> &triangles);
void bvtDump(char *fname, const std::vector<SlVector3> &meshPts, const std::vector<SlTri> &triangles);
void gtDump(char *fname, const std::vector<SlVector3> &meshPts, const std::vector<SlTri> &triangles);
void gtRead(char *fname, std::vector<SlVector3> &meshPts, std::vector<SlTri> &triangles);
void bgtDump(char *fname, const std::vector<SlVector3> &meshPts, const std::vector<SlTri> &triangles);
bool readObjFile(char *fname, std::vector<SlVector3> &pts, std::vector<SlTri> &triangles);
void objDump(char *fname, const std::vector<SlVector3> &meshPts, const std::vector<SlTri> &triangles);


inline double absmin2 (double x, double y) {
  if (fabs(x) < fabs(y)) return x;
  return y;
}

inline double absmax2 (double x, double y) {
  if (fabs(x) > fabs(y)) return x;
  return y;
}

inline int sign(double x) {
  if (x>0.0) return 1;
	else return -1;
}

inline double max3 (double x, double y, double z) {
  if (x > y && x > z) return x;
  if (y > z) return y;
  return z;
}

inline double min3 (double x, double y, double z) {
  if (x < y && x < z) return x;
  if (y < z) return y;
  return z;
}

inline int max2 (int x, int y) {
  if (x>y) return x;
  else return y;
}

inline int min2 (int x, int y) {
  if (x<y) return x;
  else return y;
}

inline double min2 (double x, double y) {
  if (x<y) return x;
  else return y;
}

inline double max2 (double x, double y) {
  if (x>y) return x;
  else return y;
}

inline double max5(double x1, double x2, double x3, double x4, double x5) {
  double m = x1;
  if (x2 > m) m = x2;
  if (x3 > m) m = x3;
  if (x4 > m) m = x4;
  if (x5 > m) m = x5;
  return m;
}

double pow2(int p);

inline SlTri &SlTri::operator=(const SlTri &that) {
	this->indices[0] = that.indices[0];
	this->indices[1] = that.indices[1];
	this->indices[2] = that.indices[2];
	return (*this);
};

struct eqUnorderedInt3 {
  bool operator()(const SlInt3 &i1, const SlInt3 &i2) const {
    return ((i1[0] == i2[0] && i1[1] == i2[1] && i1[2] == i2[2]) ||
	    (i1[1] == i2[0] && i1[0] == i2[1] && i1[2] == i2[2]) ||
	    (i1[2] == i2[0] && i1[1] == i2[1] && i1[0] == i2[2]) ||
	    (i1[0] == i2[0] && i1[2] == i2[1] && i1[1] == i2[2]) ||
	    (i1[1] == i2[0] && i1[2] == i2[1] && i1[0] == i2[2]) ||
	    (i1[2] == i2[0] && i1[0] == i2[1] && i1[1] == i2[2]));
  }
};

struct eqOrderedInt3 {
	bool operator()(const SlInt3 &i1, const SlInt3 &i2) const {
		return (i1[0] == i2[0] && i1[1] == i2[1] && i1[2] == i2[2]);
	};
};

template <typename T>
class SlArray3D {
	T *data;
	unsigned int _nx, _ny, _nz, nynz;
public:
	inline T &operator()(const unsigned int &i, const unsigned int &j, const unsigned int &k);
	inline T &operator()(const unsigned int &i, const unsigned int &j, const unsigned int &k) const;
	/*inline*/ const T &operator=(const T &d);
	inline const SlArray3D<T> &operator=(const SlArray3D<T> &d);
	inline SlArray3D();
	inline SlArray3D(const unsigned int &nx, const unsigned int &ny, const unsigned int &nz);
	inline ~SlArray3D();
	inline bool allocate(const unsigned int &nx, const unsigned int &ny, const unsigned int &nz);
	unsigned int nx() const {return _nx;};
	unsigned int ny() const {return _ny;};
	unsigned int nz() const {return _nz;};
};

template <typename T> std::istream &operator>>(std::istream &strm, SlArray3D<T> &v);
template <typename T> std::ostream &operator<<(std::ostream &strm, const SlArray3D<T> &v);

template <typename T>
inline bool SlArray3D<T>::allocate(const unsigned int &x, const unsigned int &y, const unsigned int &z) {
	if (data!=NULL && _nx==x && _ny==y && _nz==z) return true;
	_nx = x; 
	_ny = y;
	_nz = z;
	nynz = _ny*_nz;
	if (data) delete []data;
	data = new T[_nx*nynz];
	return true;
}

template<class T>
inline SlArray3D<T> operator*(const SlArray3D<T> &d, T x) {
	SlArray3D<T> y(d.nx(), d.ny(), d.nz());
	for (unsigned int i=0; i<d.nx(); i++)
		for (unsigned int j=0; j<d.ny(); j++)
			for (unsigned int k=0; k<d.nz(); k++)
				y(i,j,k) = d(i,j,k)*x;
	return y;
}

template <class T>
/*inline*/ const T &SlArray3D<T>::operator=(const T &x) {
	T *d = data;
	for (unsigned int i=0; i<_nx; i++) {
		for (unsigned int j=0; j<_ny; j++) {
			for (unsigned int k=0; k<_nz; k++, d++) {
				(*d) = x;
			}
		}
	}
	return x;
}


template <class T>
inline const SlArray3D<T> &SlArray3D<T>::operator=(const SlArray3D<T> &x) {
	T *d = data;
	T *y = x.data;
	for (unsigned int i=0; i<_nx; i++) {
		for (unsigned int j=0; j<_ny; j++) {
			for (unsigned int k=0; k<_nz; k++, d++, y++) {
				(*d) = (*y);
			}
		}
	}
	return x;
}

template <class T>
inline T &SlArray3D<T>::operator()(const unsigned int &i, const unsigned int &j, const unsigned int &k) {
	return data[i*nynz+j*_nz+k];
}

template <class T>
inline T &SlArray3D<T>::operator()(const unsigned int &i, const unsigned int &j, const unsigned int &k) const {
	return data[i*nynz+j*_nz+k];
}

template <class T>
inline SlArray3D<T>::SlArray3D() {
	data = NULL;
}

template <class T>
inline SlArray3D<T>::SlArray3D(const unsigned int &x, const unsigned int &y, const unsigned int &z) {
	data = NULL;
	allocate(x,y,z);
}

template <class T>
inline SlArray3D<T>::~SlArray3D() {
	if (data) delete [] data;
}

template <class T>
std::istream &operator>>(std::istream &strm, SlArray3D<T> &v) {
	unsigned int nx,ny,nz,i,j,k;
	std::ios::fmtflags orgFlags = strm.setf(std::ios::skipws);

	eatStr("[", strm);
	strm >> nx;
	eatStr(",", strm);
	strm >> ny;
	eatStr(",", strm);
	strm >> nz;
	eatStr("]", strm);
	
	eatStr("[", strm);
	v.allocate(nx, ny, nz);
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) {
			for (k=0; k<nz; k++) {
				strm >> v(i,j,k);
			}
		}
	}
	eatStr("]", strm);
	strm.flags(orgFlags);
	return strm;
}

template <typename T>
std::ostream &operator<<(std::ostream &strm, const SlArray3D<T> &v) {
	strm << "[";
	strm << v.nx() << "," << v.ny() << "," <<v.nz() << "]";
	unsigned int i, j, k;
	strm<<"["<<std::endl;
	for (i=0; i<v.nx(); i++) {
		for (j=0; j<v.ny(); j++) {
			for (k=0; k<v.nz(); k++) {
				strm << v(i,j,k) << " ";
			}
		}
	}
	strm<<"]\n";
	return strm;
}
#endif

