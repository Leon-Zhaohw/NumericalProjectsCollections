
#ifndef ARRAY3D_H
#define ARRAY3D_H
#include <iostream>

template <typename T>
class Array3D {
	T *data;
	unsigned int _nx, _ny, _nz, nynz;
public:
	inline T &operator()(const unsigned int &i, const unsigned int &j, const unsigned int &k);
	inline T &operator()(const unsigned int &i, const unsigned int &j, const unsigned int &k) const;
	const T &operator=(const T &d);
	inline const Array3D<T> &operator=(const Array3D<T> &d);
	inline Array3D();
	inline Array3D(const unsigned int &nx, const unsigned int &ny, const unsigned int &nz);
	inline ~Array3D();
	inline bool allocate(const unsigned int &nx, const unsigned int &ny, const unsigned int &nz);
	unsigned int nx() const {return _nx;};
	unsigned int ny() const {return _ny;};
	unsigned int nz() const {return _nz;};
};

template <typename T> std::istream &operator>>(std::istream &strm, Array3D<T> &v);
template <typename T> std::ostream &operator<<(std::ostream &strm, const Array3D<T> &v);

template <typename T>
inline bool Array3D<T>::allocate(const unsigned int &x, const unsigned int &y, const unsigned int &z) {
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
inline Array3D<T> operator*(const Array3D<T> &d, T x) {
	Array3D<T> y(d.nx(), d.ny(), d.nz());
	for (unsigned int i=0; i<d.nx(); i++)
		for (unsigned int j=0; j<d.ny(); j++)
			for (unsigned int k=0; k<d.nz(); k++)
				y(i,j,k) = d(i,j,k)*x;
	return y;
}

template<class T>
inline double dot(const Array3D<T> &x, const Array3D<T> &y) {
  double d=0.0;
  for (unsigned int i=0; i<x.nx(); i++)
	for (unsigned int j=0; j<x.ny(); j++)
	  for (unsigned int k=0; k<x.nz(); k++)
		d += x(i,j,k)*y(i,j,k);
  return d;
}

template<class T>
inline double nrm2(const Array3D<T> &x) {
  double d=0.0;
  for (unsigned int i=0; i<x.nx(); i++)
	for (unsigned int j=0; j<x.ny(); j++)
	  for (unsigned int k=0; k<x.nz(); k++)
		d += x(i,j,k)*x(i,j,k);
  return sqrt(d);
}

template<class T>
inline void daxpy(double alpha, const Array3D<T> &x, Array3D<T> &y) {
  for (unsigned int i=0; i<x.nx(); i++)
	for (unsigned int j=0; j<x.ny(); j++)
	  for (unsigned int k=0; k<x.nz(); k++)
		y(i,j,k) += alpha*x(i,j,k);
}

template <class T>
/*inline*/ const T &Array3D<T>::operator=(const T &x) {
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
inline const Array3D<T> &Array3D<T>::operator=(const Array3D<T> &x) {
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
inline T &Array3D<T>::operator()(const unsigned int &i, const unsigned int &j, const unsigned int &k) {
	return data[i*nynz+j*_nz+k];
}

template <class T>
inline T &Array3D<T>::operator()(const unsigned int &i, const unsigned int &j, const unsigned int &k) const {
	return data[i*nynz+j*_nz+k];
}

template <class T>
inline Array3D<T>::Array3D() {
	data = NULL;
}

template <class T>
inline Array3D<T>::Array3D(const unsigned int &x, const unsigned int &y, const unsigned int &z) {
	data = NULL;
	allocate(x,y,z);
}

template <class T>
inline Array3D<T>::~Array3D() {
	if (data) delete [] data;
}

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

template <class T>
std::istream &operator>>(std::istream &strm, Array3D<T> &v) {
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
std::ostream &operator<<(std::ostream &strm, const Array3D<T> &v) {
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

