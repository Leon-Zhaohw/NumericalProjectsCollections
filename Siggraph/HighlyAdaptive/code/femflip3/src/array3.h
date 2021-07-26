/*
 *	array3.h
 *	
 *	Created by Ryoichi Ando on 11/20/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <stdio.h>
#include <vector>
#include "macros.h"

#ifndef _ARRAY3_H
#define _ARRAY3_H

struct size3 {
	size3() : w(0), h(0), d(0) {}
	size3( uint w, uint h, uint d ) {
		this->w = w;
		this->h = h;
		this->d = d;
	}
	uint w;
	uint h;
	uint d;
	uint& operator [](uint idx) {
		if( idx == 0 ) return w;
		if( idx == 1 ) return h;
		if( idx == 2 ) return d;
	}
	bool operator==(const size3 &size ) {
		return w == size.w && h == size.h && d == size.d;
	}
	bool operator!=(const size3 &size ) {
		return 1 - ((*this) == size);
	}
	bool validIndex( int i, int j, int k ) {
		return i>=0 && i<w && j>=0 && j<h && k>=0 && k<d;
	}
};

template <class T> class array3 {
public:
	array3() : ptr(NULL) {}
	array3( const array3<T> &array ) : ptr(NULL) {
		copy(array);
	}
	array3( int w, int h, int d ) : ptr(NULL) {
		resize(w,h,d);
	}
	array3( size3 theSize ) : ptr(NULL) {
		resize(theSize);
	}
	~array3() {
		free3D(ptr);
	}
	bool resize( size3 theSize ) {
		if( arraySize != theSize ) {
			if( ptr ) free3D(ptr);
			arraySize = theSize;
			ptr = alloc3D(arraySize.w,arraySize.h,arraySize.d);
			return true;
		}
		return false;
	}
	bool resize( uint w, uint h, uint d ) {
		return resize(size3(w,h,d));
	}
	void copy( const array3<T> &src ) {
		if( arraySize != src.size() ) {
			resize(src.size());
		}
		FOR_EACH(arraySize.w,arraySize.h,arraySize.d) {
			ptr[i][j][k] = src[i][j][k];
		} END_FOR
	}
	void operator=( const array3<T> &src) {
		copy(src);
	}
	const T** operator[](uint idx) const {
		return (const T**)ptr[idx];
	}
	T** operator[](uint idx) {
		return ptr[idx];
	}
	const T& clampFetch(int i, int j, int k) const {
		if( i < 0 ) i = 0;
		if( i > arraySize.w-1 ) i = arraySize.w-1;
		if( j < 0 ) j = 0;
		if( j > arraySize.h-1 ) j = arraySize.h-1;
		if( k < 0 ) j = 0;
		if( k > arraySize.d-1 ) k = arraySize.d-1;
		return ptr[i][j][k];
	}
	void safeSet(int i, int j, int k, T v) {
		if( i < 0 || i > arraySize.w-1 || j < 0 || j > arraySize.h-1 || k < 0 || k > arraySize.d-1 ) return;
		ptr[i][j][k] = v;
	}
	size3 size() const {
		return arraySize;
	}
	void clear() {
		FOR_EACH(arraySize.w,arraySize.h,arraySize.d) {
			ptr[i][j][k] = T();
		} END_FOR
	}
protected:
	// Do you think we should use vector<vector<vector<T> > > instead ?
	// No, don't do it. vector class is very sensible for parallel operations
	// even though you only access in read-only fashion !
	T *** ptr;
	size3 arraySize;
private:
	T *** alloc3D( uint w, uint h, uint d ) {
		T ***ptr = new T **[w];
		for( int i=0; i<w; i++ ) {
			ptr[i] = new T *[h];
			for( int j=0; j<h; j++ ) {
				ptr[i][j] = new T[d];
				for( int k=0; k<d; k++ ) {
					ptr[i][j][k] = T();
				}
			}
		}
		return ptr;
	}
	void free3D( T ***ptr ) {
		if( ptr ) {
			for( int i=0; i<arraySize.w; i++ ) {
				for( int j=0; j<arraySize.h; j++ ) {
					delete [] ptr[i][j];
				}
				delete [] ptr [i];
			}
			delete [] ptr;
			ptr = NULL;
		}
	}
};

#endif