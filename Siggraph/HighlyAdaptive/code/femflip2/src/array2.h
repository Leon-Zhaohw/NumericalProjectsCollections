/*
 *	array2.h
 *	
 *	Created by Ryoichi Ando on 11/3/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <stdio.h>
#include <vector>
#include "macros.h"

#ifndef _ARRAY2_H
#define _ARRAY2_H

struct size2 {
	size2() : w(0), h(0) {}
	size2( uint w, uint h ) {
		this->w = w;
		this->h = h;
	}
	uint w;
	uint h;
	uint& operator [](uint idx) {
		return idx == 0 ? w : h;
	}
	bool operator==(const size2 &size ) {
		return w == size.w && h == size.h;
	}
	bool operator!=(const size2 &size ) {
		return 1 - ((*this) == size);
	}
	bool validIndex( int i, int j ) {
		return i>=0 && i<w && j>=0 && j<h;
	}
};

template <class T> class array2 {
public:
	array2() : ptr(NULL) {}
	array2( const array2<T> &array ) : ptr(NULL) {
		copy(array);
	}
	array2( int w, int h ) : ptr(NULL) {
		resize(w,h);
	}
	array2( size2 theSize ) : ptr(NULL) {
		resize(theSize);
	}
	~array2() {
		free2D(ptr);
	}
	bool resize( size2 theSize ) {
		if( arraySize != theSize ) {
			if( ptr ) free2D(ptr);
			arraySize = theSize;
			ptr = alloc2D(arraySize.w,arraySize.h);
			return true;
		}
		return false;
	}
	bool resize( uint w, uint h ) {
		return resize(size2(w,h));
	}
	void copy( const array2<T> &src ) {
		if( arraySize != src.size() ) {
			resize(src.size());
		}
		FOR_EACH(arraySize.w,arraySize.h) {
			ptr[i][j] = src[i][j];
		} END_FOR
	}
	void operator=( const array2<T> &src) {
		copy(src);
	}
	const T* operator[](uint idx) const {
		return ptr[idx];
	}
	T* operator[](uint idx) {
		return ptr[idx];
	}
	const T& clampFetch(int i, int j) const {
		if( i < 0 ) i = 0;
		if( i > arraySize.w-1 ) i = arraySize.w-1;
		if( j < 0 ) j = 0;
		if( j > arraySize.h-1 ) j = arraySize.h-1;
		return ptr[i][j];
	}
	void safeSet(int i, int j, T v) {
		if( i < 0 || i > arraySize.w-1 || j < 0 || j > arraySize.h-1 ) return;
		ptr[i][j] = v;
	}
	size2 size() const {
		return arraySize;
	}
	void clear() {
		FOR_EACH(arraySize.w,arraySize.h) {
			ptr[i][j] = T();
		} END_FOR
	}
protected:
	// Do you think we should use vector<vector<T> > instead ?
	// No, don't do it. vector class is very sensible for parallel operations
	// even though you only access in read-only fashion !
	T ** ptr;
	size2 arraySize;
private:
	T ** alloc2D( uint w, uint h ) {
		T **ptr = new T *[w];
		for( int i=0; i<w; i++ ) {
			ptr[i] = new T[h];
			for( int j=0; j<h; j++ ) ptr[i][j] = T();
		}
		return ptr;
	}
	void free2D( T **ptr ) {
		if( ptr ) {
			for( int i=0; i<arraySize.w; i++ ) delete [] ptr[i];
			delete [] ptr;
			ptr = NULL;
		}
	}
};

#endif