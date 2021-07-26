#ifndef MXHEAP_INCLUDED
#define MXHEAP_INCLUDED

/************************************************************************
 * Heap                                                                 *
 * Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.  *
 ************************************************************************/

#include <iostream>
#include <cmath>

/*******************************************************/
// stdmix.h

#ifndef _STDMIX_H_
#define _STDMIX_H_

#include <cassert>

#if SAFETY > 0
#  define AssertBound(t)  assert(t)
#  define PRECAUTION(x) x
#else
#  define AssertBound(t)
#  define PRECAUTION(x)
#endif

#endif

/*******************************************************/

#define SanityCheck(t)  assert(t)

typedef unsigned int uint;

enum MxAllocPolicy { MX_ALLOC_INPLACE };

//#define MIN(x,y) ((x)<(y))?(x):(y)

inline int MIN(int n, int m){return (n < m) ? n : m;}

inline void *operator new(size_t, void *spot, MxAllocPolicy)
{
  return spot;  // MX_ALLOC_INPLACE
}

//
// Default equality class
//
class MxEQ 
{
 public: bool operator==(const MxEQ& x) const { return this==&x; }
};

template<class T>
class MxBlock : public MxEQ
{
 private:
  uint N;
  T *block;
  
 protected:
  MxBlock() { }

  void init_block(uint n)
    {
      // Allocate memory for block
      N=n; block = (T *)malloc(sizeof(T)*n);
      
      // Initialize each entry
      for(uint i=0; i<n; i++) new((void *)&block[i], MX_ALLOC_INPLACE) T;
    }

  void resize_block(uint n)
    {
      T *old = block;
      
      // Allocate new block
      block = (T *)realloc(old, sizeof(T)*n);
      
      // Initialize all the newly allocated entries
      for(uint i=N; i<n; i++) new((void *)&block[i], MX_ALLOC_INPLACE) T;
      
      N = n;
    }
 
  void free_block()
    { free(block); }
  
  T&       raw(uint i)       { return block[i]; }
  const T& raw(uint i) const { return block[i]; }
  
 public:
  MxBlock(uint n) { init_block(n); }
  
  virtual ~MxBlock() { free_block(); }
  
  // Standard accessors
  //
  T&       ref(uint i)              { AssertBound(i<N);return block[i]; }
  const T& ref(uint i) const        { AssertBound(i<N);return block[i]; }
  T&       operator()(uint i)       { return ref(i);   }
  const T& operator()(uint i) const { return ref(i);   }

#ifdef __GNUC__
  // Most compilers can apply the casting operators and then use the
  // standard array [].  GCC (2.7.2) seems to need explicit operators.
  //
  T&       operator[](uint i)       { return block[i]; }
  const T& operator[](uint i) const { return block[i]; }
#endif
  
  operator const T*() const { return block; }
  operator       T*()       { return block; }
  uint length() const { return N; }
  
  // Primitive methods for altering the data block
  //
  void resize(uint n) { resize_block(n); }
  void bitcopy(const T *array, uint n) // copy bits directly
    { memcpy(block, array, MIN(n,N)*sizeof(T)); }
  void copy(const T *array, const uint n) // copy using assignment operator
    { for(uint i=0; i<MIN(n,N); i++) block[i] = array[i]; }
  void bitcopy(const MxBlock<T>& b) { bitcopy(b, b.length()); }
  void copy(const MxBlock<T>& b) { copy(b, b.length()); }
};

template<class T>
class MxDynBlock : public MxBlock<T>
{
 private:
  uint fill;
  
 public:
  MxDynBlock() {MxBlock<T>::init_block(8); fill=0; }
  MxDynBlock(uint n) {MxBlock<T>::init_block(n); fill=0; }

  virtual ~MxDynBlock() { }
  
  T&       operator()(uint i)       { AssertBound(i<fill); return MxBlock<T>::ref(i); }
  const T& operator()(uint i) const { AssertBound(i<fill); return MxBlock<T>::ref(i); }
  
  uint length() const { return fill; }
  uint total_space() const { return MxBlock<T>::length(); }
  T& last() { AssertBound(fill>0); return MxBlock<T>::raw(fill-1); }
  const T& last() const { AssertBound(fill>0); return MxBlock<T>::raw(fill-1); }
  
  uint add()
    {
      if( length()==total_space() ) resize(total_space() * 2);
      return fill++;
    }
  
  uint add(const T& t)
    {
      uint i=add();
      MxBlock<T>::raw(i) = t;
      return i;
    }
  
  void reset() { fill = 0; }
  T& drop() { return MxBlock<T>::raw(--fill); }
  void drop(uint d) { fill -= d; }
  bool find(const T& t, uint *index=NULL)
    {
      for(uint i=0; i<fill; i++)
	if( MxBlock<T>::raw(i) == t )
	  {
	    if( index ) *index = i;
	    return true;
	  }
      return false;
    }
  void remove(uint i)
    { AssertBound(i<fill); fill--; MxBlock<T>::raw(i) = MxBlock<T>::raw(fill); }
  void remove_inorder(uint i)
    { memmove(&MxBlock<T>::raw(i), &MxBlock<T>::raw(i+1), (--fill - i)*sizeof(T)); }
};

template<class T, int T_SIZE>
  class MxSizedDynBlock : public MxDynBlock<T>
{
 public:   
  MxSizedDynBlock(uint n=T_SIZE) : MxDynBlock<T>(n) { }
};

class MxHeapable
{
 private:
  double import;
  int token;
  
 public:
  MxHeapable() { not_in_heap(); heap_key(0.0f); }
  
  inline bool is_in_heap() { return token != -47; }
  inline void not_in_heap() { token = -47; }
  inline int get_heap_pos() { return token; }
  inline void set_heap_pos(int t) { token=t; }
  
  inline void  heap_key(double k) { import=k; }
  inline double heap_key() const  { return import; }
};

class MxHeap : private MxDynBlock<MxHeapable *>
{
 private:
  bool type_; // true - maximum | false - minimum

  inline bool cmp(double a, double b) { return (type_)? (a>b) : (a<b); }
  
  void place(MxHeapable *x, unsigned int i);
  void swap(unsigned int i, unsigned int j);
  
  unsigned int parent(unsigned int i) { return (i-1)/2; }
  unsigned int left(unsigned int i) { return 2*i+1; }
  unsigned int right(unsigned int i) { return 2*i+2; }
  
  void upheap(unsigned int i);
  void downheap(unsigned int i);
  
 public:
  MxHeap(bool t=true) : MxDynBlock<MxHeapable *>(8), type_(t) { }
  MxHeap(unsigned int n, bool t=true) : MxDynBlock<MxHeapable *>(n), type_(t) { }
  
  void insert(MxHeapable *t) { insert(t, t->heap_key()); }
  void insert(MxHeapable *, double);
  void update(MxHeapable *t) { update(t, t->heap_key()); }
  void update(MxHeapable *, double);
  
  unsigned int size() const { return length(); }
  MxHeapable       *item(uint i)       { return ref(i); }
  const MxHeapable *item(uint i) const { return ref(i); }
  MxHeapable *extract();
  MxHeapable *top() { return (length()<1 ? (MxHeapable *)NULL : raw(0)); }
  MxHeapable *remove(MxHeapable *);
};

// MXHEAP_INCLUDED
#endif
