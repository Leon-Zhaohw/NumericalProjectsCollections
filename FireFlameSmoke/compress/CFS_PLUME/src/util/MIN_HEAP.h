// MIN_HEAP.h: interface for the MIN_HEAP class.
//
//////////////////////////////////////////////////////////////////////

#ifndef MIN_HEAP_H
#define MIN_HEAP_H

#include <iostream>
#include <vector>
#include <cmath>
#include <map>

#include "SETTINGS.h"

using namespace std;

class HEAP_ENTRY {
public:
  HEAP_ENTRY() {
    distance = 0;
    index = 0;
    heapIndex = 0;
  };

  Real distance;
  int index;
  int heapIndex;
};

// the actual heap
class MIN_HEAP  
{
public:
  
	MIN_HEAP();
	virtual ~MIN_HEAP();

  // heap ops
  void insert(HEAP_ENTRY& cell);
  //void decreaseKey(HEAP_ENTRY& toChange);
  void decreaseKey(int toChange, Real newKey);
  HEAP_ENTRY popMin();
  HEAP_ENTRY heapMin();
  int size() { return _size; };

  // debugging
  void print();
 
  void clear() {
    _heap.clear();
    _heapIndex.clear();
    _size = 0;
  };

  bool empty() { return _size == 0; };
  
private:
  // tree traversal
  int parent(int i) { return i / 2; };
  int left(int i) { return i * 2; };
  int right(int i) { return i * 2 + 1; };
 
  // the heap
  vector<HEAP_ENTRY> _heap;
  
  // hash table mapping grid index to heap index
  map<int,int> _heapIndex;

  // size of current heap
  int _size;
  
  // enforce heap property 
  void heapify(int index);

  // swap two entries
  void swap(int index1, int index2);
};

#endif
