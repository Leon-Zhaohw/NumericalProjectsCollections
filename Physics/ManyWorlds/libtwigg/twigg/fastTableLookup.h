/*
  Fast table lookup
  Jernej Barbic, CMU
  2004,2005

  Avoids the slow floating point to integer cast
  See:
    http://www.cs.cmu.edu/~poladian/papers/report_ning_poladian.pdf
    http://mega-nerd.com/FPcast/
*/


#ifndef __FAST_TABLE_LOOKUP_
#define __FAST_TABLE_LOOKUP_

#include <stdio.h>
#include "fastTableLookup.h"

// implements a fast 1-D table lookup
class FastTableLookup
{
public:
  // table resolution must be a power of 2
  FastTableLookup(unsigned int resolution);

  // must call this before calling 'tableLookup'
  // call this as many times as you want
  inline void SetTableData(float * tableData) { this->tableData = tableData; }

  // a master function that does the full lookup:
  inline float tableLookup(float lookupPos) const;

  // the "classical" approach for reference:
  inline float conventionalTableLookup(float lookupPos) const;

  // this is mostly an auxiliary function to 'tableLookup'
  // input: lookupPos; assumes lookupPos is a floating point number from the interval 0 to 1
  // output: index of the bucket, coordinate within the bucket (on interval from 0 to 1)
  inline void tableLookupIndices(float lookupPos, unsigned int & bucketIndex, float & coordinateWithinBucket) const;


protected:
    // this is the number of table buckets; it must be a power of 2 for these techniques to work:
    unsigned int resolution;
    // internal parameters
    unsigned char logResolution;
    unsigned int maskIndex;
	unsigned int lowestImportantBit;

	// the table data that will be interpolated
	float * tableData;

};

inline void FastTableLookup::tableLookupIndices(float lookupPos, unsigned int & bucketIndex, float & coordinateWithinBucket) const
{

  // get bucketIndex
  // the line of code below effectively does this:
  //   bucketIndex = (unsigned int)(lookupPos * resolution);

  float pos = lookupPos + 1.0f;

  //   most significant logResolution bits in the mantissa, index mask for 32: 0x7C0000
  bucketIndex = ((unsigned int&)(pos) & maskIndex) >> lowestImportantBit;

  // assert( (bucketIndex >= 0) && (bucketIndex < resolution)
  // note: there might be a degenerate case if lookupPos is precisely 1.0

  // determine offset of lookupPos within the bucket (black magic)
  // this code does the same as:
  // coordinateWithinBucket = lookupPos * resolution - bucketIndex;
  //   barycentric mask for 32: 0xFF83FFFF
  //   to get the barycentric coordinate, has to shift bits 1..23-xe to the left by xe
  //   and then OR it with a mask that sets the bits 24..30 to 1 and 31..32 to zero
  //   exponent = 0, exponent encoding = 0 + 127 = 01111111
  //   sign = positive, sign encoding = 0
  unsigned char temp = logResolution; // won't compile if logResolution used directly below
  float temp2;

  __asm mov eax, pos  // eax = pos
  __asm mov cl, temp
  __asm shl eax, cl
  __asm or  eax, 0x03F800000 // turn on bits 24..30
  __asm and eax, 0x03FFFFFFF // switch off bits 30..31
  __asm mov temp2, eax // coordinateWithinBucket = eax

  coordinateWithinBucket = temp2 - 1.0f; // could merge this into the assembler code above
}


// a master function that does the full lookup:
inline float FastTableLookup::tableLookup(float lookupPos) const
{
  unsigned int bucketIndex;
  float coordinateWithinBucket;

  tableLookupIndices(lookupPos,bucketIndex,coordinateWithinBucket);

  return tableData[bucketIndex] * (1-coordinateWithinBucket) +
	     tableData[bucketIndex+1] * coordinateWithinBucket;
		 
}

// the "classical" slow approach for reference:
inline float FastTableLookup::conventionalTableLookup(float lookupPos) const
{
  unsigned int bucketIndex = (unsigned int)(lookupPos * resolution);
  float coordinateWithinBucket = lookupPos * resolution - bucketIndex;

  return tableData[bucketIndex] * (1-coordinateWithinBucket) +
	     tableData[bucketIndex+1] * coordinateWithinBucket;
}


#endif
