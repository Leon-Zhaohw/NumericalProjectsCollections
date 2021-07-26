/*
  Fast table lookup
  Jernej Barbic, CMU
  2004,2005

  Avoids the slow floating point to integer cast
  See:
    http://www.cs.cmu.edu/~poladian/papers/report_ning_poladian.pdf
    http://mega-nerd.com/FPcast/
*/

#include "stdafx.h"
#include "twigg/fastTableLookup.h"

// call this once before using 'tableLookup'
// before calling this routine, you must set 'resolution' to the appropriate value
FastTableLookup::FastTableLookup(unsigned int resolution)
{
  this->resolution = resolution;
  // compute logResolution, such that 2^logResolution = resolution
  logResolution = 0;
  unsigned int testVariable = resolution;
  while ( testVariable % 2 == 0 )
  {
    testVariable = testVariable >> 1;
    logResolution++;
  }
  if (testVariable != 1) // error: resolution is not a power of 2
  {
     printf("Error: resolution (%d) is not a power of 2.\n",logResolution);
  }

  lowestImportantBit = 23 - logResolution;
  maskIndex = (resolution-1) << lowestImportantBit;

}


