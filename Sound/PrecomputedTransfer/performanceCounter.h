#ifndef _PERFORMANCECOUNTER_H_
#define _PERFORMANCECOUNTER_H_

/**************** LINUX COUNTER *******************/

#ifdef __unix__

#include "stdlib.h"
#include "sys/time.h"

/*

Jernej Barbic
Performance counter to measure program execution time

Graphics Lab, CMU, 2004

Uses internal linux timers to measure program execution time

*/

class PerformanceCounter
{
  protected:

  long startCountSec,stopCountSec,startCountMicroSec,stopCountMicroSec;

  public:

  PerformanceCounter()
  {
    // also, reset the starting counter
    StartCounter();
  }

  void StartCounter(); // can reset multiple times
  void StopCounter(); // after stopping the counter, read the elapsed value by GetElapsedTime

  // returns elapsed time in seconds, with an accuracy almost up to the system clock
  // must be preceded by a call to StopCounter
  // will return the elapsed time between the last invocation of StartCounter-StopCounter pair
  double GetElapsedTime();

};

inline void PerformanceCounter::StartCounter()
{
  struct timeval tv;

  gettimeofday(&tv,NULL);

  startCountSec = tv.tv_sec;
  startCountMicroSec = tv.tv_usec;

}

inline void PerformanceCounter::StopCounter()
{
  struct timeval tv;

  gettimeofday(&tv,NULL);

  stopCountSec = tv.tv_sec;
  stopCountMicroSec = tv.tv_usec;
}


inline double PerformanceCounter::GetElapsedTime()
{
  float elapsedTime = 1.0 * (stopCountSec-startCountSec) + 1E-6 * (stopCountMicroSec - startCountMicroSec);
  return elapsedTime;
}

#else

/**************** WINDOWS COUNTER *******************/

#include <windows.h>

/*

Jernej Barbic
Performance counter to measure program execution time

Graphics Lab, CMU, 2004

Uses internal Windows timers to measure program execution time

*/

class PerformanceCounter
{
  protected:

  LARGE_INTEGER timerFrequency;
  LARGE_INTEGER startCount,stopCount;

  public:

  PerformanceCounter() 
  {
    // reset the counter frequency
    QueryPerformanceFrequency(&timerFrequency);
    // also, reset the starting counter
    StartCounter();
  }

  void StartCounter(); // can reset multiple times
  void StopCounter(); // after stopping the counter, read the elapsed value by GetElapsedTime

  // returns elapsed time in seconds, with an accuracy almost up to the system clock
  // must be preceded by a call to StopCounter
  // will return the elapsed time between the last invocation of StartCounter-StopCounter pair
  float GetElapsedTime(); 

};

inline void PerformanceCounter::StartCounter()
{
  QueryPerformanceCounter(&startCount);
}

inline void PerformanceCounter::StopCounter()
{
  QueryPerformanceCounter(&stopCount);
}


inline float PerformanceCounter::GetElapsedTime()
{
	return ((float)(stopCount.QuadPart - startCount.QuadPart))
	/((float)timerFrequency.QuadPart);
}

#endif

#endif

