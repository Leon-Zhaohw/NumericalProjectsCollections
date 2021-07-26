/*
This file is part of Cubica.
 
Cubica is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Cubica is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Cubica.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef TIMER_H
#define TIMER_H

#include <iostream>
#include <string>

//////////////////////////////////////////////////////////////////////
//  Windows counter
//  Uses QueryPerformanceCounter for microsecond accuracy.
//////////////////////////////////////////////////////////////////////

#ifdef _WIN32
#include <windows.h>
#include <winnt.h>
#else
#include <sys/time.h>
#endif

class TIMER
{
	public:
    // start the timer by default -- if a tick is called later,
    // it will just stomp it
		TIMER() {
#if _WIN32
      if (_firstCall) {
        LARGE_INTEGER ticksPerSecond;
        QueryPerformanceFrequency(&ticksPerSecond);
        _divisor = 1.0f / (double) ticksPerSecond.QuadPart;
        _firstCall = false;
      }
      QueryPerformanceCounter(&_tick);
#else
      error = gettimeofday(&_tick, 0);
#endif
    };
		~TIMER() {};

    double timing() {
#if _WIN32
      QueryPerformanceCounter(&_tock);
      return (double) (_tock.QuadPart - _tick.QuadPart) * _divisor;
#else
      if (error) return 0.0;
      error = gettimeofday(&_tock, 0);

      double beginTime = (double)_tick.tv_sec + 1e-6 * _tick.tv_usec;
      double endTime = (double)_tock.tv_sec + 1e-6 * _tock.tv_usec;
      return endTime - beginTime;
#endif
    };
    static int hours(int seconds) { return seconds / (60 * 60); };
    static int minutes(int seconds) {
     int mod = seconds % (60 * 60);
     return mod / 60;
    };
    static int seconds(int seconds) {
      int mod = seconds % (60 * 60);
      return mod % 60;
    };

	private:
#if _WIN32
		LARGE_INTEGER _tick;
		LARGE_INTEGER _tock;
		static double _divisor;
    static bool _firstCall;
#else
    timeval _tick;
    timeval _tock;
    int error;
#endif
};

#endif
