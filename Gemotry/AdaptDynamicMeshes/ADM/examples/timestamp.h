#ifndef timestamp_h
#define timestamp_h 1

#include <sys/time.h>

class Timestamp 
{
 protected:
  int S, U;

 public:
  Timestamp() 
    { now(); }

  Timestamp & operator=(Timestamp& t) 
    { S = t.S; U=t.U; return (*this); }

  double operator-(Timestamp& t) 
    {
      int msec;
      double sec;
      msec = 1000 * ( S - t.S );
      msec += (U/1000 - t.U/1000);
      sec = double(msec / 1000.0);
      return sec;
    }

  void now() 
    { 
      struct timeval tv;
      gettimeofday( &tv, 0 );
      S = (int) tv.tv_sec;
      U = (int) tv.tv_usec;
    }
};

#include <time.h>

class CPUtimer
{
 protected:
  clock_t time_;

 public:
  CPUtimer()
    { now(); }

  CPUtimer& operator=(CPUtimer& t) 
    { time_ = t.time_; return (*this); }

  double operator-(CPUtimer& t) 
    { return (time_ - t.time_) / (CLOCKS_PER_SEC); }
  
  void now() 
    { time_ = clock(); }
};

#endif
