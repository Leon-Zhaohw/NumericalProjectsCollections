#ifndef INTERP_H
#define INTERP_H

#include <cmath>

template<class T>
inline void get_barycentric(T x, int &i, T &f, int i_low, int i_high)
{
   double s=std::floor(x);
   i=(int)s;
   if(i<i_low){
      i=i_low;
      f=0;
   }else if(i>i_high-2){
      i=i_high-2;
      f=1;
   }else{
      f=(T)(x-s);
   }
}

template<class S,class T>
inline S lerp(const S &value0, const S &value1, T f)
{ return (1-f)*value0 + f*value1; }

template<class S,class T>
inline S bilerp(const S &v00, const S &v10, 
                const S &v01, const S &v11, 
                T fx, T fy)
{ 
   return lerp(    
               lerp(v00, v10, fx),
               lerp(v01, v11, fx), 
               fy);
}

template<class S,class T>
inline S trilerp(
   const S& v000, const S& v100,
   const S& v010, const S& v110,
   const S& v001, const S& v101,  
   const S& v011, const S& v111,
   T fx, T fy, T fz) 
{
   return lerp(
            bilerp(v000, v100, v010, v110, fx, fy),
            bilerp(v001, v101, v011, v111, fx, fy),
            fz);
}

template<class S,class T>
inline S quadlerp(
   const S& v0000, const S& v1000,
   const S& v0100, const S& v1100,
   const S& v0010, const S& v1010,  
   const S& v0110, const S& v1110,
   const S& v0001, const S& v1001,
   const S& v0101, const S& v1101,
   const S& v0011, const S& v1011,  
   const S& v0111, const S& v1111,
   T fx, T fy, T fz, T ft) 
{
   return lerp(
            trilerp(v0000, v1000, v0100, v1100, v0010, v1010, v0110, v1110, fx, fy, fz),
            trilerp(v0001, v1001, v0101, v1101, v0011, v1011, v0111, v1111, fx, fy, fz),
            ft);
}

#endif
