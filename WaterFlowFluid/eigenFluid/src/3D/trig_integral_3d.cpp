/*
 This file is part of SSFR (Zephyr).
 
 Zephyr is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Zephyr is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Zephyr.  If not, see <http://www.gnu.org/licenses/>.
 
 Copyright 2018 Qiaodong Cui (qiaodong@ucsb.edu)
 */

#include <math.h>

#include "3D/trig_integral_3d.h"
#include "util/util.h"

#define INTEGRAL_SS(c_, i, n) \
int c_ = 0; \
if ((i - n) == 0) { \
  c_ ++; \
} \
if ((i + n) == 0) { \
  c_ --; \
} \
if (c_ == 0) return 0;
 
#define INTEGRAL_CC(c_, i, n) \
int c_ = 0; \
if ((i + n) == 0) {\
  c_ ++; \
} \
if ((i - n) == 0) {\
    c_ ++; \
} \
if (c_ == 0) return 0;

#define INTEGRAL_SC(c_, i, n) \
double c_ = 0; \
if ((i + n) != 0 && ((i + n) % 2 != 0)) { \
  c_ += 1.0 / (static_cast<double>(i + n)); \
} \
if ((i - n) != 0 && ((i - n) % 2 != 0)) { \
  c_ += 1.0 / (static_cast<double>(i - n)); \
}

#define INTEGRAL_CS(c_, i, n) \
double c_ = 0; \
if ((i + n) != 0 && ((i + n) % 2 != 0)) { \
  c_ += 1.0 / (static_cast<double>(i + n)); \
} \
if ((i - n) != 0 && ((i - n) % 2 != 0)) { \
  c_ -= 1.0 / (static_cast<double>(i - n)); \
}

namespace {
// This function takes a xt as input, xt can only be integer, or integer offset by 0.5;
// Evaluate sin(xt*M_PI);
inline double CheapSin(double xt) {
  const int dx = static_cast<int>(xt*2.0);
  if (dx % 2 == 0) {    // integer
    return 0;
  } else {  // integer offset by 0.5
    if (dx % 4 == 1) {
      return 1;
    } else {
      return -1;
    }
  }
}

}  // namespace

double ComputeIntegralSSSCSS(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l) {
  INTEGRAL_SS(cy, i2, n)
  INTEGRAL_SS(cz, i3, l)
  double cx = 0;
  if ((i1 + m) != 0 && ((i1 + m) % 2 != 0)) {
    cx += 1.0 / (static_cast<double>(i1 + m));
  }
  if ((i1 - m) != 0 && ((i1 - m) % 2 != 0)) {
    cx += 1.0 / (static_cast<double>(i1 - m));
  }
  return cx*cy*cz*0.25*PI_SQUARE;
}

double ComputeIntegralCCSSCS(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l) {
  INTEGRAL_CC(cy, i2, n)
  INTEGRAL_SS(cz, i3, l)
  double cx = 0;
  if ((i1 + m) != 0 && ((i1 + m) % 2 != 0)) {
    cx += 1.0 / (static_cast<double>(i1 + m));
  }
  if ((i1 - m) != 0 && ((i1 - m) % 2 != 0)) {
    cx -= 1.0 / (static_cast<double>(i1 - m));
  }
  return cx*cy*cz*0.25*PI_SQUARE;
}

double ComputeIntegralCSCSSC(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l) {
  INTEGRAL_SS(cy, i2, n)
  INTEGRAL_CC(cz, i3, l)
  double cx = 0;
  if ((i1 + m) != 0 && ((i1 + m) % 2 != 0)) {
    cx += 1.0 / (static_cast<double>(i1 + m));
  }
  if ((i1 - m) != 0 && ((i1 - m) % 2 != 0)) {
    cx -= 1.0 / (static_cast<double>(i1 - m));
  }
  return cx*cy*cz*0.25*PI_SQUARE;
}

// Compute integral, SSC dot with CSS.
// \int_{0}^{\pi}sin(i1*x)cos(mx)dx *
// \int_{0}^{\pi}sin(i2*y)sin(ny)dy * 
// \int_{0}^{\pi}cos(i3*x)sin(lz)dz
double ComputeIntegralSSCCSS(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l) {
  INTEGRAL_SS(cy, i2, n)
  INTEGRAL_SC(cx, i1, m)
  INTEGRAL_CS(cz, i3, l)
  return cx*cy*cz*0.5*M_PI;
}

double ComputeIntegralCSSSSC(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l) {
  return ComputeIntegralSSCCSS(m, n, l, i1, i2, i3);
}
// Compute integral, CCC dot with SCS.
// \int_{0}^{\pi}cos(i1*x)sin(mx)dx *
// \int_{0}^{\pi}cos(i2*y)cos(ny)dy * 
// \int_{0}^{\pi}cos(i3*x)sin(lz)dz
double ComputeIntegralCCCSCS(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l) {
  INTEGRAL_CS(cx, i1, m)
  INTEGRAL_CC(cy, i2, n)
  INTEGRAL_CS(cz, i3, l)
  return cx*cy*cz*0.5*M_PI;
}
// Compute integral, SCC 1 dot with CSS 6.
// \int_{0}^{\pi}sin(i1*x)cos(mx)dx *
// \int_{0}^{\pi}cos(i2*y)sin(ny)dy * 
// \int_{0}^{\pi}cos(i3*x)sin(lz)dz
double ComputeIntegralSCCCSS(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l) {
  INTEGRAL_SC(cx, i1, m)
  INTEGRAL_CS(cy, i2, n)
  INTEGRAL_CS(cz, i3, l)
  return cx*cy*cz;
}

// Compute integral, CSC 2 dot with SCS 5.
// \int_{0}^{\pi}cos(i1*x)sin(mx)dx *
// \int_{0}^{\pi}sin(i2*y)cos(ny)dy * 
// \int_{0}^{\pi}cos(i3*x)sin(lz)dz
double ComputeIntegralCSCSCS(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l) {
  INTEGRAL_CS(cx, i1, m)
  INTEGRAL_SC(cy, i2, n)
  INTEGRAL_CS(cz, i3, l)
  return cx*cy*cz;
}

// Compute integral, CCS 3 dot with SSC 4.
// \int_{0}^{\pi}cos(i1*x)sin(mx)dx *
// \int_{0}^{\pi}cos(i2*y)sin(ny)dy * 
// \int_{0}^{\pi}sin(i3*x)cos(lz)dz
double ComputeIntegralCCSSSC(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l) {
  INTEGRAL_CS(cx, i1, m)
  INTEGRAL_CS(cy, i2, n)
  INTEGRAL_SC(cz, i3, l)
  return cx*cy*cz;
}

// \int_{0}^{\pi}cos(i1*x)cos(mx)dx * ...
// \int_{0}^{\pi}sin(i2*y)sin(ny)dy * 
// \int_{0}^{\pi}sin(i3*x)sin(lz)dz
double ComputeIntegralCSSDOUBLE(const double i1, const int i2, const int i3,
                              const double m,
                              const int n, const int l) {
  INTEGRAL_SS(cy, i2, n)
  INTEGRAL_SS(cz, i3, l)
  double cx = 0;
  if (std::abs(i1 + m) < 1e-10) {
    cx += 0.5*M_PI;
  } else {
    cx += 0.5*CheapSin((i1+m)) / (i1 + m);
  }
  if (std::abs(i1 - m) < 1e-10) {
    cx += 0.5*M_PI;
  } else {
    cx += 0.5*CheapSin((i1-m)) / (i1 - m);
  }
  return cx*cy*cz*0.25*PI_SQUARE;

}

// \int_{0}^{\pi}sin((i1)*x)sin(mx)dx * 
// \int_{0}^{\pi}cos(i2*y)cos(ny)dy * 
// \int_{0}^{\pi}sin(i3*x)sin(lz)dz
double ComputeIntegralSCSDOUBLE(const double i1, const int i2, const int i3,
                          const double m, const int n, const int l) {
  INTEGRAL_CC(cy, i2, n)
  INTEGRAL_SS(cz, i3, l)
  double cx = 0;
  if (std::abs(i1 - m) < 1e-10) {
    cx += 0.5*M_PI;
  } else {
    cx += 0.5*CheapSin((i1 - m)) / (i1 - m);
  }
  if (std::abs(i1 + m) < 1e-10) {
    cx -= 0.5*M_PI;
  } else {
    cx -= 0.5*CheapSin((i1 + m)) / (i1 + m);
  }
  return cx*cy*cz*0.25*PI_SQUARE;
}

// \int_{0}^{\pi}sin((i1)*x)sin(mx)dx * 
// \int_{0}^{\pi}sin(i2*y)sin(ny)dy * 
// \int_{0}^{\pi}cos(i3*x)cos(lz)dz
double ComputeIntegralSSCDOUBLE(const double i1, const int i2, const int i3,
                          const double m, const int n, const int l) {
  INTEGRAL_SS(cy, i2, n)
  INTEGRAL_CC(cz, i3, l)
  double cx = 0;
  if (std::abs(i1 - m) < 1e-10) {
    cx += 0.5*M_PI;
  } else {
    cx += 0.5*CheapSin((i1 - m)) / (i1 - m);
  }
  if (std::abs(i1 + m) < 1e-10) {
    cx -= 0.5*M_PI;
  } else {
    cx -= 0.5*CheapSin((i1 + m)) / (i1 + m);
  }
  return cx*cy*cz*0.25*PI_SQUARE;
}

// type 4.
double ComputeIntegralSSC(const int i1, const int i2, const int i3,
                          const int m, const int n, const int l) {
  int cx = 0;
  if ((i1 - m) == 0) {
    cx ++;
  }
  if ((i1 + m) == 0) {
    cx --;
  }
  if (cx == 0) return 0;
  
  int cy = 0;
  if ((i2 - n) == 0) {
    cy ++;
  }
  if ((i2 + n) == 0) {
    cy --;
  }
  if (cy == 0) return 0;
  
  int cz = 0;
  if ((i3 + l) == 0) {
    cz ++;
  }
  if ((i3 - l) == 0) {
    cz ++;
  }
  if (cz == 0) return 0;
  
  return PI_CUBE*0.125*cx*cy*cz;
}

// type 5
double ComputeIntegralSCS(const int i1, const int i2, const int i3,
                          const int m, const int n, const int l) {
  int cx = 0;
  if ((i1 - m) == 0) {
    cx ++;
  }
  if ((i1 + m) == 0) {
    cx --;
  }
  if (cx == 0) return 0;
  
  int cy = 0;
  if ((i2 + n) == 0) {
    cy ++;
  }
  if ((i2 - n) == 0) {
    cy ++;
  }
  if (cy == 0) return 0;
  
  int cz = 0;
  if ((i3 - l) == 0) {
    cz ++;
  }
  if ((i3 + l) == 0) {
    cz --;
  }
  if (cz == 0) return 0;
  
  return PI_CUBE*0.125*cx*cy*cz;
}

// type 6
double ComputeIntegralCSS(const int i1, const int i2, const int i3,
                          const int m, const int n, const int l) {
  int cx = 0;
  if ((i1 + m) == 0) {
    cx ++;
  }
  if ((i1 - m) == 0) {
    cx ++;
  }
  if (cx == 0) return 0;
  
  int cy = 0;
  if ((i2 - n) == 0) {
    cy ++;
  }
  if ((i2 + n) == 0) {
    cy --;
  }
  if (cy == 0) return 0;
  
  int cz = 0;
  if ((i3 - l) == 0) {
    cz ++;
  }
  if ((i3 + l) == 0) {
    cz --;
  }
  if (cz == 0) return 0;
  
  return PI_CUBE*0.125*cx*cy*cz;
}

double ComputeIntegralCSSS(const int i1, const int i2, const int m, const int n) {
  INTEGRAL_CS(cx, i1, m)
  INTEGRAL_SS(cy, i2, n)
  return cx*cy*0.5*M_PI;
}

#undef INTEGRAL_SS
#undef INTEGRAL_CC
#undef INTEGRAL_SC
#undef INTEGRAL_CS
