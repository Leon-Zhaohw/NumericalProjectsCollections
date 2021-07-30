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

#ifndef TRIG_INTEGRAL_3D_H
#define TRIG_INTEGRAL_3D_H

// Compute the integral, index 6: 
// \int_{0}^{\pi}cos(i1*x)cos(mx)dx * 
// \int_{0}^{\pi}sin(i2*y)sin(ny)dy * 
// \int_{0}^{\pi}sin(i3*x)sin(lz)dz
double ComputeIntegralCSS(const int i1, const int i2, const int i3,
                          const int m, const int n, const int l);

// Compute the integral, index 5: 
// \int_{0}^{\pi}sin(i1*x)sin(mx)dx * 
// \int_{0}^{\pi}cos(i2*y)cos(ny)dy * 
// \int_{0}^{\pi}sin(i3*x)sin(lz)dz
double ComputeIntegralSCS(const int i1, const int i2, const int i3,
                          const int m, const int n, const int l);

// Compute the integral, index 4: 
// \int_{0}^{\pi}sin(i1*x)sin(mx)dx * 
// \int_{0}^{\pi}sin(i2*y)sin(ny)dy * 
// \int_{0}^{\pi}cos(i3*x)cos(lz)dz
double ComputeIntegralSSC(const int i1, const int i2, const int i3,
                          const int m, const int n, const int l);
// Compute integral, SSS dot with CSS.
// \int_{0}^{\pi}sin(i1*x)cos(mx)dx *
// \int_{0}^{\pi}sin(i2*y)sin(ny)dy * 
// \int_{0}^{\pi}sin(i3*x)sin(lz)dz
double ComputeIntegralSSSCSS(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l);
// Compute integral, CCS dot with SCS.
// \int_{0}^{\pi}cos(i1*x)sin(mx)dx *
// \int_{0}^{\pi}cos(i2*y)cos(ny)dy * 
// \int_{0}^{\pi}sin(i3*x)sin(lz)dz
double ComputeIntegralCCSSCS(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l);
// Compute integral, CSC dot with SSC.
// \int_{0}^{\pi}cos(i1*x)sin(mx)dx *
// \int_{0}^{\pi}sin(i2*y)sin(ny)dy * 
// \int_{0}^{\pi}cos(i3*x)cos(lz)dz
double ComputeIntegralCSCSSC(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l);

// Compute integral, SSC 4 dot with CSS 6.
// \int_{0}^{\pi}sin(i1*x)cos(mx)dx *
// \int_{0}^{\pi}sin(i2*y)sin(ny)dy * 
// \int_{0}^{\pi}cos(i3*x)sin(lz)dz
double ComputeIntegralSSCCSS(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l);
// 6 with 4.
double ComputeIntegralCSSSSC(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l) ;
                             
// Compute integral, CCC dot with SCS.
// \int_{0}^{\pi}cos(i1*x)sin(mx)dx *
// \int_{0}^{\pi}cos(i2*y)cos(ny)dy * 
// \int_{0}^{\pi}cos(i3*x)sin(lz)dz
double ComputeIntegralCCCSCS(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l);

// Compute integral, SCC 1 dot with CSS 6.
// \int_{0}^{\pi}sin(i1*x)cos(mx)dx *
// \int_{0}^{\pi}cos(i2*y)sin(ny)dy * 
// \int_{0}^{\pi}cos(i3*x)sin(lz)dz
double ComputeIntegralSCCCSS(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l);

// Compute integral, CSC 2 dot with SCS 5.
// \int_{0}^{\pi}cos(i1*x)sin(mx)dx *
// \int_{0}^{\pi}sin(i2*y)cos(ny)dy * 
// \int_{0}^{\pi}cos(i3*x)sin(lz)dz
double ComputeIntegralCSCSCS(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l);

// Compute integral, CCS 3 dot with SSC 4.
// \int_{0}^{\pi}cos(i1*x)sin(mx)dx *
// \int_{0}^{\pi}cos(i2*y)sin(ny)dy * 
// \int_{0}^{\pi}sin(i3*x)cos(lz)dz
double ComputeIntegralCCSSSC(const int i1, const int i2, const int i3,
                             const int m, const int n, const int l);

// Compute the integral, index 6: 
// \int_{0}^{\pi}cos((i1)*x)cos(mx)dx * 
// \int_{0}^{\pi}sin(i2*y)sin(ny)dy * 
// \int_{0}^{\pi}sin(i3*x)sin(lz)dz
double ComputeIntegralCSSDOUBLE(const double i1, const int i2, const int i3,
                              const double m,
                              const int n, const int l);
// Compute the integral, index 5: 
// \int_{0}^{\pi}sin((i1)*x)sin(mx)dx * 
// \int_{0}^{\pi}cos(i2*y)cos(ny)dy * 
// \int_{0}^{\pi}sin(i3*x)sin(lz)dz
double ComputeIntegralSCSDOUBLE(const double i1, const int i2, const int i3,
                          const double m, const int n, const int l);
// Compute the integral, index 4: 
// \int_{0}^{\pi}sin((i1)*x)sin(mx)dx * 
// \int_{0}^{\pi}sin(i2*y)sin(ny)dy * 
// \int_{0}^{\pi}cos(i3*x)cos(lz)dz
double ComputeIntegralSSCDOUBLE(const double i1, const int i2, const int i3,
                          const double m, const int n, const int l);

// Used by 2D
// \int_{0}^{\pi}cos(i1*x)sin(mx)dx *
// \int_{0}^{\pi}sin(i2*y)sin(ny)dy
double ComputeIntegralCSSS(const int i1, const int i2, const int m, const int n);

#endif  // TRIG_INTEGRAL_3D_H
