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

#include "3D/trig_integral_3d.h"
#include "util/util.h"

int main(int argc, char ** argv) {
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  
  if (ComputeIntegralCSS(0 != 0,0,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCSS(1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCSS(-1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCSS(-1 != 1,1,1,1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCSS(-1 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralCSS(-1 != 1,-1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralCSS(1 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralCSS(1 != 1,1,1,-1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralCSS(0 != -1,1,0,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralCSS(0 != 1,1,0,1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralCSS(0 != -1,1,0,-1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralCSS(0 != -1,1,0,1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.25) << std::endl; exit(0);}
  
  if (ComputeIntegralSCS(0 != 0,0,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSCS(1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSCS(-1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSCS(-1 != 1,1,1,1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSCS(-1 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralSCS(-1 != 1,-1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralSCS(1 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralSCS(1 != 1,1,1,-1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralSCS(-1 != 0,1,1,0,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralSCS(1 != 0,1,1,0,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralSCS(-1 != 0,1,-1,0,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralSCS(-1 != 0,1,1,0,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.25) << std::endl; exit(0);}
  
  if (ComputeIntegralSSC(0 != 0,0,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSSC(1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSSC(-1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSSC(-1 != 1,1,1,1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSSC(-1 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralSSC(-1 != 1,-1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralSSC(1 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralSSC(1 != 1,1,1,-1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralSSC(-1 != 1,0,1,1,0) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralSSC(1 != 1,0,1,-1,0) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralSSC(1 != -1,0,1,-1,0) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralSSC(-1 != 1,0,1,-1,0) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.25) << std::endl; exit(0);}
  
  if (ComputeIntegralSSSCSS(1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 0) << std::endl; exit(0);}
  if (ComputeIntegralSSSCSS(1 != 1,1,2,1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 0) << std::endl; exit(0);}
  if (ComputeIntegralSSSCSS(1 != 1,1,2,2,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 0) << std::endl; exit(0);}
  if (ComputeIntegralSSSCSS(1 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 0) << std::endl; exit(0);}
  if (ComputeIntegralSSSCSS(2 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 4.0/3.0*0.25*PI_SQUARE) << std::endl; exit(0);}
  if (ComputeIntegralSSSCSS(2 != 1,1,1,1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , -4.0/3.0*0.25*PI_SQUARE) << std::endl; exit(0);}
  if (ComputeIntegralSSSCSS(2 != 1,1,1,-1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , -4.0/3.0*0.25*PI_SQUARE) << std::endl; exit(0);}
  if (ComputeIntegralSSSCSS(2 != 1,1,1,-1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 4.0/3.0*0.25*PI_SQUARE) << std::endl; exit(0);}
  if (ComputeIntegralSSSCSS(2 != 1,1,4,-1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 0) << std::endl; exit(0);}
  if (ComputeIntegralSSSCSS(2 != 1,1,8,-1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 0) << std::endl; exit(0);}
  if (ComputeIntegralSSSCSS(2 != 1,1,3,-1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , -4.0/5.0*0.25*PI_SQUARE) << std::endl; exit(0);}
  if (ComputeIntegralSSSCSS(3 != 1,1,2,-1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 6.0/5.0*0.25*PI_SQUARE) << std::endl; exit(0);}
  
  if (ComputeIntegralCCSSCS(1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 0) << std::endl; exit(0);}
  if (ComputeIntegralCCSSCS(1 != 1,1,2,1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 0) << std::endl; exit(0);}
  if (ComputeIntegralCCSSCS(1 != 1,1,2,2,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 0) << std::endl; exit(0);}
  if (ComputeIntegralCCSSCS(1 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 0) << std::endl; exit(0);}
  // std::cout <<  ComputeIntegralCCSSCS(2,1,1,1,1,1) + 2.0/3.0*0.25*PI_SQUARE << std::endl;
  // if (ComputeIntegralCCSSCS(2 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , -2.0/3.0*0.25*PI_SQUARE) << std::endl; exit(0);}
  // if (ComputeIntegralCCSSCS(2 != 1,1,1,1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 2.0/3.0*0.25*PI_SQUARE) << std::endl; exit(0);}
  // if (ComputeIntegralCCSSCS(2 != 1,1,1,-1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , -2.0/3.0*0.25*PI_SQUARE) << std::endl; exit(0);}
  // if (ComputeIntegralCCSSCS(2 != 1,1,1,-1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 2.0/3.0*0.25*PI_SQUARE) << std::endl; exit(0);}
  // if (ComputeIntegralCCSSCS(2 != 0,1,1,0,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 4.0/3.0*0.25*PI_SQUARE) << std::endl; exit(0);}

  if (ComputeIntegralCSCSSC(1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 0) << std::endl; exit(0);}
  if (ComputeIntegralCSCSSC(1 != 1,1,2,1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 0) << std::endl; exit(0);}
  if (ComputeIntegralCSCSSC(1 != 1,1,2,2,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 0) << std::endl; exit(0);}
  if (ComputeIntegralCSCSSC(1 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 0) << std::endl; exit(0);}
 // if (ComputeIntegralCSCSSC(2 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , -2.0/3.0*0.25*PI_SQUARE) << std::endl; exit(0);}
  
  if (ComputeIntegralCSSDOUBLE(0 != 0,0,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCSSDOUBLE(1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCSSDOUBLE(-1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCSSDOUBLE(-1 != 1,1,1,1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCSSDOUBLE(-1 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralCSSDOUBLE(-1 != 1,-1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralCSSDOUBLE(1 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralCSSDOUBLE(1 != 1,1,1,-1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralCSSDOUBLE(0 != -1,1,0,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralCSSDOUBLE(0 != 1,1,0,1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralCSSDOUBLE(0 != -1,1,0,-1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralCSSDOUBLE(0 != -1,1,0,1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralCSSDOUBLE(0.5 != -1,1,0.5,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralCSSDOUBLE(0.5 != 1,1,1.5,1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCSSDOUBLE(1 != -1,1,1.5,-1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 0.25*PI_SQUARE*1.2) << std::endl; exit(0);}
  if (ComputeIntegralCSSDOUBLE(1.5 != -1,1,2,1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0.25*PI_SQUARE*(1 - 1.0/7)) << std::endl; exit(0);}
  
  if (ComputeIntegralSCSDOUBLE(0 != 0,0,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSCSDOUBLE(1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSCSDOUBLE(-1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSCSDOUBLE(-1 != 1,1,1,1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSCSDOUBLE(-1 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralSCSDOUBLE(-1 != 1,-1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralSCSDOUBLE(1 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralSCSDOUBLE(1 != 1,1,1,-1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralSCSDOUBLE(-1 != 0,1,1,0,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralSCSDOUBLE(1 != 0,1,1,0,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralSCSDOUBLE(-1 != 0,1,-1,0,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralSCSDOUBLE(-1 != 0,1,1,0,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralSCSDOUBLE(0.5 != 0,1,0.5,0,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralSCSDOUBLE(1.5 != 0,1,0.5,0,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSCSDOUBLE(1 != 0,1,1.5,0,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , PI_SQUARE*0.5*0.8) << std::endl; exit(0);}
  if (ComputeIntegralSCSDOUBLE(1.5 != 0,1,2,0,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_SQUARE*0.5*(1+1.0/7.0)) << std::endl; exit(0);}
  
  if (ComputeIntegralSSCDOUBLE(0 != 0,0,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSSCDOUBLE(1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSSCDOUBLE(-1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSSCDOUBLE(-1 != 1,1,1,1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSSCDOUBLE(-1 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralSSCDOUBLE(-1 != 1,-1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralSSCDOUBLE(1 != 1,1,1,1,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralSSCDOUBLE(1 != 1,1,1,-1,-1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.125) << std::endl; exit(0);}
  if (ComputeIntegralSSCDOUBLE(-1 != 1,0,1,1,0) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralSSCDOUBLE(1 != 1,0,1,-1,0) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralSSCDOUBLE(1 != -1,0,1,-1,0) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralSSCDOUBLE(-1 != 1,0,1,-1,0) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralSSCDOUBLE(0.5 != 1,0,0.5,1,0) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_CUBE*0.25) << std::endl; exit(0);}
  if (ComputeIntegralSSCDOUBLE(1.5 != 1,0,0.5,-1,0) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSSCDOUBLE(1 != -1,0,1.5,-1,0) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,PI_SQUARE*0.5*0.8) << std::endl; exit(0);}
  if (ComputeIntegralSSCDOUBLE(1.5 != 1,0,2,-1,0) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-PI_SQUARE*0.5*(1+1.0/7.0)) << std::endl; exit(0);}
  
  if (ComputeIntegralSSCCSS(0 != 0,0,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSSCCSS(1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSSCCSS(-1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSSCCSS(-1 != 1,1,1,1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  // if (ComputeIntegralSSCCSS(2 != 1,2,1,1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0.5*M_PI*4.0/3.0*1.20) << std::endl; exit(0);}
  // if (ComputeIntegralSSCCSS(2 != 1,2,1,-1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-0.5*M_PI*4.0/3.0*1.20) << std::endl; exit(0);}
  if (ComputeIntegralSSCCSS(5 != 1,2,2,1,2) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  // if (ComputeIntegralSSCCSS(5 != 1,2,2,1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0.5*M_PI*10.0/21.0*1.20) << std::endl; exit(0);}
  
  if (ComputeIntegralCCCSCS(0 != 0,0,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCCCSCS(1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCCCSCS(-1 != 1,1,1,2,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCCCSCS(-1 != 1,1,1,1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  //if (ComputeIntegralCCCSCS(2 != 1,2,1,1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-0.5*M_PI*2.0/3.0*1.2) << std::endl; exit(0);}
  // if (ComputeIntegralCCCSCS(2 != 1,2,1,-1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-0.5*M_PI*2.0/3.0*1.2) << std::endl; exit(0);}
  // if (ComputeIntegralCCCSCS(2 != 0,2,1,0,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-M_PI*2.0/3.0*1.2) << std::endl; exit(0);}
  if (ComputeIntegralCCCSCS(5 != 1,2,1,-1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  // if (ComputeIntegralCCCSCS(5 != 1,2,2,-1,3) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,-0.5*M_PI*4.0/21.0*1.2) << std::endl; exit(0);}
  
  if (ComputeIntegralSCCCSS(1 != 1,1,1,2,2) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSCCCSS(1 != 1,1,2,1,2) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSCCCSS(1 != 1,1,2,2,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralSCCCSS(1 != 1,1,2,2,2) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , - 32.0/27.0 ) << std::endl; exit(0);}
  if (ComputeIntegralSCCCSS(3 != 5,7,2,2,2) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 1.2*16.0/21.0/45.0) << std::endl; exit(0);}
  if (ComputeIntegralSCCCSS(1 != 2,1,2,5,2) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , -80.0/9.0/21.0 ) << std::endl; exit(0);}
  
  if (ComputeIntegralCSCSCS(1 != 1,1,1,2,2) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCSCSCS(1 != 1,1,2,1,2) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCSCSCS(1 != 1,1,2,2,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCSCSCS(1 != 1,1,2,2,2) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , -32.0/27.0 ) << std::endl; exit(0);}
  // if (ComputeIntegralCSCSCS(3 != 5,7,2,2,2) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 32/21.0/45.0) << std::endl; exit(0);}
  if (ComputeIntegralCSCSCS(1 != 2,1,2,5,2) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , -64.0/9.0/21.0 ) << std::endl; exit(0);}
  
  if (ComputeIntegralCCSSSC(1 != 1,1,1,2,2) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCCSSSC(1 != 1,1,2,1,2) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  if (ComputeIntegralCCSSSC(1 != 1,1,2,2,1) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " ,0) << std::endl; exit(0);}
  // if (ComputeIntegralCCSSSC(1 != 1,1,2,2,2) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , - 32.0/27.0 ) << std::endl; exit(0);}
  //  if (ComputeIntegralCCSSSC(3 != 5,7,2,2,2) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , 3.2/21*14/45) << std::endl; exit(0);}
  // if (ComputeIntegralCCSSSC(1 != 2,1,2,5,2) {std::cout << "trig_integral_3d_test.cpp " << __LINE__ << " FATAL: " , -4.0/3.0*10.0/21.0*2.0/3.0) << std::endl; exit(0);}
  
  return 0;
}
