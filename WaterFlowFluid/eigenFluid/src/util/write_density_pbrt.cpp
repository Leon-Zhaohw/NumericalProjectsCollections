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

#include <fstream>
#include <sstream>
#include <stdlib.h> 

#include "util/write_density_pbrt.h"

void WriteDensityPBRT(const std::string& fname, const FIELD_3D& field, const float atten_fact) {
  std::ofstream out(fname);
  if (!out.is_open()) {
    std::cout << "write_density_pbrt.cpp " << __LINE__ << " FATAL: " <<  "Cannot write to file: " << fname << std::endl; exit(0);
  }
  const int xRes = field.xRes();
  const int yRes = field.yRes();
  const int zRes = field.zRes();
  const int maxRes = std::max(std::max(xRes, yRes), zRes);
  const uint total = xRes*yRes*zRes;
  
  // Header.
  //out << "Volume \"volumegrid\" \n"; // PBRT v2.
  out << "MakeNamedMedium \"smoke\" \n "; //PBRT v3
  out << "\"integer nx\" [ " << xRes << " ]\n";
  out << "\"integer ny\" [ " << yRes << " ]\n";
  out << "\"integer nz\" [ " << zRes << " ]\n";
  // Write the extent, the maximum extent is 1.
  out << "\"point p0\" [ 0.0 0.0 0.0 ] ";
  out << "\"point p1\" [ " << static_cast<float>(xRes) / maxRes << " " << static_cast<float>(yRes) / maxRes 
      << " " << static_cast<float>(zRes) / maxRes << " ] \n";
  out << "\"float density\" [\n";
  // Dump all the density values.
  for (int i = 0; i < total; i++) {
    if (std::abs(field[i]*atten_fact) < 1e-10) {
      out << "0 ";
    } else {
      out << field[i]*atten_fact << " ";
    }
    if (i % 50 == 0) {
      out << "\n";
    }
  }
  out << " ]\n";
  out << "\"string type\" [ \"heterogeneous\" ] ";  //PBRT v3
  
  out.close();
  
  std::stringstream zip_cmd;
  zip_cmd << "gzip " << fname << " &";
  // zip the frame.
  system(zip_cmd.str().c_str());
}; 
