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

#ifndef DRAWER_3D_H
#define DRAWER_3D_H

#include "Eigen"
#include <vector>

#include "3D/particle_3d.h"
#include "util/util.h"

class Drawer3D {
public:
  Drawer3D() {createDefaultHeatMapGradient();
    
  };
  ~Drawer3D() {};
  void DrawParticles(const std::vector<Particle3D>& particles, const double dx_,
                     const double dt_, const double ptl_length);
  // Visualize the coefficients as a bar graph.
  void DrawCoefficients(const Eigen::VectorXd& coefficients,
                        const double dx_, const double multi_factor);
  void DrawSpectogram(const Eigen::VectorXd& coefficients,
                      const double multi_factor);
protected:
  struct ColorPoint  // Internal class used to store colors at different points in the gradient.
  {
    float r,g,b;      // Red, green and blue values of our color.
    float val;        // Position of our color along the gradient (between 0 and 1).
    ColorPoint(float red, float green, float blue, float value)
      : r(red), g(green), b(blue), val(value) {}
  };
  std::vector<ColorPoint> color;      // An array of color points in ascending value.
  void createDefaultHeatMapGradient() {
    color.clear();
   // rainbow..
    // color.push_back(ColorPoint(0, 0, 0,   0.0f));      // Black.
    // color.push_back(ColorPoint(0, 0, 1,   0.0f));      // Blue.
    // color.push_back(ColorPoint(0, 1, 1,   0.25f));     // Cyan.
    // color.push_back(ColorPoint(0, 1, 0,   0.5f));      // Green.
    // color.push_back(ColorPoint(1, 1, 0,   0.75f));     // Yellow.
    // color.push_back(ColorPoint(1, 0, 0,   1.0f));      // Red.
    // blue to red
    color.push_back(ColorPoint(0., 0., 1.0, 0.0f)); // blue
    //color.push_back(ColorPoint(1.0, 1.0,1.0, 0.f));
    color.push_back(ColorPoint(0.48, 0.76,1.0, 0.25f)); // purple-close to blue
    color.push_back(ColorPoint(0.78, 0.78, 0.78, 0.5f)); // purple
    color.push_back(ColorPoint(0.85, 0.5, 0.5, 0.75f)); // purple - close to red
    color.push_back(ColorPoint(1.0, 0.0, 0., 1.0f)); // red
    //color.push_back(ColorPoint(0.0, 0.0, 0.0, 1.0f));
   }
  //-- Inputs a (value) between 0 and 1 and outputs the (red), (green) and (blue)
  //-- values representing that position in the gradient.
  void getColorAtValue(const float value, float &red, float &green, float &blue) {
    if(color.size()==0)
      return;
    
    float colorval = clamp(value, 0.f, 1.f  );
    
    for(int i=0; i<color.size(); i++)
    {
      ColorPoint &currC = color[i];
      if(colorval < currC.val)
      {
        ColorPoint &prevC  = color[ std::max(0,i-1) ];
        float valueDiff    = (prevC.val - currC.val);
        float fractBetween = (valueDiff==0) ? 0 : (colorval - currC.val) / valueDiff;
        red   = (prevC.r - currC.r)*fractBetween + currC.r;
        green = (prevC.g - currC.g)*fractBetween + currC.g;
        blue  = (prevC.b - currC.b)*fractBetween + currC.b;
        return;
      }
    }
    red   = color.back().r;
    green = color.back().g;
    blue  = color.back().b;
    return;
  }
  std::vector<Eigen::VectorXd> spectogramCoeff_;
  void subSampleVector(const Eigen::VectorXd& vecIn, Eigen::VectorXd& vecOut);
};

#endif  // DRAWER_3D_H
