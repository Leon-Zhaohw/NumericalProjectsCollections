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

#ifndef VFIELD2D_H
#define VFIELD2D_H

#include "Alg/VEC2.h"
#include "setting.h"
#include "2D/FIELD2D.h"
#include <cassert>

class VFIELD2D
{
public:
  VFIELD2D(){};

  VFIELD2D(int _xRes, int _yRes): xRes(_xRes), yRes(_yRes) {
    data = (VEC2F*) malloc(sizeof(VEC2F)*xRes*yRes);
    //memset(data, 0x00, xRes*yRes*sizeof(VEC2F));
    for(int i=0;i<xRes*yRes;i++)
    { data[i][0] = 0.; data[i][1] = 0.;}
    totalcell = xRes*yRes;
  }

  inline VEC2F& operator()(int x, int y) { 
    assert(  y >= 0 && y < yRes && x >= 0 && x < xRes);
    return data[ y * xRes + x]; 
  };
  inline VEC2F& operator [] (int x) {
    assert(x>=0 && x<totalcell);
    return data[x];
  }
  const VEC2F operator [] (int x) const {
    assert(x>=0 && x<totalcell);
    return data[x];
  }
  // Get the velocity at a specific position.
  VEC2 GetVelocity(const float pos_x, const float pos_y) const;
  void set_zero_border();

  void SetZeroLeft();
  void SetZeroRight();
  void SetZeroTop();
  void SetZeroBottom();
  int getyRes()const{return yRes;}
  int getxRes()const{return xRes;}
  
  void SetNeumannLeft();
  void SetNeumannRight();
  void SetNeumannTop();
  void SetNeumannBottom();
  void CopyBorderRight();
  void CopyBorderLeft();
  void CopyBorderTop();
  void CopyBorderBottom();
  void axpy(const Real alpha, const VFIELD2D&);
  void swapPointer( VFIELD2D& field);
  static void advect(const Real dt, const VFIELD2D& velocity_field, const FIELD2D& oldField, FIELD2D& newField);
  static void advect(const Real dt, const VFIELD2D& velocity_field, const VFIELD2D& oldField, VFIELD2D& newField);
   
  static void advectMacCormack(const Real dt, const VFIELD2D& velocityGrid, FIELD2D& oldField, 
                               FIELD2D& newField, FIELD2D& temp1, FIELD2D& temp2);
  static void clampExtrema(const Real dt, const VFIELD2D& velocityField, const FIELD2D& 
                                      oldField, FIELD2D& newField);
  static void clampOutsideRays(const Real dt, const VFIELD2D& velocityField, const FIELD2D&  
                                        oldField, const FIELD2D& oldAdvection, FIELD2D& newField);
  
  void clear();
  void SaveVelocity(FILE* out);
  void Draw(float magnitute);
  void Draw_X(float magnitute, int mode);
  double GetEnergy();
  double GetMaximumVelocity();
  // Do a dot product with another vector field in 2D.
  double DotProduct(const VFIELD2D& b);
private:
  int totalcell;
  VEC2F* data;
  int xRes;
  int yRes;
};

#endif
