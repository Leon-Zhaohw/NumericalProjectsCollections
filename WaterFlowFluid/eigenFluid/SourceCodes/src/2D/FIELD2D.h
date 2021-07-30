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

#ifndef FIELD2D_H
#define FIELD2D_H

#include <cassert>
#include <stdlib.h>
#include <memory>
#include <stdio.h>
#include "setting.h"

class FIELD2D
{
public:
	FIELD2D(){};

	//x minor index, y multiplied by xres
	inline Real& operator()(int x, int y) { 
		assert(  y >= 0 && y < yRes && x >= 0 && x < xRes);
		return data[ y * xRes + x]; 
	};
	inline Real& operator[](int x) { return data[x]; };
	const Real operator[](int x) const { return data[x]; };
	FIELD2D(int _xRes, int _yRes):xRes(_xRes),yRes(_yRes)
	{
		data = (Real*) malloc(sizeof(Real) * xRes* yRes);
		//memset(data, 0x00, sizeof(Real) * xRes * yRes);
		for(int i=0;i<xRes*yRes;i++)
			data[i] = 0.;
		totalsize = xRes*yRes;
	}
	void CheckNegative();
	FIELD2D& operator *=(const Real alpha);
	FIELD2D& operator +=(const Real alpha);
// 	FIELD2D& operator=(const Real& alpha);
// 	FIELD2D& operator=(const FIELD2D& A);
	FIELD2D& operator +=(const FIELD2D& input);
	FIELD2D& operator *=(const FIELD2D& input);
	void clear();
	void setZeroBorder();
    void setZeroBorder2Layers();
	Real* getdata(){return data;}
	void swapPointer(FIELD2D& field);
	int getxRes() const {return xRes;}
	int getyRes() const {return yRes;}

	void copyBorderAll();
	Real dot(const FIELD2D& in);
	void axpy(const Real& alpha, const FIELD2D& input);
	Real& get_dens(int i) {return data[i];}
	void Field2D_Save(FILE* out);
private:
	int totalsize;
	int xRes;
	int yRes;
Real* data;
};



#endif
