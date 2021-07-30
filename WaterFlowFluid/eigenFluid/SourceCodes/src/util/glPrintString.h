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

#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

void glPrintString1(float x, float y, const char *string) {
  //set the position of the text in the window using the x and y coordinates
  glRasterPos2f(x,y);
  //get the length of the string to display
  int len = (int) strlen(string);

  //loop to display character by character
  for (int i = 0; i < len; i++) 
  {
  glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, string[i]);
  }
};

