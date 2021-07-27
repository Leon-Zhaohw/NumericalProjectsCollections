/* (C)Copyright 2007                                               */
/* International Business Machines Corporation,                    */
/* All Rights Reserved.                                            */
/*                                                                 */
/* This program is made available under the terms of the           */
/* Common Public License v1.0 which accompanies this distribution. */
/* --------------------------------------------------------------- */
/* PROLOG END TAG zYx                                              */

////////////////////////////////////////////////////////////////////////////
// This project is a 2D implementation of the boiling module from:
//
//   T. Kim and M. Carlson, A Simple Boiling Module. Proceedings of
//    ACM SIGGRAPH / Eurographics Symposium on Computer Animation 2007
//
// This file contains mostly memory allocation and GL calls. The actual
// implementations are in:
// 
//    ExtendedYanagita.h (our algorithm)
//    OriginalYanagita.h (the original algorithm)
//
// Please direct any questions or comments to Ted Kim: twkim@us.ibm.com
////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <GL/glut.h>
#include <cstdio>

using namespace std;

////////////////////////////////////////////////////////////////////////////
// globals
////////////////////////////////////////////////////////////////////////////
bool animate = false;
bool useOriginal = true;
bool drawNormals = false;

int factor = 2;
int res = 64 * factor;

float* heat = NULL;
float* heat1 = NULL;
float* heat2 = NULL;
float* phase = NULL;
bool* crossings = NULL;

float* normalX = NULL;
float* normalY = NULL;

// timestep scaling factor
float dtFactor = 1.0f / ((float)factor);

// simulation constants
float eps       = 0.5f;
float alpha     = 10.0f;
float sigma     = 0.05f;
float eta       = 0.1f;
float Tc        = 10.0f;
float bottom    = 9.90f;
float top       = 9.2f;

// simulation constants, scaled by the current resolution
float diffConst    = eps * 0.25f * (float)(factor * factor) * dtFactor;
float buoyConst    = sigma * (float)factor * dtFactor;
float latentConst  = 1.0f / (float)(factor * factor) * dtFactor;
float tensionConst = 0.025f * dtFactor * (float)(factor) * 2.0f;

float zoom = 0.0f;

//////////////////////////////////////////////////////////////////////
// Implementation of our algorithm.
//////////////////////////////////////////////////////////////////////
#include "ExtendedYanagita.h"

//////////////////////////////////////////////////////////////////////
// Implementation of Yanagita' original algorithm
//////////////////////////////////////////////////////////////////////
#include "OriginalYanagita.h"

//////////////////////////////////////////////////////////////////////
// draw field to screen
//////////////////////////////////////////////////////////////////////
void drawHeat()
{
  glPushMatrix();
  
  // draw quads
  glBegin(GL_QUADS);
    int i = 0;
    for (int y = 0; y < res; y++)
      for (int x = 0; x < res; x++, i++)
      {
        float color = 5.0f * (heat[i] / bottom - 0.8f);
        glColor4f((10.0f - heat[i]) * 0.1f, (phase[i] + 1.0f) * 0.5f, color, 1.0f);

        glVertex2f(x, y);
        glVertex2f(x + 1, y);
        glVertex2f(x + 1, y + 1);
        glVertex2f(x, y + 1);
      }
  glEnd();

  glPopMatrix();
}

////////////////////////////////////////////////////////////////////////////
// window Reshape function 
////////////////////////////////////////////////////////////////////////////
void Reshape(int w, int h)
{
  if (h == 0) h = 1;
  
  glViewport(0, 0, w, h);
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(1.0f + zoom, res - 1 - zoom, 1.0f + zoom, res - 1 - zoom);
  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

////////////////////////////////////////////////////////////////////////////
// GLUT Display callback
////////////////////////////////////////////////////////////////////////////
void Display()
{
  // setup for draw
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(1.0f + zoom, res - 1 - zoom, 1.0f + zoom, res - 1 - zoom);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  
  drawHeat();

  glutSwapBuffers();
}

////////////////////////////////////////////////////////////////////////////
// GLUT Keyboard callback
////////////////////////////////////////////////////////////////////////////
void Keyboard(unsigned char key, int x, int y)
{
  switch(key) 
  {
    // quit
    case 'q': exit(0); break;
      
    // toggle simulation animation
    case ' ': animate = !animate; break;

    // toggle between original and extended Yanagita
    case 'o':
      useOriginal = !useOriginal;
      if (useOriginal)
        cout << " Using original Yanagita algorithm" << endl;
      else
        cout << " Using our extended Yanagita algorithm" << endl;
      break;
  }
  
  glutPostRedisplay();
}

////////////////////////////////////////////////////////////////////////////
// GLUT Idle callback
////////////////////////////////////////////////////////////////////////////
void Idle()
{
  static int step = 0;

  if (animate)
  {
    // set bottom to hot
    for (int x = 0; x < factor * res; x++)
    {
      float set = top /2 + step * 0.05;
      if (set > bottom) set = bottom;
      heat[x] = heat1[x] = heat2[x] = set;
    }

    // set sides to periodic
    for (int y = 1; y < res - 1; y++) {
      int i = y * res;
      heat[i]  = heat[i + res - 2];  heat[i + res - 1]  = heat[i + 1];
      heat1[i] = heat1[i + res - 2]; heat1[i + res - 1] = heat1[i + 1];
      heat2[i] = heat2[i + res - 2]; heat2[i + res - 1] = heat2[i + 1];
    }

    if (useOriginal)
      stepOriginalYanagita();
    else
      stepExtendedYanagita();

    step++;
  }
    
  glutPostRedisplay();
}

////////////////////////////////////////////////////////////////////////////
// GLUT Main 
////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{  
  cout << " Boiling Code: " << endl << endl;
  cout << "   Toggle animation: space" << endl;
  cout << "   Toggle original and extended Yanagita: o" << endl;
  cout << "   Quit:              q" << endl;
  
  srand(123456);

  heat = new float[res * res];
  heat1 = new float[res * res];
  heat2 = new float[res * res];
  phase = new float[res * res];
  crossings = new bool[res * res];
  normalX = new float[res * res];
  normalY = new float[res * res];

  // throw in some noise
  for (int x = 0; x < res * res; x++)
  {
    float jitter = bottom * 0.1f;
    heat[x] = (jitter * rand() / RAND_MAX - jitter * 0.5f) + 0.89f *bottom;
      
    heat1[x] = heat2[x] = heat[x];
    crossings[x] = false;
  }

  glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
  glutInitWindowPosition(50, 50);
  glutInitWindowSize(600, 600);
  glutInit(&argc, argv);
  glutCreateWindow("Boiling Code");

  glutDisplayFunc(Display);
  glutKeyboardFunc(Keyboard);
  glutIdleFunc(Idle);  
  Reshape(600, 600); 
  glLineWidth(2.0f);
  glPointSize(4.0f);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE);

  // Go!
  glutMainLoop();

  delete[] heat;
  delete[] heat1;
  delete[] heat2;
  delete[] phase;
  delete[] crossings;
  delete[] normalX;
  delete[] normalY;

  return 0;
}
