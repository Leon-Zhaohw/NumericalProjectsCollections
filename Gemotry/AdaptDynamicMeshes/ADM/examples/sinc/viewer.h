#ifndef _VIEWER_SINC_
#define _VIEWER_SINC_ 1

#include "drawing.h"

class Viewer : public ADMViewer
{
 public:
  Viewer(const char* label, ADM::AdaptiveMesh* m, ADM::Oracle* s, ADM::StochasticSampling* c) :
    ADMViewer(label,m,s,c)
    { }

  ~Viewer() { }

  void draw()
    {
      if (!valid())
	{ glViewport(0,0,w(),h()); }
      
      if (first_)
	{
	  SetPosition();
	  InitializeGL();
	  CreateDisplayLists();
	  first_ = false;
	}

      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      glPushMatrix(); //1
      
      glPushMatrix(); //2

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      GLfloat aspect = (GLfloat)w()/(GLfloat)h();
      gluPerspective(45.0,aspect,1.0,100.0);

      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      glTranslatef(0,0.8,0.0);
      RigidTransform();
      DrawModel();
      UpdateRot();

      glPopMatrix(); //2
      
      glPushMatrix(); //3

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glOrtho(-3,3,-3,3,-3,3);
      glMatrixMode(GL_MODELVIEW);
      
      glPushMatrix(); //4
      glTranslatef(-1.5,-1.5,0);
      glScalef(0.8,0.8*aspect,0.8);
      glRotatef(90,1,0,0);
      glCallList(displaylist_+8);
      glPopMatrix(); //4
      
      glPushMatrix(); //5
      glTranslatef(1.5,-1.5,0);
      glScalef(0.8,0.8*aspect,0.8);
      glRotatef(90,1,0,0);
      glCallList(displaylist_+9);
      glPopMatrix(); //5
      
      glPopMatrix(); //3
      
      glPopMatrix(); //1
      glFlush();

      if (animating_>0 && saving_) saveASpng();
    }
};

#endif
