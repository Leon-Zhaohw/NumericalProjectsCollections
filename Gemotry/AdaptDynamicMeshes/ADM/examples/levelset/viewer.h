#ifndef _LEVEL_SET_VIEWER_
#define _LEVEL_SET_VIEWER 1

#include "drawing.h"

class Viewer : public ADMViewer
{
 public:
  Viewer(const char* label, ADM::AdaptiveMesh* m, ADM::Oracle* s, ADM::StochasticSampling* c) :
    ADMViewer(label,m,s,c)
    { nmodes_ = 3; mode_ = 3; }

  ~Viewer() { }

  void DrawModel()
    {
      if ( wprimal_  ) glCallList(displaylist_);
      if ( mode_==3  ) glCallList(displaylist_+1);
      if ( mode_==11 ) glCallList(displaylist_+2);
    }

  void CreateDisplayLists()
    {
      displaylist_ = glGenLists(nmodes_);

      glNewList(displaylist_, GL_COMPILE);
      glPushMatrix();
      Centering();
      PrimalWireframe();
      glPopMatrix();
      glEndList();

      glNewList(displaylist_+1, GL_COMPILE);
      glPushMatrix();
      Centering();
      glEnable(GL_POLYGON_OFFSET_FILL);
      glPolygonOffset(1,1);
      Gouraud();
      glDisable(GL_POLYGON_OFFSET_FILL);
      glPopMatrix();
      glEndList();

      glNewList(displaylist_+2, GL_COMPILE);
      glPushMatrix();
      Centering();
      glEnable(GL_POLYGON_OFFSET_FILL);
      glPolygonOffset(1,1);
      Error();
      glDisable(GL_POLYGON_OFFSET_FILL);
      glPopMatrix();
      glEndList();
    }
  
};

#endif
