#include "TMViewer.h"

TMViewer::TMViewer(TML::TriMesh<Point>* mesh)
  : glViewer(512,512,"TMViewer"),
    mesh_(mesh), 
    wprimal_(false), wdual_(false), dmax_(false), dmin_(false)
{ 
  nmodes_ = 11; mode_ = 4; 
}

TMViewer::~TMViewer() { }

void TMViewer::DrawModel()
{
  if (mode_>3 && mode_<nmodes_) glCallList(displaylist_+mode_);
  if (wprimal_)  glCallList(displaylist_);
  if (wdual_)    glCallList(displaylist_+1);
  if (dmax_)     glCallList(displaylist_+2);
  if (dmin_)     glCallList(displaylist_+3);
}

void TMViewer::saveASoff()
{
  char filename[256];
  std::sprintf(filename,"%s-m%.1d-n%.3d.off",label(),mode_,counter_);
  
  mesh_->save_mesh(filename);

  std::cout << "Saving " << filename << std::endl;
}

int TMViewer::handle(int e)
{
  switch (e)
    {
    case FL_KEYBOARD:
      switch (Fl::event_key())
	{
	case FL_F+1:
	  mode_ = 4;
	  redraw();
	  break;
	case FL_F+2:
	  mode_ = 5;
	  redraw();
	  break;
	case FL_F+3:
	  mode_ = 6;
	  redraw();
	  break;
	case FL_F+4:
	  mode_ = 7;
	  redraw();
	  break;
	case FL_F+5:
	  mode_ = 8;
	  redraw();
	  break;
	case FL_F+6:
	  mode_ = 9;
	  redraw();
	  break;
	case FL_F+7:
	  mode_ = 10;
	  redraw();
	  break;
	case FL_F+8:
	  mode_ = nmodes_;
	  redraw();
	  break;
	case 'e':
	  wprimal_ = !wprimal_;
	  redraw();
	  break;
	case 'r':
	  wdual_ = !wdual_;
	  redraw();
	  break;
	case 'q':
	  dmax_ = !dmax_;
	  redraw();
	  break;
	case 'w':
	  dmin_ = !dmin_;
	  redraw();
	  break;

	case 'v':
	  std::cout << "Center " << cx_ << " " << cy_ << " " << cz_ << std::endl;
	  std::cout << "Scale " << scale_ << std::endl;
	  std::cout << "Zoom " << zoom_ << std::endl;
	  std::cout << rot_[0] << " " << rot_[4] << " " << rot_[8] << " " << rot_[12] << std::endl;
	  std::cout << rot_[1] << " " << rot_[5] << " " << rot_[9] << " " << rot_[13] << std::endl;
	  std::cout << rot_[2] << " " << rot_[6] << " " << rot_[10] << " " << rot_[14] << std::endl;
	  std::cout << rot_[3] << " " << rot_[7] << " " << rot_[11] << " " << rot_[15] << std::endl;
	  break;

	default: return glViewer::handle(e);
	};
    default: return glViewer::handle(e);
    };
  return 1;
}

void TMViewer::CreateDisplayLists()
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
  DualWireframe();
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+2, GL_COMPILE);
  glPushMatrix();
  Centering();
  MaxDirCurv();
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+3, GL_COMPILE);
  glPushMatrix();
  Centering();
  MinDirCurv();
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+4, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  Gouraud();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+5, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  MaxCurv();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+6, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  MinCurv();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+7, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  MeanCurv();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+8, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  GausCurv();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+9, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  AbsMaxCurv();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+10, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  AbsMinCurv();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();
}

void TMViewer::SetPosition()
{
  double xmax, ymax, zmax, xmin, ymin, zmin;
  mesh_->bounding_box(xmax,ymax,zmax,xmin,ymin,zmin);
  cx_ = float(0.5*(xmax+xmin));
  cy_ = float(0.5*(ymax+ymin));
  cz_ = float(0.5*(zmax+zmin));
  scale_ = std::max( (xmax-xmin), std::max( (ymax-ymin), (zmax-zmin) ) );
  scale_ = 2. / scale_;
}
