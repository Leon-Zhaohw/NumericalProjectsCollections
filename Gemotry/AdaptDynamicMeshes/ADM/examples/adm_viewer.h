#ifndef _ADM_VIEWER_
#define _ADM_VIEWER_ 1

#include <adm.h>
#include "glViewer.h"

class ADMViewer : public glViewer
{
 public:
  typedef ADM::R3 R3;
  typedef ADM::AdaptiveMesh::VertexIter VertexIter;
  typedef ADM::AdaptiveMesh::EdgeIter   EdgeIter;
  typedef ADM::AdaptiveMesh::FacetIter  FacetIter;
  
 protected:
  ADM::AdaptiveMesh*       mesh_;
  ADM::Oracle*             surf_;
  ADM::StochasticSampling* criteria_;
  
  bool wprimal_, dmax_, dmin_;
  
 public:
  ADMViewer(const char* label, ADM::AdaptiveMesh* m, ADM::Oracle* s, ADM::StochasticSampling* c);
  virtual ~ADMViewer();
  
  virtual void DrawModel();
  void saveASoff();
  int  handle(int e);

  void KeyboardOptions()
    {
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
      std::cout << "\t q - max principal direction" << std::endl;
      std::cout << "\t w - min principal direction" << std::endl;
      std::cout << "\t e - primal mesh" << std::endl;
      std::cout << "\t F1 - Gouraud" << std::endl;
      std::cout << "\t F2 - max curvature [min - 0 - max]" << std::endl;
      std::cout << "\t F3 - min curvature [min - 0 - max]" << std::endl;
      std::cout << "\t F4 - mean curvature    [min - 0 - max]" << std::endl;
      std::cout << "\t F5 - gaussian curvature [min - 0 - max]" << std::endl;
      std::cout << "\t F6 - max Abs curvature [min - max]" << std::endl;
      std::cout << "\t F7 - min Abs curvature [min - max]" << std::endl;
      std::cout << "\t F8 - samples" << std::endl;
      std::cout << "\t F9 - error" << std::endl;
      std::cout << "\t F10 - checker" << std::endl;
      std::cout << "\t F11 - none" << std::endl;
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
      std::cout << "\t x - deform and adapt"          << std::endl;
      std::cout << "\t a - adapt"                     << std::endl;
      std::cout << "\t d - deform"                    << std::endl;
      std::cout << "\t t - strutural phase"           << std::endl;
      std::cout << "\t r - refine"                    << std::endl;
      std::cout << "\t s - simplify"                  << std::endl;
      std::cout << "\t g - smoothing phase"           << std::endl;
      std::cout << "\t n - refine_step"               << std::endl;
      std::cout << "\t m - simplify_step"             << std::endl;
      std::cout << "\t / - play refinement animation" << std::endl;
      glViewer::KeyboardOptions();
    }
  
  friend void animation_func(void* o);
  
 protected:
  virtual void CreateDisplayLists();
  virtual void SetPosition();

  void PrimalWireframe(); // 0
  void MaxDirCurv();      // 1
  void MinDirCurv();      // 2
  void Gouraud();         // 3
  void MaxCurv();         // 4
  void MinCurv();         // 5
  void MeanCurv();        // 6
  void GausCurv();        // 7
  void AbsMaxCurv();      // 8
  void AbsMinCurv();      // 9 
  void Samples();         // 10
  void Error();           // 11
  void Checker();         // 12

  void rec_samples(int level, R3 v0, R3 v1, R3 v2);
  double triangle_error(int level, R3 v0, R3 v1, R3 v2);

  inline double AbsMax(double a, double b)
    { return std::max( std::fabs(a), std::fabs(b) ); }
  inline double AbsMin(double a, double b)
    { return std::min( std::fabs(a), std::fabs(b) ); }
};

/*************************************************************************************/

ADMViewer::ADMViewer(const char* label, ADM::AdaptiveMesh* m, 
		     ADM::Oracle* s, ADM::StochasticSampling* c)
  : glViewer(512,512,label),
     mesh_(m), surf_(s), criteria_(c),
     wprimal_(false), dmax_(false), dmin_(false)
{ 
#ifdef CHECKER
  nmodes_ = 13;
#else
  nmodes_ = 12; 
#endif
  mode_ = 3; 
  Fl::add_timeout(0.5,animation_func,this); 
}

ADMViewer::~ADMViewer() { }

void ADMViewer::DrawModel()
{
  if ( mode_>2 && mode_<nmodes_) glCallList(displaylist_+mode_);
  if (wprimal_)  glCallList(displaylist_);
  if (dmax_)     glCallList(displaylist_+1);
  if (dmin_)     glCallList(displaylist_+2);
}

void ADMViewer::saveASoff()
{
  char filename[256];
  std::sprintf(filename,"%s-m%.1d-n%.3d.off",label(),mode_,counter_);

  mesh_->save_mesh(filename);

  /*
  std::ofstream output(filename);
  output << "OFF" << std::endl;
  output << mesh_->num_verts() << " "
         << mesh_->num_facets() << " "
         << "0" << std::endl;
  
  VertexIter vi;
  for (vi=mesh_->verts_begin(); vi!=mesh_->verts_end(); vi++)
    output << (*vi)->pos() << std::endl;

  EdgeIter ei;
  for (ei=mesh_->edges_begin(); ei!=mesh_->edges_end(); ei++)
    {
      if ( !ADM::edge_is_subdiv2( (*ei)->hedge(0) ) ) continue;
      output << "3 " 
	     << (*ei)->hedge(0)->org()->id() << " "
	     << (*ei)->hedge(0)->dst()->id() << " "
	     << (*ei)->hedge(0)->next()->dst()->id() << std::endl;
      output << "3 " 
	     << (*ei)->hedge(1)->org()->id() << " "
	     << (*ei)->hedge(1)->dst()->id() << " "
	     << (*ei)->hedge(1)->next()->dst()->id() << std::endl;
    }
  output.close();
  */

  std::cout << "Saving " << filename << std::endl;
}

int ADMViewer::handle(int e)
{
  switch (e)
    {
    case FL_KEYBOARD:
      switch (Fl::event_key())
        {
	  // VIEW MODES
        case FL_F+1: // gouraud
          mode_ = 3;
          redraw();
          break;
        case FL_F+2:
          mode_ = 4; // max curv
          redraw();
          break;
        case FL_F+3:
          mode_ = 5; // min curv
          redraw();
          break;
        case FL_F+4:
          mode_ = 6; // mean curv
          redraw();
          break;
        case FL_F+5:
          mode_ = 7; // gaus curv
          redraw();
          break;
        case FL_F+6:
          mode_ = 8; // abs max curv
          redraw();
          break;
        case FL_F+7:
          mode_ = 9; // abs min curv
          redraw();
          break;
        case FL_F+8:
          mode_ = 10; // samples
          redraw();
          break;
        case FL_F+9:
          mode_ = 11; // error
          redraw();
          break;
        case FL_F+10:
          mode_ = 12; // checker
          redraw();
          break;
        case FL_F+11:
          mode_ = nmodes_; // none
          redraw();
          break;
        case 'e':
          wprimal_ = !wprimal_;
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

	  // ADAPTATION
	case 'x':
          glDeleteLists(displaylist_,nmodes_);
          mesh_->deform();
          mesh_->adapt();
          CreateDisplayLists();
          redraw();
          break;
        case 'a':
          glDeleteLists(displaylist_,nmodes_);
          mesh_->adapt();
          CreateDisplayLists();
          redraw();
          break;
        case 'd':
          glDeleteLists(displaylist_,nmodes_);
          mesh_->deform();
          CreateDisplayLists();
          redraw();
          break;
        case 't':
          glDeleteLists(displaylist_,nmodes_);
          mesh_->strutural_phase();
          CreateDisplayLists();
          redraw();
          break;
        case 'r':
          glDeleteLists(displaylist_,nmodes_);
          mesh_->refine_phase();
          CreateDisplayLists();
          redraw();
          break;
        case 's':
          glDeleteLists(displaylist_,nmodes_);
          mesh_->simplify_phase();
          CreateDisplayLists();
          redraw();
          break;
        case 'g':
          glDeleteLists(displaylist_,nmodes_);
          mesh_->smoothing_phase();
          CreateDisplayLists();
          redraw();
          break;
        case 'n':
          glDeleteLists(displaylist_,nmodes_);
          mesh_->refine_one_step();
          CreateDisplayLists();
          redraw();
          break;
	case 'm':
          glDeleteLists(displaylist_,nmodes_);
          mesh_->simplify_one_step();
          CreateDisplayLists();
          redraw();
          break;

          //ANIMATION BUTTONS
        case '/':
          animating_ = 2;
          break;

	default: return glViewer::handle(e);
        };
    default: return glViewer::handle(e);
    };
  return 1;
}

void ADMViewer::CreateDisplayLists()
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
  MaxDirCurv();
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+2, GL_COMPILE);
  glPushMatrix();
  Centering();
  MinDirCurv();
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+3, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  Gouraud();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+4, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  MaxCurv();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+5, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  MinCurv();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+6, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  MeanCurv();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+7, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  GausCurv();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+8, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  AbsMaxCurv();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+9, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  AbsMinCurv();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+10, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  Samples();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();

  glNewList(displaylist_+11, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  Error();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();

#ifdef CHECKER
  glNewList(displaylist_+12, GL_COMPILE);
  glPushMatrix();
  Centering();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  Checker();
  glDisable(GL_POLYGON_OFFSET_FILL);
  glPopMatrix();
  glEndList();
#endif
}

void ADMViewer::SetPosition()
{
  double xmax, ymax, zmax, xmin, ymin, zmin;
  mesh_->bounding_box(xmax,ymax,zmax,xmin,ymin,zmin);
  cx_ = float(0.5*(xmax+xmin));
  cy_ = float(0.5*(ymax+ymin));
  cz_ = float(0.5*(zmax+zmin));
  scale_ = std::max( (xmax-xmin), std::max( (ymax-ymin), (zmax-zmin) ) );
  scale_ = 2. / scale_;
}

/*************************************************************************************/

void animation_func(void* o)
{
  ADMViewer* v = (ADMViewer*)o;
  switch (v->animating_)
    {
    case 1:
      glDeleteLists( v->displaylist_, v->nmodes_ );
      v->mesh_->deform();
      v->mesh_->adapt();
      v->CreateDisplayLists();
      v->redraw();
      break;
    case 2:
      glDeleteLists( v->displaylist_, v->nmodes_ );
      v->mesh_->refine_one_step();
      v->CreateDisplayLists();
      v->redraw();
      break;
    default: break;
    };
  Fl::repeat_timeout(0.5,animation_func,v);
}

#endif
