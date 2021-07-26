#include "glViewer.h"
#include <tml.h>
#include <cassert>
#include <iostream>
#include <stack>

int compute_num_components(TML::TriMesh<>& mesh)
{
  int num = -1;
  TML::Vertex* w = NULL;
  TML::Halfedge* h;
  TML::TriMesh<>::VertexIter vi;
  std::stack<TML::Vertex*> Q;
  for (vi=mesh.verts_begin(); vi!=mesh.verts_end(); vi++)
    {
      if ( (*vi)->mark()!=-1 ) continue;
      num++;
      Q.push( *vi );
      while (!Q.empty())
	{
	  w = Q.top(); Q.pop();
	  w->set_mark(num);
	  for (h=w->star_first(); h!=NULL; h=w->star_next(h))
	    {
	      if (h->org()->mark()!=-1) continue;
	      Q.push( h->org() );
	    }
	}
    }
  return (num+1);
}

//----------------------------------------------------------//

class Viewer : public glViewer
{
protected:
  TML::TriMesh<> mesh_;
  bool wireframe_active_;

public:
  Viewer(const char* filename)
    : glViewer(720,480,"Viewer"),
      mesh_(filename), wireframe_active_(false)
  {
    nmodes_ = 2;
    mode_ = 1;

    // size
    std::cout << "(V,E,F) = ("
	      << mesh_.num_verts()  << ", "
	      << mesh_.num_edges()  << ", "
	      << mesh_.num_facets() << ")" << std::endl;
    
    // checking boundary
    bool has_bdry = false;
    TML::TriMesh<>::VertexIter vi;
    for (vi=mesh_.verts_begin(); vi!=mesh_.verts_end(); vi++)
      if ( (*vi)->is_bdry() ) { has_bdry = true; break; }
    if (has_bdry) std::cout << "YES, I do have boundary(ies)." << std::endl;
    else          std::cout << "NO, I do not have boundary." << std::endl;

    // couting components
    std::cout << "I have " << compute_num_components(mesh_) << " component(s)." << std::endl;
  }

  ~Viewer() { }
  
  void DrawModel()
  {
    if (wireframe_active_) glCallList(displaylist_);
    glCallList(displaylist_+1);
  }

  int handle(int e)
  {
    switch (e)
      {
      case FL_KEYBOARD:
	switch (Fl::event_key())
	  {
	  case 'w':
	    wireframe_active_ = !wireframe_active_;
	    redraw();
	    break;
	  default: return glViewer::handle(e);
          };
      default: return glViewer::handle(e);
      };
    return 1;
  }
  
  void saveASoff() { }

  void KeyboardOptions()
  {
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    std::cout << "\t w - (des)active wireframe" << std::endl;
    glViewer::KeyboardOptions();
  }

  void CreateDisplayLists()
  {
    displaylist_ = glGenLists(nmodes_);

    glNewList(displaylist_, GL_COMPILE);
    glPushMatrix();
    Centering();
    Wireframe();
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
  }

  void SetPosition()
  {
    double xmax, ymax, zmax, xmin, ymin, zmin;
    mesh_.bounding_box(xmax,ymax,zmax,xmin,ymin,zmin);
    cx_ = float(0.5*(xmax+xmin));
    cy_ = float(0.5*(ymax+ymin));
    cz_ = float(0.5*(zmax+zmin));
    scale_ = std::max( (xmax-xmin), std::max( (ymax-ymin), (zmax-zmin) ) );
    scale_ = 2. / scale_;
  }

  void Wireframe()
  {
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    
    TML::R3 p;
    TML::TriMesh<>::FacetIter fi;
    
    glColor3f(0.,0.,0.);
    glLineWidth(1);
    for (fi=mesh_.facets_begin(); fi!=mesh_.facets_end(); fi++)
      {
	glBegin(GL_TRIANGLES);
	for (int i=0; i<3; i++)
	  {
	    p = (*fi)->vertex(i)->pos();
	    glVertex3d( p.x, p.y, p.z );
	  }
	glEnd();
      }
  }
  
  void Gouraud()
  {
    glEnable(GL_CULL_FACE);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_RESCALE_NORMAL);
    glEnable(GL_COLOR_MATERIAL);

    glPolygonMode(GL_BACK,GL_LINE);
    glPolygonMode(GL_FRONT,GL_FILL);
    glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
    
    glLineWidth(1);

    TML::R3 p, n;
    TML::TriMesh<>::FacetIter fi;

    glColor3d(0.8,0.8,0.8);
    for (fi=mesh_.facets_begin(); fi!=mesh_.facets_end(); fi++)
      {
        glBegin(GL_TRIANGLES);
        for (int i=0; i<3; i++)
          {
            p = (*fi)->vertex(i)->pos();
            n = (*fi)->vertex(i)->normal();
            glNormal3d( n.x, n.y, n.z );
            glVertex3d( p.x, p.y, p.z );
          }
        glEnd();
      }
    
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_RESCALE_NORMAL);
    glDisable(GL_LIGHT0);
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
  }
};

//----------------------------------------------------------//

int main(int argc, char** argv)
{
  if (argc!=2)
    throw TML::Error("USAGE: ./viewer <smf | off | ply | ifs>");

  Viewer win(argv[1]);
  //win.KeyboardOptions();
  win.show();
  Fl::run();
  
  return 1;
}
