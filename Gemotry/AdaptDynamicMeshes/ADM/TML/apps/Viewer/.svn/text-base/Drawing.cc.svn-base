#include "TMViewer.h"

void TMViewer::PrimalWireframe()
{ 
  glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
  
  R3 p;
  FacetIter fi;
  
  glColor3f(0.,0.,0.);
  glLineWidth(1);
  for (fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
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

void TMViewer::DualWireframe()
{ 
  glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

  R3 p;
  VertexIter vi;
  TML::Halfedge* hv;

  glColor3f(0.,0.7,0.);
  glLineWidth(1);
  for (vi=mesh_->verts_begin(); vi!=mesh_->verts_end(); vi++)
    {
      glBegin(GL_POLYGON);
      for (hv=(*vi)->star_first(); hv!=NULL; hv=(*vi)->star_next(hv))
	{
	  if (hv->is_bdry()) continue;
	  p = hv->facet()->barycenter();
	  glVertex3d( p.x, p.y, p.z );
	}
      glEnd();
    }
}

void TMViewer::Gouraud()
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

  R3 p, n;
  FacetIter fi;
  
  glColor3d(0.8,0.8,0.8);
  for (fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
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

void TMViewer::MaxCurv()
{ 
  glEnable(GL_CULL_FACE);

  glPolygonMode(GL_BACK,GL_LINE);
  glPolygonMode(GL_FRONT,GL_FILL);

  glLineWidth(1);

  ///////////////////////////////////////////////////
  // Computing MaxCurv //////////////////////////////
  ///////////////////////////////////////////////////
  VertexIter vi = mesh_->verts_begin();
  double max = (*vi)->max_curvature();
  double min = 0.;

  vi++;
  double val;
  for ( ; vi!=mesh_->verts_end(); vi++)
    {
      val = (*vi)->max_curvature();
      max = std::max( max, val );
      min = std::min( min, val );
    }
  //////////////////////////////////////////////////

  int i;
  R3 p;
  Point* vp;
  FacetIter fi;
  for (fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
    {
      glBegin(GL_TRIANGLES);
      for (i=0; i<3; i++)
	{
	  vp = (Point*)( (*fi)->vertex(i) );
	  Color( vp->max_curvature(), min, 0., max );
	  p = vp->pos();
	  glVertex3d( p.x, p.y, p.z );
	}
      glEnd();
    }

  glDisable(GL_CULL_FACE);
}

void TMViewer::MinCurv()
{
  glEnable(GL_CULL_FACE);
  
  glPolygonMode(GL_BACK,GL_LINE);
  glPolygonMode(GL_FRONT,GL_FILL);
  
  glLineWidth(1);

  ///////////////////////////////////////////////////
  // Computing MinCurv //////////////////////////////
  ///////////////////////////////////////////////////
  VertexIter vi = mesh_->verts_begin();
  double max = (*vi)->min_curvature();
  double min = 0.;

  vi++;
  double val;
  for ( ; vi!=mesh_->verts_end(); vi++)
    {
      val = (*vi)->min_curvature();
      max = std::max( max, val );
      min = std::min( min, val );
    }
  //////////////////////////////////////////////////

  int i;
  R3 p;
  Point* vp;
  FacetIter fi;
  for (fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
    {
      glBegin(GL_TRIANGLES);
      for (i=0; i<3; i++)
	{
	  vp = (Point*)( (*fi)->vertex(i) );
	  Color( vp->min_curvature(), min, 0., max );
	  p = vp->pos();
	  glVertex3d( p.x, p.y, p.z );
	}
      glEnd();
    }
  
  glDisable(GL_CULL_FACE);
}

void TMViewer::MeanCurv()
{
  glEnable(GL_CULL_FACE);

  glPolygonMode(GL_BACK,GL_LINE);
  glPolygonMode(GL_FRONT,GL_FILL);

  glLineWidth(1);

  ///////////////////////////////////////////////////
  // Computing MeanCurv /////////////////////////////
  ///////////////////////////////////////////////////
  VertexIter vi = mesh_->verts_begin();
  double max = (*vi)->mean_curvature();
  double min = 0.;

  vi++;
  double kmean;
  for ( ; vi!=mesh_->verts_end(); vi++)
    {
      kmean = (*vi)->mean_curvature();
      max = std::max( max, kmean );
      min = std::min( min, kmean );
    }
  //////////////////////////////////////////////////

  int i;
  R3 p;
  FacetIter fi; 
  for (fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
    {
      glBegin(GL_TRIANGLES);
      for (i=0; i<3; i++)
	{
	  Color( (*fi)->vertex(i)->mean_curvature(), min, 0., max );
	  p = (*fi)->vertex(i)->pos();
	  glVertex3d( p.x, p.y, p.z );
	}
      glEnd();
    }
  
  glDisable(GL_CULL_FACE);
}

void TMViewer::GausCurv()
{ 
  glEnable(GL_CULL_FACE);
  
  glPolygonMode(GL_BACK,GL_LINE);
  glPolygonMode(GL_FRONT,GL_FILL);
  
  glLineWidth(1);

  ///////////////////////////////////////////////////
  // Computing GausCurv /////////////////////////////
  ///////////////////////////////////////////////////
  VertexIter vi = mesh_->verts_begin();
  double max = (*vi)->gaus_curvature();
  double min = 0.;

  vi++;
  double kgaus;
  for ( ; vi!=mesh_->verts_end(); vi++)
    {
      kgaus = (*vi)->gaus_curvature();
      max = std::max( max, kgaus );
      min = std::min( min, kgaus );
    }
  //////////////////////////////////////////////////

  int i;
  R3 p;
  FacetIter fi;
  for (fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
    {
      glBegin(GL_TRIANGLES);
      for (i=0; i<3; i++)
	{
	  Color( (*fi)->vertex(i)->gaus_curvature(), min, 0., max );
	  p = (*fi)->vertex(i)->pos();
	  glVertex3d( p.x, p.y, p.z );
	}
      glEnd();
    }

  glDisable(GL_CULL_FACE);
}

void TMViewer::AbsMaxCurv()
{
  glEnable(GL_CULL_FACE);
  
  glPolygonMode(GL_BACK,GL_LINE);
  glPolygonMode(GL_FRONT,GL_FILL);
  
  glLineWidth(1);

  ///////////////////////////////////////////////////
  // Computing AbsMaxCurv ///////////////////////////
  ///////////////////////////////////////////////////
  VertexIter vi = mesh_->verts_begin();
  double max = AbsMax( (*vi)->max_curvature(), (*vi)->min_curvature() );
  double min = max;

  vi++;
  double val;
  for ( ; vi!=mesh_->verts_end(); vi++)
    {
      val = AbsMax( (*vi)->max_curvature(), (*vi)->min_curvature() );
      max = std::max( max, val );
      min = std::min( min, val );
    }
  
  double zero = 0.5*(max+min);
  //////////////////////////////////////////////////

  int i;
  R3 p;
  Point* vp;
  FacetIter fi;
  for (fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
    {
      glBegin(GL_TRIANGLES);
      for (i=0; i<3; i++)
	{
	  vp = (Point*)( (*fi)->vertex(i) );
	  Color( AbsMax( vp->max_curvature(), vp->min_curvature() ), min, zero, max );
	  p = vp->pos();
	  glVertex3d( p.x, p.y, p.z );
	}
      glEnd();
    }
  
  glDisable(GL_CULL_FACE);
}

void TMViewer::AbsMinCurv()
{ 
  glEnable(GL_CULL_FACE);

  glPolygonMode(GL_BACK,GL_LINE);
  glPolygonMode(GL_FRONT,GL_FILL);

  glLineWidth(1);

  ///////////////////////////////////////////////////
  // Computing AbsMinCurv ///////////////////////////
  ///////////////////////////////////////////////////
  VertexIter vi = mesh_->verts_begin();
  double max = AbsMin( (*vi)->max_curvature(), (*vi)->min_curvature() );
  double min = max;

  vi++;
  double val;
  for ( ; vi!=mesh_->verts_end(); vi++)
    {
      val = AbsMin( (*vi)->max_curvature(), (*vi)->min_curvature() );
      max = std::max( max, val );
      min = std::min( min, val );
    }

  double zero = 0.5*(max+min);
  //////////////////////////////////////////////////

  int i;
  R3 p;
  Point* vp;
  FacetIter fi;
  for (fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
    {
      glBegin(GL_TRIANGLES);
      for (i=0; i<3; i++)
	{
	  vp = (Point*)( (*fi)->vertex(i) );
	  Color( AbsMin( vp->max_curvature(), vp->min_curvature() ), min, zero, max );
	  p = vp->pos();
	  glVertex3d( p.x, p.y, p.z );
	}
      glEnd();
    }
  
  glDisable(GL_CULL_FACE);
}

void TMViewer::MaxDirCurv()
{
  R3 p, d;
  double s, x1, y1, z1, x2, y2, z2;
  
  glColor3f(1.,0.,0.);
  glLineWidth(2);
  
  VertexIter vi;
  for (vi=mesh_->verts_begin(); vi!=mesh_->verts_end(); vi++)
    {
      p = (*vi)->pos();
      d = (*vi)->max_direction();
      
      s  = 0.75*(*vi)->min_edge_length();
      x1 =  s*d.x;
      y1 =  s*d.y;
      z1 =  s*d.z;
      x2 = -s*d.x;
      y2 = -s*d.y;
      z2 = -s*d.z;
  
      glPushMatrix();
      glTranslatef(p.x,p.y,p.z);
      glBegin(GL_LINES);
      glVertex3d(x1,y1,z1);
      glVertex3d(x2,y2,z2);
      glEnd();
      glPopMatrix();
    }
}

void TMViewer::MinDirCurv()
{
  R3 p, d;
  double s, x1, y1, z1, x2, y2, z2;

  glColor3f(0.,0.,1.);
  glLineWidth(2);

  VertexIter vi;
  for (vi=mesh_->verts_begin(); vi!=mesh_->verts_end(); vi++)
    {
      p = (*vi)->pos();
      d = (*vi)->min_direction();
      
      s  = 0.75*(*vi)->min_edge_length();
      x1 =  s*d.x;
      y1 =  s*d.y;
      z1 =  s*d.z;
      x2 = -s*d.x;
      y2 = -s*d.y;
      z2 = -s*d.z;
  
      glPushMatrix();
      glTranslatef(p.x,p.y,p.z);
      glBegin(GL_LINES);
      glVertex3d(x1,y1,z1);
      glVertex3d(x2,y2,z2);
      glEnd();
      glPopMatrix();
    }
}
