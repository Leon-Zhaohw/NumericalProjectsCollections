#include "adm_viewer.h"

void ADMViewer::PrimalWireframe()
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

  std::cout << "SIZE: " << mesh_->num_verts() << " " << mesh_->num_facets() << std::endl;
}

void ADMViewer::Gouraud()
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
	  n = ADM::der_cast( (*fi)->vertex(i) )->normal();
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

void ADMViewer::MaxCurv()
{ 
  glEnable(GL_CULL_FACE);

  glPolygonMode(GL_BACK,GL_LINE);
  glPolygonMode(GL_FRONT,GL_FILL);

  glLineWidth(1);

  ///////////////////////////////////////////////////
  // Computing MaxCurv //////////////////////////////
  ///////////////////////////////////////////////////
  VertexIter vi = mesh_->verts_begin();
  double max = (*vi)->max_curv();
  double min = 0.;

  vi++;
  double val;
  for ( ; vi!=mesh_->verts_end(); vi++)
    {
      val = (*vi)->max_curv();
      max = std::max( max, val );
      min = std::min( min, val );
    }
  //////////////////////////////////////////////////

  int i;
  R3 p;
  ADM::A48Vertex* vp;
  FacetIter fi;
  for (fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
    {
      glBegin(GL_TRIANGLES);
      for (i=0; i<3; i++)
	{
	  vp = ADM::der_cast( (*fi)->vertex(i) );
	  Color( vp->max_curv(), min, 0., max );
	  p = vp->pos();
	  glVertex3d( p.x, p.y, p.z );
	}
      glEnd();
    }
  
  glDisable(GL_CULL_FACE);
}

void ADMViewer::MinCurv()
{
  glEnable(GL_CULL_FACE);
  
  glPolygonMode(GL_BACK,GL_LINE);
  glPolygonMode(GL_FRONT,GL_FILL);
  
  glLineWidth(1);

  ///////////////////////////////////////////////////
  // Computing MinCurv //////////////////////////////
  ///////////////////////////////////////////////////
  VertexIter vi = mesh_->verts_begin();
  double max = (*vi)->min_curv();
  double min = 0.;

  vi++;
  double val;
  for ( ; vi!=mesh_->verts_end(); vi++)
    {
      val = (*vi)->min_curv();
      max = std::max( max, val );
      min = std::min( min, val );
    }
  //////////////////////////////////////////////////

  int i;
  R3 p;
  ADM::A48Vertex* vp;
  FacetIter fi;
  for (fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
    {
      glBegin(GL_TRIANGLES);
      for (i=0; i<3; i++)
	{
	  vp = ADM::der_cast( (*fi)->vertex(i) );
	  Color( vp->min_curv(), min, 0., max );
	  p = vp->pos();
	  glVertex3d( p.x, p.y, p.z );
	}
      glEnd();
    }
  
  glDisable(GL_CULL_FACE);
}

void ADMViewer::MeanCurv()
{
  glEnable(GL_CULL_FACE);

  glPolygonMode(GL_BACK,GL_LINE);
  glPolygonMode(GL_FRONT,GL_FILL);

  glLineWidth(1);

  ///////////////////////////////////////////////////
  // Computing MeanCurv /////////////////////////////
  ///////////////////////////////////////////////////
  VertexIter vi = mesh_->verts_begin();
  double max = (*vi)->mean_curv();
  double min = 0.;

  vi++;
  double kmean;
  for ( ; vi!=mesh_->verts_end(); vi++)
    {
      kmean = (*vi)->mean_curv();
      max = std::max( max, kmean );
      min = std::min( min, kmean );
    }
  //////////////////////////////////////////////////

  int i;
  R3 p;
  ADM::A48Vertex* vp;
  FacetIter fi; 
  for (fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
    {
      glBegin(GL_TRIANGLES);
      for (i=0; i<3; i++)
	{
	  vp = ADM::der_cast( (*fi)->vertex(i) );
	  Color( vp->mean_curv(), min, 0., max );
	  p = (*fi)->vertex(i)->pos();
	  glVertex3d( p.x, p.y, p.z );
	}
      glEnd();
    }
  
  glDisable(GL_CULL_FACE);
}

void ADMViewer::GausCurv()
{ 
  glEnable(GL_CULL_FACE);
  
  glPolygonMode(GL_BACK,GL_LINE);
  glPolygonMode(GL_FRONT,GL_FILL);
  
  glLineWidth(1);

  ///////////////////////////////////////////////////
  // Computing GausCurv /////////////////////////////
  ///////////////////////////////////////////////////
  VertexIter vi = mesh_->verts_begin();
  double max = (*vi)->gaus_curv();
  double min = 0.;

  vi++;
  double kgaus;
  for ( ; vi!=mesh_->verts_end(); vi++)
    {
      kgaus = (*vi)->gaus_curv();
      max = std::max( max, kgaus );
      min = std::min( min, kgaus );
    }
  //////////////////////////////////////////////////

  int i;
  R3 p;
  ADM::A48Vertex* vp;
  FacetIter fi;
  for (fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
    {
      glBegin(GL_TRIANGLES);
      for (i=0; i<3; i++)
	{
	  vp = ADM::der_cast( (*fi)->vertex(i) );
	  Color( vp->gaus_curv(), min, 0., max );
	  p = (*fi)->vertex(i)->pos();
	  glVertex3d( p.x, p.y, p.z );
	}
      glEnd();
    }
  
  glDisable(GL_CULL_FACE);
}

void ADMViewer::AbsMaxCurv()
{
  glEnable(GL_CULL_FACE);
  
  glPolygonMode(GL_BACK,GL_LINE);
  glPolygonMode(GL_FRONT,GL_FILL);
  
  glLineWidth(1);

  ///////////////////////////////////////////////////
  // Computing AbsMaxCurv ///////////////////////////
  ///////////////////////////////////////////////////
  VertexIter vi = mesh_->verts_begin();
  double max = AbsMax( (*vi)->max_curv(), (*vi)->min_curv() );
  double min = max;

  vi++;
  double val;
  for ( ; vi!=mesh_->verts_end(); vi++)
    {
      val = AbsMax( (*vi)->max_curv(), (*vi)->min_curv() );
      max = std::max( max, val );
      min = std::min( min, val );
    }
  
  double zero = 0.5*(max+min);
  //////////////////////////////////////////////////

  int i;
  R3 p;
  ADM::A48Vertex* vp;
  FacetIter fi;
  for (fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
    {
      glBegin(GL_TRIANGLES);
      for (i=0; i<3; i++)
	{
	  vp = ADM::der_cast( (*fi)->vertex(i) );
	  Color( AbsMax( vp->max_curv(), vp->min_curv() ), min, zero, max );
	  p = vp->pos();
	  glVertex3d( p.x, p.y, p.z );
	}
      glEnd();
    }
  
  glDisable(GL_CULL_FACE);
}

void ADMViewer::AbsMinCurv()
{ 
  glEnable(GL_CULL_FACE);

  glPolygonMode(GL_BACK,GL_LINE);
  glPolygonMode(GL_FRONT,GL_FILL);

  glLineWidth(1);

  ///////////////////////////////////////////////////
  // Computing AbsMinCurv ///////////////////////////
  ///////////////////////////////////////////////////
  VertexIter vi = mesh_->verts_begin();
  double max = AbsMin( (*vi)->max_curv(), (*vi)->min_curv() );
  double min = max;

  vi++;
  double val;
  for ( ; vi!=mesh_->verts_end(); vi++)
    {
      val = AbsMin( (*vi)->max_curv(), (*vi)->min_curv() );
      max = std::max( max, val );
      min = std::min( min, val );
    }

  double zero = 0.5*(max+min);
  //////////////////////////////////////////////////

  int i;
  R3 p;
  ADM::A48Vertex* vp;
  FacetIter fi;
  for (fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
    {
      glBegin(GL_TRIANGLES);
      for (i=0; i<3; i++)
	{
	  vp = ADM::der_cast( (*fi)->vertex(i) );
	  Color( AbsMin( vp->max_curv(), vp->min_curv() ), min, zero, max );
	  p = vp->pos();
	  glVertex3d( p.x, p.y, p.z );
	}
      glEnd();
    }
  
  glDisable(GL_CULL_FACE);
}

void ADMViewer::MaxDirCurv()
{
  R3 p, d;
  double s, x1, y1, z1, x2, y2, z2;
  
  glColor3f(1.,0.,0.);
  glLineWidth(2);
  
  VertexIter vi;
  for (vi=mesh_->verts_begin(); vi!=mesh_->verts_end(); vi++)
    {
      if ( (*vi)->is_bdry() ) continue;

      p = (*vi)->pos();
      d = (*vi)->max_direc();
      
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

void ADMViewer::MinDirCurv()
{
  R3 p, d;
  double s, x1, y1, z1, x2, y2, z2;

  glColor3f(0.,0.,1.);
  glLineWidth(2);

  VertexIter vi;
  for (vi=mesh_->verts_begin(); vi!=mesh_->verts_end(); vi++)
    {
      if ( (*vi)->is_bdry() ) continue;

      p = (*vi)->pos();
      d = (*vi)->min_direc();
      
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

void ADMViewer::Samples()
{
  double lv;
  glPointSize(4);
  glColor3f(0,0,0);

  for (FacetIter fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
    {
      lv = std::floor( (*fi)->area()* criteria_->density() );
      glBegin(GL_POINTS);
      rec_samples( int(lv), (*fi)->vertex(0)->pos(), (*fi)->vertex(1)->pos(), (*fi)->vertex(2)->pos() );
      glEnd();
    }
}

void ADMViewer::Error()
{
  glEnable(GL_CULL_FACE);

  glPolygonMode(GL_BACK,GL_LINE);
  glPolygonMode(GL_FRONT,GL_FILL);
  
  R3 v;
  double lv, ns, error;
  for (FacetIter fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
    {
      lv = std::floor( (*fi)->area() * criteria_->density() );
      ns = std::pow(4.,lv);
      error = triangle_error( int(lv), 
			      (*fi)->vertex(0)->pos(), (*fi)->vertex(1)->pos(), (*fi)->vertex(2)->pos() );
      error /= ns;

      Color( error, 0., mesh_->threshold(), error );
      glBegin(GL_TRIANGLES);
      for (int i=0; i<3; i++)
        {
          v = ADM::der_cast( (*fi)->vertex(i) )->pos();
          glVertex3f(v.x, v.y, v.z);
        }
      glEnd();
    }

  glDisable(GL_CULL_FACE);
}

/*****************************************************************/

void ADMViewer::rec_samples(int level, R3 v0, R3 v1, R3 v2)
{
  if (level==0)
    {
      R3 sample = (v0 + v1 + v2) / 3;
      glVertex3d(sample.x,sample.y,sample.z);
    } else {
      R3 subd[3];
      subd[0] = (v1 + v2)/2;
      subd[1] = (v2 + v0)/2;
      subd[2] = (v0 + v1)/2;

      level--;
      rec_samples(level,v0,subd[2],subd[1]);
      rec_samples(level,v1,subd[0],subd[2]);
      rec_samples(level,v2,subd[1],subd[0]);
      rec_samples(level,subd[0],subd[1],subd[2]);
    }
}

double ADMViewer::triangle_error(int level, R3 v0, R3 v1, R3 v2)
{
  if (level==0) { R3 b = (v0+v1+v2)/3.; return surf_->error(b); }

  R3 subd[3];
  subd[0] = 0.5*(v1 + v2);
  subd[1] = 0.5*(v2 + v0);
  subd[2] = 0.5*(v0 + v1);
  level--;
  return ( triangle_error( level, v0,      subd[2], subd[1]) +
           triangle_error( level, v1,      subd[0], subd[2]) +
           triangle_error( level, v2,      subd[1], subd[0]) +
           triangle_error( level, subd[0], subd[1], subd[2] ) );
}

/*****************************************************************/

void ADMViewer::Checker()
{
  /*
  std::cout << "CHECKER\n";

  EdgeIter ei;
  for (ei=mesh_->edges_begin(); ei!=mesh_->edges_end(); ei++)
    std::cout << "(" << (*ei)->hedge(0)->mark() << "," 
	      << (*ei)->hedge(1)->mark() << ")" << std::endl;
  */
  glEnable(GL_CULL_FACE);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_RESCALE_NORMAL);
  glEnable(GL_COLOR_MATERIAL);

  glPolygonMode(GL_BACK,GL_LINE);
  glPolygonMode(GL_FRONT,GL_FILL);
  glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);

  glLineWidth(1);
  
  int val;
  R3 p, n;
  FacetIter fi;
  
  for (fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
    {
      val = (*fi)->subd_edge()->mark();
      if (val==-1) glColor3d(0,0,0);
      else if (val==0) glColor3d(1,1,1);
      else glColor3d(1,0,0);
      
      /*
	std::cout << (*fi)->subd_edge()->org()->pos() << " -> "
	<< (*fi)->subd_edge()->dst()->pos() << " | "
	<< (*fi)->subd_edge()->mark() << std::endl;
      */
      
      glBegin(GL_TRIANGLES);
      for (int i=0; i<3; i++)
	{
	  p = (*fi)->vertex(i)->pos();
	  n = ADM::der_cast( (*fi)->vertex(i) )->normal();
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
