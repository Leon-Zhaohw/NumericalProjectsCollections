#include "tml.h"

using namespace TML;

double Vertex::min_edge_length()
{
  Halfedge* h = star_first();
  double val = h->edge()->length(); 
  for ( ; h!=NULL; h=star_next(h) )
    val = std::min( val, h->edge()->length() );
  return val;
}

double Vertex::avg_edge_length()
{
  double val = 0.; int n; 
  Halfedge* h;
  for (h=star_first(), n=0; h!=NULL; h=star_next(h), n++)
    val += h->edge()->length();
  return ( val / double(n) );
}

R3 Vertex::star_center()
{
  R3 c = R3(); int n; 
  Halfedge* h;
  for (h=star_first(), n=0; h!=NULL; h=star_next(h), n++)
    c += h->org()->pos();
  return ( c / double(n) );
}

double Vertex::compute_area()
{
  /************************************
   * Mixed Area -> Meyer et. al. 2002 *
   ************************************
   double A = 0.;
   for (Halfedge* h=star_first(); h!=NULL; h=star_next(h))
   {
   if (h->is_bdry()) continue;
      
   if ( h->next()->is_obtuse() ) 
   A += 0.5*h->facet()->area();
   else if ( h->is_obtuse() || h->prev()->is_obtuse() )
   A += 0.25*h->facet()->area();
   else
   A += ( h->edge()->length2()        *h->cotangent()        + 
   h->next()->edge()->length2()*h->next()->cotangent() ) / 8.;
   }
   return A;
  ************************************/

  double A = 0.;
  for (Halfedge* h=star_first(); h!=NULL; h=star_next(h))
    {
      if (h->is_bdry()) continue;
      A += h->facet()->area();
    }
  A /= 3.;
  return A;
}

/**************************
 * Average Faces' Normals *
 **************************/
R3 Vertex::compute_normal()
{
  R3 n = R3();
  for (Halfedge* h=star_first(); h!=NULL; h=star_next(h))
    {
      if ( h->is_bdry() ) continue;
      n += h->facet()->normal();
    }
  n.normalize();
  return n;
}

/*************************************************
 * Integrated Mean Curvature -> Meyer et. al. 02 *
 *************************************************/
double Vertex::compute_mean_curvature()
{
  R3 v, K = R3();
  for (Halfedge* h=star_first(); h!=NULL; h=star_next(h))
    {
      v = pos() - h->org()->pos();
      K += h->edge()->cotangent()*v;
    }

  K /= (4*area());  // K = mean*n

  double abs_mean = K.length();
  if ( TML::zero(abs_mean,0.) ) return 0.;
  else return (dot(K,normal())>=0.)? abs_mean : -abs_mean;
}

/*****************************************************
 * Integrated Gaussian Curvature -> Meyer et. al. 02 *
 *****************************************************/
double Vertex::compute_gaus_curvature()
{
  double kg = 2.*M_PI;
  for (Halfedge* h=star_first(); h!=NULL; h=star_next(h))
    {
      if ( h->is_bdry() ) continue;
      kg -= h->next()->angle();
    }
  kg /= area();
  
  if ( TML::zero(kg,0.) ) return 0.;
  return kg;
}

/****************************************************************************
 * Principal Curvature Values and Directions -> Cohen-Steiner and Morvan 03 *
 *               based on Alliez's implementation.                          *
 ****************************************************************************/
// k1, d1 -> max
// k2, d2 -> min
void Vertex::curvature(double& k1, R3& d1, double& k2, R3& d2)
{
  // initialize
  Matrix T(3,3);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      T[i][j] = 0.;

  // Compute curvature tensor 3x3
  double beta, len_edge;
  R3 e;
  for (Halfedge* h=star_first(); h!=NULL; h=star_next(h))
    {
      if ( h->edge()->is_bdry() ) continue;

      len_edge = h->edge()->length();
      if ( TML::zero( len_edge, 0. ) ) continue;

      beta = h->edge()->dihedral();
      e    = h->direction();
      add_tensor( T, beta*len_edge, e );
    }

  // Extract eigenvalues and eigenvectors
  Array  E(3);
  Matrix V(3,3);
  Eigen eig(T);
  eig.getRealEigenvalues(E);
  eig.getV(V);

  // Reorder them
  int indices[3] = {0,1,2};
  reorder( E, V, indices );

  d1 = R3( V[0][indices[1]], V[1][indices[1]], V[2][indices[1]] );
  d2 = R3( V[0][indices[2]], V[1][indices[2]], V[2][indices[2]] );

  //k1 = E[indices[2]];
  //k2 = E[indices[1]];

  // Setting principal curvatures values throw kmean and kgaus
  double mean  = mean_curvature();
  double delta = std::pow(mean,2) - gaus_curvature();
  if (delta<0.) delta = 0.;
  k1 = mean + std::sqrt(delta);
  k2 = mean - std::sqrt(delta);
}

void Vertex::add_tensor(Matrix& T, const double& coeff, const R3& e)
{
  T[0][0] += coeff * e.x*e.x;   T[0][1] += coeff * e.x*e.y;   T[0][2] += coeff * e.x*e.z;
  T[1][0] += coeff * e.y*e.x;   T[1][1] += coeff * e.y*e.y;   T[1][2] += coeff * e.y*e.z;
  T[2][0] += coeff * e.z*e.x;   T[2][1] += coeff * e.z*e.y;   T[2][2] += coeff * e.z*e.z;
}

// Ordering eigenvalues: [0] - normal, [1] - min curv, [2] - max curv
void Vertex::reorder(Array& E, Matrix& V, int* indices)
{
  R3 n = normal();

  double abs_dot[3];
  for (int i=0; i<3; i++)
    abs_dot[i] = std::fabs( dot( n, R3(V[0][i], V[1][i], V[2][i]) ) );
  
  if ( abs_dot[0]>abs_dot[1] && abs_dot[0]>abs_dot[2] )
    {
      indices[0] = 0;
      if (E[1]<E[2])
        { indices[1] = 1; indices[2] = 2; }
      else
        { indices[1] = 2; indices[2] = 1; }
    }
  else if ( abs_dot[1]>abs_dot[0] && abs_dot[1]>abs_dot[2] )
    {
      indices[0] = 1;
      if (E[0]<E[2])
        { indices[1] = 0; indices[2] = 2; }
      else
        { indices[1] = 2; indices[2] = 0; }
    } else {
      indices[0] = 2;
      if (E[0]<E[1])
        { indices[1] = 0; indices[2] = 1; }
      else
        { indices[1] = 1; indices[2] = 0; }
    }
}
