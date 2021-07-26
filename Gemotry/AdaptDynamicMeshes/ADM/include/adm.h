#ifndef _ADM_
#define _ADM_ 1

#include <tml.h>

namespace ADM
{
  // Classes
  class AdaptiveMesh;
  class A48Vertex;
  class A48Facet;
  class A48Edge;

  class AdaptCriteria;
  class StochasticSampling;

  class Oracle;

  // Typedefs
  typedef TML::R3       R3;
  typedef TML::Halfedge Halfedge;
  typedef TML::Edge     Edge;
  typedef TML::Vertex   Vertex;
  typedef TML::Facet    Facet;

  typedef TML::Matrix Matrix;
  typedef TML::Array  Array;
  typedef TML::LU     LU;
  typedef TML::Eigen  Eigen;
}

#include "edge.h"
#include "facet.h"
#include "vertex.h"
#include "oracle.h"
#include "criteria.h"
#include "mesh.h"

#endif
