#ifndef _TML_
#define _TML_ 1

// STD
#include <iostream>
#include <exception>
#include <fstream>
#include <string>
#include <cmath>
#include <set>
#include <map>
#include <vector>
#include <stack>
#include <iterator>
#include <algorithm>

// TNT
#include <tnt.h>
#include <jama_lu.h>
#include <jama_eig.h>

// EXTRA
#include "heap.h"
#include "rply.h"

namespace TML
{
  // Classes
  template <class V, class E, class F> class TriMesh;
  class Facet;
  class Vertex;
  class Edge;
  class Halfedge;

  class TriMeshIO;
  class Polygon;
  class R3;

  class Markable;
  class Error;

  // Typedefs
  typedef TNT::Array2D<double>     Matrix;
  typedef TNT::Array1D<double>     Array;
  typedef JAMA::LU<double>         LU;
  typedef JAMA::Eigenvalue<double> Eigen;

  typedef std::pair<int,int> Ipair;
  
  /*------------------------------------------------------------*/

  // class Markable
  class Markable 
    {
    protected:
      int mark_; 
    public:
      Markable() : mark_(-1) { };
      void set_mark(int m) { mark_ = m; };
      int  mark() const { return mark_; };
    };
  
  // class Error
  class Error 
    {
    public:
      Error(char *s="")
	{ std::cerr << "TML::Error: " << s << std::endl; exit(0); }
    };
}

#include "r3.h"
#include "io.h"
#include "edge.h"
#include "vertex.h"
#include "facet.h"
#include "trimesh.h"

#endif
