#ifndef _TML_VIEWER_
#define _TML_VIEWER_ 1

#include <tml.h>
#include "glViewer.h"

class Point : public TML::Vertex
{
 protected:
  double  area_;
  double  kmax_, kmin_, kmean_, kgaus_;
  TML::R3 dmax_, dmin_, normal_;

 public:
  Point(TML::R3 p=TML::R3()) 
    : TML::Vertex(p), 
    area_(0.), kmax_(0.), kmin_(0.), kmean_(0.), kgaus_(0.),
    dmax_(), dmin_(), normal_() 
    { } 

  void set_attr(int i)
    {
      set_id(i);
      normal_ = compute_normal();
      area_   = compute_area();
      kmean_  = compute_mean_curvature();
      kgaus_  = compute_gaus_curvature();
      curvature( kmax_, dmax_, kmin_, dmin_ );
    }

  inline double  area()   { return area_; }
  inline TML::R3 normal() { return normal_; }
  inline double mean_curvature() { return kmean_; }
  inline double gaus_curvature() { return kgaus_; }
  inline double max_curvature()  { return kmax_;  }
  inline double min_curvature()  { return kmin_;  }
  inline TML::R3 max_direction() { return dmax_;  }
  inline TML::R3 min_direction() { return dmin_;  }
};

class TMViewer : public glViewer
{
 public:
  typedef TML::R3 R3;
  typedef TML::TriMesh<Point>::VertexIter VertexIter;
  typedef TML::TriMesh<Point>::EdgeIter   EdgeIter;
  typedef TML::TriMesh<Point>::FacetIter  FacetIter;
  
 private:
  TML::TriMesh<Point>* mesh_;
  
  bool wprimal_, wdual_, dmax_, dmin_;

 public:
  TMViewer(TML::TriMesh<Point>* mesh);
  ~TMViewer();
  
  void DrawModel();
  void saveASoff();
  int handle(int e);

  void KeyboardOptions()
    {
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
      std::cout << "\t q - max principal direction" << std::endl;
      std::cout << "\t w - min principal direction" << std::endl;
      std::cout << "\t e - primal mesh" << std::endl;
      std::cout << "\t r - dual mesh" << std::endl;
      std::cout << "\t F1 - Gouraud" << std::endl;
      std::cout << "\t F2 - max curvature [min - 0 - max]" << std::endl;
      std::cout << "\t F3 - min curvature [min - 0 - max]" << std::endl;
      std::cout << "\t F4 - mean curvature    [min - 0 - max]" << std::endl;
      std::cout << "\t F5 - gaussian curvature [min - 0 - max]" << std::endl;
      std::cout << "\t F6 - Max Abs curvature [min - max]" << std::endl;
      std::cout << "\t F7 - Min Abs curvature [min - max]" << std::endl;
      std::cout << "\t F8 - none" << std::endl;
      glViewer::KeyboardOptions();
    }
  
 protected:
  void CreateDisplayLists();
  void SetPosition();

  void PrimalWireframe(); // 0
  void DualWireframe();   // 1
  void Gouraud();         // 4
  void MaxCurv();         // 5
  void MinCurv();         // 6
  void MeanCurv();        // 7
  void GausCurv();        // 8
  void AbsMaxCurv();      // 9
  void AbsMinCurv();      // 10
  void MaxDirCurv();      // 2
  void MinDirCurv();      // 3

  inline double AbsMax(double a, double b)
    { return std::max( std::fabs(a), std::fabs(b) ); }
  inline double AbsMin(double a, double b)
    { return std::min( std::fabs(a), std::fabs(b) ); }
};

#endif
