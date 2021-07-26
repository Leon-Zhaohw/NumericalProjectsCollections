#ifndef _TYPES_H_
#define _TYPES_H_ 1

// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// local
#include "primitives.h"
#include "dt2.h"

#undef min
#undef max

// Kernel
// typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// Basic types
typedef Kernel::FT         FT;
typedef Kernel::Point_2    Point;
typedef Kernel::Vector_2   Vector;
typedef Kernel::Ray_2      Ray;
typedef Kernel::Line_2     Line;
typedef Kernel::Segment_2  Segment;
typedef Kernel::Triangle_2 Triangle;

// Triangulation
typedef CGAL::Triangulation_vertex_base_2<Kernel> Vb;
typedef My_vertex_base<Kernel,Vb> MVb;
typedef CGAL::Triangulation_face_base_2<Kernel> Fb;
typedef My_face_base<Kernel,Fb> MFb;
typedef CGAL::Triangulation_data_structure_2<MVb,MFb> TDS;
typedef CGAL::Delaunay_triangulation_2<Kernel,TDS> MDT;
typedef Triangulation<MDT> DT;

typedef DT::Vertex                   Vertex;
typedef DT::Vertex_handle            Vertex_handle;
typedef DT::Vertex_iterator          Vertex_iterator;
typedef DT::Vertex_circulator        Vertex_circulator;
typedef DT::Finite_vertices_iterator Finite_vertices_iterator;

typedef DT::Edge                  Edge;
typedef DT::Edge_iterator         Edge_iterator;
typedef DT::Edge_circulator       Edge_circulator;
typedef DT::Finite_edges_iterator Finite_edges_iterator;

typedef DT::Face                  Face;
typedef DT::Face_handle           Face_handle;
typedef DT::Face_iterator         Face_iterator;
typedef DT::Face_circulator       Face_circulator;
typedef DT::Finite_faces_iterator Finite_faces_iterator;

typedef DT::Vertex_handle_map Vertex_handle_map;
typedef DT::Face_handle_map Face_handle_map;

typedef DT::Vertex_handle_set Vertex_handle_set;
typedef DT::Edge_set Edge_set;

typedef DT::Edge_list Edge_list;

typedef DT::Cost Cost;
typedef DT::Sample Sample;
typedef DT::Sample_list Sample_list;
typedef DT::Sample_list_const_iterator Sample_list_const_iterator;

typedef DT::Point_list Point_list;
typedef DT::Point_list_const_iterator Point_list_const_iterator;

typedef DT::PSample PSample;
typedef DT::SQueue  SQueue;

typedef DT::PEdge  PEdge;
typedef DT::PQueue PQueue;

#endif
