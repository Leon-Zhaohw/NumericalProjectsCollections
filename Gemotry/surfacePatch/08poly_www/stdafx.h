#pragma once

// Disable warnings about deprecated string functions
#pragma warning( disable:4995 )

#include <dxut.h>
#include <d3dx10.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <math.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

using namespace std;

#include "PolygonMesh.h"

//typedef Kernel::Vector_3                                     Vector;
typedef PolygonMesh                                          Polyhedron;
typedef Polyhedron::Vertex                                   Vertex;
typedef Polyhedron::Halfedge                                 Halfedge;
typedef Polyhedron::Facet                                    Facet;
typedef Polyhedron::Vertex_handle                            Vertex_handle;
typedef Polyhedron::Halfedge_handle                          Halfedge_handle;
typedef Polyhedron::Facet_handle                             Facet_handle;
typedef Polyhedron::Vertex_iterator                          Vertex_iterator;
typedef Polyhedron::Edge_iterator                            Edge_iterator;
typedef Polyhedron::Facet_iterator                           Facet_iterator;
typedef Polyhedron::Halfedge_around_vertex_circulator        HV_circulator;
typedef Polyhedron::Halfedge_around_facet_circulator         HF_circulator;

// Global variables
string name;

// Set packing to 1 so the compiler doesn't generate any surprise padding
//#pragma pack(push,1) // never use #pragma pack without a corresponding pop

void RenderText();
