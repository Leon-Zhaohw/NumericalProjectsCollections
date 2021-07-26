#ifndef _VIEWER_
#define _VIEWER_

#include "TMViewer.h"

int main(int argc, char** argv)
{
  if (argc!=2)
    { std::cout << "Usage: ./viewer <off | ply | smf | ifs>\n"; std::exit(1); }

  TML::TriMesh<Point> mesh(argv[1]);

  std::cout << mesh.num_verts() << " "
            << mesh.num_edges()    << " "
            << mesh.num_facets()   << std::endl;
  
  TMViewer win( &mesh );
  //win.KeyboardOptions();
  win.show();
  Fl::run();

  return 1;
}

#endif
