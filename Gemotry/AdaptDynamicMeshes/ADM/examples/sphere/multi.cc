#include "sphere.h"
#include "multi_criteria.h"
#include "drawing.h"

int main(int argc, char** argv)
{
  Sphere surf(1.);

  MultiCriteria criteria( &surf, 0.2 );

  ADM::AdaptiveMesh mesh( "../base_mesh/cube.off", &criteria );
  mesh.set_threshold( 0.05 );

  ADMViewer win( "MultiSphere", &mesh, &surf, &criteria );
  //win.KeyboardOptions();
  win.show();
  Fl::run();

  return 1;
}
