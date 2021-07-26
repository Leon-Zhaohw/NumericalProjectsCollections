#ifndef _HEIGHT_SIN_
#define _HEIGHT_SIN_ 1

#include "sinc.h"
//#include "drawing.h"
#include "viewer.h"
#include "timestamp.h"
//#include "save_pov.h"

int main(int argc, char** argv)
{
  HeightSinc surf( 0., 3*M_PI, 0.9, 0., 0. );
  
  ADM::StochasticSampling criteria( &surf );

  ADM::AdaptiveMesh mesh( "../base_mesh/plane.off", &criteria );
  //mesh.set_threshold(0.0035);
  //mesh.set_hysteresis(5e-4);

  //ADMViewer win( "Height_sinc", &mesh, &surf, &criteria );
  Viewer win( "Height_sinc", &mesh, &surf, &criteria );
  //win.KeyboardOptions();
  win.show();
  Fl::run();
  
  /*
  // init
  Timestamp init_tstart; CPUtimer  init_cstart;
  mesh.adapt();
  CPUtimer  init_cend; Timestamp init_tend;
  std::cout << "INIT: \n"
  << "\t CPU " << (init_cend - init_cstart) 
  << "\n \t TOTAL " << (init_tend - init_tstart)
  << std::endl;
  
  // loop
  int i, num_iters = 30;
  double cpu_deform, cpu_adapt, total_deform, total_adapt;

  Timestamp tstart, tend;
  CPUtimer  cstart, cend;

  cpu_deform = cpu_adapt = total_deform = total_adapt = 0.;
  for (i=0; i<num_iters; i++)
  {
  tstart.now(); cstart.now();
  mesh.deform();
  cend.now(); tend.now(); 
  
  cpu_deform   += (cend - cstart);
  total_deform += (tend - tstart);
  
  tstart.now(); cstart.now();
  mesh.adapt();
  cend.now(); tend.now();

  cpu_adapt   += (cend - cstart);
  total_adapt += (tend - tstart);
  }
  
  double n = double(num_iters);
  std::cout << "STATS: \n" 
  << "       \t CPU  \t TOTAL \n"
  << "deform \t " << cpu_deform/n << " \t " << total_deform/n << std::endl
  << "adapt  \t " << cpu_adapt/n << " \t "  << total_adapt/n  << std::endl;
  */

  /*
    int i, num_iters = 81;
    char filename[100];
    for (i=0; i<num_iters; i++)
    {
    if (i>0) mesh.deform();
    mesh.adapt();
    
    sprintf( filename, "./pov/mesh/height%.3d.inc", i );
    save_pov( mesh, filename );
    }
  */

  return 1;
}

#endif
