#include "lp3d.h"
#include "drawing.h"
#include "timestamp.h"
//#include "off2pov.h"

int main(int argc, char** argv)
{
  Lp3D surf( 2., 1. );

  ADM::StochasticSampling criteria( &surf );

  ADM::AdaptiveMesh mesh( "../base_mesh/cube.off", &criteria );
  //mesh.set_threshold(0.005);
  
  ADMViewer win( "CubeLp", &mesh, &surf, &criteria );
  //win.KeyboardOptions();
  win.show();
  Fl::run();

  /*
    int i, num_iters = 30;
    char filename[100];
    
    for (i=0; i<num_iters; i++)
    {
    if (i>0) mesh.deform();
    mesh.adapt();
    
    //sprintf( filename, "./mesh/cubelp%.3d.off", i );
    //mesh.save_mesh( filename );
    
    sprintf( filename, "./pov/mesh/cubelp0%.3d.inc", i );
    save_pov( mesh, filename, 0 );
    sprintf( filename, "./pov/mesh/cubelp1%.3d.inc", i );
    save_pov( mesh, filename, 1 );
    }
  */
  
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
  int i, num_iters = 10;
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

  return 1;
}
