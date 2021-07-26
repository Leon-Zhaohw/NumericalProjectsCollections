#include "sinoidal.h"
#include "drawing.h"
#include "timestamp.h"

int main(int argc, char** argv)
{
  ParamSin surf( 0.8, 0.5*M_PI, 0.5*M_PI, 0.);

  ADM::StochasticSampling criteria( &surf );

  ADM::AdaptiveMesh mesh( "../base_mesh/plane.off", &criteria );
  
  ADMViewer win( "Sinoidal", &mesh, &surf, &criteria );
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
