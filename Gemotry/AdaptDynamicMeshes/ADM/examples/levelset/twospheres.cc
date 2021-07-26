#include "surfaces/levelset.h"
#include "drawing.h"
#include "timestamp.h"

class MinusNormal : public CPLS3D::Velocity3D
{
protected:
  CPLS3D::Grid3D grid;

public:
  MinusNormal(int nx, int ny, int nz) : Velocity3D(), grid(nx,ny,nz,HH) { }

  virtual ~MinusNormal() { }

  void set_grid(const CPLS3D::Grid3D& g)
  { grid = g; }

  void GetVelocity(const CPLS3D::Vec3D& p, CPLS3D::Vec3D& u)
  {
    grid.Normal(p,u);
    u *= -1.;
  }
};

/*-------------------------------------------------------------------------------*/

void TwoSpheres(CPLS3D::Grid3D& grid, const CPLS3D::Vec3D& c1, const CPLS3D::Vec3D& c2, double radius)
{
  int nx = grid.GetNx();
  int ny = grid.GetNy();
  int nz = grid.GetNz();
  double val1, val2;

  for(int k=0; k<(nz+2); k++)
    for(int j=0; j<(ny+2); j++)
      for(int i=0; i<(nx+2); i++)
        {
          val1 = ( c1 - CPLS3D::Vec3D(i,j,k) ).Length() - radius;
          val2 = ( c2 - CPLS3D::Vec3D(i,j,k) ).Length() - radius;
          grid(i,j,k) = std::min( val1, val2 );
        }
}

/*-------------------------------------------------------------------------------*/

int main(int argc, char** argv)
{
  int nx, ny, nz;
  nx = 100;
  ny = 100;
  nz = 100;

  double radius = 0.12*double(ny);
  double dt = 1.; //0.5;
  
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  
  CPLS3D::Grid3D grid(nx,ny,nz,HH);
  TwoSpheres( grid, CPLS3D::Vec3D(0.4*nx,0.5*ny,0.5*nz), CPLS3D::Vec3D(0.6*nx,0.5*ny,0.5*nz), radius );
  
  MinusNormal vel(nx,ny,nz);
  
  CPLS3D::PLS3D::RESEED_FLAG = 5;
  CPLS3D::PLS3D pls( nx, ny, nz, HH, grid, &(vel), dt );
  
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  
  LevelSet surf( &pls );

  ADM::StochasticSampling criteria( &surf );
  criteria.set_density( criteria.density() / std::pow(radius,2) );
  
  ADM::AdaptiveMesh mesh( "./twospheres100.off", &criteria );
  mesh.set_threshold( radius*0.01 );
  mesh.set_stiffness( -1 );
  mesh.set_smooth_steps( 5 );

  ADMViewer win( "TwoSpheres", &mesh, &surf, &criteria );
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
  int i, num_iters = 4;
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
