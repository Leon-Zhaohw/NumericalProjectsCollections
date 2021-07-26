#include "surfaces/levelset.h"
#include "drawing.h"
#include "timestamp.h"
//#include "off2pov.h"

class Enright : public CPLS3D::Velocity3D
{
protected:
  double s[3];
  
public:
  Enright(int nx, int ny, int nz) 
  { 
    s[0] = double(nx);
    s[1] = double(ny);
    s[2] = double(nz);
  }
  
  virtual ~Enright() { }
  
  void GetVelocity(const CPLS3D::Vec3D& p, CPLS3D::Vec3D& u)
  {
    // dt <= 4.9 / (2*nx+ny+nz)
    
    double x, y, z;
    x = p[0]/s[0];
    y = p[1]/s[1];
    z = p[2]/s[2];
    
    u = CPLS3D::Vec3D( s[0]*2*std::pow(std::sin(M_PI*x),2)*std::sin(2*M_PI*y)*std::sin(2*M_PI*z),
		       - s[1]*std::sin(2*M_PI*x)*std::pow(std::sin(M_PI*y),2)*std::sin(2*M_PI*z),
		       - s[2]*std::sin(2*M_PI*x)*std::sin(2*M_PI*y)*std::pow(std::sin(M_PI*z),2) );
    
    u *= std::cos( M_PI*time() / 1.5 );
  }
};

CPLS3D::Grid3D Sphere(int nx, int ny, int nz, const CPLS3D::Vec3D& c, double radius)
{
  CPLS3D::Grid3D grid(nx,ny,nz,HH);
  for(int k=0; k<(nz+2); k++)
    for(int j=0; j<(ny+2); j++)
      for(int i=0; i<(nx+2); i++)
	grid(i,j,k) = ( c - CPLS3D::Vec3D(i,j,k) ).Length() - radius;
  return grid;
}

/*-------------------------------------------------------------------------------*/

int main(int argc, char** argv)
{
  int nx, ny, nz;
  nx = NX;
  ny = NY;
  nz = NZ;

  double radius = 0.15*double(ny);
  double dt = 4.8 / (2*double(nx) + double(ny) + double(nz));
  
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  
  Enright vel(nx,ny,nz);
  
  CPLS3D::PLS3D::RESEED_FLAG = 5;
  CPLS3D::PLS3D pls( nx, ny, nz, HH, 
		     Sphere( nx, ny, nz, CPLS3D::Vec3D(0.35*nx,0.35*ny,0.35*nz), radius ),
		     &(vel), dt );
  
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  LevelSet surf( &pls );

  ADM::StochasticSampling criteria( &surf );
  criteria.set_density( criteria.density() / std::pow(radius,2) );
  
  ADM::AdaptiveMesh mesh( "./enright100b.off", &criteria );
  mesh.set_threshold( radius*0.01 );
  mesh.set_stiffness( -1 );
  mesh.set_smooth_steps( 10 );

  ADMViewer win( "Enright", &mesh, &surf, &criteria );
  //win.KeyboardOptions();
  win.show();
  Fl::run();
  
  /*
    char file[100];
    sprintf( file, "m_enright0.off" );
    mesh.adapt();
    mesh.save_mesh( file );
    while ( vel.time() <= 1.5 )
    {
    mesh.deform();
    mesh.adapt();
    sprintf( file, "m_enright/off/m_enright%.3f.off", vel.time() );
    mesh.save_mesh( file );
    //sprintf( file, "pov/enright/enright%.3f.inc", vel.time() );
    //save_pov(mesh, file);
    std::cout << "Save " << file << std::endl;
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

  return 1;
}
