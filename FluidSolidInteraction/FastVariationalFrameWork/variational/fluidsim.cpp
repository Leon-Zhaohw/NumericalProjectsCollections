#include "fluidsim.h"

#include "array2_utils.h"

#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"

float fraction_inside(float phi_left, float phi_right);
void extrapolate(Array2f& grid, const Array2f& grid_weight,  Array2c& valid, Array2c old_valid);

void FluidSim::initialize(float width, int ni_, int nj_) {
   ni = ni_;
   nj = nj_;
   dx = width / (float)ni;
   u.resize(ni+1,nj); temp_u.resize(ni+1,nj); u_weights.resize(ni+1,nj);
   v.resize(ni,nj+1); temp_v.resize(ni,nj+1); v_weights.resize(ni,nj+1);
   u.set_zero();
   v.set_zero();
   nodal_solid_phi.resize(ni+1,nj+1);
   valid.resize(ni+1, nj+1);
   old_valid.resize(ni+1, nj+1);
}

//Initialize the grid-based signed distance field that dictates the position of the solid boundary
void FluidSim::set_boundary(float (*phi)(const Vec2f&)) {
   
   for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni+1; ++i) {
      Vec2f pos(i*dx,j*dx);
      nodal_solid_phi(i,j) = phi(pos);
   }

}

//The main fluid simulation step
void FluidSim::advance(float dt) {
   
   //Passively advect particles
   advect_particles(dt);

   //Advance the velocity
   advect(dt);
   add_force(dt);
   project(dt); 
   
   //Pressure projection only produces valid velocities in faces with non-zero associated face area.
   //Because the advection step may interpolate from these invalid faces, 
   //we must extrapolate velocities from the fluid domain into these zero-area faces.
   extrapolate(u, u_weights, valid, old_valid);
   extrapolate(v, v_weights, valid, old_valid);
   
   //For extrapolated velocities, replace the normal component with
   //that of the object.
   constrain_velocity();
   
}

void FluidSim::add_force(float dt) {
   int imid = ni/2;
   int jmid = nj/2;
   v(imid,jmid) = 10;
   v(imid+1,jmid) = 10;
}

//For extrapolated points, replace the normal component
//of velocity with the object velocity (in this case zero).
void FluidSim::constrain_velocity() {
   temp_u = u;
   temp_v = v;
   
   //(At lower grid resolutions, the normal estimate from the signed
   //distance function is poor, so it doesn't work quite as well.
   //An exact normal would do better.)
   
   //constrain u
   for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i) {
      if(u_weights(i,j) == 0) {
         //apply constraint
         Vec2f pos(i*dx, (j+0.5)*dx);
         Vec2f vel = get_velocity(pos);
         Vec2f normal(0,0);
         interpolate_gradient(normal, pos/dx, nodal_solid_phi); 
         normalize(normal);
         float perp_component = dot(vel, normal);
         vel -= perp_component*normal;
         temp_u(i,j) = vel[0];
      }
   }
   
   //constrain v
   for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i) {
      if(v_weights(i,j) == 0) {
         //apply constraint
         Vec2f pos((i+0.5)*dx, j*dx);
         Vec2f vel = get_velocity(pos);
         Vec2f normal(0,0);
         interpolate_gradient(normal, pos/dx, nodal_solid_phi); 
         normalize(normal);
         float perp_component = dot(vel, normal);
         vel -= perp_component*normal;
         temp_v(i,j) = vel[1];
      }
   }
   
   //update
   u = temp_u;
   v = temp_v;

}

//Add a tracer particle for visualization
void FluidSim::add_particle(const Vec2f& position) {
   particles.push_back(position);
}

//Basic first order semi-Lagrangian advection of velocities
void FluidSim::advect(float dt) {
   
   //semi-Lagrangian advection on u-component of velocity
   for(int j = 0; j < nj; ++j) for(int i = 0; i < ni+1; ++i) {
      Vec2f pos(i*dx, (j+0.5)*dx);
      pos = trace_rk2(pos, -dt);
      temp_u(i,j) = get_velocity(pos)[0];  
   }

   //semi-Lagrangian advection on v-component of velocity
   for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni; ++i) {
      Vec2f pos((i+0.5)*dx, j*dx);
      pos = trace_rk2(pos, -dt);
      temp_v(i,j) = get_velocity(pos)[1];
   }

   //move update velocities into u/v vectors
   u = temp_u;
   v = temp_v;
}

//Perform 2nd order Runge Kutta to move the particles in the fluid
void FluidSim::advect_particles(float dt) {
   
   for(int p = 0; p < particles.size(); ++p) {
      particles[p] = trace_rk2(particles[p], dt);
      
      //Particles can still occasionally leave the domain due to truncation errors, 
      //interpolation error, or large timesteps, so we project them back in for good measure.
      
      //Try commenting this section out to see the degree of accumulated error.
      float phi_value = interpolate_value(particles[p]/dx, nodal_solid_phi);
      if(phi_value < 0) {
         Vec2f normal;
         interpolate_gradient(normal, particles[p]/dx, nodal_solid_phi);
         normalize(normal);
         particles[p] -= phi_value*normal;
      }
   }

}

void FluidSim::project(float dt) {

   //Compute finite-volume type face area weight for each velocity sample.
   compute_weights();
   
   //Set up and solve the variational pressure solve.
   solve_pressure(dt);

}

//Apply RK2 to advect a point in the domain.
Vec2f FluidSim::trace_rk2(const Vec2f& position, float dt) {
   Vec2f velocity = get_velocity(position);
   velocity = get_velocity(position + 0.5f*dt*velocity);
   return position + dt*velocity;
}

//Interpolate velocity from the MAC grid.
Vec2f FluidSim::get_velocity(const Vec2f& position) {
   
   //Interpolate the velocity from the u and v grids
   float u_value = interpolate_value(position / dx - Vec2f(0, 0.5f), u);
   float v_value = interpolate_value(position / dx - Vec2f(0.5f, 0), v);
   
   return Vec2f(u_value, v_value);
}


//Given two signed distance values, determine what fraction of a connecting segment is "inside"
float fraction_inside(float phi_left, float phi_right) {
   if(phi_left < 0 && phi_right < 0)
      return 1;
   if (phi_left < 0 && phi_right >= 0)
      return phi_left / (phi_left - phi_right);
   if(phi_left >= 0 && phi_right < 0)
      return phi_right / (phi_right - phi_left);
   else
      return 0;
}

//Compute finite-volume style face-weights for fluid from nodal signed distances
void FluidSim::compute_weights() {

   for(int j = 0; j < u_weights.nj; ++j) for(int i = 0; i < u_weights.ni; ++i) {
      u_weights(i,j) = 1 - fraction_inside(nodal_solid_phi(i,j+1), nodal_solid_phi(i,j));
   }
   for(int j = 0; j < v_weights.nj; ++j) for(int i = 0; i < v_weights.ni; ++i) {
      v_weights(i,j) = 1 - fraction_inside(nodal_solid_phi(i+1,j), nodal_solid_phi(i,j));
   }

}

//An implementation of the variational pressure projection solve for static geometry
void FluidSim::solve_pressure(float dt) {
   
   //This linear system could be simplified, but I've left it as is for clarity 
   //and consistency with the standard naive discretization
   
   int ni = v.ni;
   int nj = u.nj;
   int system_size = ni*nj;
   if(rhs.size() != system_size) {
      rhs.resize(system_size);
      pressure.resize(system_size);
      matrix.resize(system_size);
   }
   matrix.zero();

   //Build the linear system for pressure
   for(int j = 1; j < nj-1; ++j) {
      for(int i = 1; i < ni-1; ++i) {
         int index = i + ni*j;
         
         rhs[index] = 0;

         //right neighbour
         float term = u_weights(i+1,j) * dt / sqr(dx);
         matrix.add_to_element(index, index, term);
         matrix.add_to_element(index, index + 1, -term);
         rhs[index] -= u_weights(i+1,j)*u(i+1,j) / dx;
         
         //left neighbour
         term = u_weights(i,j) * dt / sqr(dx);
         matrix.add_to_element(index, index, term);
         matrix.add_to_element(index, index - 1, -term);
         rhs[index] += u_weights(i,j)*u(i,j) / dx;
         
         //top neighbour
         term = v_weights(i,j+1) * dt / sqr(dx);
         matrix.add_to_element(index, index, term);
         matrix.add_to_element(index, index + ni, -term);
         rhs[index] -= v_weights(i,j+1)*v(i,j+1) / dx;
         
         //bottom neighbour
         term = v_weights(i,j) * dt / sqr(dx);
         matrix.add_to_element(index, index, term);
         matrix.add_to_element(index, index - ni, -term);
         rhs[index] += v_weights(i,j)*v(i,j) / dx;
      }
   }

   //Solve the system using Robert Bridson's incomplete Cholesky PCG solver
   
   double tolerance;
   int iterations;
   bool success = solver.solve(matrix, rhs, pressure, tolerance, iterations);
   if(!success)
      printf("WARNING: Pressure solve failed!\n");
   
   //Apply the velocity update
   for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i) {
      int index = i + j*ni;
      if(u_weights(i,j) > 0)
         u(i,j) -= dt  * (pressure[index] - pressure[index-1]) / dx; 
      else
         u(i,j) = 0;
   }
   for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i) {
      int index = i + j*ni;
      if(v_weights(i,j) > 0) 
         v(i,j) -= dt  * (pressure[index] - pressure[index-ni]) / dx; 
      else
         v(i,j) = 0;
   }

}


//Apply several iterations of a very simple "Jacobi"-style propagation of valid velocity data in all directions
void extrapolate(Array2f& grid, const Array2f& grid_weight, Array2c& valid, Array2c old_valid) {

   //Initialize the list of valid cells
   for(int j = 0; j < grid.nj; ++j) for(int i = 0; i < grid.ni; ++i)
      valid(i,j) = grid_weight(i,j) > 0;
   
   for(int layers = 0; layers < 5; ++layers) {
      old_valid = valid;
      for(int j = 1; j < grid.nj-1; ++j) for(int i = 1; i < grid.ni-1; ++i) {
         float sum = 0;
         int count = 0;
         
         if(!old_valid(i,j)) {
            
            if(old_valid(i+1,j)) {
               sum += grid(i+1,j);\
               ++count;
            }
            if(old_valid(i-1,j)) {
               sum += grid(i-1,j);\
               ++count;
            }
            if(old_valid(i,j+1)) {
               sum += grid(i,j+1);\
               ++count;
            }
            if(old_valid(i,j-1)) {
               sum += grid(i,j-1);\
               ++count;
            }

            //If any of neighbour cells were valid, 
            //assign the cell their average value and tag it as valid
            if(count > 0) {
               grid(i,j) = sum /(float)count;
               valid(i,j) = 1;
            }

         }
      }

   }

}
