
#include <fstream>

#include "dynamicsurface.h"
#include "surftrack.h"
#include "pressure_vor2d.h"
#include "cgsolver.h"



std::vector<double> pressure_solve_voronoi(TriMesh2D& mesh, 
                                           eltopo2d::SurfTrack& surface,
                                          std::vector<float>& face_velocities, 
                                          const std::vector<float>& solid_weights, 
                                          const std::vector<float>& liquid_phi,
                                          const std::vector<float>& wall_velocities) 
{


      
   float thetaClamp = 0.001f;
   //Clamping it slightly away from zero (0.001f) avoids division by zero
   
   SparseMatrixd matrix(mesh.vertices.size());
   std::vector<double> rhs(mesh.vertices.size(), 0);
   std::vector<double> pressure(mesh.vertices.size(), 0);
   
   std::cout << "pressure_solve_voronoi: Building matrix" << std::endl;
   
   for(unsigned int i = 0; i < mesh.vertices.size(); ++i) {
      //consider an element, and set up its row
      float self_phi = liquid_phi[i];
      if(self_phi >= 0) continue;

      bool empty = true;
      int nbr_count = mesh.vert_to_edge_map[i].size();
      for(int j = 0; j < nbr_count; ++j) {
         //consider each face of an element, and set up its column
         int edge = mesh.vert_to_edge_map[i][j];
         if(solid_weights[edge] > 1e-8)
            empty = false;
      }
      
      if(empty) continue;
      
      for(int j = 0; j < nbr_count; ++j) {
         //consider an edge of a tri, and set up its column
         int edge = mesh.vert_to_edge_map[i][j];
         
         int nbr_vert = mesh.edges[edge][0]==i? mesh.edges[edge][1] : mesh.edges[edge][0];
         
         if(nbr_vert == -1)
            continue;

         float dist = mesh.edge_lengths[edge];
         
         if(dist < 1e-7f) {
            dist = 1e-7f;
         }

         float nbr_phi = liquid_phi[nbr_vert];

         if(nbr_phi < 0) {
            matrix.add_to_element(i,i, solid_weights[edge] * mesh.dual_edge_lengths[edge] / dist);
            matrix.add_to_element(i,nbr_vert, -solid_weights[edge] * mesh.dual_edge_lengths[edge] / dist);
         }
         else  {
            float theta = max(thetaClamp, self_phi / (self_phi - nbr_phi));
            matrix.add_to_element(i,i, solid_weights[edge] * mesh.dual_edge_lengths[edge] / dist / theta);
         }

         rhs[i] -= solid_weights[edge] * (face_velocities[edge] * mesh.dual_edge_signs[i][j]) * mesh.dual_edge_lengths[edge];
         
      }
   }
   
   std::cout << "pressure_solve_voronoi: Solving.." << std::endl;
   
   //solve
   CGSolver<double> solver;
   solver.tolerance_factor = 1e-12;
   solver.max_iterations = 1000;
 
   std::vector<double> sln(rhs.size()); 

   double residual;
   int iterations;

   solver.solve(matrix, rhs, pressure, residual, iterations);
      
   for(unsigned int i = 0; i < mesh.edges.size(); ++i) {
      Vec2ui verts = mesh.edges[i]; //already ordered correctly for the next step
      
      if(solid_weights[i] > 0) {
         float phi0 = liquid_phi[verts[0]];
         float phi1 = liquid_phi[verts[1]];
         
         if(phi0 >= 0 && phi1 >= 0) {
            face_velocities[i] = 0.0f;
            continue; //skip air faces
         }
         
         double p0 = pressure[verts[0]];
         double p1 = pressure[verts[1]];
         float theta = 1;
         if(phi0 >= 0) {
            p0 = 0;
            theta = max(thetaClamp, phi1 / (phi1 - phi0));
         }
         if(phi1 >= 0) {
            p1 = 0;
            theta = max(thetaClamp, phi0 / (phi0 - phi1));
         }
         float dist = mesh.edge_lengths[i];
         
         if(dist < 1e-7f) {
            dist = 1e-7f;
         }
         face_velocities[i] -= (float)(p1 - p0) / dist /  theta; 
         
      }
      else {
         face_velocities[i] = 0.0f;
      }
   }
   
   return pressure;
}