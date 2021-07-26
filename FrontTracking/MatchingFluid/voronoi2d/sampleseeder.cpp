/*
 *  sampleseeder.cpp
 *  eltopo2d_project
 *
 *  Created by tyson on 17/11/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "sampleseeder.h"
#include "dynamicsurface.h"


void SampleSeeder::generate_samples( const eltopo2d::DynamicSurface& surface, 
                                     double desired_dx,
                                     std::vector<Vec2f>& samples, 
                                     std::vector<int>& expected_signs)
{
   // for each vertex on the surface
   
   for ( unsigned int i = 0; i < surface.m_positions.size(); ++i )
   {
      if ( surface.m_mesh.vtxedge[i].size() < 2 ) { continue; }
      
      const Vec2d& ray_origin = surface.m_positions[i];
      const Vec2d normal = surface.get_vertex_normal(i);

      if ( mag(normal) < 0.2 ) { continue; }
      
      //
      // fire a ray in the positive and negative normal direction
      //
      
      // ignore incident edges
      const std::vector<unsigned int>& incident_edges = surface.m_mesh.vtxedge[i];
      assert( incident_edges.size() == 2 );
      
      for ( int sign = -1; sign < 2; sign += 2 )
      {      
         Vec2d ray_end = ray_origin + (double)(sign) * desired_dx * normal;
         
         std::vector<double> hit_ss;
         std::vector<unsigned int> hit_edges; 
         surface.get_segment_collisions( ray_origin, ray_end, hit_ss, hit_edges );
                  
         // get first ray intersection time (bounded by desired_dx)
         double min_hit = desired_dx;
         for ( unsigned int j = 0; j < hit_ss.size(); ++j )
         {
            if ( hit_edges[j] == incident_edges[0] || hit_edges[j] == incident_edges[1] ) { continue; }
            
            min_hit = min( min_hit, desired_dx * hit_ss[j] );
         }
      
         //std::cout << "desired_dx: " << desired_dx << ", min_hit: " << min_hit << std::endl;
         
         Vec2f new_sample;
         if ( sign < 0.0 )
         {
            new_sample = Vec2f( ray_origin - 0.5 * min_hit * normal ); //outside
         }
         else
         {
            new_sample = Vec2f( ray_origin + 0.25 * min_hit * normal ); //inside
         }
         
         std::vector<double> test_ss;
         std::vector<unsigned int> test_edges; 
         
         surface.get_segment_collisions( ray_origin, Vec2d(new_sample), test_ss, test_edges );

         for ( unsigned int j = 0; j < test_ss.size(); ++j )
         {
            if ( test_edges[j] == incident_edges[0] || test_edges[j] == incident_edges[1] ) { continue; }
            std::cout << "desired_dx = " << desired_dx << std::endl;
            std::cout << "min_hit = " << min_hit << std::endl;
            std::cout << "dist( ray_origin, new_sample ) = " << dist( ray_origin, Vec2d(new_sample) ) << std::endl;
            std::cout << "test_ss * desired_dx = " << test_ss[j] * desired_dx << std::endl;
            assert(0);
         }            
         
         assert( new_sample[0] == new_sample[0] );
         assert( new_sample[1] == new_sample[1] );
         
         bool already_exists = false;
         for ( unsigned int i = 0; i < samples.size(); ++i )
         {
            if ( dist( new_sample, samples[i] ) < 1e-6 ) 
            {
               already_exists = true;
               break;
            }               
         }
         
         if ( !already_exists )
         {
            //only add a point if either it's inside, or it's outside and the next surface is further than some threshold.
            if(sign > 0 || min_hit > 0.2*desired_dx) {
               samples.push_back( new_sample );
               if(sign > 0)
                  expected_signs.push_back(-1);
               else
                  expected_signs.push_back(0); //don't bother, since we don't care about this case (thin bit of air in a fluid misidentified)
            }
            
         }
         
      }
      
   }
   
}



void SampleSeeder::generate_regular_samples( const Vec2d& domain_min,
                                             const Vec2d& domain_max,
                                             double dx,
                                             std::vector<Vec2f>& samples,
                                             std::vector<int>& expected_signs  )
{

   float dy =  (float)(0.5f * dx * tan(60*M_PI/180));

   int cols = (int)((domain_max[0] - domain_min[0]) / dx);
   int rows = (int)((domain_max[1] - domain_min[1]) / dy);
   
   for (int row = 0; row <= rows; ++row) 
   {
      for (int col = 0; col <= cols; ++col) 
      {
         int rparity = row % 2;
         if(rparity == 0) 
         {
            float i = (float)col;
            float j = (float)row;
            samples.push_back( Vec2f(domain_min) + Vec2f( i*(float)dx, j*dy) );
            expected_signs.push_back(0); //no particular sign is expected.
         }
         else 
         {
            float i = (float)col + 0.5f;
            float j = (float)row;
            samples.push_back( Vec2f(domain_min) + Vec2f( i*(float)dx, j*dy) );
            expected_signs.push_back(0); //no particular sign is expected.
         }
         
         
      }
      
   }
  
   
}



