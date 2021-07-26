
#include "eltopo.h"
#include "surftrack.h"


void el_topo_static_operations( const int in_num_vertices,                   // input
                                const double* in_vertex_locations,           
                                const double* in_vertex_masses,                 
                                const int in_num_edges,                   
                                const int* in_edges, 
                                const ElTopoGeneralOptions& general_options,       // options
                                const ElTopoStaticOperationsOptions& options, 
                                int *out_num_vertices,                       // output
                                double** out_vertex_locations,  
                                double** out_vertex_masses,       
                                int *out_num_edges,
                                int** out_edges ) 
{
   //
   // data wrangling
   //
   
   std::vector<Vec2d> vs;
   std::vector<double> masses;
   for (  int i = 0; i < in_num_vertices; ++i )
   {
      vs.push_back( Vec2d( in_vertex_locations[2*i], in_vertex_locations[2*i + 1] ) );
      masses.push_back( in_vertex_masses[i] );
   }

   std::vector<Vec2ui> es;
   for (  int i = 0; i < in_num_edges; ++i )
   {
      es.push_back( Vec2ui( in_edges[2*i], in_edges[2*i + 1] ) );
   }
   

   // =================================================================================
   
   //
   // do the actual operations
   //
   
   // build a SurfTrack
   eltopo2d::SurfTrackInitializationParameters construction_parameters;
   
   construction_parameters.m_verbose = general_options.m_verbose;
   construction_parameters.m_collision_safety = general_options.m_collision_safety;
   construction_parameters.m_max_area_change = options.m_max_area_change;
   construction_parameters.m_min_edge_length = options.m_min_edge_length;
   construction_parameters.m_max_edge_length = options.m_max_edge_length;
   construction_parameters.m_merge_proximity_epsilon = options.m_merge_proximity_epsilon;
   construction_parameters.m_allow_topology_changes = options.m_allow_topology_changes;
   construction_parameters.m_perform_improvement = options.m_perform_improvement;
   
   eltopo2d::SurfTrack surface_tracker( vs, es, masses, construction_parameters );

   // perform mesh improvement
   surface_tracker.improve_mesh();
   
   // do merging
   surface_tracker.topology_changes();
   
   // =================================================================================

   //
   // data wrangling
   //
   
   *out_num_vertices = surface_tracker.m_positions.size();
   *out_vertex_locations = (double*) malloc( 2 * (*out_num_vertices) * sizeof(double) );
   *out_vertex_masses = (double*) malloc( (*out_num_vertices) * sizeof(double) );
   
   for (  int i = 0; i < (*out_num_vertices); ++i )
   {
      (*out_vertex_locations)[2*i] = surface_tracker.m_positions[i][0];
      (*out_vertex_locations)[2*i + 1] = surface_tracker.m_positions[i][1];      
      (*out_vertex_masses)[i] = surface_tracker.m_masses[i];
   }
   
   *out_num_edges = surface_tracker.m_mesh.edges.size();
   *out_edges = (int*) malloc( 2 * (*out_num_edges) * sizeof(int) );
   for (  int i = 0; i < (*out_num_edges); ++i )
   {
      (*out_edges)[2*i] = surface_tracker.m_mesh.edges[i][0];
      (*out_edges)[2*i + 1] = surface_tracker.m_mesh.edges[i][1]; 
   }
   
}


void el_topo_free_static_operations_results( double* out_vertex_locations, double* out_vertex_masses, int* out_edges )
{
   free( out_vertex_locations );
   free( out_vertex_masses );
   free( out_edges );
}
   
void el_topo_integrate( const int num_vertices, 
                const double *in_vertex_locations, 
                const double *in_vertex_new_locations, 
                const double* in_vertex_masses,
                const int num_edges, 
                const int *edges, 
                const ElTopoGeneralOptions& general_options,
                const ElTopoIntegrationOptions& options,
                double **out_vertex_locations )
{
   //
   // data wrangling
   //
   
   std::vector<Vec2d> vs;
   std::vector<double> masses;
   for ( int i = 0; i < num_vertices; ++i )
   {
      vs.push_back( Vec2d( in_vertex_locations[2*i], in_vertex_locations[2*i + 1] ) );
      
      masses.push_back( in_vertex_masses[i] );
   }
   
   std::vector<Vec2ui> es;
   for (  int i = 0; i < num_edges; ++i )
   {
      es.push_back( Vec2ui( edges[2*i], edges[2*i + 1] ) );
   }
   
   
   // =================================================================================
   
   //
   // do the integration
   //
   
   // build a DynamicSurface
   eltopo2d::DynamicSurface dynamic_surface( vs, es, masses, options.m_proximity_epsilon, general_options.m_collision_safety, general_options.m_verbose );
   
   // set velocities
   dynamic_surface.m_velocities.resize( num_vertices );
   for (  int i = 0; i < num_vertices; ++i )
   {
      dynamic_surface.m_newpositions[i] = Vec2d( in_vertex_new_locations[2*i], in_vertex_new_locations[2*i + 1] );
   }
   
   // advance by dt
   dynamic_surface.integrate( options.m_dt );
   
   
   // =================================================================================
   
   
   //
   // data wrangling
   //
   
   *out_vertex_locations = (double*) malloc( 2 * num_vertices * sizeof(double) );
   for (  int i = 0; i < num_vertices; ++i )
   {
      (*out_vertex_locations)[2*i] = dynamic_surface.m_positions[i][0];
      (*out_vertex_locations)[2*i + 1] = dynamic_surface.m_positions[i][1];      
   }
      
}

void el_topo_free_integrate_results( double* out_vertex_locations )
{
   free( out_vertex_locations );
}

   