
#ifndef ELTOPO2D_H
#define ELTOPO2D_H

#ifdef __cplusplus
extern "C" {
#endif
   
   
struct ElTopoGeneralOptions
{
   bool m_verbose;
   bool m_collision_safety;
};
   
struct ElTopoStaticOperationsOptions
{
   double m_max_area_change;   
   double m_min_edge_length;   
   double m_max_edge_length;   
   double m_merge_proximity_epsilon;
   bool m_perform_improvement;   
   bool m_allow_topology_changes;
};

struct ElTopoIntegrationOptions
{
   double m_proximity_epsilon;
   double m_dt;
};
   
 
void el_topo_static_operations( const int in_num_vertices,                   
                                const double* in_vertex_locations,           
                                const double* in_vertex_masses,                                                
                                const int num_edges,                   
                                const int* in_edges, 
                                const ElTopoGeneralOptions& general_options,       
                                const ElTopoStaticOperationsOptions& options, 
                                int *out_num_vertices,                       
                                double** out_vertex_locations,  
                                double** out_vertex_masses,
                                int *out_num_edges,
                                int** out_edges );
 
void el_topo_free_static_operations_results( double* out_vertex_locations, double* out_vertex_masses, int* out_edges );
   
void el_topo_integrate( const int num_vertices, 
                        const double *in_vertex_locations, 
                        const double *in_vertex_new_locations, 
                        const double* in_vertex_masses,
                        const int num_edges, 
                        const int *edges, 
                        const ElTopoGeneralOptions& general_options,
                        const ElTopoIntegrationOptions& options,
                        double **out_vertex_locations );
   

void el_topo_free_integrate_results( double* out_vertex_locations );
   

#ifdef __cplusplus
}
#endif

#endif
