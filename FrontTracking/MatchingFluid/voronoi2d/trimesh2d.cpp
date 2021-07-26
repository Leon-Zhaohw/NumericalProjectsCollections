#include "trimesh2d.h"

#include <cfloat>
#include "wallclocktime.h"

TriMesh2D::TriMesh2D() {}

void TriMesh2D::setup_mesh(const std::vector<Vec2f>& v, const std::vector<Vec3ui>& t) {
   vertices = v;
   tris = t;

   std::cout << "Setting up mesh data structures." << std::endl;
   
   compute_tris();
   compute_edges();
   compute_nbrs();
   compute_accelerator();
}


Vec2f circumcentre(Vec2f v0, Vec2f v1, Vec2f v2) {
   //Due to Dave Watson, from Eppstein's geometry junkyard

   float D = (v0[0] - v2[0]) * (v1[1] - v2[1]) - (v1[0] - v2[0]) * (v0[1] - v2[1]);
   
   float p_0 = (((v0[0] - v2[0]) * (v0[0] + v2[0]) + (v0[1] - v2[1]) * (v0[1] + v2[1])) / 2 * (v1[1] - v2[1]) 
    -  ((v1[0] - v2[0]) * (v1[0] + v2[0]) + (v1[1] - v2[1]) * (v1[1] + v2[1])) / 2 * (v0[1] - v2[1])) 
    / D;

   float p_1 = (((v1[0] - v2[0]) * (v1[0] + v2[0]) + (v1[1] - v2[1]) * (v1[1] + v2[1])) / 2 * (v0[0] - v2[0])
    -  ((v0[0] - v2[0]) * (v0[0] + v2[0]) + (v0[1] - v2[1]) * (v0[1] + v2[1])) / 2 * (v1[0] - v2[0]))
    / D;

   return Vec2f(p_0,p_1);

}

Vec2f barycentre(Vec2f v0, Vec2f v1, Vec2f v2) {
   return (v0 + v1 + v2) / 3.0f;
}

void TriMesh2D::compute_tris() {
   
   tri_circumcentres.clear();
   tri_barycentres.clear();
   
   for(unsigned int i = 0; i < tris.size(); ++i) {
      //compute area
      Vec2f v0 = vertices[tris[i][0]];
      Vec2f v1 = vertices[tris[i][1]];
      Vec2f v2 = vertices[tris[i][2]];
      
      float area = 0.5f*cross(v1-v0,v2-v0);
      
      //flip vertex order if necessary for consistency
      if(area < 0) {
         int t = tris[i][1];
         tris[i][1] = tris[i][2];
         tris[i][2] = t;
         area = -area;
      }
      tri_areas.push_back(area);
      tri_circumcentres.push_back(circumcentre(v0,v1,v2));
      tri_barycentres.push_back(barycentre(v0,v1,v2));      
   }
   
}

bool vec2ui_cmp(Vec2ui e0, Vec2ui e1) {
   if(e0[0] == e1[0])
      return e0[1] < e1[1];
   else 
      return e0[0] < e1[0];
}

struct EdgeWithTris {
   Vec2ui edge;
   Vec2i tris;
   EdgeWithTris(int e0, int e1, int t):edge(e0,e1), tris(t,-1) {};
   
   bool operator==(const EdgeWithTris& other) const {
      return edge == other.edge && tris == other.tris;
   }
};

bool ewt_cmp(EdgeWithTris e0, EdgeWithTris e1) {
   if(e0.edge[0] == e1.edge[0])
      return e0.edge[1] < e1.edge[1];
   else 
      return e0.edge[0] < e1.edge[0];
}

void TriMesh2D::compute_edges() {
   
   
   std::vector<EdgeWithTris> edges_temp;
   for(unsigned int f = 0; f < tris.size(); ++f) {
      //compute edges
      EdgeWithTris e0(tris[f][0], tris[f][1], f);
      EdgeWithTris e1(tris[f][1], tris[f][2], f);
      EdgeWithTris e2(tris[f][2], tris[f][0], f);

      //ensure consistency, by sorting indices
      if(e0.edge[1]<e0.edge[0])
         swap(e0.edge[1],e0.edge[0]);
      if(e1.edge[1]<e1.edge[0])
         swap(e1.edge[1],e1.edge[0]);
      if(e2.edge[1]<e2.edge[0])
         swap(e2.edge[1],e2.edge[0]);

      //add them
      edges_temp.push_back(e0);
      edges_temp.push_back(e1);
      edges_temp.push_back(e2);
   
   }
   
   sort(edges_temp.begin(), edges_temp.end(), ewt_cmp);
   
   //go through pairs, and copy tri info over.
   unsigned int pos = 0;
   while(pos < edges_temp.size()-1) {
      if(edges_temp[pos].edge == edges_temp[pos+1].edge) {
         //determine the triangle pair
         Vec2i tri(edges_temp[pos].tris[0], edges_temp[pos+1].tris[0]);
         //copy into both         
         edges_temp[pos].tris = tri;
         edges_temp[pos+1].tris = tri;
         ++pos; //can skip ahead one more, since we've handled this pair.
      }
      ++pos;
   }

   edges_temp.erase( unique( edges_temp.begin(), edges_temp.end() ), edges_temp.end() );
   //we should now have a collection of winged edges, in some sorted order
   
   //copy the edges into the edge list...
   
   edges.resize(edges_temp.size());
   for(unsigned int i = 0; i < edges_temp.size(); ++i) {
      edges[i] = edges_temp[i].edge;
   }

   //now, construct the tri_to_edge_map
   tri_to_edge_map.resize(tris.size());
   std::vector<int> tri_edges_so_far(tris.size(), 0);
   for(unsigned int i = 0; i < edges_temp.size(); ++i) {
      EdgeWithTris cur = edges_temp[i];
      int tri0 = cur.tris[0];
      tri_to_edge_map[tri0][ tri_edges_so_far[tri0] ] = i;
      ++tri_edges_so_far[tri0];
      
      int tri1 = cur.tris[1];
      if(tri1 != -1) { //if there is a second tri bordering this edge...
         tri_to_edge_map[tri1][ tri_edges_so_far[tri1] ] = i;
         ++tri_edges_so_far[tri1];
      }
   }

   for(unsigned int i = 0; i < tri_to_edge_map.size(); ++i) {
      //reorder the edges in the map... the code apparently depends on the
      //fact that the ordering of the edges lines up with the ordering of the vertices
      Vec3ui new_order;

      Vec2ui e0(tris[i][0], tris[i][1]);
      if(e0[0] > e0[1]) swap(e0[0],e0[1]);
      new_order[0] = (e0 == edges[tri_to_edge_map[i][0]]? 
                           tri_to_edge_map[i][0] :
                           (e0 == edges[tri_to_edge_map[i][1]]? 
                                 tri_to_edge_map[i][1] : 
                                 tri_to_edge_map[i][2])
                      );

      Vec2ui e1(tris[i][1], tris[i][2]);
      if(e1[0] > e1[1]) swap(e1[0],e1[1]);
      new_order[1] = (e1 == edges[tri_to_edge_map[i][0]]? 
                           tri_to_edge_map[i][0] :
                           (e1 == edges[tri_to_edge_map[i][1]]? 
                                 tri_to_edge_map[i][1] : 
                                 tri_to_edge_map[i][2])
                      );

      Vec2ui e2(tris[i][2], tris[i][0]);
      if(e2[0] > e2[1]) swap(e2[0],e2[1]);
      new_order[2] = (e2 == edges[tri_to_edge_map[i][0]]? 
                           tri_to_edge_map[i][0] :
                           (e2 == edges[tri_to_edge_map[i][1]]? 
                                 tri_to_edge_map[i][1] : 
                                 tri_to_edge_map[i][2])
                      );
      tri_to_edge_map[i] = new_order;
   }  

   //compute midpoints and normals
   for(unsigned int e = 0; e < edges.size(); ++e) {
    
      Vec2f v0 = vertices[edges[e][0]];
      Vec2f v1 = vertices[edges[e][1]];
      Vec2f edge_vec = v1-v0;
      Vec2f normal = perp(edge_vec);
      float length = mag(edge_vec);
      edge_lengths.push_back(length);
      normal /= length;
      edge_normals.push_back(normal);
      edge_midpoints.push_back(0.5f*(v0+v1));
   }

}


void TriMesh2D::compute_accelerator() {
   Vec2f minPt = vertices[0], maxPt = vertices[0];
   for(unsigned int i = 0; i < vertices.size(); ++i) {
      minPt = min_union(vertices[i], minPt);
      maxPt = max_union(vertices[i], maxPt);
   }
   
   float avglen = 0;
   for(unsigned int i = 0; i < edge_lengths.size(); ++i) { 
      avglen += edge_lengths[i];
   }
   avglen /= (float)edge_lengths.size();
   accel_dx = avglen;
   inv_accel_dx = 1.0f/accel_dx;
   accel_ni = (int)((maxPt[0] - minPt[0]) / accel_dx)+4;
   accel_nj = (int)((maxPt[1] - minPt[1]) / accel_dx)+4;
   accel_origin = minPt - Vec2f(accel_dx,accel_dx);

   accel_grid.resize(accel_ni*accel_nj);
   vert_accel_grid.resize(accel_ni*accel_nj);

   for(unsigned int t = 0; t < tris.size(); ++t) {
      //for each tri, compute its bound box
      Vec2f a = vertices[tris[t][0]];
      Vec2f b = vertices[tris[t][1]];
      Vec2f c = vertices[tris[t][2]];
      Vec2f minPt, maxPt;
      get_bound_box(a,b,c,minPt,maxPt);
      
      //use the bound box to determine a range of cells to look at
      int i_lower = (int)((minPt[0] - accel_origin[0])*inv_accel_dx-1);
      int j_lower = (int)((minPt[1] - accel_origin[1])*inv_accel_dx-1);
      int i_upper = (int)((maxPt[0] - accel_origin[0])*inv_accel_dx+1);
      int j_upper = (int)((maxPt[1] - accel_origin[1])*inv_accel_dx+1);

      //visit those cells to see if there are intersections
      for(int j = j_lower; j < j_upper; ++j) {
         for(int i = i_lower; i < i_upper; ++i) {
            int grid_ind = i+j*accel_ni;
            accel_grid[grid_ind].push_back(t);
         }
      }
   }


   for(unsigned int v = 0; v < vertices.size(); ++v) {
      //for each vertex (voronoi seed), compute bound box of adjacent tri circumcentres
      //(ie. of the voronoi vertices)

      //Get bounding box for the voronoi cell
      Vec2f minPt(1e30f,1e30f), 
         maxPt = Vec2f(-1e30f,-1e30f);
      for(unsigned int t = 0; t < vert_to_tri_map[v].size(); ++t) {
         int tri_ind = vert_to_tri_map[v][t];
         minPt = min_union(minPt, tri_circumcentres[tri_ind]);
         maxPt = max_union(maxPt, tri_circumcentres[tri_ind]);
      }
      
      //use the bound box to determine the range of cells 
      int i_lower = (int)((minPt[0] - accel_origin[0])*inv_accel_dx-1);
      int j_lower = (int)((minPt[1] - accel_origin[1])*inv_accel_dx-1);
      int i_upper = (int)((maxPt[0] - accel_origin[0])*inv_accel_dx+1);
      int j_upper = (int)((maxPt[1] - accel_origin[1])*inv_accel_dx+1);

      //store those cells as intersections
      for(int j = j_lower; j < j_upper; ++j) {
         for(int i = i_lower; i < i_upper; ++i) {
            int grid_ind = i+j*accel_ni;
            if(i < 0 || i >= accel_ni || j < 0 || j > accel_nj) continue;
            vert_accel_grid[grid_ind].push_back(v);
         }
      }
   }

}

void TriMesh2D::compute_nbrs() {
   
   //tri_edge_counts.resize(tris.size());
  /*    
   //old brute force code; this is now done smarter in compute_edges
   std::vector<Vec3ui> tri_to_edge_old(tris.size());
   tri_to_edge_old.resize(tris.size());
   for(unsigned int i = 0; i < tris.size(); ++i) {
      Vec2ui e0(tris[i][0], tris[i][1]);
      if(e0[0] > e0[1]) swap(e0[0],e0[1]);
      int ind = 0;
      while(edges[ind] != e0) ++ind;
      tri_to_edge_old[i][0] = ind;
      
      Vec2ui e1(tris[i][1], tris[i][2]);
      if(e1[0] > e1[1]) swap(e1[0],e1[1]);
      ind = 0;
      while(edges[ind] != e1) ++ind;
      tri_to_edge_old[i][1] = ind;
      
      Vec2ui e2(tris[i][2], tris[i][0]);
      if(e2[0] > e2[1]) swap(e2[0],e2[1]);
      ind = 0;
      while(edges[ind] != e2) ++ind;
      tri_to_edge_old[i][2] = ind;
   }*/

   //create edge_to_tri_map
   std::vector<int> tri_count(edges.size());
   edge_to_tri_map.resize(edges.size(), Vec2i(-1,-1));
   for(unsigned int i = 0; i < tri_to_edge_map.size(); ++i) {
      int edge0 = tri_to_edge_map[i][0];
      edge_to_tri_map[edge0][tri_count[edge0]++] = i;

      int edge1 = tri_to_edge_map[i][1];
      edge_to_tri_map[edge1][tri_count[edge1]++] = i;

      int edge2 = tri_to_edge_map[i][2];
      edge_to_tri_map[edge2][tri_count[edge2]++] = i;
   }

   //determine the signs of the faces relative to the relevant triangle
   tri_edge_signs.resize(tris.size());
   for(unsigned int i = 0; i < tris.size(); ++i) {
      Vec3ui tri_vert_ids = tris[i];

      //get the edges for this triangle
      Vec3ui edge_ids = tri_to_edge_map[i]; 
      for(unsigned int e = 0; e < 3; ++e) {
         //pick an edge      
         int edge_id = edge_ids[e];
         Vec2ui edge_vert_ids = edges[edge_id];
         
         //find the index of the first vertex in the edge
         int ind = 0;
         while(tri_vert_ids[ind] != edge_vert_ids[0]) ++ind;

         //check if the next vertex in the tri is the next vertex in the edge
         //if so, they are orientated the same way, no need to flip the normal
         tri_edge_signs[i][e] = (tri_vert_ids[(ind+1)%3] == edge_vert_ids[1] ? 1 : -1);
      }
   }

   //swap the order of the edge_to_tri_map so that it corresponds with the
   //direction of the edge normal.
   //ie. the triangle that the normal points at should be second
   //this corresponds to the first tri having a positive sign for the edge, the second a negative sign
   for(unsigned int i = 0; i < edge_to_tri_map.size(); ++i) {
      Vec2i tris = edge_to_tri_map[i];
      int edge_ind = 0;
      while(tri_to_edge_map[tris[0]][edge_ind] != i) ++edge_ind;

      if(tri_edge_signs[tris[0]][edge_ind] == -1)
         swap(edge_to_tri_map[i][0], edge_to_tri_map[i][1]);
   }

   //create vert_to_tri_map (ideally each should be ordered, but I don't think it matters)
   vert_to_tri_map.resize(vertices.size());
   for(unsigned int i = 0; i < tris.size(); ++i) {
      Vec3ui vert_ids = tris[i];
      vert_to_tri_map[vert_ids[0]].push_back(i);
      vert_to_tri_map[vert_ids[1]].push_back(i);
      vert_to_tri_map[vert_ids[2]].push_back(i);
   }

   //create vert_to_edge_map
   vert_to_edge_map.resize(vertices.size());
   for(unsigned int i = 0; i < edges.size(); ++i) {
      //use the edge list (ie. edge_to_vert) to build the
      //vert to edge list
      Vec2ui edge = edges[i];
      vert_to_edge_map[edge[0]].push_back(i);
      vert_to_edge_map[edge[1]].push_back(i);
   }

   //compute dual edges
   dual_edge_lengths.resize(edges.size());
   dual_edge_midpoints.resize(edges.size());
   dual_edge_normals.resize(edges.size());
   for(unsigned int i = 0; i < edges.size(); ++i) {
      Vec2i triCentreIds = edge_to_tri_map[i];
      if(triCentreIds[0] == -1 || triCentreIds[1] == -1) {
         dual_edge_lengths[i] = 0;
         dual_edge_midpoints[i] = edge_midpoints[i];
         dual_edge_normals[i] = perp(edge_normals[i]);
         continue;
      }

      Vec2f c0 = tri_circumcentres[triCentreIds[0]];
      Vec2f c1 = tri_circumcentres[triCentreIds[1]];
      dual_edge_lengths[i] = dist(c0,c1);
      dual_edge_midpoints[i] = 0.5f*(c0+c1);
      dual_edge_normals[i] = perp(edge_normals[i]); //safer, because sometimes c0 and c1 are co-located, so no normal to be found.
   }

   //set the signs of dual edges
   dual_edge_signs.resize(vertices.size());
   for(unsigned int v = 0; v < vertices.size(); ++v) {
      
      dual_edge_signs[v].resize(vert_to_edge_map[v].size());
      for(unsigned int e = 0; e < vert_to_edge_map[v].size();++e) {
         int edgeId = vert_to_edge_map[v][e];
         Vec2ui edgeEnds = edges[edgeId];
         //positive one implies the edge is oriented from the voronoi centre outwards
         dual_edge_signs[v][e] = edgeEnds[0] == v? 1 : -1;
      }
   }


}

int TriMesh2D::get_containing_tri(const Vec2f& point) const {
   
   Vec2f cell = (point - accel_origin)*inv_accel_dx;
   unsigned int offset = (int)cell[0] + accel_ni * (int)cell[1];
   if(offset < accel_grid.size()) {
      for(unsigned int i = 0; i < accel_grid[offset].size(); ++i)
         if(point_in_tri(point, accel_grid[offset][i]))
            return accel_grid[offset][i];
   }
   
   return -1;
   
}



void TriMesh2D::get_overlapping_cells( const Vec2f& centre, float radius, std::vector<int>& cells ) const 
{
   float min_x = centre[0] - radius;
   int min_i = (int)( (min_x - accel_origin[0])*inv_accel_dx );
   min_i = max( 0, min_i );
   
   float max_x = centre[0] + radius;
   int max_i = (int)( (max_x - accel_origin[0])*inv_accel_dx );
   max_i = min( accel_ni, max_i );
   
   float min_y = centre[1] - radius;
   int min_j = (int)( (min_y - accel_origin[1])*inv_accel_dx );
   min_j = max( 0, min_j );
   
   float max_y = centre[1] + radius;
   int max_j = (int) ((max_y - accel_origin[1])*inv_accel_dx );
   max_j = min( accel_nj, max_j );
   
//   std::cout << "i: " << min_i << ", " << max_i << std::endl;
//   std::cout << "j: " << min_j << ", " << max_j << std::endl;
   
   for ( int i = min_i; i <= max_i; ++i )
   {
      for ( int j = min_j; j <= max_j; ++j )
      {
         cells.push_back( i + accel_ni * j );
      }
   }
   
}

int TriMesh2D::get_nearest_vertex(const Vec2f& point) const 
{
      
   /*  
   //the brute force version, for verifying.
   int check_min_index = -1;
   float check_min_distance = 1e30f;
   for ( unsigned int i = 0; i < vertices.size(); ++i ) 
   {
      if ( vert_to_tri_map[i].size() > 0 )
      {
         float dist2val = dist2( vertices[i], point );
         if(dist2val < check_min_distance) {
            check_min_distance = dist2val;
            check_min_index = i;
         }
      }
   }*/
   
   
   int i_ind = (int)((point[0] - accel_origin[0])*inv_accel_dx);
   int j_ind = (int)((point[1] - accel_origin[1])*inv_accel_dx);
   
   i_ind = clamp(i_ind,0, accel_ni-1);
   j_ind = clamp(j_ind,0, accel_nj-1);
   int accel_ind = i_ind + j_ind*accel_ni;

   int check_min_index = -1;
   float check_min_distance = 1e30f;

   for(unsigned int p = 0; p < vert_accel_grid[accel_ind].size(); ++p) {
      int v_ind = vert_accel_grid[accel_ind][p];
      if ( vert_to_tri_map[v_ind].size() > 0 )
      {
         float dist2val = dist2( vertices[v_ind], point );
         if(dist2val < check_min_distance) {
            check_min_distance = dist2val;
            check_min_index = v_ind;
         }
      }
   }
   
  
   return check_min_index;   
   
}


////walking point location
//int TriMesh2D::get_containing_tri(const Vec2f& point, int start_tri) const {
//   //int answer = get_containing_tri(point);
//
//   int tri = start_tri;
//   while(tri != -1 && !point_in_tri(point, tri)) {
//      //determine neighbours
//      Vec3ui nbr_edges = tri_to_edge_map[tri];
//      
//      //look at each edge, figure out which side the point is on
//      for(int i = 0; i < 3; ++i) {     
//         int edge = nbr_edges[i];
//         Vec2ui vert_indices = edges[edge];
//         int sign = tri_edge_signs[tri][i];
//         float value = 0;
//         Vec2f v0 = vertices[vert_indices[0]];
//         Vec2f v1 = vertices[vert_indices[1]];
//         if(sign == 1)
//            value = cross(v1 - v0, point - v0);
//         else if(sign == -1)
//            value = cross(v0 - v1, point - v1);
//         if(value <= 0) {
//            int choice0 = edge_to_tri_map[edge][0];
//            int choice1 = edge_to_tri_map[edge][1];
//            if(choice0 == tri)
//               tri = choice1;
//            else if(choice1 == tri)
//               tri = choice0;
//         }
//         if(tri == -1)
//            break;
//      }
//
//   }
//   
//   return tri;
//}



bool TriMesh2D::point_in_tri(const Vec2f& point, int tri_id) const {
   Vec2f p0 = vertices[tris[tri_id][0]];
   Vec2f p1 = vertices[tris[tri_id][1]];
   Vec2f p2 = vertices[tris[tri_id][2]];

   Vec2f v0 = p2 - p0;
   Vec2f v1 = p1 - p0;
   Vec2f v2 = point - p0;

   // Compute dot products
   double dot00 = dot(v0, v0);
   double dot01 = dot(v0, v1);
   double dot02 = dot(v0, v2);
   double dot11 = dot(v1, v1);
   double dot12 = dot(v1, v2);

   // Compute barycentric coordinates
   double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
   double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
   double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

   // Check if point is in triangle (TODO Make this more consistent)
   return (u >= -1e-5) && (v >= -1e-5) && (u + v <= 1+2e-5);
}

/*
void TriMesh2D::get_interpolation_weights(const Vec2f& point, int& closestPoint, std::vector<float>& weights, int tri) const {
   
   //TODO Use acceleration structure
   float min_dist = FLT_MAX;
   int i = 0;

   if(tri == -1 || !point_in_tri(point, tri)) {
      //determine closest point, that's tells us which voronoi cell to use
      for(unsigned int j = 0; j < vertices.size(); ++j) {
         float d = dist(vertices[j], point);
         if(d < min_dist && vert_to_tri_map[j].size() > 0) {
            min_dist = d;
            i = j;
         }
      }
   }
   else { //just check the three triangle verts
      Vec3ui verts_of_tri = tris[tri];
      for(unsigned int j = 0; j < 3; ++j) {
         float d = dist(vertices[verts_of_tri[j]], point);
         if(d < min_dist) {
            min_dist = d;
            i = verts_of_tri[j];
         }
      }
   }
   
   Vec2f vert_location = point;

   //need to know all the tris surrounding a vertex
   const std::vector<int>& surrounding_tris = vert_to_tri_map[i];
   
   float weight_sum = 0;
   Vec2f velocity(0,0);
   weights.clear();
   weights.resize(surrounding_tris.size());
   for(unsigned int j = 0; j < surrounding_tris.size(); ++j) {
      Vec3ui tri_edge_ids = tri_to_edge_map[surrounding_tris[j]];

      //compute the vector joining the interpolation point to the tri's circumcentre
      Vec2f centre = tri_circumcentres[surrounding_tris[j]];
      Vec2f centre_to_vert = centre - vert_location;
      
      //Test this distance: if less than some threshold, just return the tri_centre velocity, since otherwise a denominator
      //will be zero, and it'll break things.
      //if(mag(centre_to_vert) < 1e-10) {
      //   return velocity_data[surrounding_tris[j]];
      //}

      //determine an edge containing this vertex, and compute denominator terms for it
      int edge_ind0 = 0;
      while(edges[tri_edge_ids[edge_ind0]][0] != i && edges[tri_edge_ids[edge_ind0]][1] != i) ++edge_ind0;
      int edge0_id = tri_edge_ids[edge_ind0];
      Vec2f edge0_v0 = vertices[ edges[edge0_id][0] ];
      Vec2f edge0_v1 = vertices[ edges[edge0_id][1] ];
      float den0 = dot(edge0_v1 - edge0_v0, centre_to_vert) * (i == edges[edge0_id][0]?1.0f:-1.0f);

      //ditto for other edge
      int edge_ind1 = edge_ind0+1;
      while(edges[tri_edge_ids[edge_ind1]][0] != i && edges[tri_edge_ids[edge_ind1]][1] != i) ++edge_ind1;
      int edge1_id = tri_edge_ids[edge_ind1];
      Vec2f edge1_v0 = vertices[ edges[edge1_id][0] ];
      Vec2f edge1_v1 = vertices[ edges[edge1_id][1] ];
      float den1 = dot(edge1_v1 - edge1_v0, centre_to_vert) * (i == edges[edge1_id][0]?1.0f:-1.0f);


      //compute the final weight
      double denTerm = den0*den1;
      if(fabs(denTerm) < 1e-10)
         denTerm = 1e-10;
      float weight = 2.0f * tri_areas[j] / (float)denTerm;
      weight_sum += weight;
   
      weights[j] = weight;
   }

   for(unsigned int j = 0; j < surrounding_tris.size(); ++j)
      weights[j] /= weight_sum;

   closestPoint = i;
   
}
*/

void TriMesh2D::get_barycentric_weights_tri(int containing_tri, const Vec2f& pt, Vec3f& weights) const {
   Vec2f p0 = vertices[tris[containing_tri][0]];
   Vec2f p1 = vertices[tris[containing_tri][1]];
   Vec2f p2 = vertices[tris[containing_tri][2]];

   // Compute twice area of triangle ABC
   float AreaABC = cross(p1-p0,p2-p0);
   float inv = 1/AreaABC;
   
   // Compute a
   float AreaPBC = cross(p1-pt,p2-pt);
   float a = AreaPBC * inv;

   // Compute b
   float AreaPCA = cross(p2-pt,p0-pt);
   float b = AreaPCA * inv;

   // Compute c
   float c = 1.0f - a - b;

   weights = Vec3f(a,b,c);
}

void TriMesh2D::get_barycentric_weights_voronoi(int closest_vertex, const Vec2f& point, std::vector<float>& weights) const{
   
   unsigned int i = closest_vertex;
   
   Vec2f vert_location = point;

   //need to know all the tris surrounding a vertex
   const std::vector<int>& surrounding_tris = vert_to_tri_map[i];
   
   float weight_sum = 0;
   Vec2f velocity(0,0);
   weights.clear();
   weights.resize(surrounding_tris.size());
   for(unsigned int j = 0; j < surrounding_tris.size(); ++j) {
      Vec3ui tri_edge_ids = tri_to_edge_map[surrounding_tris[j]];

      //compute the vector joining the interpolation point to the tri's circumcentre
      Vec2f centre = tri_circumcentres[surrounding_tris[j]];
      Vec2f centre_to_vert = centre - vert_location;
      
      //Test this distance: if less than some threshold, just return the tri_centre velocity, since otherwise a denominator
      //will be zero, and it'll break things.
      //if(mag(centre_to_vert) < 1e-10) {
      //   return velocity_data[surrounding_tris[j]];
      //}

      //determine an edge containing this vertex, and compute denominator terms for it
      int edge_ind0 = 0;
      while(edges[tri_edge_ids[edge_ind0]][0] != i && edges[tri_edge_ids[edge_ind0]][1] != i) ++edge_ind0;
      int edge0_id = tri_edge_ids[edge_ind0];
      Vec2f edge0_v0 = vertices[ edges[edge0_id][0] ];
      Vec2f edge0_v1 = vertices[ edges[edge0_id][1] ];
      float den0 = dot(edge0_v1 - edge0_v0, centre_to_vert) * (i == edges[edge0_id][0]?1.0f:-1.0f);

      //ditto for other edge
      int edge_ind1 = edge_ind0+1;
      while(edges[tri_edge_ids[edge_ind1]][0] != i && edges[tri_edge_ids[edge_ind1]][1] != i) ++edge_ind1;
      int edge1_id = tri_edge_ids[edge_ind1];
      Vec2f edge1_v0 = vertices[ edges[edge1_id][0] ];
      Vec2f edge1_v1 = vertices[ edges[edge1_id][1] ];
      float den1 = dot(edge1_v1 - edge1_v0, centre_to_vert) * (i == edges[edge1_id][0]?1.0f:-1.0f);


      //compute the final weight
      double denTerm = den0*den1;
      if(fabs(denTerm) < 1e-10)
         denTerm = 1e-10;
      float weight = 2.0f * tri_areas[j] / (float)denTerm;
      weight_sum += weight;
   
      weights[j] = weight;
   }

   for(unsigned int j = 0; j < surrounding_tris.size(); ++j)
      weights[j] /= weight_sum;

}



int TriMesh2D::get_edge_index( const Vec2ui& edge_query ) const
{
   const std::vector<int>& incident_edges_a = vert_to_edge_map[ edge_query[0] ];
   const std::vector<int>& incident_edges_b = vert_to_edge_map[ edge_query[1] ];
   
   for ( unsigned int i = 0; i < incident_edges_a.size(); ++i )
   {
      for ( unsigned int j = 0; j < incident_edges_b.size(); ++j )
      {
         if ( incident_edges_a[i] == incident_edges_b[j] )
         {
            return incident_edges_a[i];
         }
      }      
   }
   
   return -1;
   
}


int TriMesh2D::get_shared_edge_index( unsigned int tri_a, unsigned int tri_b ) const
{
   const Vec3ui& incident_edges_a = tri_to_edge_map[tri_a];
   const Vec3ui& incident_edges_b = tri_to_edge_map[tri_b];
   
   for ( unsigned int i = 0; i < 3; ++i ) for ( unsigned int j = 0; j < 3; ++j )
   {
      if ( incident_edges_a[i] == incident_edges_b[j] )
      {
         return incident_edges_a[i];
      }
   }
   
   return -1;
   
}



