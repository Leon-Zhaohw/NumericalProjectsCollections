#ifndef PRESSURE_VOR2D_H
#define PRESSURE_VOR2D_H

#include "trimesh2d.h"
#include "sparse_matrix.h"

#include "surftrack.h"

std::vector<double> pressure_solve_voronoi(TriMesh2D& mesh, 
                            eltopo2d::SurfTrack& surface,
                            std::vector<float>& face_velocities, 
                            const std::vector<float>& solid_weights,
                            const std::vector<float>& liquid_weights,
                            const std::vector<float>& wall_velocities);


#endif

