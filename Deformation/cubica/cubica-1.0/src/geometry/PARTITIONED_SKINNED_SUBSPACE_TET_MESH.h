/*
This file is part of Cubica.
 
Cubica is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Cubica is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Cubica.  If not, see <http://www.gnu.org/licenses/>.
*/
// PARTITIONED_SKINNED_SUBSPACE_TET_MESH.h: interface for the PARTITIONED_SKINNED_SUBSPACE_TET_MESH class.
//
//////////////////////////////////////////////////////////////////////

#ifndef PARTITIONED_SKINNED_SUBSPACE_TET_MESH_H
#define PARTITIONED_SKINNED_SUBSPACE_TET_MESH_H

#include <PARTITIONED_SUBSPACE_TET_MESH.h>
#include <SKELETON.h>

class PARTITIONED_SKINNED_SUBSPACE_TET_MESH : public PARTITIONED_SUBSPACE_TET_MESH {

public:
  PARTITIONED_SKINNED_SUBSPACE_TET_MESH(SKELETON* skeleton, const char* filename, MATERIAL** materials, int totalMaterials, int partitions, Real springConst = 2.0, bool simulate = false, bool loadOriginal = true, string partitionPath = string(""));
  virtual ~PARTITIONED_SKINNED_SUBSPACE_TET_MESH();

  // load skinned simulation data snapshot, including the mocap data
  void loadSimulationFrame(string dataPath, int frame, VECTOR& position);
 
  // load skinned Ode data snapshot
  void loadOdeFrame(string dataPath, int frame, VECTOR& position);

  // load skeleton data only -- simulation should take of everything else
  void loadSkeletonFrame(string dataPath, int frame);

  // load ODE skeleton data only -- simulation should take of everything else
  void loadOdeSkeletonFrame(string dataPath, int frame);

  // load skeleton data only -- simulation should take of everything else
  void loadInterpolatedSkeletonFrame(string dataPath, float frame);

  // load a mocap frame into the skeleton, and propagate the bone transforms
  // to the partitions
  void loadMocapFrame(int frame);

  // load a mocap frame into the skeleton, and propagate the bone transforms
  // to the partitions
  void loadInterpolatedMocapFrame(float frame);

  // scatter the constrained node positions to the submeshes -- in the case of a skinned
  // animation, these positions will actually change along with the skeleton
  void scatterConstrainedNodes();

  // undo the translation in a snapshot for a given partition
  void subtractTranslation(int partition, VECTOR& snapshot);

  // undo the rotation in a snapshot for a given partition
  void subtractRotation(int partition, VECTOR& snapshot);

  // hack to determine if this is a PARTITIONED_SKINNED_SUBSPACE_TET_MESH or not
  virtual bool hasSkeleton() { return _skeleton; };

  // set the rest poses of the meshes according to the first ODE frame
  void setOdeRestRigids(string dataPath);

  // scatter the current skeleton bone transforms to the meshes
  void scatterSkeletonTransforms();

protected:
  SKELETON* _skeleton;

  // propagate bone transforms to partitions
  void scatterBoneTransforms();
};

#endif
