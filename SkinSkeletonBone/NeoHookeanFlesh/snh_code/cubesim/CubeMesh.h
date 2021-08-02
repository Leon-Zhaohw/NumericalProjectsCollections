///////////////////////////////////////////////////////////////////////////////////////////////////
// I. LICENSE CONDITIONS
//
// Copyright (c) 2018 by Disney-Pixar
//
// Permission is hereby granted to use this software solely for non-commercial applications
// and purposes including academic or industrial research, evaluation and not-for-profit media
// production. All other rights are retained by Pixar. For use for or in connection with
// commercial applications and purposes, including without limitation in or in connection with
// software products offered for sale or for-profit media production, please contact Pixar at
// tech-licensing@pixar.com.
//
// THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
// NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PIXAR OR ITS AFFILIATES BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef CUBE_MESH_H
#define CUBE_MESH_H

#include "Cube.h"
#include "LocalHessianLocation.h"
#include "Settings.h"

#include <vector>

namespace CubeSim
{

struct InitialMesh;
class Material;

struct DynamicState final
{
    // Deformed positions of the hex lattice
    std::vector<Vector3> vertices;

    // TODO: Break this out of here
    // Temporary storage caches
    mutable VectorX elementEnergyCache; // Used in ComputeEnergy
    mutable std::vector<Eigen::Matrix<Scalar,24,1>> elementForceCache; // Used in ComputeForces
    mutable std::vector<Eigen::Matrix<Scalar,24,24>> localKs; // Used in _ComputeElasticStiffnessMatrix
};

class CubeMesh final
{

public:

    CubeMesh();

    explicit CubeMesh(const InitialMesh& initial_mesh);

    void CreateInitialDynamicState(const InitialMesh& initial_mesh, DynamicState& dynamicState);

    void CreateInitialSparseMatrix(SparseMatrixX& K);

    CubeMesh(CubeMesh&&) = default;
    CubeMesh& operator=(CubeMesh&&) = default;

    CubeMesh(const CubeMesh&) = delete;
    CubeMesh& operator=(const CubeMesh&) = delete;

    // Integrates the strain energy density over the entire mesh.
    Scalar ComputeEnergy(const DynamicState& dynamicState, const Material& material) const;

    // Computes the global force vector.
    void ComputeForces(const DynamicState& dynamicState, const Material& material, VectorX& f) const;

    // Computes the global sitffness matrix.
    void ComputeStiffnessMatrixSparse(const Material& material, DynamicState& dynamicState, SparseMatrixX& K) const;

    // True if a vertex is kinematically scripted.
    bool IsVertexFixed(const int vrtIdx) const;

    // Return the number of simulatd degrees of freedom in the mesh.
    int GetNumSimulatedDOFs() const;

    // Scatter the displacements to the vertices.
    void ScatterDisplacements(const VectorX& u, DynamicState& dynamicState) const;

    // Computes the displacement at each vertex.
    void ComputeDisplacements(const DynamicState& dynamicState, VectorX& u) const;

private:

    void _ComputeInternalLocalStiffnessMatrices(const Material& material, DynamicState& dynamicState) const;

    // Rest positions of the hex lattice
    std::vector<Vector3> _restVertices;
    // Hexahedral elements
    std::vector<Cube> _cubes;
    // For each vertex, maps from the global index to the start of the vertex in the simulated DoFs.
    // -1 for kinematic vertices. Multiple of 3 for simulated vertices.
    std::vector<int> _vertexIdxToSimulatedDof;
    // Base indices of simulated cubes in the stiffness matrix
    std::vector<long> _cubeBaseIndices;
    // The number of entries in the stiffness matrix (as triplets) from the interal forces
    long _internalStiffnessTripletCount;
    // The number of degrees of freedom in the system (3 * simulated_vertex_count)
    int _DOFs;
    // For accelerating the internal elastic stiffness matrix construction
    std::vector<std::vector<LocalHessianLocation>> _compressedIndexToLocalHessianMap;
    // Map from off-diagonal entries to their corresponding entries
    std::vector<std::pair<int,int>> _symmetricEntryMap;

};

}

#endif
