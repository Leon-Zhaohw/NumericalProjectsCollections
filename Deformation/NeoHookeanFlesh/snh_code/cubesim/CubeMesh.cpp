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

#include "CubeMesh.h"

#include "InitialMesh.h"
#include "Material.h"

#include <set>

#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace CubeSim
{

static std::vector<Vector3> EigenVertsToSTDVerts(const Matrix3X& verts)
{
    std::vector<Vector3> stdVerts(static_cast<unsigned long>(verts.cols()));
    #pragma omp parallel for
    for (int vrtIdx = 0; vrtIdx < verts.cols(); vrtIdx++)
    {
        stdVerts[static_cast<unsigned long>(vrtIdx)] = verts.col(vrtIdx);
    }
    return stdVerts;
}

static std::vector<Cube> GenerateCubes(const InitialMesh& initialMesh, const std::vector<Vector3>& vertices, const std::vector<Vector3>& restVertices)
{
    assert(int(initialMesh.kinematic.size()) == initialMesh.vertsInit.cols());
    std::vector<Cube> cubes(static_cast<unsigned long>(initialMesh.cubes.cols()));
    #pragma omp parallel for
    for (int cidx = 0; cidx < initialMesh.cubes.cols(); cidx++)
    {
        const std::array<int,8> vidx {
            initialMesh.cubes(0, cidx),
            initialMesh.cubes(1, cidx),
            initialMesh.cubes(2, cidx),
            initialMesh.cubes(3, cidx),
            initialMesh.cubes(4, cidx),
            initialMesh.cubes(5, cidx),
            initialMesh.cubes(6, cidx),
            initialMesh.cubes(7, cidx)
        };
        const bool cubeFixed = std::all_of(vidx.begin(), vidx.end(), [&initialMesh](const int idx){return initialMesh.kinematic[static_cast<unsigned long>(idx)];});
        cubes[static_cast<unsigned long>(cidx)] = Cube(cubeFixed, vidx, vertices, restVertices);
    }
    return cubes;
}

static std::vector<int> ComputeVertexToIndexTable(const std::vector<bool>& vertexIsKinematic)
{
    std::vector<int> vertexIdxToSimulatedDof(vertexIsKinematic.size());
    int curDof = 0;
    for (int vidx = 0; vidx < int(vertexIsKinematic.size()); vidx++)
    {
        if (vertexIsKinematic[static_cast<unsigned long>(vidx)])
        {
            vertexIdxToSimulatedDof[static_cast<unsigned long>(vidx)] = -1;
        }
        else
        {
            vertexIdxToSimulatedDof[static_cast<unsigned long>(vidx)] = curDof;
            curDof += 3;
        }
    }
    return vertexIdxToSimulatedDof;
}

static std::vector<long> ComputeCubeBaseIndices(const std::vector<Cube>& cubes, const std::vector<int>& vertexIdxToSimulatedDof)
{
    std::vector<long> cubeBaseIndices(cubes.size());

    long internalStiffnessTripletCount = 0;
    for (std::vector<Cube>::size_type cubeIdx = 0; cubeIdx < cubes.size(); cubeIdx++)
    {
        if (cubes[cubeIdx].AllCornersFixed())
        {
            cubeBaseIndices[cubeIdx] = -1;
            continue;
        }

        cubeBaseIndices[cubeIdx] = internalStiffnessTripletCount;

        // Map the cube's verts to DoFs in the global stiffness matrix
        std::array<int,8> indices;
        for (int y = 0; y < 8; y++)
        {
            const int vidx = cubes[cubeIdx].VertexIndex(y);
            indices[static_cast<unsigned long>(y)] = vertexIdxToSimulatedDof[static_cast<unsigned long>(vidx)];
        }

        // Add the local stiffness matrix into the global stiffness matrix
        for (int j = 0; j < 8; j++)
        {
            for (int i = 0; i < 8; i++)
            {
                if (indices[static_cast<unsigned long>(i)] == -1 || indices[static_cast<unsigned long>(j)] == -1)
                {
                    continue;
                }
                internalStiffnessTripletCount += 9;
            }
        }
    }

    return cubeBaseIndices;
}

// TODO: This recomputes a bunch of stuff from ComputeCubeBaseIndices, instead just grab the total count from there
static long ComputeInternalStiffnessTripletCount(const std::vector<Cube>& cubes, const std::vector<int>& vertexIdxToSimulatedDof)
{
    long internalStiffnessTripletCount = 0;
    for (std::vector<Cube>::size_type cubeIdx = 0; cubeIdx < cubes.size(); cubeIdx++)
    {
        if (cubes[cubeIdx].AllCornersFixed())
        {
            continue;
        }

        // Map the cube's verts to DoFs in the global stiffness matrix
        std::array<int,8> indices;
        for (int y = 0; y < 8; y++)
        {
            const int vidx = cubes[cubeIdx].VertexIndex(y);
            indices[static_cast<unsigned long>(y)] = vertexIdxToSimulatedDof[static_cast<unsigned long>(vidx)];
        }

        // Add the local stiffness matrix into the global stiffness matrix
        for (int j = 0; j < 8; j++)
        {
            for (int i = 0; i < 8; i++)
            {
                if (indices[static_cast<unsigned long>(i)] == -1 || indices[static_cast<unsigned long>(j)] == -1)
                {
                    continue;
                }
                internalStiffnessTripletCount += 9;
            }
        }
    }

    return internalStiffnessTripletCount;
}

static int ComputeSimulatedDOFCount(const std::vector<int>& vertexIdxToSimulatedDof)
{
    return 3 * int(std::count_if(vertexIdxToSimulatedDof.begin(), vertexIdxToSimulatedDof.end(), [](const int idx){return idx != -1;}));
}

static void CreateMatrixWithPlaceholders(const long internalStiffnessTripletCount, const std::vector<Cube>& cubes, const std::vector<int>& vertexIdxToSimulatedDof, const std::vector<long>& cubeBaseIndices, SparseMatrixX& KInternal)
{
    // Construct the sparse matrix (for the sparsity pattern only)
    std::vector<Eigen::Triplet<Scalar>> triplets(static_cast<unsigned long>(internalStiffnessTripletCount));
    #pragma omp parallel for
    for (std::vector<Cube>::size_type cubeIdx = 0; cubeIdx < cubes.size(); cubeIdx++)
    {
        if (cubes[cubeIdx].AllCornersFixed())
        {
            continue;
        }

        const Eigen::Matrix<Scalar,24,24> localK = Eigen::Matrix<Scalar,24,24>::Ones();

        // Map the cube's verts to DoFs in the global stiffness matrix
        std::array<int,8> indices;
        for (int y = 0; y < 8; y++)
        {
            const int vrtIdx = cubes[cubeIdx].VertexIndex(y);
            indices[static_cast<unsigned long>(y)] = vertexIdxToSimulatedDof[static_cast<unsigned long>(vrtIdx)];
        }

        long current_index = cubeBaseIndices[cubeIdx];

        // Add the local stiffness matrix into the global stiffness matrix
        for (int j = 0; j < 8; j++) for (int i = 0; i < 8; i++)
        {
            if (indices[static_cast<unsigned long>(i)] == -1 || indices[static_cast<unsigned long>(j)] == -1)
            {
                continue;
            }

            const int rowStart = indices[static_cast<unsigned long>(i)];
            const int colStart = indices[static_cast<unsigned long>(j)];

            const int i3 = 3 * i;
            const int j3 = 3 * j;

            for (int m = 0; m < 3; m++) for (int n = 0; n < 3; n++)
            {
                triplets[static_cast<unsigned long>(current_index++)] = Eigen::Triplet<Scalar>(rowStart + m, colStart + n, localK(i3 + m, j3 + n));
            }
        }
    }

    KInternal.setFromTriplets(std::begin(triplets), std::end(triplets));
    assert(KInternal.innerNonZeroPtr() == nullptr); // Ensure that K is in compressed mode
}

// NB: This only builds the map for the upper diagonal part
static std::vector<std::vector<LocalHessianLocation>> PrecomputeCompressedToLocalHessianMap(const std::vector<Cube>& cubes, const std::vector<int>& vertexIdxToSimulatedDof, const SparseMatrixX& KInternal)
{
    // Build a map from row and col indices to compressed coefficient locations
    std::map<std::pair<int,int>,int> rowColToIndexMap;
    for (int k = 0; k < KInternal.outerSize(); k++)
    {
        for (SparseMatrixX::InnerIterator it(KInternal, k); it; ++it)
        {
            const int index = int(&it.value() - KInternal.valuePtr());

            if (it.row() < it.col())
            {
                continue;
            }

            std::pair<int,int> rowCol(it.row(), it.col());
            rowColToIndexMap.insert(std::pair<std::pair<int,int>,int>(rowCol, index));
        }
    }

    std::vector<std::vector<LocalHessianLocation>> compressedIndexToLocalHessianMap(static_cast<unsigned long>(KInternal.nonZeros()));

    // For each cube
    for (std::vector<Cube>::size_type cubeIdx = 0; cubeIdx < cubes.size(); cubeIdx++)
    {
        if (cubes[cubeIdx].AllCornersFixed())
        {
            continue;
        }

        // Map the cube's verts to DoFs in the global stiffness matrix
        std::array<int,8> indices;
        for (int y = 0; y < 8; y++)
        {
            const int vrtIdx = cubes[cubeIdx].VertexIndex(y);
            indices[static_cast<unsigned long>(y)] = vertexIdxToSimulatedDof[static_cast<unsigned long>(vrtIdx)];
        }

        // Add the local stiffness matrix into the global stiffness matrix
        for (int j = 0; j < 8; j++) for (int i = 0; i < 8; i++)
        {
            if (indices[static_cast<unsigned long>(i)] == -1 || indices[static_cast<unsigned long>(j)] == -1)
            {
                continue;
            }

            const int rowStart = indices[static_cast<unsigned long>(i)];
            const int colStart = indices[static_cast<unsigned long>(j)];

            const int i3 = 3 * i;
            const int j3 = 3 * j;

            for (int m = 0; m < 3; m++) for (int n = 0; n < 3; n++)
            {
                const int globalRow = rowStart + m;
                const int globalCol = colStart + n;

                if (globalRow < globalCol)
                {
                    continue;
                }

                // Get the compressed index of this coefficient
                const auto itr = rowColToIndexMap.find(std::make_pair(globalRow, globalCol));
                assert(itr != rowColToIndexMap.end());
                const int compressedIndex = itr->second;
                assert(compressedIndex >= 0);
                compressedIndexToLocalHessianMap[static_cast<unsigned long>(compressedIndex)].emplace_back(cubeIdx, i3 + m, j3 + n);
            }
        }
    }

    return compressedIndexToLocalHessianMap;
}

static std::vector<std::pair<int,int>> PrecomputeSymmetricHessianMap(const SparseMatrixX& KInternal)
{
    // Identify all entries above the diagonal
    std::map<std::pair<int,int>,int> upperEntries;
    for (int k = 0; k < KInternal.outerSize(); k++)
    {
        for (SparseMatrixX::InnerIterator it(KInternal, k); it; ++it)
        {
            if (it.col() > it.row())
            {
                const int destination = int(&it.value() - KInternal.valuePtr());
                std::pair<int,int> rowCol(it.row(), it.col());
                upperEntries.insert(std::pair<std::pair<int,int>,int>(rowCol, destination));
            }
        }
    }

    std::vector<std::pair<int,int>> correspondences;
    correspondences.reserve(upperEntries.size());

    // Identify the symmetric entries
    for (int k = 0; k < KInternal.outerSize(); k++)
    {
        for (SparseMatrixX::InnerIterator it(KInternal, k); it; ++it)
        {
            if (it.col() < it.row())
            {
                std::pair<int,int> colRow(it.col(), it.row());
                std::map<std::pair<int,int>,int>::iterator itr = upperEntries.find(colRow);
                assert(itr != upperEntries.end());

                const int destination = itr->second;
                const int source = int(&it.value() - KInternal.valuePtr());
                correspondences.emplace_back(source, destination);
            }
        }
    }

    return correspondences;
}

CubeMesh::CubeMesh()
: _internalStiffnessTripletCount(0)
, _DOFs(0)
{}

CubeMesh::CubeMesh(const InitialMesh& initialMesh)
: _restVertices(EigenVertsToSTDVerts(initialMesh.vertsRest))
, _cubes(GenerateCubes(initialMesh, EigenVertsToSTDVerts(initialMesh.vertsInit), _restVertices))
, _vertexIdxToSimulatedDof(ComputeVertexToIndexTable(initialMesh.kinematic))
, _cubeBaseIndices(ComputeCubeBaseIndices(_cubes, _vertexIdxToSimulatedDof))
, _internalStiffnessTripletCount(ComputeInternalStiffnessTripletCount(_cubes, _vertexIdxToSimulatedDof))
, _DOFs(ComputeSimulatedDOFCount(_vertexIdxToSimulatedDof))
, _compressedIndexToLocalHessianMap()
, _symmetricEntryMap()
{
    // Construct the sparse matrix (for the sparsity pattern only)
    SparseMatrixX KInternal(_DOFs, _DOFs);
    CreateMatrixWithPlaceholders(_internalStiffnessTripletCount, _cubes, _vertexIdxToSimulatedDof, _cubeBaseIndices, KInternal);

    _compressedIndexToLocalHessianMap = PrecomputeCompressedToLocalHessianMap(_cubes, _vertexIdxToSimulatedDof, KInternal);
    _symmetricEntryMap = PrecomputeSymmetricHessianMap(KInternal);
}

void CubeMesh::CreateInitialDynamicState(const InitialMesh& initialMesh, DynamicState& dynamicState)
{
    dynamicState.vertices = EigenVertsToSTDVerts(initialMesh.vertsInit);
    dynamicState.elementEnergyCache.resize(static_cast<long>(_cubes.size()));
    dynamicState.elementForceCache.resize(_cubes.size());
    dynamicState.localKs.resize(_cubes.size());
}

void CubeMesh::CreateInitialSparseMatrix(SparseMatrixX& K)
{
    K.resize(_DOFs, _DOFs);
    CreateMatrixWithPlaceholders(_internalStiffnessTripletCount, _cubes, _vertexIdxToSimulatedDof, _cubeBaseIndices, K);
}

void CubeMesh::ComputeForces(const DynamicState& dynamicState, const Material& material, VectorX& f) const
{
    assert(f.size() == _DOFs);
    assert(dynamicState.elementForceCache.size() == _cubes.size());

    f.setZero();

    // Compute the per-element force contributions
    #pragma omp parallel for
    for (std::vector<Cube>::size_type cidx = 0; cidx < _cubes.size(); cidx++)
    {
        if (_cubes[cidx].AllCornersFixed())
        {
            continue;
        }

        dynamicState.elementForceCache[cidx] = _cubes[cidx].ComputeForces(material, dynamicState.vertices);
    }

    // Accumulate the per-element contributions
    for (std::vector<Cube>::size_type cidx = 0; cidx < _cubes.size(); cidx++)
    {
        if (_cubes[cidx].AllCornersFixed())
        {
            continue;
        }

        // Scatter the forces
        for (int y = 0; y < 8; y++)
        {
            const int vidx = _cubes[cidx].VertexIndex(y);

            // No forces on constrained vertices
            if (_vertexIdxToSimulatedDof[static_cast<unsigned long>(vidx)] == -1)
            {
                continue;
            }

            // Copy the local force vector to the global force vector
            f.segment<3>(_vertexIdxToSimulatedDof[static_cast<unsigned long>(vidx)]) += dynamicState.elementForceCache[cidx].segment<3>(3 * y);
        }
    }
}

bool CubeMesh::IsVertexFixed(const int vidx) const
{
    assert(vidx >= 0);
    assert(vidx < int(_vertexIdxToSimulatedDof.size()));
    return _vertexIdxToSimulatedDof[static_cast<unsigned long>(vidx)] == -1;
}

int CubeMesh::GetNumSimulatedDOFs() const
{
    return _DOFs;
}

void CubeMesh::_ComputeInternalLocalStiffnessMatrices(const Material& material, DynamicState& dynamicState) const
{
    assert(dynamicState.localKs.size() == _cubes.size());

    #pragma omp parallel for
    for (std::vector<Cube>::size_type cubeIdx = 0; cubeIdx < _cubes.size(); cubeIdx++)
    {
        if (_cubes[cubeIdx].AllCornersFixed())
        {
            continue;
        }
        dynamicState.localKs[cubeIdx] = _cubes[cubeIdx].ComputeForceJacobian(material, dynamicState.vertices);
    }
}

void CubeMesh::ComputeStiffnessMatrixSparse(const Material& material, DynamicState& dynamicState, SparseMatrixX& K) const
{
    _ComputeInternalLocalStiffnessMatrices(material, dynamicState);

    const int numEntries = int(K.nonZeros());

    #pragma omp parallel for
    for (int entryNum = 0; entryNum < numEntries; entryNum++)
    {
        Scalar& value = K.valuePtr()[entryNum];
        value = 0.0;
        assert(entryNum < int(_compressedIndexToLocalHessianMap.size()));
        for (const LocalHessianLocation& location : _compressedIndexToLocalHessianMap[static_cast<unsigned long>(entryNum)])
        {
            assert(location.elementIdx >= 0);
            assert(location.elementIdx < int(dynamicState.localKs.size()));
            value += dynamicState.localKs[static_cast<unsigned long>(location.elementIdx)](location.localRow, location.localCol);
        }
    }

    // Copy the symmetric half over
    #pragma omp parallel for
    for (int entryNum = 0; entryNum < int(_symmetricEntryMap.size()); entryNum++)
    {
        const std::pair<int,int>& symmetricPair = _symmetricEntryMap[static_cast<unsigned long>(entryNum)];
        const int source = symmetricPair.first;
        const int destination = symmetricPair.second;
        K.valuePtr()[destination] = K.valuePtr()[source];
    }
}

Scalar CubeMesh::ComputeEnergy(const DynamicState& dynamicState, const Material& material) const
{
    assert(dynamicState.elementEnergyCache.size() == int(_cubes.size()));

    #pragma omp parallel for
    for (std::vector<Cube>::size_type cidx = 0; cidx < _cubes.size(); cidx++)
    {
        if (!_cubes[cidx].AllCornersFixed())
        {
            dynamicState.elementEnergyCache[static_cast<long>(cidx)] = _cubes[cidx].ComputeInternalEnergy(material, dynamicState.vertices);
        }
        else
        {
            dynamicState.elementEnergyCache[static_cast<long>(cidx)] = 0.0;
        }
    }

    return dynamicState.elementEnergyCache.sum();
}

void CubeMesh::ScatterDisplacements(const VectorX& u, DynamicState& dynamicState) const
{
    // For each vertex
    for (std::vector<int>::size_type vidx = 0; vidx < _vertexIdxToSimulatedDof.size(); vidx++)
    {
        // If the vertex is not fixed
        assert(vidx < _vertexIdxToSimulatedDof.size());
        const int simDofIndex = _vertexIdxToSimulatedDof[vidx];
        if (simDofIndex != -1)
        {
            assert(simDofIndex % 3 == 0);
            assert(vidx < dynamicState.vertices.size());
            assert(vidx < _restVertices.size());
            dynamicState.vertices[vidx] = _restVertices[vidx] + u.segment<3>(simDofIndex);
        }
    }
}

void CubeMesh::ComputeDisplacements(const DynamicState& dynamicState, VectorX& u) const
{
    assert(u.size() == _DOFs);

    // For each vertex
    for (std::vector<int>::size_type vidx = 0; vidx < _vertexIdxToSimulatedDof.size(); vidx++)
    {
        // If the vertex is not fixed
        const int simDofIndex = _vertexIdxToSimulatedDof[vidx];
        if (simDofIndex != -1)
        {
            assert(simDofIndex % 3 == 0);
            u.segment<3>(simDofIndex) = dynamicState.vertices[vidx] - _restVertices[vidx];
        }
    }
}

}
