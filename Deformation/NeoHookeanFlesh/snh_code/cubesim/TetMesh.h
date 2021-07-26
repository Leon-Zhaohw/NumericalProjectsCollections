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

#ifndef TET_MESH_H
#define TET_MESH_H

#include "Settings.h"

namespace CubeSim
{

class Material;

struct TetMeshCachedState final
{
    // The deformation gradient
    std::vector<Matrix3> F;
    // The singular value decomposition
    std::vector<Matrix3> U;
    std::vector<Matrix3> V;
    std::vector<Vector3> sigma;
};

struct ForceContributor final
{
    ForceContributor(int e, int r, int c)
    : element(e)
    , row(r)
    , col(c)
    {}

    int element;
    int row;
    int col;
};

class TetMesh final
{

public:

    TetMesh() = default;
    TetMesh(const std::vector<Vector3>& restVerts, const std::vector<Vector3>& verts, const std::vector<Vector4i>& tets, const std::vector<bool>& vertIsKinematic);

    int GetNumVerts() const;

    // TODO: Shorten the names of these
    const std::vector<Vector3>& GetVertices() const;
    const Vector3& GetVertex(const int idx) const;
    void SetVertex(const int idx, const Vector3& vert);

    void PerturbSimulatedVertices(const VectorX& dx);

    bool IsVertexFixed(const int vrtIdx);

    int GetNumSimulatedVertices() const;

    int GetNumTets() const;

    Matrix3 ComputeF(const int tetIdx) const;
    Matrix3 ComputeF(const int tetIdx, const VectorX& perturbation) const;

    Scalar ComputePotentialEnergy(const Material& material) const;
    Scalar ComputePerturbedPotentialEnergy(const Material& material) const;

    // Computes forces on *simulated* degrees of freedom
    void ComputeForce(const Material& material, VectorX& force) const;

    // Computes the Hessian on *simulated* degrees of freedom
    void ComputeHessian(const Material& material, SparseMatrixX& H) const;

    // Builds a Hessian with the correct sparsity pattern, but placeholder one values in lieu of the actual values
    void BuildHessianWithPlaceholderValues(SparseMatrixX& H) const;

    // TODO: Remove const
    void UpdateCurrentCachedState() const;

    // TODO: Remove const
    void UpdatePerturbedCachedState(const VectorX& perturbation) const;

    void SwapStateCaches();

private:

    #ifndef NDEBUG
    void _VerifyTetOrientation() const;
    #endif

    // True if all four verts of a tet are fixed
    bool _IsTetFixed(const int tetIdx) const;

    std::vector<std::vector<ForceContributor>> _BuildForceContributorMap() const;
    std::vector<std::vector<ForceContributor>> _BuildHessianContributorMap() const;

    std::vector<Vector3> _vertices;
    std::vector<Vector3> _restVertices; // Constant
    std::vector<Vector4i> _tets; // Constant
    std::vector<Scalar> _restVolumes; // Precomputed, constant
    std::vector<Matrix3> _DmInv; // Precomputed, constant
    std::vector<Matrix34> _Bm; // Precomputed, constant
    std::vector<Matrix912> _pFpU; // Precomputed, constant
    std::vector<bool> _vertIsFixed; // Precomputed, constant
    int _numSimulatedVerts;
    // For each vertex, maps from the global index to the simulated index (start dof is 3 * simIdx). -1 for kinematic vertices.
    std::vector<int> _vIdxToSimIdx;

    // Precomputed DoF maps for accelerating force and Hessian construction
    std::vector<std::vector<ForceContributor>> _forceContributorMap;
    std::vector<std::vector<ForceContributor>> _hessianContributorMap;
    // Cache space to aid in computing energy, forces, and the Hessian
    mutable VectorX _energySumCache;
    mutable std::vector<Matrix34> _elementForceCache;
    mutable std::vector<Matrix1212> _elementHessianCache;

    // TODO: Remove mutable
    mutable TetMeshCachedState _cachedCurrentState;
    mutable TetMeshCachedState _cachedPerturbedState;

};

}

#endif
