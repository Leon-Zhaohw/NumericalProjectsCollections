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

#include "TetMesh.h"

#include "Material.h"

#include <cassert>
#include <iostream>
#ifndef NDEBUG
#include <set>
#endif

namespace CubeSim
{

//  3           2
//    o_______o
//    |      /|
//    |`    / |
//    | `  /  |
//    |   `   |
//    |  / `  |
//    | /   ` |
//    |/     `|
//    o-------o
//  0           1

static void RotationVariantSVD(const Matrix3& F, Matrix3& U, Matrix3& V, Vector3& S)
{
    const Eigen::JacobiSVD<Matrix3,Eigen::NoQRPreconditioner> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    V = svd.matrixV();
    S = svd.singularValues();
    if (U.determinant() < 0.0)
    {
        U.col(0) *= -1.0;
        S(0) *= -1.0;
    }
    if (V.determinant() < 0.0)
    {
        V.col(0) *= -1.0;
        S(0) *= -1.0;
    }
}

#ifndef NDEBUG
void TetMesh::_VerifyTetOrientation() const
{
    for (const Vector4i& tet : _tets)
    {
        const Vector3& v0 = _restVertices[static_cast<unsigned long>(tet(0))];
        const Vector3& v1 = _restVertices[static_cast<unsigned long>(tet(1))];
        const Vector3& v2 = _restVertices[static_cast<unsigned long>(tet(2))];
        const Vector3& v3 = _restVertices[static_cast<unsigned long>(tet(3))];
        const Scalar scaledVolume = (v3 - v0).dot((v1 - v0).cross(v2 - v0));
        assert(scaledVolume > 0.0);
    }
}
#endif

bool TetMesh::_IsTetFixed(const int tetIdx) const
{
    assert(tetIdx >= 0);
    assert(tetIdx < int(_tets.size()));
    const Vector4i& t = _tets[static_cast<unsigned long>(tetIdx)];
    assert((t.array() >= 0).all());
    assert((t.array() < int(_vertIsFixed.size())).all());
    return _vertIsFixed[static_cast<unsigned long>(t[0])] && _vertIsFixed[static_cast<unsigned long>(t[1])] && _vertIsFixed[static_cast<unsigned long>(t[2])] && _vertIsFixed[static_cast<unsigned long>(t[3])];
}

static int ComputeNumSimulatedVerts(const std::vector<bool>& vertIsFixed)
{
    return int(std::count(vertIsFixed.cbegin(), vertIsFixed.cend(), false));
}

static std::vector<int> ComputeGlobalToSimIdx(const std::vector<bool>& vertIsFixed)
{
    std::vector<int> globalToSimMap(vertIsFixed.size());
    int simIdx = 0;
    for (int globalIdx = 0; globalIdx < int(vertIsFixed.size()); globalIdx++)
    {
        if (!vertIsFixed[static_cast<unsigned long>(globalIdx)])
        {
            globalToSimMap[static_cast<unsigned long>(globalIdx)] = simIdx++;
        }
        else
        {
            globalToSimMap[static_cast<unsigned long>(globalIdx)] = -1;
        }
    }
    return globalToSimMap;
}

static Scalar ComputeTetVolume(const Vector3& v0, const Vector3& v1, const Vector3& v2, const Vector3& v3)
{
    return (v3 - v0).dot((v1 - v0).cross(v2 - v0)) / 6.0;
}

static std::vector<Scalar> ComputeRestVolumes(const std::vector<Vector4i>& tets, const std::vector<Vector3>& vertices)
{
    std::vector<Scalar> restVolumes(tets.size());
    for (int tetIdx = 0; tetIdx < int(tets.size()); tetIdx++)
    {
        const Vector4i& t = tets[static_cast<unsigned long>(tetIdx)];
        const Vector3& v0 = vertices[static_cast<unsigned long>(t[0])];
        const Vector3& v1 = vertices[static_cast<unsigned long>(t[1])];
        const Vector3& v2 = vertices[static_cast<unsigned long>(t[2])];
        const Vector3& v3 = vertices[static_cast<unsigned long>(t[3])];
        restVolumes[static_cast<unsigned long>(tetIdx)] = ComputeTetVolume(v0, v1, v2, v3);
        assert(restVolumes[static_cast<unsigned long>(tetIdx)] > 0.0);
    }
    return restVolumes;
}

static std::vector<Matrix3> ComputeDmInv(const std::vector<Vector4i>& tets, const std::vector<Vector3>& vertices)
{
    std::vector<Matrix3> DmInv(tets.size());
    for (int tetIdx = 0; tetIdx < int(tets.size()); tetIdx++)
    {
        const Vector4i& t = tets[static_cast<unsigned long>(tetIdx)];
        Matrix3 Dm;
        Dm.col(0) = vertices[static_cast<unsigned long>(t[1])] - vertices[static_cast<unsigned long>(t[0])];
        Dm.col(1) = vertices[static_cast<unsigned long>(t[2])] - vertices[static_cast<unsigned long>(t[0])];
        Dm.col(2) = vertices[static_cast<unsigned long>(t[3])] - vertices[static_cast<unsigned long>(t[0])];
        DmInv[static_cast<unsigned long>(tetIdx)] = Dm.fullPivLu().inverse();
        assert((Dm * DmInv[static_cast<unsigned long>(tetIdx)] - Matrix3::Identity()).lpNorm<Eigen::Infinity>() <= 1.0e-6);
    }
    return DmInv;
}

static std::vector<Matrix34> ComputeBm(const std::vector<Vector4i>& tets, const std::vector<Vector3>& vertices)
{
    // TODO: Eliminate this class, it is just used to get the area normal
    class Triangle final
    {
        public:
            Triangle(const Vector3& v0, const Vector3& v1, const Vector3& v2)
            : _x0(v0)
            , _x1(v1)
            , _x2(v2)
            {}

            // TODO: Make this a constructor
            static Triangle GenerateFace(const int faceNum, const Vector4i& tet, const std::vector<Vector3>& verts)
            {
                assert(faceNum >= 0); assert(faceNum < 4);
                if (faceNum == 0)
                {
                    return Triangle(verts[static_cast<unsigned long>(tet[0])], verts[static_cast<unsigned long>(tet[1])], verts[static_cast<unsigned long>(tet[3])]);
                }
                else if (faceNum == 1)
                {
                    return Triangle(verts[static_cast<unsigned long>(tet[0])], verts[static_cast<unsigned long>(tet[2])], verts[static_cast<unsigned long>(tet[1])]);
                }
                else if (faceNum == 2)
                {
                    return Triangle(verts[static_cast<unsigned long>(tet[3])], verts[static_cast<unsigned long>(tet[2])], verts[static_cast<unsigned long>(tet[0])]);
                }
                else if (faceNum == 3)
                {
                    return Triangle(verts[static_cast<unsigned long>(tet[1])], verts[static_cast<unsigned long>(tet[2])], verts[static_cast<unsigned long>(tet[3])]);
                }
                else
                {
                    std::cerr << "Error, impossible code path hit." << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            }

            // TODO: Roll these into one computation
            Vector3 normal() const
            {
                return ((_x1 - _x0).cross(_x2 - _x0)).normalized();
            }
            Scalar area() const
            {
                return 0.5 * ((_x1 - _x0).cross(_x2 - _x0)).norm();
            }
        private:
            const Vector3& _x0;
            const Vector3& _x1;
            const Vector3& _x2;
    };

    std::vector<Matrix34> Bm(tets.size());
    for (int tetIdx = 0; tetIdx < int(tets.size()); tetIdx++)
    {
        const Triangle tri0 = Triangle::GenerateFace(0, tets[static_cast<unsigned long>(tetIdx)], vertices);
        const Triangle tri1 = Triangle::GenerateFace(1, tets[static_cast<unsigned long>(tetIdx)], vertices);
        const Triangle tri2 = Triangle::GenerateFace(2, tets[static_cast<unsigned long>(tetIdx)], vertices);
        const Triangle tri3 = Triangle::GenerateFace(3, tets[static_cast<unsigned long>(tetIdx)], vertices);

        #ifndef NDEBUG
        const Vector3 resid = tri0.area() * tri0.normal() + tri1.area() * tri1.normal() + tri2.area() * tri2.normal() + tri3.area() * tri3.normal();
        assert( resid.lpNorm<Eigen::Infinity>() <= 1.0e-6 );
        #endif

        // Calculate the area vectors
        // v0 is incident on faces (0,1,2)
        Bm[static_cast<unsigned long>(tetIdx)].col(0) = tri0.normal() * tri0.area() + tri1.normal() * tri1.area() + tri2.normal() * tri2.area();
        // v1 is incident on faces (0,1,3)
        Bm[static_cast<unsigned long>(tetIdx)].col(1) = tri0.normal() * tri0.area() + tri1.normal() * tri1.area() + tri3.normal() * tri3.area();
        // v2 is incident on faces (1,2,3)
        Bm[static_cast<unsigned long>(tetIdx)].col(2) = tri1.normal() * tri1.area() + tri2.normal() * tri2.area() + tri3.normal() * tri3.area();
        // v3 is incident on faces (0,2,3)
        Bm[static_cast<unsigned long>(tetIdx)].col(3) = tri0.normal() * tri0.area() + tri2.normal() * tri2.area() + tri3.normal() * tri3.area();
        Bm[static_cast<unsigned long>(tetIdx)] /= -3.0;
    }

    return Bm;
}

static Matrix912 ComputePFPuEntry(const Matrix3& matInv)
{
    const Scalar m = matInv(0, 0);
    const Scalar n = matInv(0, 1);
    const Scalar o = matInv(0, 2);
    const Scalar p = matInv(1, 0);
    const Scalar q = matInv(1, 1);
    const Scalar r = matInv(1, 2);
    const Scalar s = matInv(2, 0);
    const Scalar t = matInv(2, 1);
    const Scalar u = matInv(2, 2);

    const Scalar t1 = - m - p - s;
    const Scalar t2 = - n - q - t;
    const Scalar t3 = - o - r - u;

    Matrix912 PFPu;
    PFPu(0, 0) = t1;
    PFPu(0, 1) = 0.0;
    PFPu(0, 2) = 0.0;
    PFPu(0, 3) = m;
    PFPu(0, 4) = 0.0;
    PFPu(0, 5) = 0.0;
    PFPu(0, 6) = p;
    PFPu(0, 7) = 0.0;
    PFPu(0, 8) = 0.0;
    PFPu(0, 9) = s;
    PFPu(0, 10) = 0.0;
    PFPu(0, 11) = 0.0;
    PFPu(1, 0) = 0.0;
    PFPu(1, 1) = t1;
    PFPu(1, 2) = 0.0;
    PFPu(1, 3) = 0.0;
    PFPu(1, 4) = m;
    PFPu(1, 5) = 0.0;
    PFPu(1, 6) = 0.0;
    PFPu(1, 7) = p;
    PFPu(1, 8) = 0.0;
    PFPu(1, 9) = 0.0;
    PFPu(1, 10) = s;
    PFPu(1, 11) = 0.0;
    PFPu(2, 0) = 0.0;
    PFPu(2, 1) = 0.0;
    PFPu(2, 2) = t1;
    PFPu(2, 3) = 0.0;
    PFPu(2, 4) = 0.0;
    PFPu(2, 5) = m;
    PFPu(2, 6) = 0.0;
    PFPu(2, 7) = 0.0;
    PFPu(2, 8) = p;
    PFPu(2, 9) = 0.0;
    PFPu(2, 10) = 0.0;
    PFPu(2, 11) = s;
    PFPu(3, 0) = t2;
    PFPu(3, 1) = 0.0;
    PFPu(3, 2) = 0.0;
    PFPu(3, 3) = n;
    PFPu(3, 4) = 0.0;
    PFPu(3, 5) = 0.0;
    PFPu(3, 6) = q;
    PFPu(3, 7) = 0.0;
    PFPu(3, 8) = 0.0;
    PFPu(3, 9) = t;
    PFPu(3, 10) = 0.0;
    PFPu(3, 11) = 0.0;
    PFPu(4, 0) = 0.0;
    PFPu(4, 1) = t2;
    PFPu(4, 2) = 0.0;
    PFPu(4, 3) = 0.0;
    PFPu(4, 4) = n;
    PFPu(4, 5) = 0.0;
    PFPu(4, 6) = 0.0;
    PFPu(4, 7) = q;
    PFPu(4, 8) = 0.0;
    PFPu(4, 9) = 0.0;
    PFPu(4, 10) = t;
    PFPu(4, 11) = 0.0;
    PFPu(5, 0) = 0.0;
    PFPu(5, 1) = 0.0;
    PFPu(5, 2) = t2;
    PFPu(5, 3) = 0.0;
    PFPu(5, 4) = 0.0;
    PFPu(5, 5) = n;
    PFPu(5, 6) = 0.0;
    PFPu(5, 7) = 0.0;
    PFPu(5, 8) = q;
    PFPu(5, 9) = 0.0;
    PFPu(5, 10) = 0.0;
    PFPu(5, 11) = t;
    PFPu(6, 0) = t3;
    PFPu(6, 1) = 0.0;
    PFPu(6, 2) = 0.0;
    PFPu(6, 3) = o;
    PFPu(6, 4) = 0.0;
    PFPu(6, 5) = 0.0;
    PFPu(6, 6) = r;
    PFPu(6, 7) = 0.0;
    PFPu(6, 8) = 0.0;
    PFPu(6, 9) = u;
    PFPu(6, 10) = 0.0;
    PFPu(6, 11) = 0.0;
    PFPu(7, 0) = 0.0;
    PFPu(7, 1) = t3;
    PFPu(7, 2) = 0.0;
    PFPu(7, 3) = 0.0;
    PFPu(7, 4) = o;
    PFPu(7, 5) = 0.0;
    PFPu(7, 6) = 0.0;
    PFPu(7, 7) = r;
    PFPu(7, 8) = 0.0;
    PFPu(7, 9) = 0.0;
    PFPu(7, 10) = u;
    PFPu(7, 11) = 0.0;
    PFPu(8, 0) = 0.0;
    PFPu(8, 1) = 0.0;
    PFPu(8, 2) = t3;
    PFPu(8, 3) = 0.0;
    PFPu(8, 4) = 0.0;
    PFPu(8, 5) = o;
    PFPu(8, 6) = 0.0;
    PFPu(8, 7) = 0.0;
    PFPu(8, 8) = r;
    PFPu(8, 9) = 0.0;
    PFPu(8, 10) = 0.0;
    PFPu(8, 11) = u;

    return PFPu;
}

static std::vector<Matrix912> ComputepFpU(const std::vector<Vector4i>& tets, const std::vector<Matrix3>& DmInv)
{
    std::vector<Matrix912> pFpU(tets.size());
    for (int tetIdx = 0; tetIdx < int(tets.size()); tetIdx++)
    {
        pFpU[static_cast<unsigned long>(tetIdx)] = ComputePFPuEntry(DmInv[static_cast<unsigned long>(tetIdx)]);
    }
    return pFpU;
}

#ifndef NDEBUG
static bool AllVerticesReferenced(const std::vector<Vector4i>& tets, const int numVerts)
{
    std::vector<bool> referenced(static_cast<unsigned long>(numVerts), false);
    for (const Vector4i& tet : tets)
    {
        for (int i = 0; i < tet.size(); i++)
        {
            assert(tet[i] >= 0);
            assert(tet[i] < int(referenced.size()));
            referenced[static_cast<unsigned long>(tet[i])] = true;
        }
    }
    return std::all_of(referenced.cbegin(), referenced.cend(), [](bool b){return b;});
}
#endif

#ifndef NDEBUG
// Check for duplicated tets
static bool AllTetsUnique(const std::vector<Vector4i>& tets)
{
    std::set<std::tuple<int,int,int,int>> tetSet;
    for (const Vector4i& tet : tets)
    {
        Vector4i sortedTet = tet;
        std::sort(&sortedTet[0], &sortedTet[0] + 4);
        assert(sortedTet[0] <= sortedTet[1]);
        assert(sortedTet[1] <= sortedTet[2]);
        assert(sortedTet[2] <= sortedTet[3]);
        auto setInsert = tetSet.insert(std::tie(sortedTet[0], sortedTet[1], sortedTet[2], sortedTet[3]));
        if (!setInsert.second)
        {
            return false;
        }
    }
    assert(tets.size() == tetSet.size());
    return true;
}
#endif

static std::vector<Vector4i> FixTetOrientations(const std::vector<Vector4i>& tets, const std::vector<Vector3>& verts)
{
    std::vector<Vector4i> fixedTets = tets;
    for (int tetIdx = 0; tetIdx < int(tets.size()); tetIdx++)
    {
        assert((fixedTets[static_cast<unsigned long>(tetIdx)].array() >= 0).all());
        assert((fixedTets[static_cast<unsigned long>(tetIdx)].array() < int(verts.size())).all());
        const Vector3& v0 = verts[static_cast<unsigned long>(fixedTets[static_cast<unsigned long>(tetIdx)](0))];
        const Vector3& v1 = verts[static_cast<unsigned long>(fixedTets[static_cast<unsigned long>(tetIdx)](1))];
        const Vector3& v2 = verts[static_cast<unsigned long>(fixedTets[static_cast<unsigned long>(tetIdx)](2))];
        const Vector3& v3 = verts[static_cast<unsigned long>(fixedTets[static_cast<unsigned long>(tetIdx)](3))];
        const Scalar scaledVolume = (v3 - v0).dot((v1 - v0).cross(v2 - v0));
        if (scaledVolume < 0.0)
        {
            std::swap(fixedTets[static_cast<unsigned long>(tetIdx)](2), fixedTets[static_cast<unsigned long>(tetIdx)](3));
        }
    }
    #ifndef NDEBUG
    for (int tetIdx = 0; tetIdx < int(tets.size()); tetIdx++)
    {
        const Vector3& v0 = verts[static_cast<unsigned long>(fixedTets[static_cast<unsigned long>(tetIdx)](0))];
        const Vector3& v1 = verts[static_cast<unsigned long>(fixedTets[static_cast<unsigned long>(tetIdx)](1))];
        const Vector3& v2 = verts[static_cast<unsigned long>(fixedTets[static_cast<unsigned long>(tetIdx)](2))];
        const Vector3& v3 = verts[static_cast<unsigned long>(fixedTets[static_cast<unsigned long>(tetIdx)](3))];
        const Scalar scaledVolume = (v3 - v0).dot((v1 - v0).cross(v2 - v0));
        assert(scaledVolume > 0.0);
    }
    #endif
    return fixedTets;
}

std::vector<std::vector<ForceContributor>> TetMesh::_BuildForceContributorMap() const
{
    std::vector<std::vector<ForceContributor>> forceContribMap(static_cast<unsigned long>(3 * _numSimulatedVerts));

    for (int tetIdx = 0; tetIdx < int(_tets.size()); tetIdx++)
    {
        if (_IsTetFixed(tetIdx))
        {
            continue;
        }

        if (!_vertIsFixed[static_cast<unsigned long>(_tets[static_cast<unsigned long>(tetIdx)][0])])
        {
            assert(_vIdxToSimIdx[_tets[static_cast<unsigned long>(tetIdx)][0] != -1]);
            const int baseIdx = 3 * _vIdxToSimIdx[static_cast<unsigned long>(_tets[static_cast<unsigned long>(tetIdx)][0])];
            forceContribMap[static_cast<unsigned long>(baseIdx + 0)].emplace_back(tetIdx, 0, 0);
            forceContribMap[static_cast<unsigned long>(baseIdx + 1)].emplace_back(tetIdx, 1, 0);
            forceContribMap[static_cast<unsigned long>(baseIdx + 2)].emplace_back(tetIdx, 2, 0);
        }
        if (!_vertIsFixed[static_cast<unsigned long>(_tets[static_cast<unsigned long>(tetIdx)][1])])
        {
            assert(_vIdxToSimIdx[_tets[static_cast<unsigned long>(tetIdx)][1] != -1]);
            const int baseIdx = 3 * _vIdxToSimIdx[static_cast<unsigned long>(_tets[static_cast<unsigned long>(tetIdx)][1])];
            forceContribMap[static_cast<unsigned long>(baseIdx + 0)].emplace_back(tetIdx, 0, 1);
            forceContribMap[static_cast<unsigned long>(baseIdx + 1)].emplace_back(tetIdx, 1, 1);
            forceContribMap[static_cast<unsigned long>(baseIdx + 2)].emplace_back(tetIdx, 2, 1);
        }
        if (!_vertIsFixed[static_cast<unsigned long>(_tets[static_cast<unsigned long>(tetIdx)][2])])
        {
            assert(_vIdxToSimIdx[_tets[static_cast<unsigned long>(tetIdx)][2] != -1]);
            const int baseIdx = 3 * _vIdxToSimIdx[static_cast<unsigned long>(_tets[static_cast<unsigned long>(tetIdx)][2])];
            forceContribMap[static_cast<unsigned long>(baseIdx + 0)].emplace_back(tetIdx, 0, 2);
            forceContribMap[static_cast<unsigned long>(baseIdx + 1)].emplace_back(tetIdx, 1, 2);
            forceContribMap[static_cast<unsigned long>(baseIdx + 2)].emplace_back(tetIdx, 2, 2);
        }
        if (!_vertIsFixed[static_cast<unsigned long>(_tets[static_cast<unsigned long>(tetIdx)][3])])
        {
            assert(_vIdxToSimIdx[_tets[static_cast<unsigned long>(tetIdx)][3] != -1]);
            const int baseIdx = 3 * _vIdxToSimIdx[static_cast<unsigned long>(_tets[static_cast<unsigned long>(tetIdx)][3])];
            forceContribMap[static_cast<unsigned long>(baseIdx + 0)].emplace_back(tetIdx, 0, 3);
            forceContribMap[static_cast<unsigned long>(baseIdx + 1)].emplace_back(tetIdx, 1, 3);
            forceContribMap[static_cast<unsigned long>(baseIdx + 2)].emplace_back(tetIdx, 2, 3);
        }
    }

    return forceContribMap;
}

void TetMesh::BuildHessianWithPlaceholderValues(SparseMatrixX& H) const
{
    H.resize(3 * _numSimulatedVerts, 3 * _numSimulatedVerts);
    std::vector<Eigen::Triplet<Scalar>> triplets;
    for (int tetIdx = 0; tetIdx < int(_tets.size()); tetIdx++)
    {
        if (_IsTetFixed(tetIdx))
        {
            continue;
        }

        const Vector4i& t = _tets[static_cast<unsigned long>(tetIdx)];

        for (int i = 0; i < 4; i++)
        {
            if (_vertIsFixed[static_cast<unsigned long>(t(i))])
            {
                continue;
            }
            const int simIdx0 = _vIdxToSimIdx[static_cast<unsigned long>(t(i))];
            assert(simIdx0 != -1);
            for (int j = 0; j < 4; j++)
            {
                if (_vertIsFixed[static_cast<unsigned long>(t(j))])
                {
                    continue;
                }
                const int simIdx1 = _vIdxToSimIdx[static_cast<unsigned long>(t(j))];
                assert(simIdx1 != -1);
                for (int a = 0; a < 3; a++)
                {
                    for (int b = 0; b < 3; b++)
                    {
                        triplets.emplace_back(3 * simIdx0 + a, 3 * simIdx1 + b, 1.0);
                    }
                }
            }
        }
    }
    H.setFromTriplets(triplets.cbegin(), triplets.cend());
    assert(H.innerNonZeroPtr() == nullptr); // Double check that H is in compressed mode
    // TODO: Assert symmetry
}

std::vector<std::vector<ForceContributor>> TetMesh::_BuildHessianContributorMap() const
{
    // Construct the sparse Hessian to get the non-zero structure
    SparseMatrixX H;
    BuildHessianWithPlaceholderValues(H);

    // Map from *global* row and column pairs to global value locations
    std::map<std::pair<int,int>,int> rowColToIndexMap;
    for (int k = 0; k < H.outerSize(); k++)
    {
        for (SparseMatrixX::InnerIterator it(H, k); it; ++it)
        {
            const int index = int(&it.value() - H.valuePtr());
            std::pair<int,int> rowCol(it.row(), it.col());
            rowColToIndexMap.insert(std::pair<std::pair<int,int>,int>(rowCol, index));
        }
    }
    assert(int(rowColToIndexMap.size()) == H.nonZeros());

    // Build a map from compressed coefficient locations to entries in the local Hessians
    std::vector<std::vector<ForceContributor>> hessianContribMap(static_cast<unsigned long>(H.nonZeros()));
    for (int tetIdx = 0; tetIdx < int(_tets.size()); tetIdx++)
    {
        if (_IsTetFixed(tetIdx))
        {
            continue;
        }
        const Vector4i& t = _tets[static_cast<unsigned long>(tetIdx)];
        for (int i = 0; i < 4; i++)
        {
            if (_vertIsFixed[static_cast<unsigned long>(t(i))])
            {
                continue;
            }
            const int simIdx0 = _vIdxToSimIdx[static_cast<unsigned long>(t(i))];
            assert(simIdx0 != -1);
            for (int j = 0; j < 4; j++)
            {
                if (_vertIsFixed[static_cast<unsigned long>(t(j))])
                {
                    continue;
                }
                const int simIdx1 = _vIdxToSimIdx[static_cast<unsigned long>(t(j))];
                assert(simIdx1 != -1);
                for (int a = 0; a < 3; a++)
                {
                    for (int b = 0; b < 3; b++)
                    {
                        const int globalRow = 3 * simIdx0 + a;
                        const int globalCol = 3 * simIdx1 + b;
                        // Get the compressed index of this coefficient
                        const auto itr = rowColToIndexMap.find(std::make_pair(globalRow, globalCol));
                        assert(itr != rowColToIndexMap.end());
                        const int compressedIndex = itr->second;
                        assert(compressedIndex >= 0);
                        hessianContribMap[static_cast<unsigned long>(compressedIndex)].emplace_back(tetIdx, 3 * i + a, 3 * j + b);
                    }
                }
            }
        }
    }

    return hessianContribMap;
}

TetMesh::TetMesh(const std::vector<Vector3>& restVerts, const std::vector<Vector3>& verts, const std::vector<Vector4i>& tetsIn, const std::vector<bool>& vertIsKinematic)
: _vertices(verts)
, _restVertices(restVerts)
, _tets(FixTetOrientations(tetsIn, restVerts))
, _restVolumes(ComputeRestVolumes(_tets, restVerts))
, _DmInv(ComputeDmInv(_tets, restVerts))
, _Bm(ComputeBm(_tets, restVerts))
, _pFpU(ComputepFpU(_tets, _DmInv))
, _vertIsFixed(vertIsKinematic)
, _numSimulatedVerts(ComputeNumSimulatedVerts(_vertIsFixed))
, _vIdxToSimIdx(ComputeGlobalToSimIdx(_vertIsFixed))
, _forceContributorMap(_BuildForceContributorMap())
, _hessianContributorMap(_BuildHessianContributorMap())
, _energySumCache(_tets.size())
, _elementForceCache(_tets.size())
, _elementHessianCache(_tets.size())
{
    assert(_vertices.size() == _restVertices.size());
    assert(std::all_of(_tets.cbegin(), _tets.cend(), [](const Vector4i& t){return (t.array() >= 0).all();}));
    assert(std::all_of(_tets.cbegin(), _tets.cend(), [this](const Vector4i& t){return (t.array() < int(this->_vertices.size())).all();}));
    assert(AllTetsUnique(_tets));
    #ifndef NDEBUG
    _VerifyTetOrientation();
    #endif
    assert(AllVerticesReferenced(_tets, int(_vertices.size())));
    // TODO: Check that surface triangles have valid indices
    assert(_restVolumes.size() == _tets.size());
    // TODO: Check that rest volumes are positive
    assert(_DmInv.size() == _tets.size());
    // TODO: Check Dm * DmInv
    assert(_Bm.size() == _tets.size());
    assert(_pFpU.size() == _tets.size());
    assert(_vertIsFixed.size() == _vertices.size());
    assert(_vIdxToSimIdx.size() == _vertices.size());
    assert(std::all_of(_vIdxToSimIdx.cbegin(), _vIdxToSimIdx.cend(), [](const int i){return i >= -1;}));
    assert(std::all_of(_vIdxToSimIdx.cbegin(), _vIdxToSimIdx.cend(), [this](const int i){return i < this->_numSimulatedVerts;}));

    UpdateCurrentCachedState();
    _cachedPerturbedState = _cachedCurrentState;
}

int TetMesh::GetNumVerts() const
{
    assert(_vertices.size() == _restVertices.size());
    return int(_vertices.size());
}

const std::vector<Vector3>& TetMesh::GetVertices() const
{
    return _vertices;
}

const Vector3& TetMesh::GetVertex(const int idx) const
{
    assert(idx >= 0);
    assert(idx < GetNumVerts());
    return _vertices[static_cast<unsigned long>(idx)];
}

void TetMesh::SetVertex(const int idx, const Vector3& vert)
{
    _vertices[static_cast<unsigned long>(idx)] = vert;
}

void TetMesh::PerturbSimulatedVertices(const VectorX& dx)
{
    assert(dx.size() == 3 * _numSimulatedVerts);
    int curSimIdx = 0;
    for (int vidx = 0; vidx < GetNumVerts(); vidx++)
    {
        if (!_vertIsFixed[static_cast<unsigned long>(vidx)])
        {
            _vertices[static_cast<unsigned long>(vidx)] += dx.segment<3>(3 * curSimIdx);
            curSimIdx++;
        }
    }
}

bool TetMesh::IsVertexFixed(const int vrtIdx)
{
    assert(vrtIdx >= 0);
    assert(vrtIdx < int(_vertIsFixed.size()));
    return _vertIsFixed[static_cast<unsigned long>(vrtIdx)];
}

int TetMesh::GetNumSimulatedVertices() const
{
    return _numSimulatedVerts;
}

int TetMesh::GetNumTets() const
{
    return int(_tets.size());
}

Matrix3 TetMesh::ComputeF(const int tetIdx) const
{
    assert(tetIdx >= 0);
    assert(tetIdx < GetNumTets());
    const Vector4i& t = _tets[static_cast<unsigned long>(tetIdx)];
    Matrix3 Ds;
    Eigen::Map<Vector3>(Ds.data()) = _vertices[static_cast<unsigned long>(t[1])] - _vertices[static_cast<unsigned long>(t[0])];
    Eigen::Map<Vector3>(Ds.data() + 3) = _vertices[static_cast<unsigned long>(t[2])] - _vertices[static_cast<unsigned long>(t[0])];
    Eigen::Map<Vector3>(Ds.data() + 6) = _vertices[static_cast<unsigned long>(t[3])] - _vertices[static_cast<unsigned long>(t[0])];
    return Ds * _DmInv[static_cast<unsigned long>(tetIdx)];
}

Matrix3 TetMesh::ComputeF(const int tetIdx, const VectorX& perturbation) const
{
    assert(tetIdx >= 0);
    assert(tetIdx < GetNumTets());
    assert(perturbation.size() % 3 == 0);

    const Vector4i& t = _tets[static_cast<unsigned long>(tetIdx)];

    Vector3 x0 = _vertices[static_cast<unsigned long>(t[0])];
    if (_vIdxToSimIdx[static_cast<unsigned long>(t[0])] != -1)
    {
        assert(_vIdxToSimIdx[static_cast<unsigned long>(t[0])] < perturbation.size() / 3);
        x0 += perturbation.segment<3>(3 * _vIdxToSimIdx[static_cast<unsigned long>(t[0])]);
    }

    Vector3 x1 = _vertices[static_cast<unsigned long>(t[1])];
    if (_vIdxToSimIdx[static_cast<unsigned long>(t[1])] != -1)
    {
        assert(_vIdxToSimIdx[static_cast<unsigned long>(t[1])] < perturbation.size() / 3);
        x1 += perturbation.segment<3>(3 * _vIdxToSimIdx[static_cast<unsigned long>(t[1])]);
    }

    Vector3 x2 = _vertices[static_cast<unsigned long>(t[2])];
    if (_vIdxToSimIdx[static_cast<unsigned long>(t[2])] != -1)
    {
        assert(_vIdxToSimIdx[static_cast<unsigned long>(t[2])] < perturbation.size() / 3);
        x2 += perturbation.segment<3>(3 * _vIdxToSimIdx[static_cast<unsigned long>(t[2])]);
    }

    Vector3 x3 = _vertices[static_cast<unsigned long>(t[3])];
    if (_vIdxToSimIdx[static_cast<unsigned long>(t[3])] != -1)
    {
        assert(_vIdxToSimIdx[static_cast<unsigned long>(t[3])] < perturbation.size() / 3);
        x3 += perturbation.segment<3>(3 * _vIdxToSimIdx[static_cast<unsigned long>(t[3])]);
    }

    Matrix3 Ds;
    Eigen::Map<Vector3>(Ds.data()) = x1 - x0;
    Eigen::Map<Vector3>(Ds.data() + 3) = x2 - x0;
    Eigen::Map<Vector3>(Ds.data() + 6) = x3 - x0;
    return Ds * _DmInv[static_cast<unsigned long>(tetIdx)];
}

Scalar TetMesh::ComputePotentialEnergy(const Material& material) const
{
    assert(_tets.size() == _restVolumes.size());
    assert(int(_tets.size()) == _energySumCache.size());

    #pragma omp parallel for
    for (int tetIdx = 0; tetIdx < int(_tets.size()); tetIdx++)
    {
        if (_IsTetFixed(tetIdx))
        {
            _energySumCache[tetIdx] = 0.0;
            continue;
        }
        const Matrix3& F = _cachedCurrentState.F[static_cast<unsigned long>(tetIdx)];
        // Ensure the cache is up to date
        #ifndef NDEBUG
        {
            const Matrix3 Ftest = ComputeF(tetIdx);
            assert((Ftest - F).lpNorm<Eigen::Infinity>() <= 1.0e-9);

            Matrix3 UTest;
            Matrix3 VTest;
            Vector3 sigmaTest;
            RotationVariantSVD(_cachedCurrentState.F[static_cast<unsigned long>(tetIdx)], UTest, VTest, sigmaTest);
            const Scalar Uresid = (UTest - _cachedCurrentState.U[static_cast<unsigned long>(tetIdx)]).lpNorm<Eigen::Infinity>();
            assert(Uresid <= 1.0e-9);
            const Scalar Vresid = (VTest - _cachedCurrentState.V[static_cast<unsigned long>(tetIdx)]).lpNorm<Eigen::Infinity>();
            assert(Vresid <= 1.0e-9);
            const Scalar Sresid = (sigmaTest - _cachedCurrentState.sigma[static_cast<unsigned long>(tetIdx)]).lpNorm<Eigen::Infinity>();
            assert(Sresid <= 1.0e-9);
        }
        #endif
        const Matrix3& U = _cachedCurrentState.U[static_cast<unsigned long>(tetIdx)];
        const Matrix3& V = _cachedCurrentState.V[static_cast<unsigned long>(tetIdx)];
        const Vector3& sigma = _cachedCurrentState.sigma[static_cast<unsigned long>(tetIdx)];
        _energySumCache[tetIdx] = _restVolumes[static_cast<unsigned long>(tetIdx)] * material.Psi(F, U, V, sigma);
    }

    return _energySumCache.sum();
}

Scalar TetMesh::ComputePerturbedPotentialEnergy(const Material& material) const
{
    assert(_tets.size() == _restVolumes.size());
    assert(int(_tets.size()) == _energySumCache.size());

    #pragma omp parallel for
    for (int tetIdx = 0; tetIdx < int(_tets.size()); tetIdx++)
    {
        if (_IsTetFixed(tetIdx))
        {
            _energySumCache[tetIdx] = 0.0;
            continue;
        }
        const Matrix3& F = _cachedPerturbedState.F[static_cast<unsigned long>(tetIdx)];
        const Matrix3& U = _cachedPerturbedState.U[static_cast<unsigned long>(tetIdx)];
        const Matrix3& V = _cachedPerturbedState.V[static_cast<unsigned long>(tetIdx)];
        const Vector3& sigma = _cachedPerturbedState.sigma[static_cast<unsigned long>(tetIdx)];
        _energySumCache[tetIdx] = _restVolumes[static_cast<unsigned long>(tetIdx)] * material.Psi(F, U, V, sigma);
    }

    return _energySumCache.sum();
}

void TetMesh::ComputeForce(const Material& material, VectorX& force) const
{
    assert(_tets.size() == _elementForceCache.size());
    assert(_tets.size() == _Bm.size());
    #pragma omp parallel for
    for (int tetIdx = 0; tetIdx < int(_tets.size()); tetIdx++)
    {
        if (_IsTetFixed(tetIdx))
        {
            continue;
        }
        const Matrix3& F = _cachedCurrentState.F[static_cast<unsigned long>(tetIdx)];
        // Ensure the cache is up to date
        #ifndef NDEBUG
        {
            const Matrix3 Ftest = ComputeF(tetIdx);
            assert((Ftest - F).lpNorm<Eigen::Infinity>() <= 1.0e-9);

            Matrix3 UTest;
            Matrix3 VTest;
            Vector3 sigmaTest;
            RotationVariantSVD(_cachedCurrentState.F[static_cast<unsigned long>(tetIdx)], UTest, VTest, sigmaTest);
            const Scalar Uresid = (UTest - _cachedCurrentState.U[static_cast<unsigned long>(tetIdx)]).lpNorm<Eigen::Infinity>();
            assert(Uresid <= 1.0e-9);
            const Scalar Vresid = (VTest - _cachedCurrentState.V[static_cast<unsigned long>(tetIdx)]).lpNorm<Eigen::Infinity>();
            assert(Vresid <= 1.0e-9);
            const Scalar Sresid = (sigmaTest - _cachedCurrentState.sigma[static_cast<unsigned long>(tetIdx)]).lpNorm<Eigen::Infinity>();
            assert(Sresid <= 1.0e-9);
        }
        #endif
        const Matrix3& U = _cachedCurrentState.U[static_cast<unsigned long>(tetIdx)];
        const Matrix3& V = _cachedCurrentState.V[static_cast<unsigned long>(tetIdx)];
        const Vector3& sigma = _cachedCurrentState.sigma[static_cast<unsigned long>(tetIdx)];
        _elementForceCache[static_cast<unsigned long>(tetIdx)] = material.PK1(F, U, V, sigma) * _Bm[static_cast<unsigned long>(tetIdx)];
    }

    assert(force.size() == 3 * _numSimulatedVerts);
    assert(int(_forceContributorMap.size()) == force.size());
    #pragma omp parallel for
    for (int dof = 0; dof < force.size(); dof++)
    {
        force[dof] = 0.0;
        for (const ForceContributor& c : _forceContributorMap[static_cast<unsigned long>(dof)])
        {
            force[static_cast<long>(dof)] += _elementForceCache[static_cast<unsigned long>(c.element)](c.row, c.col);
        }
    }
}

void TetMesh::ComputeHessian(const Material& material, SparseMatrixX& H) const
{
    assert(_elementHessianCache.size() == _tets.size());

    #pragma omp parallel for
    for (int tetIdx = 0; tetIdx < int(_tets.size()); tetIdx++)
    {
        if (_IsTetFixed(tetIdx))
        {
            continue;
        }
        // Ensure the cache is up to date
        const Matrix3& F = _cachedCurrentState.F[static_cast<unsigned long>(tetIdx)];
        #ifndef NDEBUG
        {
            const Matrix3 Ftest = ComputeF(tetIdx);
            assert((Ftest - F).lpNorm<Eigen::Infinity>() <= 1.0e-9);

            Matrix3 UTest;
            Matrix3 VTest;
            Vector3 sigmaTest;
            RotationVariantSVD(_cachedCurrentState.F[static_cast<unsigned long>(tetIdx)], UTest, VTest, sigmaTest);
            const Scalar Uresid = (UTest - _cachedCurrentState.U[static_cast<unsigned long>(tetIdx)]).lpNorm<Eigen::Infinity>();
            assert(Uresid <= 1.0e-9);
            const Scalar Vresid = (VTest - _cachedCurrentState.V[static_cast<unsigned long>(tetIdx)]).lpNorm<Eigen::Infinity>();
            assert(Vresid <= 1.0e-9);
            const Scalar Sresid = (sigmaTest - _cachedCurrentState.sigma[static_cast<unsigned long>(tetIdx)]).lpNorm<Eigen::Infinity>();
            assert(Sresid <= 1.0e-9);
        }
        #endif

        const Matrix3& U = _cachedCurrentState.U[static_cast<unsigned long>(tetIdx)];
        const Matrix3& V = _cachedCurrentState.V[static_cast<unsigned long>(tetIdx)];
        const Vector3& sigma = _cachedCurrentState.sigma[static_cast<unsigned long>(tetIdx)];
        const Matrix9 pPpF = _restVolumes[static_cast<unsigned long>(tetIdx)] * material.ClampedPartialPpartialF(F, U, V, sigma);
        _elementHessianCache[static_cast<unsigned long>(tetIdx)] = (_pFpU[static_cast<unsigned long>(tetIdx)].transpose() * pPpF) * _pFpU[static_cast<unsigned long>(tetIdx)];
    }

    assert(H.nonZeros() == int(_hessianContributorMap.size()));
    #pragma omp parallel for
    for (int dofIdx = 0; dofIdx < int(_hessianContributorMap.size()); dofIdx++)
    {
        H.valuePtr()[dofIdx] = 0.0;
        for (const ForceContributor& c : _hessianContributorMap[static_cast<unsigned long>(dofIdx)])
        {
            H.valuePtr()[dofIdx] += _elementHessianCache[static_cast<unsigned long>(c.element)](c.row, c.col);
        }
    }
}

void TetMesh::UpdateCurrentCachedState() const
{
    const int ntets = GetNumTets();
    _cachedCurrentState.F.resize(static_cast<unsigned long>(ntets));
    _cachedCurrentState.U.resize(static_cast<unsigned long>(ntets));
    _cachedCurrentState.V.resize(static_cast<unsigned long>(ntets));
    _cachedCurrentState.sigma.resize(static_cast<unsigned long>(ntets));

    #pragma omp parallel for
    for (int tetIdx = 0; tetIdx < ntets; tetIdx++)
    {
        _cachedCurrentState.F[static_cast<unsigned long>(tetIdx)] = ComputeF(tetIdx);
    }

    #pragma omp parallel for
    for (int tetIdx = 0; tetIdx < ntets; tetIdx++)
    {
        RotationVariantSVD(_cachedCurrentState.F[static_cast<unsigned long>(tetIdx)], _cachedCurrentState.U[static_cast<unsigned long>(tetIdx)], _cachedCurrentState.V[static_cast<unsigned long>(tetIdx)], _cachedCurrentState.sigma[static_cast<unsigned long>(tetIdx)]);
    }
}

void TetMesh::UpdatePerturbedCachedState(const VectorX& perturbation) const
{
    const int ntets = GetNumTets();
    _cachedPerturbedState.F.resize(static_cast<unsigned long>(ntets));
    _cachedPerturbedState.U.resize(static_cast<unsigned long>(ntets));
    _cachedPerturbedState.V.resize(static_cast<unsigned long>(ntets));
    _cachedPerturbedState.sigma.resize(static_cast<unsigned long>(ntets));

    #pragma omp parallel for
    for (int tetIdx = 0; tetIdx < ntets; tetIdx++)
    {
        _cachedPerturbedState.F[static_cast<unsigned long>(tetIdx)] = ComputeF(tetIdx, perturbation);
    }

    #pragma omp parallel for
    for (int tetIdx = 0; tetIdx < ntets; tetIdx++)
    {
        RotationVariantSVD(_cachedPerturbedState.F[static_cast<unsigned long>(tetIdx)], _cachedPerturbedState.U[static_cast<unsigned long>(tetIdx)], _cachedPerturbedState.V[static_cast<unsigned long>(tetIdx)], _cachedPerturbedState.sigma[static_cast<unsigned long>(tetIdx)]);
    }
}

void TetMesh::SwapStateCaches()
{
    std::swap(_cachedCurrentState, _cachedPerturbedState);
}

}
