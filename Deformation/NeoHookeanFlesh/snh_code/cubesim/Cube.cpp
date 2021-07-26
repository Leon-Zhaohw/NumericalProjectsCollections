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

#include "Cube.h"
#include "Material.h"

namespace CubeSim
{

static const Scalar invSqrt3 = 1.0 / std::sqrt(3.0);
static const std::array<Vector3,8> restGaussPoints = {invSqrt3 * Vector3(-1, -1, -1),
                                                      invSqrt3 * Vector3( 1, -1, -1),
                                                      invSqrt3 * Vector3(-1,  1, -1),
                                                      invSqrt3 * Vector3( 1,  1, -1),
                                                      invSqrt3 * Vector3(-1, -1,  1),
                                                      invSqrt3 * Vector3( 1, -1,  1),
                                                      invSqrt3 * Vector3(-1,  1,  1),
                                                      invSqrt3 * Vector3( 1,  1,  1)};

static std::array<Vector8,8> ComputeGaussPointWeights()
{
    std::array<Vector8,8> gpWeights;
    for (int gpIdx = 0; gpIdx < 8; gpIdx++)
    {
        const Array3 bc = 0.5 * (1.0 + restGaussPoints[static_cast<unsigned long>(gpIdx)].array());
        assert((bc >= 0.0).all());
        assert((bc <= 1.0).all());
        const Array3 bcm = 1.0 - bc;
        gpWeights[static_cast<unsigned long>(gpIdx)](0) = bcm.x() * bcm.y() * bcm.z();
        gpWeights[static_cast<unsigned long>(gpIdx)](1) =  bc.x() * bcm.y() * bcm.z();
        gpWeights[static_cast<unsigned long>(gpIdx)](2) = bcm.x() *  bc.y() * bcm.z();
        gpWeights[static_cast<unsigned long>(gpIdx)](3) =  bc.x() *  bc.y() * bcm.z();
        gpWeights[static_cast<unsigned long>(gpIdx)](4) = bcm.x() * bcm.y() *  bc.z();
        gpWeights[static_cast<unsigned long>(gpIdx)](5) =  bc.x() * bcm.y() *  bc.z();
        gpWeights[static_cast<unsigned long>(gpIdx)](6) = bcm.x() *  bc.y() *  bc.z();
        gpWeights[static_cast<unsigned long>(gpIdx)](7) =  bc.x() *  bc.y() *  bc.z();
        assert(fabs(gpWeights[static_cast<unsigned long>(gpIdx)].sum() - 1.0) <= 1.0e-6);
    }
    return gpWeights;
}

static const std::array<Vector8,8> gaussPointWeights = ComputeGaussPointWeights();

static Vector9 Flatten(const Matrix3& A)
{
    Vector9 flattened;
    int index = 0;
    for (int y = 0; y < 3; y++)
    {
        for (int x = 0; x < 3; x++, index++)
        {
            flattened[index] = A(x, y);
        }
    }
    return flattened;
}

static Vector24 Flatten(const Matrix38& A )
{
    Vector24 flattened;
    int index = 0;
    for (int y = 0; y < 8; y++)
    {
        for (int x = 0; x < 3; x++, index++)
        {
            flattened[index] = A(x, y);
        }
    }
    return flattened;
}

static Matrix38 VertexMatrix(const std::vector<Vector3>& verts, const std::array<int,8>& vertexIndices)
{
    Matrix38 output;
    for (int vidx = 0; vidx < 8; vidx++)
    {
        output.col(static_cast<long>(vidx)) = verts[static_cast<unsigned long>(vertexIndices[static_cast<unsigned long>(vidx)])];
    }
    return output;
}

// Shape function derivative matrix
static Matrix83 ComputeH(const Vector3& Xi)
{
    const Scalar Xi1Minus = (1.0 - Xi[0]);
    const Scalar Xi1Plus  = (1.0 + Xi[0]);
    const Scalar Xi2Minus = (1.0 - Xi[1]);
    const Scalar Xi2Plus  = (1.0 + Xi[1]);
    const Scalar Xi3Minus = (1.0 - Xi[2]);
    const Scalar Xi3Plus  = (1.0 + Xi[2]);

    Matrix83 H;

    H(0,0) = -Xi2Minus * Xi3Minus;
    H(0,1) = -Xi1Minus * Xi3Minus;
    H(0,2) = -Xi1Minus * Xi2Minus;

    H(1,0) =  Xi2Minus * Xi3Minus;
    H(1,1) = -Xi1Plus * Xi3Minus;
    H(1,2) = -Xi1Plus * Xi2Minus;

    H(2,0) = -Xi2Plus * Xi3Minus;
    H(2,1) =  Xi1Minus * Xi3Minus;
    H(2,2) = -Xi1Minus * Xi2Plus;

    H(3,0) =  Xi2Plus * Xi3Minus;
    H(3,1) =  Xi1Plus * Xi3Minus;
    H(3,2) = -Xi1Plus * Xi2Plus;

    H(4,0) = -Xi2Minus * Xi3Plus;
    H(4,1) = -Xi1Minus * Xi3Plus;
    H(4,2) =  Xi1Minus * Xi2Minus;

    H(5,0) =  Xi2Minus * Xi3Plus;
    H(5,1) = -Xi1Plus * Xi3Plus;
    H(5,2) =  Xi1Plus * Xi2Minus;

    H(6,0) = -Xi2Plus * Xi3Plus;
    H(6,1) =  Xi1Minus * Xi3Plus;
    H(6,2) =  Xi1Minus * Xi2Plus;

    H(7,0) =  Xi2Plus * Xi3Plus;
    H(7,1) =  Xi1Plus * Xi3Plus;
    H(7,2) =  Xi1Plus * Xi2Plus;

    return H / 8.0;
}

static Scalar ComputeRestVolume(const std::array<int,8>& vertexIndices, const std::vector<Vector3>& restVerts)
{
    const Vector3& x0 = restVerts[static_cast<unsigned long>(vertexIndices[0])];
    const Vector3& x1 = restVerts[static_cast<unsigned long>(vertexIndices[1])];
    const Vector3& x2 = restVerts[static_cast<unsigned long>(vertexIndices[2])];
    const Vector3& x4 = restVerts[static_cast<unsigned long>(vertexIndices[4])];
    return (x1.x() - x0.x()) * (x2.y() - x0.y()) * (x4.z() - x0.z());
}

static std::array<Matrix38,8> ComputeBmg(const Matrix38& Dm, const Scalar& restVolume)
{
    std::array<Matrix38,8> Bmg;

    for (int x = 0; x < 8; x++)
    {
        const Matrix83 H = ComputeH(restGaussPoints[static_cast<unsigned long>(x)]);
        const Matrix3 DmHg = Dm * H;
        const Matrix3 DmHgInverse = DmHg.inverse();

        Bmg[static_cast<unsigned long>(x)] = DmHgInverse.transpose() * H.transpose() * restVolume;
    }

    return Bmg;
}

static std::array<Matrix83,8> ComputeHDmHInv(const Matrix38& Dm)
{
    std::array<Matrix83,8> HDmHInv;

    for (int x = 0; x < 8; x++)
    {
        const Matrix83 H = ComputeH(restGaussPoints[static_cast<unsigned long>(x)]);
        const Matrix3 DmHg = Dm * H;
        const Matrix3 DmHgInverse = DmHg.inverse();

        HDmHInv[static_cast<unsigned long>(x)] = H * DmHgInverse;
    }

    return HDmHInv;
}

static Matrix3 ComputePartialFpartialu(const Matrix38& Dm, const Vector3& Xi, const int index)
{
    const Matrix83 H = ComputeH(Xi);

    // Everything except the u-varying component;
    const Matrix83 A = H * (Dm * H).inverse();

    // Probe each entry with a 1 to build each column
    Matrix38 delta = Matrix38::Zero();
    delta(index % 3, index / 3) = 1.0;

    return delta * A;
}

static std::array<Matrix924,8> ComputepFpuMatrix(const Matrix38& Dm)
{
    std::array<Matrix924,8> pFpuMatrix;

    for (int x = 0; x < 8; x++)
    {
        for (int i = 0; i < 24; i++)
        {
            const Matrix3 pFpu = ComputePartialFpartialu(Dm, restGaussPoints[static_cast<unsigned long>(x)], i);
            const Vector9 pFpuFlat = Flatten(pFpu);
            pFpuMatrix[static_cast<unsigned long>(x)].col(i) = pFpuFlat;
        }
    }

    return pFpuMatrix;
}

Cube::Cube(const bool cubeFixed, const std::array<int,8>& vertexIndices, const std::vector<Vector3>& verts, const std::vector<Vector3>& restVerts)
: _cubeFixed(cubeFixed)
, _vertexIndices(vertexIndices)
, _Dm(VertexMatrix(verts, vertexIndices))
, _restVolume(ComputeRestVolume(vertexIndices, restVerts))
, _Bmg(ComputeBmg(_Dm, _restVolume))
, _HDmHInv(ComputeHDmHInv(_Dm))
, _pFpuMatrix(ComputepFpuMatrix(_Dm))
{}

bool Cube::AllCornersFixed() const
{
    return _cubeFixed;
}

Matrix3 Cube::ComputeF(const std::vector<Vector3>& verts, const int gpIdx) const
{
    const Matrix38 Ds = VertexMatrix(verts, _vertexIndices);
    return Ds * _HDmHInv[static_cast<unsigned long>(gpIdx)];
}

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

Scalar Cube::ComputeInternalEnergy(const Material& material, const std::vector<Vector3>& verts) const
{
    Scalar totalEnergy = 0.0;
    for (int gpIdx = 0; gpIdx < 8; gpIdx++)
    {
        const Matrix3 F = ComputeF(verts, gpIdx);
        Matrix3 U;
        Matrix3 V;
        Vector3 S;
        if (material.EnergyNeedsSVD())
        {
            RotationVariantSVD(F, U, V, S);
        }
        totalEnergy += material.Psi(F, U, V, S);
    }
    return _restVolume * totalEnergy;
}

Vector24 Cube::_ComputeForces(const Material& material, const std::vector<Vector3>& verts, const int gpIdx) const
{
    const Matrix3 F = ComputeF(verts, gpIdx);
    Matrix3 U;
    Matrix3 V;
    Vector3 S;
    if (material.PK1NeedsSVD())
    {
        RotationVariantSVD(F, U, V, S);
    }

    const Matrix3 P = material.PK1(F, U, V, S);

    const Matrix38 G = - P * _Bmg[static_cast<unsigned long>(gpIdx)];

    const Vector24 force = Flatten(G);

    return force;
}

Vector24 Cube::ComputeForces(const Material& material, const std::vector<Vector3>& verts) const
{
    Vector24 forceVector = _ComputeForces(material, verts, 0);
    for (int gpIdx = 1; gpIdx < 8; gpIdx++)
    {
        forceVector += _ComputeForces(material, verts, gpIdx);
    }
    return forceVector;
}

Matrix2424 Cube::_ComputeForceJacobian(const Material& material, const std::vector<Vector3>& verts, const int gaussPoint) const
{
    const Matrix3 F = ComputeF(verts, gaussPoint);

    Matrix3 U;
    Matrix3 V;
    Vector3 S;
    RotationVariantSVD(F, U, V, S);

    return _pFpuMatrix[static_cast<unsigned long>(gaussPoint)].transpose() * (-_restVolume * material.ClampedPartialPpartialF(F, U, V, S)) * _pFpuMatrix[static_cast<unsigned long>(gaussPoint)];
}

Matrix2424 Cube::ComputeForceJacobian(const Material& material, const std::vector<Vector3>& verts) const
{
    Matrix2424 K = _ComputeForceJacobian(material, verts, 0);
    for (int x = 1; x < 8; x++)
    {
        K += _ComputeForceJacobian(material, verts, x);
    }
    return K;
}

int Cube::VertexIndex(const int vertexNumber) const
{
    assert(vertexNumber < 8);
    return _vertexIndices[static_cast<unsigned long>(vertexNumber)];
}

}
