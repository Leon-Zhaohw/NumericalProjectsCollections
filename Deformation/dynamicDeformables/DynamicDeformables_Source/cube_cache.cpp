// quadrature point locations along the inside
// of a canonical cube
static const Scalar invSqrt3 = 1.0 / std::sqrt(3.0);
static const std::array<Vector3,8> restGaussPoints = 
{
    invSqrt3 * Vector3(-1, -1, -1),
    invSqrt3 * Vector3( 1, -1, -1),
    invSqrt3 * Vector3(-1,  1, -1),
    invSqrt3 * Vector3( 1,  1, -1),
    invSqrt3 * Vector3(-1, -1,  1),
    invSqrt3 * Vector3( 1, -1,  1),
    invSqrt3 * Vector3(-1,  1,  1),
    invSqrt3 * Vector3( 1,  1,  1)
};

// cache DmInv at startup
std::array<Matrix8x3,8> ComputeHDmInv(const Matrix3x8& Xbar)
{
    std::array<Matrix8x3,8> HDmInv;

    for (int x = 0; x < 8; x++)
    {
        // compute Dm
        const Matrix8x3 H = ComputeH(restGaussPoints[x]);
        const Matrix3 Dm = Xbar * H;
        const Matrix3 DmInverse = Dm.inverse();

        // cache it left-multiplied by H
        HDmInv[x] = H * DmInverse;
    }

    return HDmInv;
}
