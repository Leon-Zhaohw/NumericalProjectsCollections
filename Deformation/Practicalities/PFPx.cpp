static Matrix9x12 ComputePFPx(const Matrix3& DmInv)
{
    const Scalar m = DmInv(0, 0);
    const Scalar n = DmInv(0, 1);
    const Scalar o = DmInv(0, 2);
    const Scalar p = DmInv(1, 0);
    const Scalar q = DmInv(1, 1);
    const Scalar r = DmInv(1, 2);
    const Scalar s = DmInv(2, 0);
    const Scalar t = DmInv(2, 1);
    const Scalar u = DmInv(2, 2);

    const Scalar t1 = - m - p - s;
    const Scalar t2 = - n - q - t;
    const Scalar t3 = - o - r - u;

    Matrix9x12 PFPx;
    PFPx.setZero();
    PFPx(0, 0) = t1;
    PFPx(0, 3) = m;
    PFPx(0, 6) = p;
    PFPx(0, 9) = s;
    PFPx(1, 1) = t1;
    PFPx(1, 4) = m;
    PFPx(1, 7) = p;
    PFPx(1, 10) = s;
    PFPx(2, 2) = t1;
    PFPx(2, 5) = m;
    PFPx(2, 8) = p;
    PFPx(2, 11) = s;
    PFPx(3, 0) = t2;
    PFPx(3, 3) = n;
    PFPx(3, 6) = q;
    PFPx(3, 9) = t;
    PFPx(4, 1) = t2;
    PFPx(4, 4) = n;
    PFPx(4, 7) = q;
    PFPx(4, 10) = t;
    PFPx(5, 2) = t2;
    PFPx(5, 5) = n;
    PFPx(5, 8) = q;
    PFPx(5, 11) = t;
    PFPx(6, 0) = t3;
    PFPx(6, 3) = o;
    PFPx(6, 6) = r;
    PFPx(6, 9) = u;
    PFPx(7, 1) = t3;
    PFPx(7, 4) = o;
    PFPx(7, 7) = r;
    PFPx(7, 10) = u;
    PFPx(8, 2) = t3;
    PFPx(8, 5) = o;
    PFPx(8, 8) = r;
    PFPx(8, 11) = u;

    return PFPx;
}
