static Vector9 flatten(const Matrix3x3& A)
{
    Vector9 flattened;
    int index = 0;
    for (int y = 0; y < 3; y++)
        for (int x = 0; x < 3; x++, index++)
            flattened[index] = A(x, y);
    return flattened;
}
