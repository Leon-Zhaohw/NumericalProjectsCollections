#include "stdafx.h"
#include "twigg/util.h"

const double maxValue = 1e3;

// Debugging stuff
template <>
bool isBad(const double& number)
{
	return !(number < maxValue && number > -maxValue);
}

template <>
bool isBad(const float& number)
{
	return !(number < maxValue && number > -maxValue);
}

template <>
bool isBad(const vl::Vec4d& vec)
{
	return isBad(vec[0]) || isBad(vec[1]) || isBad(vec[2]) || isBad(vec[3]);
}

template <>
bool isBad(const vl::Vec4f& vec)
{
	return isBad(vec[0]) || isBad(vec[1]) || isBad(vec[2]) || isBad(vec[3]);
}

template <>
bool isBad(const vl::Vec3d& vec)
{
	return isBad(vec[0]) || isBad(vec[1]) || isBad(vec[2]);
}

template <>
bool isBad(const vl::Vec3f& vec)
{
	return isBad(vec[0]) || isBad(vec[1]) || isBad(vec[2]);
}

template <>
bool isBad(const vl::Vec2d& vec)
{
	return isBad(vec[0]) || isBad(vec[1]);
}

template <>
bool isBad(const vl::Vec2f& vec)
{
	return isBad(vec[0]) || isBad(vec[1]);
}

template <>
bool isBad(const vl::Mat3d& mat)
{
	return isBad(mat[0]) || isBad(mat[1]) || isBad(mat[2]);
}

template <>
bool isBad(const vl::Mat3f& mat)
{
	return isBad(mat[0]) || isBad(mat[1]) || isBad(mat[2]);
}

template <>
bool isBad(const vl::Mat4d& mat)
{
	return isBad(mat[0]) || isBad(mat[1]) || isBad(mat[2]) || isBad(mat[3]);
}

template <>
bool isBad(const vl::Mat4f& mat)
{
	return isBad(mat[0]) || isBad(mat[1]) || isBad(mat[2]) || isBad(mat[3]);
}
