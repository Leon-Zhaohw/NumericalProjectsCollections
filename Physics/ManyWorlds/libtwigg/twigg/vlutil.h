#ifndef __VLUTIL_H__
#define __VLUTIL_H__

#include <vl/VL.h>
#include <vl/VLf.h>
#include <vl/VLd.h>

#include <algorithm>

#if 0
#define CL_FLOAT

namespace vl
{
  typedef Vec2f Vec2;
  typedef Vec3f Vec3;
  typedef Vec4f Vec4;
  typedef Mat2f Mat2;
  typedef Mat3f Mat3;
  typedef Mat4f Mat4;
}
#else
namespace vl
{
  typedef Vec2d Vec2;
  typedef Vec3d Vec3;
  typedef Vec4d Vec4;
  typedef Mat2d Mat2;
  typedef Mat3d Mat3;
  typedef Mat4d Mat4;
}
#endif

namespace vl
{

	
inline vl::Mat3d outer_product( const vl::Vec3d& left, const vl::Vec3d& right )
{
	return vl::Mat3d( 
		left[0]*right[0], left[0]*right[1], left[0]*right[2],
		left[1]*right[0], left[1]*right[1], left[1]*right[2],
		left[2]*right[0], left[2]*right[1], left[2]*right[2] );
}

inline vl::Vec3d min( const vl::Vec3d& left, const vl::Vec3d& right )
{
	return vl::Vec3d( 
		std::min<double>(left[0], right[0]), 
		std::min<double>(left[1], right[1]), 
		std::min<double>(left[2], right[2]) );
}

inline vl::Vec3d max( const vl::Vec3d& left, const vl::Vec3d& right )
{
	return vl::Vec3d( 
		std::max<double>(left[0], right[0]), 
		std::max<double>(left[1], right[1]), 
		std::max<double>(left[2], right[2]) );
}

inline vl::Vec3f min( const vl::Vec3f& left, const vl::Vec3f& right )
{
	return vl::Vec3f( 
		std::min<float>(left[0], right[0]), 
		std::min<float>(left[1], right[1]), 
		std::min<float>(left[2], right[2]) );
}

inline vl::Vec3f max( const vl::Vec3f& left, const vl::Vec3f& right )
{
	return vl::Vec3f( 
		std::max<float>(left[0], right[0]), 
		std::max<float>(left[1], right[1]), 
		std::max<float>(left[2], right[2]) );
}

inline vl::Vec2f strip( const vl::Vec3f& v )
{
	return vl::Vec2f( v[0], v[1] );
}

inline vl::Vec2d strip( const vl::Vec3d& v )
{
	return vl::Vec2d( v[0], v[1] );
}

inline vl::Vec3d strip( const vl::Vec4d& v )
{
	return vl::Vec3d( v[0], v[1], v[2] );
}

inline vl::Vec3f strip( const vl::Vec4f& v )
{
	return vl::Vec3f( v[0], v[1], v[2] );
}

inline vl::Vec2f toVec2f( const vl::Vec2d& v )
{
	return vl::Vec2f( v[0], v[1] );
}

inline vl::Vec2f toVec2f( const vl::Vec2f& v )
{
	return v;
}

inline vl::Vec2d toVec2d( const vl::Vec2f& v )
{
	return vl::Vec2d( v[0], v[1] );
}

inline vl::Vec2d toVec2d( const vl::Vec2d& v )
{
	return v;
}

inline vl::Vec3f toVec3f( const vl::Vec3d& v )
{
	return vl::Vec3f( v[0], v[1], v[2] );
}

inline vl::Vec3f toVec3f( const vl::Vec3f& v )
{
	return v;
}

inline vl::Vec3d toVec3d( const vl::Vec3f& v )
{
	return vl::Vec3d( v[0], v[1], v[2] );
}

inline vl::Vec3d toVec3d( const vl::Vec3d& v )
{
	return v;
}

inline vl::Vec3f toVec3f( const double* v )
{
	return vl::Vec3f( v[0], v[1], v[2] );
}

inline vl::Vec3f toVec3f( const float* v )
{
	return vl::Vec3f( v[0], v[1], v[2] );
}

inline vl::Vec3d toVec3d( const double* v )
{
	return vl::Vec3d( v[0], v[1], v[2] );
}

inline vl::Vec3d toVec3d( const float* v )
{
	return vl::Vec3d( v[0], v[1], v[2] );
}

inline vl::Vec4f toVec4f( const vl::Vec4d& v )
{
	return vl::Vec4f( v[0], v[1], v[2], v[3] );
}

inline vl::Vec4f toVec4f( const vl::Vec4f& v )
{
	return v;
}

inline vl::Vec4d toVec4d( const vl::Vec4f& v )
{
	return vl::Vec4d( v[0], v[1], v[2], v[3] );
}

inline vl::Vec4d toVec4d( const vl::Vec4d& v )
{
	return v;
}

inline vl::Mat4f toMat4f( const vl::Mat4d& v )
{
	vl::Mat4f result;
	for( vl::Int i = 0; i < 4; ++i )
		for( vl::Int j = 0; j < 4; ++j )
			result[i][j] = v[i][j];

	return result;
}

inline vl::Mat3f toMat3f( const vl::Mat3f& v )
{
	return v;
}

inline vl::Mat3f toMat3f( const vl::Mat4d& v )
{
	vl::Mat3f result;
	for( vl::Int i = 0; i < 3; ++i )
		for( vl::Int j = 0; j < 3; ++j )
			result[i][j] = v[i][j];

	return result;
}

inline vl::Mat3f toMat3f( const vl::Mat4f& v )
{
	vl::Mat3f result;
	for( vl::Int i = 0; i < 3; ++i )
		for( vl::Int j = 0; j < 3; ++j )
			result[i][j] = v[i][j];

	return result;
}

inline vl::Mat4d toMat4d( const vl::Mat4f& v )
{
	vl::Mat4d result;
	for( vl::Int i = 0; i < 4; ++i )
		for( vl::Int j = 0; j < 4; ++j )
			result[i][j] = v[i][j];

	return result;
}

inline vl::Mat4d toMat4d( const vl::Mat4d& v )
{
	return v;
}

inline vl::Mat3d toMat3d( const vl::Mat4d& v )
{
	vl::Mat3d result;
	for( vl::Int i = 0; i < 3; ++i )
		for( vl::Int j = 0; j < 3; ++j )
			result[i][j] = v[i][j];

	return result;
}

inline vl::Mat3d toMat3d( const vl::Mat3d& v )
{
	return v;
}

inline vl::Mat3d toMat3d( const vl::Mat4f& v )
{
	vl::Mat3d result;
	for( vl::Int i = 0; i < 3; ++i )
		for( vl::Int j = 0; j < 3; ++j )
			result[i][j] = v[i][j];

	return result;
}

}

#endif
