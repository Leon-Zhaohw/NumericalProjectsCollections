#ifndef __BULLETWRAPPERS_H__
#define __BULLETWRAPPERS_H__

#include "physicsFwd.h"

#ifdef USE_BULLET
inline btVector3 toBtVector3( const vl::Vec3f& v )
{
	return btVector3(v[0], v[1], v[2]);
}

inline btVector3 toBtVector3( const vl::Vec3d& v )
{
	return btVector3(v[0], v[1], v[2]);
}

inline btMatrix3x3 toBtMatrix3x3( const vl::Mat3f& mat )
{
	return btMatrix3x3( 
		mat[0][0], mat[0][1], mat[0][2],
		mat[1][0], mat[1][1], mat[1][2],
		mat[2][0], mat[2][1], mat[2][2] );
}

inline btMatrix3x3 toBtMatrix3x3( const vl::Mat3d& mat )
{
	return btMatrix3x3( 
		mat[0][0], mat[0][1], mat[0][2],
		mat[1][0], mat[1][1], mat[1][2],
		mat[2][0], mat[2][1], mat[2][2] );
}

inline btMatrix3x3 toBtMatrix3x3( const vl::Mat4f& mat )
{
	return btMatrix3x3( 
		mat[0][0], mat[0][1], mat[0][2],
		mat[1][0], mat[1][1], mat[1][2],
		mat[2][0], mat[2][1], mat[2][2] );
}

inline btMatrix3x3 toBtMatrix3x3( const vl::Mat4d& mat )
{
	return btMatrix3x3( 
		mat[0][0], mat[0][1], mat[0][2],
		mat[1][0], mat[1][1], mat[1][2],
		mat[2][0], mat[2][1], mat[2][2] );
}

inline btTransform toBtTransform( const vl::Mat4f& mat )
{
	return btTransform( toBtMatrix3x3(mat), 
			btVector3(mat[0][3], mat[1][3], mat[2][3]) );
}

inline btTransform toBtTransform( const vl::Mat4d& mat )
{
	return btTransform( toBtMatrix3x3(mat), 
			btVector3(mat[0][3], mat[1][3], mat[2][3]) );
}
#endif

#endif
