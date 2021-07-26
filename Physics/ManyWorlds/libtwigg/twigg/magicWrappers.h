#ifndef __MAGICWRAPPERS_H__

#define __MAGICWRAPPERS_H__

#include "twigg/util.h"
#include "Wm4Triangle3.h"
#include "Wm4Segment3.h"
#include "Wm4Segment2.h"

inline Wm4::Vector3d toMagicVec( const vl::Vec3d& v )
{
	return Wm4::Vector3d( v[0], v[1], v[2] );
}

inline Wm4::Vector2d toMagicVec( const vl::Vec2d& v )
{
	return Wm4::Vector2d( v[0], v[1] );
}

inline Wm4::Vector3f toMagicVec( const vl::Vec3f& v )
{
	return Wm4::Vector3f( v[0], v[1], v[2] );
}

inline Wm4::Vector2f toMagicVec( const vl::Vec2f& v )
{
	return Wm4::Vector2f( v[0], v[1] );
}


#endif

