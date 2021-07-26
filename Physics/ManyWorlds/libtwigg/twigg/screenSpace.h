#ifndef __SCREENSPACE_H__
#define __SCREENSPACE_H__

#include "twigg/vlutil.h"

class ScreenSpaceConverter
{
public:
	virtual ~ScreenSpaceConverter() {}
	virtual vl::Vec3 toScreenSpace( const vl::Vec3& worldSpace ) const = 0;
	virtual vl::Vec3 toWorldSpace( const vl::Vec3& screenSpace ) const = 0;
};

#endif
