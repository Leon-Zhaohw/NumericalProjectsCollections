#ifndef __IMAGELIB_H__

#define __IMAGELIB_H__

#include "twigg/util.h"
#include "twigg/vlutil.h"

/*
#include <IL/il.h>
#include <IL/ilu.h>


class ILImageWrapper
{
public:
	ILImageWrapper();
	~ILImageWrapper();

	ILuint imageId() const
	{
		return image_;
	}

private:
	ILuint image_;
};
*/

class Image
{
public:
	Image( const std::string& filename );

    // Retrieve the value stored in a physical location
    // (with integer coordinates) in the bitmap.
    // Should be called from getMappedValue in order to
    // do bilinear interpolation.
    TinyVec<unsigned char, 3> getPixelAt( int x, int y ) const;

	const unsigned char* data() const	{ return &data_[0]; }
	unsigned int width() const			{ return width_; }
	unsigned int height() const			{ return height_; }
	std::string filename() const		{ return filename_; }
	TinyVec<unsigned char, 3> averageColor() const
	{
		return avgColor_;
	}

private:
    std::string filename_;
    int width_;
    int height_;
	std::vector<unsigned char> data_;

	TinyVec<unsigned char, 3> avgColor_;
};

class TextureMap
	: public Image
{
public:
	TextureMap( const std::string& filename );

	vl::Vec3 normalizedPixel( int x, int y ) const;

    // Return the mapped value; here the coordinate
    // is assumed to be within the parametrization space:
    // [0, 1] x [0, 1]
    // (i.e., {(u, v): 0 <= u <= 1 and 0 <= v <= 1}
	vl::Vec3 getMappedValue( const vl::Vec2& coord ) const;
};

inline double toLuminance( const vl::Vec3& value )
{
	return (0.299 * value[0]) + (0.587 * value[1]) + (0.114 * value[2]);
}

#endif


