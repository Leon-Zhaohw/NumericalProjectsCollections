#include "stdafx.h"
#include "twigg/imagelib.h"

#include <boost/scoped_array.hpp>
#include <sstream>

#ifdef USE_MAGICK
#include <Magick++.h>
#else
#include <wx/image.h>
#endif

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


Image::Image( const std::string& filename )
	: filename_(filename)
{
/*
	ILImageWrapper ilImageWrapper;

	ilBindImage( ilImageWrapper.imageId() );
	std::cout << "Loading image '" << filename << "'...  " << std::flush;
	boost::scoped_array<char> fn( new char[filename.size() + 1] );
	std::copy( filename.begin(), filename.end(), fn.get() );
	fn[filename.size()] = 0;
	ILboolean result = ilLoadImage( fn.get() );

	ILenum err = ilGetError();
	if( err != IL_NO_ERROR )
	{
		std::ostringstream errMsg;
		errMsg << "Unable to load image '" << filename << "': ";
		do
		{
			errMsg << "     " << iluErrorString( err ) << std::endl;
			err = ilGetError();
		} while( err != IL_NO_ERROR );

		throw CreationException( errMsg.str() );
	}
	if( !result )
	{
		std::ostringstream errMsg;
		errMsg << "Unable to load image '" << filename << "': Unknown error.";
		throw CreationException( errMsg.str() );
	}

	width_ = ilGetInteger( IL_IMAGE_WIDTH );
	height_ = ilGetInteger( IL_IMAGE_HEIGHT );
	data_.resize( 3*width_*height_ );
	ilCopyPixels( 
			0,				// x offset
			0,				// y offset
			0,				// z offset
			width_,		
			height_,
			1,				// image depth
			IL_RGB,
			IL_UNSIGNED_BYTE,
			&data_[0] );

	std::cout << "Done." << std::endl;
*/
#ifdef USE_MAGICK
	try
	{
		std::cout << "Loading image '" << filename << "'...  " << std::flush;
		Magick::Image image( filename.c_str() );
		image.flip();

		Magick::Geometry geom = image.size();
		width_ = geom.width();
		height_ = geom.height();
		data_.resize( 3*width_*height_ );
		image.write( 0, 0, geom.width(), geom.height(), "RGB", Magick::CharPixel, &data_[0] );
		std::cout << "Done." << std::endl;

		// pull out the average color:
		TinyVec<unsigned long long, 3> avgColor(0, 0, 0);
		for( unsigned int i = 0; i < data_.size(); i += 3 )
			for( unsigned int iColor = 0; iColor < 3; ++iColor )
				avgColor[iColor] += data_[i + iColor];

		for( unsigned int i = 0; i < 3; ++i )
			avgColor[i] /= (data_.size()/3);

		avgColor_ = TinyVec<unsigned char, 3>( avgColor[0], avgColor[1], avgColor[2] );
	}
	catch( Magick::Exception& e )
	{
		throw CreationException( std::string("Error reading image file '") + filename + ":" + e.what() );
	}
#else
	wxImage image;
	if( !image.LoadFile( wxT(filename.c_str()) ) )
		throw CreationException( "Unable to load image '" + filename + "'." );

	image = image.Mirror(false);

	unsigned char* imageData = image.GetData();

	this->width_ = image.GetWidth();
	this->height_ = image.GetHeight();

	this->data_.resize( 3*this->width_*this->height_ );
	std::copy( imageData, imageData + data_.size(), 
		data_.begin() );
#endif
}


TinyVec<unsigned char, 3> Image::getPixelAt( int x, int y ) const
{
	x = std::min<int>( std::max<int>( x, 0 ), width() - 1 );
	y = std::min<int>( std::max<int>( y, 0 ), height() - 1 );

    // Find the position in the big data array...
    int pos = (y * width() + x) * 3;
	return TinyVec<unsigned char, 3>( data_[pos], data_[pos+1], data_[pos+2] );
}

TextureMap::TextureMap( const std::string& filename )
	: Image(filename)
{
}

vl::Vec3 TextureMap::getMappedValue( const vl::Vec2& coord ) const
{
	// YOUR CODE HERE

    // In order to add texture mapping support to the 
    // raytracer, you need to implement this function.
    // What this function should do is convert from
    // parametric space which is the unit square
    // [0, 1] x [0, 1] in 2-space to bitmap coordinates,
    // and use these to perform bilinear interpolation
    // of the values.

    float xCoord = (coord[0] * static_cast<double>(width()));
    float yCoord = (coord[1] * static_cast<double>(height()));

    float xDiff = xCoord - floor(xCoord);
    float yDiff = yCoord - floor(yCoord);

    int x = static_cast<int>( xCoord );
    int y = static_cast<int>( yCoord );

    // Should be bilinear interpolation here:
    return (xDiff * yDiff * normalizedPixel(x+1, y+1))
      + ((1-xDiff) * yDiff * normalizedPixel(x, y+1))
      + (xDiff * (1-yDiff) * normalizedPixel(x+1, y))
      + ((1-xDiff) * (1-yDiff) * normalizedPixel(x, y));
}

vl::Vec3 TextureMap::normalizedPixel( int x, int y ) const
{
	TinyVec<unsigned char, 3> val = getPixelAt(x, y);
	return vl::Vec3( 
		static_cast<double>(val[0]) / 255.0,
		static_cast<double>(val[1]) / 255.0,
		static_cast<double>(val[2]) / 255.0 );
}


/*
ILImageWrapper::ILImageWrapper()
{
	{
		static bool il_initialized = false;
		if( !il_initialized )
		{
			ilInit();
			iluInit();
			il_initialized = true;
		}
	}

	ILuint images[1];
	ilGenImages( 1, images );
	image_ = images[0];
}

ILImageWrapper::~ILImageWrapper()
{
	ILuint images[1] = { image_ };
	ilDeleteImages( 1, images );
}
*/
