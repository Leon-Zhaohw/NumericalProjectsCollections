#include "stdafx.h"
#include "renderLines.h"

#include "twigg/ioutil.h"
#include "twigg/linalg.h"

#include <cairo.h>

#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/utility.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <iomanip>
#include <deque>
#include <sstream>

namespace planning {

class PathWriter
	: boost::noncopyable
{
public:
	PathWriter( boost::shared_ptr<cairo_surface_t> surface, const vl::Mat4f& transform )
		: transform_(transform), surface_( surface )
	{
		width_ = cairo_image_surface_get_width( surface_.get() );
		height_ = cairo_image_surface_get_height( surface_.get() );

		cr_ = cairo_create( surface_.get() );

		cairo_set_source_rgba(cr_, 0.0, 0.0, 0.0, 1.0);
		cairo_paint( cr_ );

		/*
		cairo_set_font_size (cr_, 10.0);

		drawCoordinateText( 10.0, 10.0 );
		drawCoordinateText( 10.0, height_ - 10.0 );
		drawCoordinateText( width_ - 50.0, height_ - 10.0 );
		drawCoordinateText( width_ - 50.0, 10.0 );
		drawCoordinateText( 0.5*width_, 0.5*height_ );
		*/

		assert( cairo_status(cr_) == CAIRO_STATUS_SUCCESS );
	}

	~PathWriter()
	{
		cairo_destroy(cr_);
	}

	void writePath( const std::vector<vl::Vec3f>& path, const vl::Vec4f& color, float thickness )
	{
		assert( thickness > 0.0f );

		if( path.empty() )
			return;

		cairo_new_path( cr_ );

		vl::Vec2d firstPoint = toLocalCoords( path.front() );
		cairo_move_to (cr_, firstPoint[0], firstPoint[1]);
		for( std::vector<vl::Vec3f>::const_iterator iter = path.begin() + 1;
			iter != path.end(); ++iter )
		{
			vl::Vec2d curPoint = toLocalCoords( *iter );
			cairo_line_to (cr_, curPoint[0], curPoint[1]);
		}

		cairo_set_line_width( cr_, thickness );
		cairo_set_source_rgba(cr_, color[0], color[1], color[2], color[3]);

		cairo_stroke (cr_);
		assert( cairo_status(cr_) == CAIRO_STATUS_SUCCESS );
	}

private:
	void drawCoordinateText( double x, double y )
	{
		std::ostringstream oss;
		oss << "(" << x << ", " << y << ")";
		std::string str = oss.str();

		cairo_move_to (cr_, x, y);
		cairo_show_text (cr_, str.c_str());
	}

	vl::Vec2d toLocalCoords( const vl::Vec3f& point )
	{
		vl::Vec3f transformed = vl::xform( transform_, point );
		vl::Vec2d result( (0.5*transformed[0] + 0.5)*static_cast<double>(width_),
			(0.5 - 0.5*transformed[1])*static_cast<double>(height_) );
		return result;
	}
	
	vl::Mat4f transform_;

	int width_;
	int height_;

	boost::shared_ptr<cairo_surface_t> surface_;
	cairo_t* cr_;
};

void renderLines( 
	const std::string& linesFilename,
	const vl::Mat4f& transformMatrix,
	const std::string outFilePrefix,
	size_t width,
	size_t height )
{
	std::ifstream ifs( linesFilename.c_str(), std::ios::binary );
	if( !ifs )
		throw IOException( std::string("Unable to open file '") + linesFilename + "' for reading." );

	std::ifstream::pos_type fileSize = filelen( ifs );

	typedef std::vector<vl::Vec3f> Path;
	typedef std::deque< Path > PathList;

	cairo_format_t format = CAIRO_FORMAT_ARGB32;

	typedef std::pair<PathList, size_t> PathWithTime;
	std::list< PathWithTime > activePaths;

	size_t iFrame = 0;

	size_t lengthBetweenPaths = 3;
	size_t lengthOfFade = 10;
	double fadePerFrame = 0.01;

	float thicknesses[] =  { 0.5, 1.0, 2.0, 3.0, 2.0, 1.5, 1.0, 1.0, 1.0, 1.0 };
	float colorWeights[] = { 0.8, 0.9, 1.0, 1.0, 1.0, 0.9, 0.7, 0.5, 0.3, 0.1 };
	assert( sizeof(thicknesses) == lengthOfFade * sizeof(float) );
	assert( sizeof(colorWeights) == lengthOfFade * sizeof(float) );

	vl::Vec3f initialColor( 1.0f, 1.0f, 0.5f );
	vl::Vec3f fadedColor( 0.9f, 0.9f, 0.9f );

	while( !activePaths.empty() || ifs.tellg() < fileSize )
	{
		if( (activePaths.empty() || activePaths.back().second > lengthBetweenPaths)
			&& ifs.tellg() < fileSize )
		{
			int numPaths;
			ifs.read( reinterpret_cast<char*>( &numPaths ), sizeof(int) );

			// introduce a new path
			activePaths.push_back( PathWithTime(PathList(), 0) );
			assert( numPaths > 0 );
			PathList& currentPathList = activePaths.back().first;

			for( int i = 0; i < numPaths; ++i )
			{
				currentPathList.push_back( std::vector<vl::Vec3f>() );
				readArray( currentPathList.back(), ifs );
			}
		}

		// start with background plate
		boost::shared_ptr<cairo_surface_t> currentSurface(
			cairo_image_surface_create(format, width, height),
			&cairo_surface_destroy );

		{
			PathWriter pathWriter( currentSurface, transformMatrix );

			for( std::list< PathWithTime >::iterator currentPathItr = activePaths.begin();
				currentPathItr != activePaths.end(); )
			{
				if( currentPathItr->second >= lengthOfFade )
				{
					size_t numSteps = currentPathItr->second - lengthOfFade;
					double alpha = 1.0 - static_cast<double>( numSteps )*fadePerFrame;
					if( alpha <= 0.0 )
					{
						currentPathItr = activePaths.erase( currentPathItr );
					}
					else
					{
						float thickness = thicknesses[lengthOfFade - 1];
						const PathList& paths = currentPathItr->first;
						std::for_each( paths.begin(), paths.end(), 
							boost::bind( &PathWriter::writePath, boost::ref(pathWriter), _1, vl::Vec4f(fadedColor, alpha), thickness ) );

						++currentPathItr->second;
						++currentPathItr;
					}
				}
				else
				{
					float colorWeight = colorWeights[currentPathItr->second];
					float thickness = thicknesses[currentPathItr->second];
					vl::Vec3f color = colorWeight*initialColor + (1.0f - colorWeight)*fadedColor;
					const PathList& paths = currentPathItr->first;
					std::for_each( paths.begin(), paths.end(), 
						boost::bind( &PathWriter::writePath, boost::ref(pathWriter), _1, vl::Vec4f(color, 1.0), thickness ) );

					++currentPathItr->second;
					++currentPathItr;
				}
			}
		}

		{
			std::ostringstream oss;
			oss << outFilePrefix << "." << std::setw(4) << std::setfill('0') << iFrame << ".png";
			std::string outFilename = oss.str();
			cairo_surface_write_to_png( currentSurface.get(), outFilename.c_str() );
		}

		++iFrame;
	}
}

} // namespace planning
