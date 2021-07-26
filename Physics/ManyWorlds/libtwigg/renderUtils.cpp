#include "stdafx.h"

#ifndef CONDOR
#include "twigg/renderUtils.h"
#include "twigg/imagelib.h"

#ifdef __APPLE__
#include <OpenGL/glu.h>
#include <OpenGL/glut.h>
#else
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#ifdef USE_FTGL
#ifdef __linux__
	#define FONT_FILE "/nfs/hn00/cdtwigg/share/fonts/arial.ttf"
#endif
#ifdef __APPLE_CC__
	#define FONT_FILE "Something reasonable should go here"
#endif
#ifdef WIN32
	#define FONT_FILE "C:\\Windows\\Fonts\\arial.ttf"
#endif

#include "FTGL/FTGLPixmapFont.h"
#endif

#include <boost/numeric/conversion/bounds.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <sstream>
#include <iostream>
#include <iomanip>


#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#include "twigg/imagelib.h"

class BitmapFont
{
public:
	BitmapFont( void* font )
		: font_(font)
	{
		if( font_ == GLUT_BITMAP_8_BY_13 )
			height_ = 13;
		else if( font_ == GLUT_BITMAP_9_BY_15 )
			height_ = 15;
		else if( font_ == GLUT_BITMAP_TIMES_ROMAN_10 )
			height_ = 10;
		else if( font_ == GLUT_BITMAP_TIMES_ROMAN_24 )
			height_ = 24;
		else if( font_ == GLUT_BITMAP_HELVETICA_10 )
			height_ = 10;
		else if( font_ == GLUT_BITMAP_HELVETICA_12 )
			height_ = 12;
		else if( font_ == GLUT_BITMAP_HELVETICA_18 )
			height_ = 18;
		else
			assert( false );
	}

	void render( const char* s )
	{
		while( *s )
			glutBitmapCharacter( font_, *s++ );
	}


	vl::Vec2d size( const char* s )
	{
		vl::Vec2d result( 0.0, boost::numeric_cast<double>(height_) );
		while( *s )
			result[0] += glutBitmapWidth( font_, *s++);

		return result;
	}

private:
	void* font_;
	unsigned int height_;
};


/*
class BitmapFont
{
public:
	BitmapFont()
	{
		base_ = glGenLists(256);							// Storage For all 256 Characters

		Image image( "g:/users/cdtwigg/systemFont.bmp" );

		charWidth_ = image.width()/16;
		charHeight_ = image.height()/16;

		unsigned int totalSize = charWidth_*charHeight_;
		for (unsigned int iChar = 0; iChar < 256; ++iChar )
		{
			unsigned int cx = iChar/16;
			unsigned int cy = iChar%16;

			std::vector<unsigned char> bitmap( totalSize/8 + 1, 0 );
			for( unsigned int i = 0; i < charHeight_; ++i )
			{
				for( unsigned int j = 0; j < charWidth_; ++j )
				{
					unsigned int pos = i*charWidth_ + j;
					if( image.getPixelAt( cx*charHeight_+i, cy*charWidth_+j )[0] )
						bitmap[ pos/8 ] |= (1 << pos%8);
				}
			}

			glNewList(base_+iChar, GL_COMPILE);	// Start Building A List
			glBitmap(charWidth_, charHeight_, 0.0, 2.0, 10.0, 0.0, &bitmap[0]);
			glEndList();								// Done Building The Display List
		}
	}

	~BitmapFont()
	{
		glDeleteLists(base_, 256);							// Delete All 96 Characters
	}

	void render( const char* s )
	{
		glPushAttrib(GL_LIST_BIT);
		glListBase(base_);
		glCallLists(strlen(s), GL_UNSIGNED_BYTE, (GLubyte *) s);
		glPopAttrib ();
	}

	vl::Vec2d size( const char* s )
	{
		return vl::Vec2d( strlen(s)*charWidth_, charHeight_ );
	}

private:
	GLuint base_;

	unsigned int charWidth_;
	unsigned int charHeight_;
};
*/

/*
GLuint makeRasterFont(void)
{
   GLuint i, j;
   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

   GLuint fontOffset = glGenLists (128);
   for (i = 0,j = 'A'; i < 26; i++,j++) {
      glNewList(fontOffset + j, GL_COMPILE);
      glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, letters[i]);
      glEndList();
   }

   for (i = 0,j = 'a'; i < 26; i++,j++) {
      glNewList(fontOffset + j, GL_COMPILE);
      glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, letters[i]);
      glEndList();
   }

   glNewList(fontOffset + ' ', GL_COMPILE);
   glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, space);
   glEndList();

   return fontOffset;
}
*/

void checkGLError(const std::string& location)
{
	return;
#ifndef _DEBUG
	return;
#else
    GLenum err = glGetError();
    if (err != GL_NO_ERROR)
	{
		std::ostringstream oss;
		oss << "Caught OpenGL error " << gluErrorString(err) << " at location " << location;
		throw OpenGLException( oss.str() );
	}
#endif
}

BoundingBox3d printString(const char *s, vl::Vec3 position, HAlign horiz, VAlign vert, bool output)
{
	if( s == 0 )
		return BoundingBox3d(position, position);

	size_t sLen = strlen(s);

#ifdef USE_FTGL
	static FTGLPixmapFont font(FONT_FILE);
	static bool initialized = false;
	if( font.Error() )
	{
		std::cerr << "Failed to load font file '" << FONT_FILE << "'." << std::endl;
		return BoundingBox3d(position, position);
	}
	
	if( !initialized )
	{
		font.CharMap(ft_encoding_unicode);
		int pointSize = 12;
		if (!font.FaceSize(pointSize))
			std::cerr << "ERROR: Unable to set font face size " << pointSize << "\n";
		initialized = true;
	}
#else
	static BitmapFont font( GLUT_BITMAP_HELVETICA_10);
#endif

	OpenGLScreenSpaceConverter converter;
	vl::Vec3 screenPos = converter.toScreenSpace( position );

#ifdef USE_FTGL
	vl::Vec2 size( font.Advance( s ), font.Ascender() + font.Descender() );
#else
	// width of the default font
	vl::Vec2 size = font.size(s);
#endif
	{
		// horizontal alignment
		if( ALIGN_HORIZ_RIGHT == horiz )
			screenPos[0] -= size[0];
		else if( ALIGN_HORIZ_CENTER == horiz )
			screenPos[0] -= 0.5*size[0];
	}

#ifdef USE_FTGL
	{
		if( ALIGN_VERT_TOP == vert )
			screenPos[1] -= font.Ascender();
		else if( ALIGN_VERT_CENTER == vert )
			screenPos[1] -= 0.5*font.Ascender();
		else if( ALIGN_VERT_BOTTOM == vert )
			screenPos[1] += font.Descender();
	}
#else
	{
		if( ALIGN_VERT_TOP == vert )
			screenPos[1] -= size[1];
		else if( ALIGN_VERT_CENTER == vert )
			screenPos[1] -= 0.5*size[1];
		else if( ALIGN_VERT_BOTTOM == vert )
			screenPos[1] += 0;
	}
#endif

	position = converter.toWorldSpace( screenPos );

	if( output )
	{
		glRasterPos3d(position[0], position[1], position[2]);
		{
			glDisable( GL_LIGHTING );
			glShadeModel(GL_SMOOTH);
			GLDisableHandler depthTest( GL_DEPTH_TEST);
#ifdef USE_FTGL
			font.Render( s );
#else
			font.render( s );
#endif
		}
	}

	BoundingBox3d screenBox;
	screenBox.expand( screenPos );
	screenBox.expand( screenPos + vl::Vec3(size[0], 0.0, 0.0) );
#ifdef USE_FTGL
	screenBox.expand( screenPos - vl::Vec3(0.0, font.Descender(), 0.0) );
	screenBox.expand( screenPos + vl::Vec3(0.0, font.Ascender(), 0.0) );
#else
	screenBox.expand( screenPos - vl::Vec3(0.0, 0.0, 0.0) );
	screenBox.expand( screenPos + vl::Vec3(0.0, size[1], 0.0) );
#endif

	BoundingBox3d result;
	result.expand( converter.toWorldSpace(screenBox.minimum()) );
	result.expand( converter.toWorldSpace(screenBox.maximum()) );
	return result;
}

void printStrings( const vl::Vec3& upperLeft, const std::deque< std::pair<std::string, std::string> >& items )
{
	if( items.empty() )
		return;

	typedef std::pair<std::string, std::string> Item;

	GLOverlayHandler overlay;

	glDisable(GL_LIGHTING);
	GLDisableHandler depthTest( GL_DEPTH_TEST );
	glColor4ub(255, 255, 255, 255);

	double left = upperLeft[0];
	double top = upperLeft[1];

	vl::Vec3 maxLeftBox( vl::vl_0 );
	vl::Vec3 maxRightBox( vl::vl_0 );
	double ySkip;

	// first, figure out metrics:
	for( std::deque<Item>::const_iterator itemItr = items.begin();
		itemItr != items.end();
		++itemItr )
	{
		BoundingBox3d leftBox = printString( itemItr->first.c_str(), vl::Vec3(left, top, 0.0),
			ALIGN_HORIZ_LEFT, ALIGN_VERT_TOP, false );
		BoundingBox3d rightBox = printString( itemItr->second.c_str(), vl::Vec3(left + maxLeftBox[0], top, 0.0),
			ALIGN_HORIZ_RIGHT, ALIGN_VERT_TOP, false );

		maxLeftBox = vl::max( maxLeftBox, leftBox.maximum() - leftBox.minimum() );
		maxRightBox = vl::max( maxRightBox, rightBox.maximum() - rightBox.minimum() );
		ySkip = 1.2 * std::max( maxLeftBox[1], maxRightBox[1] );
		top -= ySkip;
	}

	// now, actually render:
	top = upperLeft[1];
	double right = left + 1.3*(maxLeftBox[0] + maxRightBox[0]);
	double bottom = top - ySkip*items.size();

	for( std::deque<Item>::const_iterator itemItr = items.begin();
		itemItr != items.end();
		++itemItr )
	{
		printString( itemItr->first.c_str(), vl::Vec3(left, top, 0.0),
			ALIGN_HORIZ_LEFT, ALIGN_VERT_TOP );
		printString( itemItr->second.c_str(), vl::Vec3(right, top, 0.0),
			ALIGN_HORIZ_RIGHT, ALIGN_VERT_TOP );
		top -= ySkip;
	}
}

GLMatrixStackHandler::GLMatrixStackHandler()
{
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
}

GLMatrixStackHandler::~GLMatrixStackHandler()
{
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

GLActionHandler::GLActionHandler( GLenum mode )
{
	glBegin( mode );
}

GLActionHandler::~GLActionHandler()
{
	glEnd();
}

GLEnableHandler::GLEnableHandler(GLenum capability)
	: capability_( capability )
{
	glEnable(capability_);
}

GLEnableHandler::~GLEnableHandler()
{
	glDisable(capability_);
}

GLDisableHandler::GLDisableHandler(GLenum capability)
	: capability_( capability )
{
	glDisable(capability_);
}

GLDisableHandler::~GLDisableHandler()
{
	glEnable(capability_);
}

GLTexture::GLTexture()
{
	glGenTextures(1, &glTextureName_);
}

GLTexture::~GLTexture()
{
	glDeleteTextures(1, &glTextureName_);
}

GLuint GLTexture::textureName() const
{
	return glTextureName_;
}

GLBindTextureHandler::GLBindTextureHandler(const GLTexture& texture)
{
	const GLfloat Kd[] = { 1.0, 1.0, 1.0, 1.0f };
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, Kd);
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, texture.textureName());
}

GLBindTextureHandler::~GLBindTextureHandler()
{
	glDisable(GL_TEXTURE_2D);
}

GLEnableClientStateHandler::GLEnableClientStateHandler(GLenum array)
	: array_(array)
{
	glEnableClientState(array_);
}

GLEnableClientStateHandler::~GLEnableClientStateHandler()
{
	glDisableClientState( array_ );
}

GLBillboardHandler::GLBillboardHandler( const vl::Vec3& position )
{
	float modelview[16];
	glTranslated( position[0], position[1], position[2] );
	glGetFloatv(GL_MODELVIEW_MATRIX, modelview);

	// undo all rotations/scales
	for( unsigned int i=0; i<3; i++ )
		for( unsigned int j=0; j<3; j++ ) {
			if ( i==j )
				modelview[i*4+j] = 1.0;
			else
				modelview[i*4+j] = 0.0;
		}

	glLoadMatrixf(modelview);
}

template <typename RealType, typename Vec3Type, typename Mat3Type>
Mat3Type buildBasisImp( const Vec3Type& up )
{
	Mat3Type result;

	result[1] = vl::norm(up);
	double maxCross = 0;
    for( vl::Int iAxis = 0; iAxis < 3; ++iAxis )
	{
		Vec3Type axisVec = vl::vl_0;
		axisVec[iAxis] = 1.0;
		Vec3Type crossProd = cross( result[1], axisVec );
		RealType crossProdLen = len( crossProd );
		if( crossProdLen > maxCross )
		{
			maxCross = crossProdLen;
			result[0] = crossProd;
			result[0] /= crossProdLen;
		}
	}
	result[2] = vl::norm( cross( result[0], result[1] ) );
	result[0] = vl::norm( cross( result[1], result[2] ) );

	return result;
}

vl::Mat3d buildBasis( const vl::Vec3d& up ) { return buildBasisImp<double, vl::Vec3d, vl::Mat3d>( up ); }
vl::Mat3f buildBasis( const vl::Vec3f& up ) { return buildBasisImp<float,  vl::Vec3f, vl::Mat3f>( up ); }


void drawShadedArrow( const vl::Vec3d& position, 
					 const vl::Vec3d& direction, 
					 double scale, 
					 const vl::Vec4f& color,
					 bool shadedColumn )
{
    drawShadedArrow( toVec3f( position ), toVec3f( direction ), scale, color, shadedColumn );
}

void drawShadedArrow( const vl::Vec3f& position, 
					 const vl::Vec3f& direction, 
					 float scale, 
					 const vl::Vec4f& color,
					 bool shadedColumn )
{
	vl::Vec3f end( position + scale*direction );

	if( len(direction) < 1e-7 )
		return;

	vl::Mat3f basis = buildBasis( direction );


	glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, color.Ref() );

	{
		const int triCount = 40;
		const float length = 0.25 * scale;
		const float radius = 0.08 * scale;

		vl::Vec3f capPos = end - length*basis[1];

		if( shadedColumn )
		{
			glEnable( GL_LIGHTING );
			float innerRadius = 0.03 * scale;
			GLActionHandler handler( GL_TRIANGLE_STRIP );
			{
                float angleStep = (2.0f * 3.1415926535f) / static_cast<float>(triCount);
                float angle = 0.0f;
				for( int iCircle = 0; iCircle <= triCount; iCircle += 1, angle += angleStep )
				{
					vl::Vec3f dir = cos(angle)*basis[0] + sin(angle)*basis[2];

					glNormal3fv( norm(dir).Ref() );
					glVertex3fv( (position + innerRadius*dir).Ref() );
					glVertex3fv( (capPos + innerRadius*dir).Ref() );
				}				
			}
		}
		else
		{
			glDisable( GL_LIGHTING );
			glColor4fv( color.Ref() );
			GLActionHandler lines( GL_LINES );
			glVertex3fv( position.Ref() );
			glVertex3fv( capPos.Ref() );
		}

		glEnable( GL_LIGHTING );

		
		vl::Vec3f tip = end;
		{
			GLActionHandler handler( GL_TRIANGLES );
			vl::Vec3f prevDir = basis[0];
            float angle = 0.0f;
            float angleStep = (2.0f*3.14159265f) / static_cast<float>(triCount);
			for( int iCircle = 0; iCircle <= triCount; iCircle += 1, angle += angleStep )
			{
				vl::Vec3f dir = cos(angle)*basis[0] + sin(angle)*basis[2];

				glNormal3fv( norm(prevDir).Ref() );
				glVertex3fv( (capPos + radius*prevDir).Ref() );

				glNormal3fv( norm(0.5*(prevDir + dir)).Ref() );
				glVertex3fv( tip.Ref() );

				glNormal3fv( norm(dir).Ref() );
				glVertex3fv( (capPos + radius*dir).Ref() );

				prevDir = dir;
			}
		}

		{
			glNormal3fv( (-basis[1]).Ref() );
			GLActionHandler handler( GL_TRIANGLE_FAN );
			glVertex3fv( capPos.Ref() );

            float angle = 0.0f;
            float angleStep = (2.0f*3.14159265f) / static_cast<float>( triCount );
			for( int iCircle = 0; iCircle <= triCount; iCircle += 1, angle += angleStep )
			{
				vl::Vec3f dir = cos(angle)*basis[0] + sin(angle)*basis[2];
				glVertex3fv( (capPos + radius*dir).Ref() );
			}
		}
	}

}


void drawArrow( const vl::Vec3& position, const vl::Vec3& direction, double scale, const vl::Vec4f& color )
{
	vl::Vec3 end( position + scale*direction );

	{
		glDisable( GL_LIGHTING );
		GLActionHandler handler( GL_LINES );
		glColor4fv( color.Ref() );
		glVertex3dv( position.Ref() );
		glVertex3dv( end.Ref() );
	}

	if( len(direction) > 1e-7 )
	{
		vl::Vec3 basis1 = vl::norm(direction);
		vl::Vec3 basis2;
		double maxCross = 0;
		for( int iAxis = 0; iAxis < 3; ++iAxis )
		{
			vl::Vec3 axisVec = vl::vl_0;
			axisVec[iAxis] = 1.0;
			if( len( cross( basis1, axisVec ) ) > maxCross )
			{
				maxCross = len( cross( basis1, axisVec ) );
				basis2 = vl::norm( cross( basis1, axisVec ) );
			}
		}
		vl::Vec3 basis3 = vl::norm( cross( basis1, basis2 ) );
		basis2 = vl::norm( cross( basis1, basis3 ) );

		/*
		vl::Mat3 rotMat(
			basis1[0], basis1[1], basis1[2],
			basis2[0], basis2[1], basis2[2],
			basis3[0], basis3[1], basis3[2] );
		rotMat = trans(rotMat);
		Real matrix[] = 
			{	rotMat[0][0], rotMat[1][0], rotMat[2][0], 0.0,
				rotMat[0][1], rotMat[1][1], rotMat[2][1], 0.0,
				rotMat[0][2], rotMat[1][2], rotMat[2][2], 0.0,
				0.0, 0.0, 0.0, 1.0 };
		glMultMatrix( matrix );
		*/

		{
			const int triCount = 10;
			const double length = 0.35 * scale;
			const double radius = 0.1 * scale;

			vl::Vec3 capPos = end;
			{
				GLActionHandler handler( GL_TRIANGLE_FAN );
				glNormal3dv( direction.Ref() );
				glVertex3dv( (end + length*basis1).Ref() );
				for( int iCircle = 0; iCircle <= triCount; ++iCircle )
				{
					double angle = 2.0 * 3.1415926535 * 
						static_cast<double>(iCircle)/static_cast<double>(triCount);
					vl::Vec3 dir = cos(angle)*basis2 + sin(angle)*basis3;
					glVertex3dv( (capPos + radius*dir).Ref() );
				}
			}

			{
				glNormal3dv( (-basis1).Ref() );
				GLActionHandler handler( GL_TRIANGLE_FAN );
				glVertex3dv( capPos.Ref() );
				for( int iCircle = 0; iCircle <= triCount; ++iCircle )
				{
					double angle = 2.0 * 3.1415926535 * 
						static_cast<double>(iCircle)/static_cast<double>(triCount);
					vl::Vec3 dir = cos(angle)*basis2 + sin(angle)*basis3;
					glVertex3dv( (capPos + radius*dir).Ref() );
				}
			}
		}
	}
}

vl::Vec3f jetColor(float x)
{
	if(x < 0.f) 
		return vl::Vec3f(0.f,0.f,0.f);
	else if (x < 0.125f) {
	    float a = x/0.125f;
		return vl::Vec3f(0.f, 0.f, 0.5f+0.5f*a);
	}
	else if (x < 0.375f) {
	    float a = (x - 0.125f)/0.25f;
	    return vl::Vec3f(0.f, a, 1.f);
	}
	else if (x < 0.625f) {
	    float a = (x - 0.375f)/0.25f;
	    return vl::Vec3f(a, 1.f, 1.f-a);
	}
	else if (x < 0.875f) {
	    float a = (x - 0.625f)/0.25f;
	    return vl::Vec3f(1.f, 1.f-a, 0.f);
	}
	else if (x <= 1.0f) {
	    float a = (x - 0.875f)/0.125f;
	    return vl::Vec3f(1.f-0.5f*a, 0.f, 0.f);
	}
	else {
	    return vl::Vec3f(1.f,1.f,1.f);
	}
}

void cube()
{
    glBegin( GL_QUADS );
    
    glNormal3d( 0.0, 0.0, -1.0 );
    glVertex3d( 0.0, 0.0, 0.0 ); glVertex3d( 0.0, 1.0, 0.0 );
    glVertex3d( 1.0, 1.0, 0.0 ); glVertex3d( 1.0, 0.0, 0.0 );
    
    glNormal3d( 0.0, -1.0, 0.0 );
    glVertex3d( 0.0, 0.0, 0.0 ); glVertex3d( 1.0, 0.0, 0.0 );
    glVertex3d( 1.0, 0.0, 1.0 ); glVertex3d( 0.0, 0.0, 1.0 );
    
    glNormal3d( -1.0, 0.0, 0.0 );
    glVertex3d( 0.0, 0.0, 0.0 ); glVertex3d( 0.0, 0.0, 1.0 );
    glVertex3d( 0.0, 1.0, 1.0 ); glVertex3d( 0.0, 1.0, 0.0 );
    
    glNormal3d( 0.0, 0.0, 1.0 );
    glVertex3d( 0.0, 0.0, 1.0 ); glVertex3d( 1.0, 0.0, 1.0 );
    glVertex3d( 1.0, 1.0, 1.0 ); glVertex3d( 0.0, 1.0, 1.0 );
    
    glNormal3d( 0.0, 1.0, 0.0 );
    glVertex3d( 0.0, 1.0, 0.0 ); glVertex3d( 0.0, 1.0, 1.0 );
    glVertex3d( 1.0, 1.0, 1.0 ); glVertex3d( 1.0, 1.0, 0.0 );
    
    glNormal3d( 1.0, 0.0, 0.0 );
    glVertex3d( 1.0, 0.0, 0.0 ); glVertex3d( 1.0, 1.0, 0.0 );
    glVertex3d( 1.0, 1.0, 1.0 ); glVertex3d( 1.0, 0.0, 1.0 );
    
    glEnd();
}

GLOverlayHandler::GLOverlayHandler()
{
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(-100.0, 100.0, -100.0, 100.0, -1.0, 1.0);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
}

GLOverlayHandler::~GLOverlayHandler()
{
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}


void legend(double min, double max)
{
	GLOverlayHandler overlay;

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	GLDisableHandler( GL_DEPTH_TEST );
	glDisable(GL_LIGHTING);

	glShadeModel(GL_SMOOTH);

	{
		GLActionHandler triangleStrip( GL_TRIANGLE_STRIP );
		for( int i = 0; i < 40; ++i )
		{
			glColor3fv( jetColor(static_cast<double>(i) / static_cast<double>(40)).Ref() );
			glVertex2d( -95.0, (95.0 - 40.0) + static_cast<double>(i) );
			glVertex2d( -90.0, (95.0 - 40.0) + static_cast<double>(i) );
		}
	}

	glColor4ub( 255, 255, 255, 255 );
	{
		std::ostringstream oss;
		oss << max;
		printString(oss.str().c_str(), vl::Vec3(-50, 95.0, 0.0));
	}

	{
		std::ostringstream oss;
		oss << min;
		printString(oss.str().c_str(), vl::Vec3(-50.0, 55.0, 0.0));
	}
}

OpenGLScreenSpaceConverter::OpenGLScreenSpaceConverter( const double mvmatrix[16], const double projmatrix[16], const int viewport[4] )
{
	// In case we need to store these beforehand.
	std::copy( mvmatrix, mvmatrix + 16, modelViewMatrix_.begin() );
	std::copy( projmatrix, projmatrix + 16, projectionMatrix_.begin() );
	std::copy( viewport, viewport + 4, viewPort_.begin() );
}

OpenGLScreenSpaceConverter::OpenGLScreenSpaceConverter()
{
	// Pull all the quantities directly from OpenGL
	glGetIntegerv(GL_VIEWPORT, &viewPort_[0]);
	glGetDoublev(GL_MODELVIEW_MATRIX, &modelViewMatrix_[0]);
	glGetDoublev(GL_PROJECTION_MATRIX, &projectionMatrix_[0]);
}

vl::Vec3 OpenGLScreenSpaceConverter::toScreenSpace( const vl::Vec3& worldSpace ) const
{
	vl::Vec3d projPoint;
	gluProject( worldSpace[0], worldSpace[1], worldSpace[2], 
		&modelViewMatrix_[0], &projectionMatrix_[0], &viewPort_[0], 
		&projPoint[0], &projPoint[1], &projPoint[2] );
	return projPoint;
}

vl::Vec3 OpenGLScreenSpaceConverter::toWorldSpace( const vl::Vec3& screenSpace ) const
{
	vl::Vec3 result;
	gluUnProject (screenSpace[0], screenSpace[1], screenSpace[2],
		&modelViewMatrix_[0], &projectionMatrix_[0], &viewPort_[0], 
		&(result[0]), &(result[1]), &(result[2]));
	return result;
}


ARBVertexProgram::BindProgram::BindProgram(const ARBVertexProgram& program)
	: program_(program)
{
	glEnable( GL_VERTEX_PROGRAM_ARB );
	glBindProgramARB( GL_VERTEX_PROGRAM_ARB, program.progID_ );
}

ARBVertexProgram::BindProgram::~BindProgram()
{
	glDisable( GL_VERTEX_PROGRAM_ARB );
}

ARBVertexProgram::ARBVertexProgram( const std::string& prog )
	: program_(prog)
{
	assert( glGetError() == GL_NO_ERROR );

	if( GLEW_ARB_vertex_program == 0 )
		throw OpenGLException( "Extension ARB_vertex_program not supported on this card." );

	// split the vertex program into lines
	std::string strippedProgram;
	strippedProgram.reserve( prog.size() );
	std::string::const_iterator iter = prog.begin();
	while( iter != prog.end() )
	{
		bool comment = (*iter == '#');
		while( iter != prog.end() )
		{
			if( !comment )
				strippedProgram.push_back( *iter );

			if( *iter == '\n' )
			{
				++iter;
				break;
			}
			else
				++iter;
		}
	}

	strippedProgram_ = strippedProgram;

	glGenProgramsARB( 1, &progID_ );

	{
		BindProgram bindProgram( *this );
        
		glProgramStringARB( GL_VERTEX_PROGRAM_ARB,
			GL_PROGRAM_FORMAT_ASCII_ARB,
			boost::numeric_cast<GLsizei>(strippedProgram.size()), strippedProgram.c_str() );

		if( GL_INVALID_OPERATION == glGetError() )
		{
			// Fix the error position
			GLint errPos;
			glGetIntegerv( GL_PROGRAM_ERROR_POSITION_ARB, &errPos );

			// Print implementation-dependent program errors and warnings string
			const GLubyte* errString = glGetString( GL_PROGRAM_ERROR_STRING_ARB );

			std::ostringstream oss;
			if( errPos >= 0 )
				oss << "Caught error at position: " << errPos << ": '" << strippedProgram.substr( errPos, 
					std::min<size_t>(errPos+10, prog.size()) ) << "'" << std::endl;

			oss << errString;
			throw OpenGLException( oss.str() );
		}
	}
}

ARBVertexProgram::~ARBVertexProgram()
{
	glDeleteProgramsARB( 1, &progID_ );
}

GLBindBuffer::~GLBindBuffer()
{
}

GLBuffer::~GLBuffer()
{
}

GLMemoryManager::~GLMemoryManager()
{
}

ARBVertexBufferObject::ARBVertexBufferObject(const char* data, unsigned int size, GLMemoryManager::BufferType type)
{
	if( GLEW_ARB_vertex_buffer_object == 0 )
		throw OpenGLException( "ARB_vertex_buffer_object not supported on this card." );

	if( type == GLMemoryManager::ARRAY )
		this->bufferType_ = GL_ARRAY_BUFFER_ARB;
	else
		this->bufferType_ = GL_ELEMENT_ARRAY_BUFFER_ARB;

	glGenBuffersARB( 1, &bufferId_ );
		checkGLError( "glGenBuffersARB" );
	BindBuffer bindBuffer( *this, bufferType_ );
		checkGLError( "bindBuffer" );
	glBufferDataARB( bufferType_, size, data, GL_STATIC_DRAW_ARB );
		checkGLError( "glBufferDataARB" );
}

#define BUFFER_OFFSET(i) ((char *)NULL + (i))

const char* ARBVertexBufferObject::buffer() const
{
	return BUFFER_OFFSET(0);
}

ARBVertexBufferObject::~ARBVertexBufferObject()
{
	glDeleteBuffersARB( 1, &bufferId_ );
}

std::auto_ptr<GLBindBuffer> ARBVertexBufferObject::bind() const
{
	return std::auto_ptr<GLBindBuffer>( 
		new ARBVertexBufferObject::BindBuffer( *this, this->bufferType_ ) );
}

ARBVertexBufferObject::BindBuffer::BindBuffer( const ARBVertexBufferObject& bufferObject, GLenum target )
	: target_(target)
{
	glBindBufferARB( target, bufferObject.bufferId_ );
}

ARBVertexBufferObject::BindBuffer::~BindBuffer()
{
	// This is apparently the correct way to disable the buffer.
	glBindBufferARB( target_, 0 );
}

GLBasicBuffer::GLBasicBuffer(const char* data, unsigned int size)
{
	data_.resize(size);
	std::copy( data, data+size, data_.begin() );
}

GLBasicBuffer::~GLBasicBuffer()
{
}

std::auto_ptr<GLBindBuffer> GLBasicBuffer::bind() const
{
	// Don't really need to bind these buiffers, so we return
	//   a placeholder.
	return std::auto_ptr<GLBindBuffer>( new GLBindBuffer() );
}

const char* GLBasicBuffer::buffer() const
{
	return &data_[0];
}

boost::shared_ptr<GLBuffer> VBOMemoryManager::buffer( const char* data, unsigned int size, GLMemoryManager::BufferType type )
{
	return boost::shared_ptr<GLBuffer>( new ARBVertexBufferObject(data, size, type) );
}

VBOMemoryManager::~VBOMemoryManager()
{
}

boost::shared_ptr<GLBuffer> GLBasicMemoryManager::buffer( const char* data, unsigned int size, GLMemoryManager::BufferType type )
{
	return boost::shared_ptr<GLBuffer>( new GLBasicBuffer(data, size) );
}


GLBasicMemoryManager::~GLBasicMemoryManager()
{
}

GLDisplayList::GLDisplayList()
{
	id_ = glGenLists(1);
}

GLDisplayList::~GLDisplayList()
{
	glDeleteLists( id_, 1 );
}

void GLDisplayList::render()
{
	glCallList(id_);
}

GLDisplayList::CreateList::CreateList(const GLDisplayList& list)
	: list_(list)
{
	glNewList(list.id_, GL_COMPILE);
}

GLDisplayList::CreateList::~CreateList()
{
	glEndList();
}

ViewFrustumCuller::ViewFrustumCuller()
{
	vl::Mat4f modelView;
	glGetFloatv( GL_MODELVIEW_MATRIX, modelView.Ref() );

	vl::Mat4f projectionMatrix;
	glGetFloatv( GL_PROJECTION_MATRIX, projectionMatrix.Ref() );

	// The "trans" is to get it in row order.
	// Someday I'd like to have a library where I can just declare this :)
	transform_ = trans(projectionMatrix)*trans(modelView);
}

bool ViewFrustumCuller::isCulled( const BoundingBox3f& box ) const
{
	boost::array<bool, 3> culledLower = {{true, true, true}};
	boost::array<bool, 3> culledUpper = {{true, true, true}};

	const vl::Vec3f min = box.minimum();
	const vl::Vec3f max = box.maximum();
	for( int i = 0; i < 2; ++i )
		for( int j = 0; j < 2; ++j )
			for( int k = 0; k < 2; ++k )
			{
				vl::Vec3f point( 
					i == 0 ? min[0] : max[0],
					j == 0 ? min[1] : max[1],
					k == 0 ? min[2] : max[2] );
				point = xform(transform_, point);

				// for this box to be culled, all of its vertices
				// must lie on one side or the other of a particular
				// clipping plane
				for( unsigned int coord = 0; coord < 3; ++coord )
				{
					culledLower[coord] = culledLower[coord] && (point[coord] < -1);
					culledUpper[coord] = culledUpper[coord] && (point[coord] > 1);
				}
			}

	for( unsigned int i = 0; i < 3; ++i )
	{
		if( culledLower[i] )
			return true;
		if( culledUpper[i] )
			return true;
	}
	return false;
}

BoundingBox2d ViewFrustumCuller::projection( const BoundingBox3f& box ) const
{
	BoundingBox2d result;

	const vl::Vec3f min = box.minimum();
	const vl::Vec3f max = box.maximum();
	for( int i = 0; i < 2; ++i )
		for( int j = 0; j < 2; ++j )
			for( int k = 0; k < 2; ++k )
			{
				vl::Vec3f point( 
					i == 0 ? min[0] : max[0],
					j == 0 ? min[1] : max[1],
					k == 0 ? min[2] : max[2] );
				point = xform(transform_, point);
				result.expand( vl::Vec2d(point[0], point[1]) );
			}

	return result;
}

GLDisableDepthMaskHandler::GLDisableDepthMaskHandler()
{
	glDepthMask(GL_FALSE);
}

GLDisableDepthMaskHandler::~GLDisableDepthMaskHandler()
{
	glDepthMask(GL_TRUE);
}


std::string prettyPrint(unsigned int value)
{
	std::ostringstream oss;
	if( value > 1000000000 )
		oss << std::setiosflags(std::ios::fixed) << std::setprecision(1) << 
			(static_cast<double>(value)/1000000000.) << "bn";
	else if( value > 1000000 )
		oss << std::setiosflags(std::ios::fixed) << std::setprecision(1) << 
			(static_cast<double>(value)/1000000.) << "m";
	else if( value > 1000 )
		oss << std::setiosflags(std::ios::fixed) << std::setprecision(1) << 
			(static_cast<double>(value)/1000.) << "k";
	else 
		oss << value;

	return oss.str();
}

TileRendererWrapper::TileRendererWrapper()
	: context_( trNew() )
{
}

TileRendererWrapper::~TileRendererWrapper()
{
	trDelete( context_ );
}

TRcontext* TileRendererWrapper::get()
{
	return context_;
}

const TRcontext* TileRendererWrapper::get() const
{
	return context_;
}


TextureStore::TextureStore()
{
}

TextureStore::~TextureStore()
{
}

// may not be thread-safe
TextureStore& TextureStore::instance()
{
	static TextureStore result;
	return result;
}

TextureStore::TexturePtr TextureStore::getTexture( const std::string& filename )
{
	TMap::iterator iter = textureMap_.find( filename );
	if( iter == textureMap_.end() )
	{
		Image image( filename );

		TexturePtr texture( new GLTexture() );
		GLBindTextureHandler bindTexture( *texture );

		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

		gluBuild2DMipmaps( GL_TEXTURE_2D, 3, image.width(), image.height(), GL_RGB,
			GL_UNSIGNED_BYTE, image.data() );

		iter = textureMap_.insert( std::make_pair(filename, texture) ).first;
	}

	return iter->second;
}

GLUQuadricWrapper::GLUQuadricWrapper()
{
	quadric_ = gluNewQuadric();
	gluQuadricOrientation( quadric_, GLU_OUTSIDE );
}

GLUQuadricWrapper::~GLUQuadricWrapper()
{
	gluDeleteQuadric( quadric_ );
}

GLUQuadricWrapper::operator GLUquadric*()
{
	return quadric_;
}


#endif // CONDOR

