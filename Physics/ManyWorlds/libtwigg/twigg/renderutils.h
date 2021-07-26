#ifndef __RENDERUTILS_H__
#define __RENDERUTILS_H__

#ifndef CONDOR
#include "twigg/util.h"
#include "twigg/boundingbox.h"
#include "twigg/tr.h"
#include "twigg/screenSpace.h"

#include <boost/shared_ptr.hpp>
#include <map>
#include <string>

void checkGLError(const std::string& location);

// GL utility function
enum HAlign
{
	ALIGN_HORIZ_LEFT,
	ALIGN_HORIZ_RIGHT,
	ALIGN_HORIZ_CENTER
};

enum VAlign
{
	ALIGN_VERT_TOP,
	ALIGN_VERT_BOTTOM,
	ALIGN_VERT_CENTER,
	ALIGN_VERT_BASELINE,
};


class OpenGLException
	: public Exception
{
public:
	OpenGLException( const std::string& message )
		: Exception( message ) {}
};

BoundingBox3d printString(const char *s, vl::Vec3 position, HAlign horiz = ALIGN_HORIZ_LEFT, VAlign vert = ALIGN_VERT_BASELINE, bool output = true);
void printStrings( const vl::Vec3& upperLeft, const std::deque< std::pair<std::string, std::string> >& items );
std::string prettyPrint(unsigned int value);

void drawArrow( const vl::Vec3& position, const vl::Vec3& direction, double scale = 1.0, const vl::Vec4f& color = vl::Vec4f(vl::vl_1) );
void drawShadedArrow( const vl::Vec3& position, const vl::Vec3& direction, double scale = 1.0, const vl::Vec4f& color = vl::Vec4f(vl::vl_1), bool shadedColumn = true );
void drawShadedArrow( const vl::Vec3f& position, const vl::Vec3f& direction, float scale = 1.0, const vl::Vec4f& color = vl::Vec4f(vl::vl_1), bool shadedColumn = true );
vl::Vec3f jetColor(float x);
void cube();
void legend(double min, double max);

class OpenGLScreenSpaceConverter
	: public ScreenSpaceConverter
{
public:
	OpenGLScreenSpaceConverter( const double mvmatrix[16], const double projmatrix[16], const int viewport[4] );
	OpenGLScreenSpaceConverter();

	vl::Vec3 toScreenSpace( const vl::Vec3& worldSpace ) const;
	vl::Vec3 toWorldSpace( const vl::Vec3& screenSpace ) const;

private:
	boost::array<double, 16> modelViewMatrix_;
	boost::array<double, 16> projectionMatrix_;
	boost::array<int, 4> viewPort_;
};

class GLMatrixStackHandler
{
public:
	GLMatrixStackHandler();
	~GLMatrixStackHandler();
};

class GLActionHandler
{
public:
	GLActionHandler( GLenum mode );
	~GLActionHandler();
};

class GLEnableHandler
{
public:
	GLEnableHandler( GLenum type );
	~GLEnableHandler();

private:
	GLenum capability_;
};

class GLDisableHandler
{
public:
	GLDisableHandler( GLenum type );
	~GLDisableHandler();

private:
	GLenum capability_;
};

class GLEnableClientStateHandler
{
public:
	GLEnableClientStateHandler(GLenum array);
	~GLEnableClientStateHandler();

private:
	GLenum array_;
};

class GLTexture
{
public:
	GLTexture();
	~GLTexture();

	friend class GLBindTextureHandler;

private:
	GLuint textureName() const;

private:
	GLuint glTextureName_;
};

// this class indexes textures by filename so that each is only loaded once
class TextureStore
{
private:
	TextureStore();
	~TextureStore();

public:
	static TextureStore& instance();

	typedef boost::shared_ptr<GLTexture> TexturePtr;
	TexturePtr getTexture( const std::string& filename );

private:
	typedef std::map< std::string, TexturePtr > TMap;
	TMap textureMap_;
};

class GLBindTextureHandler
{
public:
	GLBindTextureHandler(const GLTexture& texture);
	~GLBindTextureHandler();
};

class GLDisableDepthMaskHandler
{
public:
	GLDisableDepthMaskHandler();
	~GLDisableDepthMaskHandler();
};

class GLBillboardHandler
{
public:
	GLBillboardHandler( const vl::Vec3& position );

private:
	GLMatrixStackHandler handler_;
};

class GLOverlayHandler
{
public:
	GLOverlayHandler();
	~GLOverlayHandler();
};


class ARBVertexProgram
{
	GLuint progID_;
	std::string program_;
	std::string strippedProgram_;

public:
	ARBVertexProgram( const std::string& prog );
	~ARBVertexProgram();

	std::string program() const { return program_; }
	std::string strippedProgram() const { return strippedProgram_; }

	friend class BindProgram;
	class BindProgram
	{
	public:
		BindProgram( const ARBVertexProgram& program );
		~BindProgram();

	private:
		const ARBVertexProgram& program_;
	};
};

class GLBindBuffer
{
public:
	virtual ~GLBindBuffer();
};

class GLBuffer
{
public:
	virtual ~GLBuffer();
	virtual std::auto_ptr<GLBindBuffer> bind() const = 0;
	virtual const char* buffer() const = 0;
};

class GLMemoryManager
{
public:
	enum BufferType
	{
		INDEXED,
		ARRAY
	};

	template <typename T>
		boost::shared_ptr<GLBuffer> allocateBuffer( const std::vector<T>& data, GLMemoryManager::BufferType type )
	{
		return buffer( reinterpret_cast<const char*>(&data[0]), sizeof(T)*data.size(), type );
	}

	virtual ~GLMemoryManager();

protected:
	virtual boost::shared_ptr<GLBuffer> buffer( 
		const char* data, 
		unsigned int size, 
		GLMemoryManager::BufferType type ) = 0;
};

class ARBVertexBufferObject
	: public GLBuffer
{
public:
	ARBVertexBufferObject(const char* data, unsigned int size, GLMemoryManager::BufferType type);
	virtual ~ARBVertexBufferObject();

	friend class BindBuffer;
	class BindBuffer
		: public GLBindBuffer
	{
	public:
		BindBuffer( const ARBVertexBufferObject& bufferObject, GLenum target );
		~BindBuffer();

	private:
		GLenum target_;
	};

	virtual std::auto_ptr<GLBindBuffer> bind() const;
	const char* buffer() const;

private:
	GLenum bufferType_;
	GLuint bufferId_;
};

class VBOMemoryManager
	: public GLMemoryManager
{
public:
	virtual ~VBOMemoryManager();

protected:
	boost::shared_ptr<GLBuffer> buffer( 
		const char* data, 
		unsigned int size, 
		GLMemoryManager::BufferType type );
};

class GLBasicBuffer
	: public GLBuffer
{
public:
	GLBasicBuffer(const char* data, unsigned int size);
	virtual ~GLBasicBuffer();

	virtual std::auto_ptr<GLBindBuffer> bind() const;
	const char* buffer() const;

private:
	std::vector<char> data_;
};

class GLBasicMemoryManager
	: public GLMemoryManager
{
public:
	virtual ~GLBasicMemoryManager();

protected:
	boost::shared_ptr<GLBuffer> buffer( 
		const char* data, 
		unsigned int size, 
		GLMemoryManager::BufferType type );
};

class GLDisplayList
{
public:
	GLDisplayList();
	~GLDisplayList();
	void render();

	friend class CreateList;
	class CreateList
	{
	public:
		CreateList( const GLDisplayList& list );
		~CreateList();

	private:
		const GLDisplayList& list_;
	};

private:
	GLuint id_;
};

class ViewFrustumCuller
{
public:
	ViewFrustumCuller();
	bool isCulled( const BoundingBox3f& box ) const;
	BoundingBox2d projection( const BoundingBox3f& box ) const;

private:
	vl::Mat4f transform_;
};

class TileRendererWrapper
{
public:
	TileRendererWrapper();
	~TileRendererWrapper();

	TRcontext* get();
	const TRcontext* get() const;

private:
	TRcontext* context_;
};

class GLUquadric;

class GLUQuadricWrapper
{
public:
	GLUQuadricWrapper();
	~GLUQuadricWrapper();

	operator GLUquadric*();

private:
	GLUquadric* quadric_;
};

#endif // CONDOR

#endif

