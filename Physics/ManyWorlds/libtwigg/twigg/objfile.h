#ifndef __OBJFILE_H__

#define __OBJFILE_H__

#include "twigg/util.h"
#include "twigg/vlutil.h"

#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include <string>
#include <vector>
#include <deque>
#include <iosfwd>

#ifndef CONDOR
class GLTexture;
#endif

class ObjMaterial
	: public Named, boost::noncopyable
{
public:
	ObjMaterial( std::string name );
	~ObjMaterial();

	vl::Vec3 Ka() const					{ return Ka_; }
	vl::Vec3 Kd() const					{ return Kd_; }
	vl::Vec3 Ks() const					{ return Ks_; }
	double Tr() const					{ return Tr_; }
	double Ns() const					{ return Ns_; }
	double Ni() const					{ return Ni_; }
	unsigned int illum() const			{ return illuminationModel_; }

	vl::Vec3 textureAvg() const
	{
#ifndef CONDOR
		glTextureName();
#endif
		return textureAvg_;
	}

	void Ka( vl::Vec3 k )				{ Ka_ = k; }
	void Kd( vl::Vec3 k )				{ Kd_ = k; }
	void Ks( vl::Vec3 k )				{ Ks_ = k; }
	void Tr( double k )					{ Tr_ = k; }
	void Ns( double k )					{ Ns_ = k; }
	void Ni( double k )					{ Ni_ = k; }
	void illum( unsigned int k )		{ illuminationModel_ = k; }

	void setTexture( const std::string& filename );
	std::string getTexture() const { return textureFilename_; }
	bool hasTexture() const;
#ifndef CONDOR
	boost::shared_ptr<GLTexture> glTextureName() const;
#endif

	void textureOffset( const vl::Vec3d& offset );
	void textureRepeat( const vl::Vec3d& repeat );
	vl::Vec3d textureOffset() const { return textureOffset_; }
	vl::Vec3d textureRepeat() const { return textureRepeat_; }

private:
	void loadTexture();

	vl::Vec3 Ka_;
	vl::Vec3 Kd_;
	vl::Vec3 Ks_;
	double Tr_;
	double Ns_;
	double Ni_;
	unsigned int illuminationModel_;

	std::string textureFilename_;
	vl::Vec3 textureOffset_;
	vl::Vec3 textureRepeat_;

#ifndef CONDOR
	mutable boost::shared_ptr<GLTexture> glTexture_;
#endif
	mutable vl::Vec3 textureAvg_;
};


class ObjMaterialLibrary
	: public Named
{
public:
	ObjMaterialLibrary( const std::string& filename, bool readFile = true );
	void add( boost::shared_ptr<const ObjMaterial> material );
	boost::shared_ptr<const ObjMaterial> material( const std::string& name ) const;
	void write( std::ostream& out ) const;
	void write( const std::string& filename ) const;

private:
	typedef boost::shared_ptr<const ObjMaterial> ObjMaterialPtr;
	typedef std::deque<ObjMaterialPtr> ObjMaterialPtrList;
	ObjMaterialPtrList materials_;
};

// Use this as a wrapper for an .obj file
class ObjFile
	: public Named
{
public:

	template <typename T>
	class Vertex
	{
	public:
		Vertex( const T& pos )
			:	position_(pos),
				texture_( std::make_pair( false, T() ) ),
				normal_( std::make_pair( false, T() ) ) {}

		Vertex( const T& pos, const T& texture )
			:	position_(pos),
				texture_( std::make_pair( true, texture ) ),
				normal_( std::make_pair( false, T() ) ) {}
				
		Vertex( const T& pos, const T& texture, const T& normal )
			:	position_(pos),
				texture_( std::make_pair( true, texture ) ),
				normal_( std::make_pair( true, normal ) ) {}

		Vertex( const T pos, std::pair<bool, T> texture, std::pair<bool, T> normal )
			:	position_(pos),
				texture_(texture),
				normal_(normal) {}

		T position() const		{ return position_; }
		T normal() const		{ assert( normal_.first ); return normal_.second; }
		T texturePos() const	{ assert( texture_.first ); return texture_.second; }

		bool hasNormal() const	{ return normal_.first; }
		bool hasTextureCoordinate() const
								{ return texture_.first; }

	private:
		T position_; 
		std::pair< bool, T > texture_;
		std::pair< bool, T > normal_;	
	};

	template <typename T>
	class Face
	{
	public:
		Face() { vertices_.reserve( 3 ); }
		Face( const Vertex<T>& v1, const Vertex<T>& v2, const Vertex<T>& v3 )
		{
			vertices_.reserve(3);
			vertices_.push_back(v1);
			vertices_.push_back(v2);
			vertices_.push_back(v3);
		}

		void addVertex( const Vertex<T>& v )
		{
			vertices_.push_back( v );
		}

		size_t vertexCount() const { return vertices_.size(); }
		Vertex<T> vertex(unsigned int vert) const { assert( vert < vertexCount() ); return vertices_[vert]; }

	private:
		std::vector< Vertex<T> > vertices_;
	};

	class Group
		: public Named
	{
	public:
		explicit Group( std::string name )
			: Named(name)
		{
		}

		void addFace( const Face<unsigned int>& face )
		{
			faces_.push_back( face );
		}

		void setMaterial( boost::shared_ptr<const ObjMaterial> material )
		{
			material_ = material;
		}

		size_t faceCount() const { return faces_.size(); }
		Face<unsigned int> face(unsigned int face) const { return faces_[face]; }
		boost::shared_ptr<const ObjMaterial> material() const { return material_; }

	private:
		std::vector< Face<unsigned int> > faces_;
		boost::shared_ptr<const ObjMaterial> material_;
	};

	typedef std::deque<ObjMaterialLibrary> MaterialLibraryList;

	ObjFile();
	ObjFile( const std::string& filename );
	ObjFile( const std::string& filename, const std::string& name );

	ObjFile( const Group& group, 
		std::vector<vl::Vec3> vertexPositions, 
		std::vector<vl::Vec3> vertexNormals, 
		std::vector<vl::Vec3> textureCoords );

	template <typename GroupListType>
	ObjFile( 
		const std::string& name,
		const GroupListType& groups, 
		std::vector<vl::Vec3> vertexPositions,
		std::vector<vl::Vec3> vertexNormals,
		std::vector<vl::Vec3> textureCoords,
		const MaterialLibraryList& materialLibraries )
		:	Named( name ),
			groups_( groups.begin(), groups.end() ),
			vertexPositions_(vertexPositions),
			textureCoords_( textureCoords ),
			normals_( vertexNormals ),
			materialLibraries_( materialLibraries )
	{
	}

	void addGroups( const std::deque<ObjFile::Group>& groups, 
		const std::vector<vl::Vec3>& vertexPositions, 
		const std::vector<vl::Vec3>& textureCoordinates,
		const std::vector<vl::Vec3>& normals );
	void append( const ObjFile& other );

	Group group( std::string name ) const;
	size_t numGroups() const;
	Group group( size_t id ) const;

	// obj files are 1-indexed
	vl::Vec3 vertexPosition( unsigned int iPos ) const
	{
		assert( (iPos > 0) && (iPos <= vertexPositions_.size()) );
		return vertexPositions_[iPos - 1];
	}

	size_t vertexCount() const { return vertexPositions_.size(); }
	size_t normalCount() const { return normals_.size(); }
	size_t textureCoordinateCount() const { return textureCoords_.size(); }

	vl::Vec3 textureCoordinate( unsigned int iTex ) const
	{
		assert( (iTex > 0) && (iTex <= textureCoords_.size()) );
		return textureCoords_[iTex - 1];
	}

	vl::Vec3 normal( unsigned int iNorm ) const
	{
		assert( (iNorm > 0) && (iNorm <= normals_.size()) );
		return normals_[iNorm - 1];
	}

#ifndef CONDOR
	void render( const std::string& groupName ) const;
#endif // CONDOR
	std::string filename() const { return filename_; }

	void setFilename( const std::string& filename )
	{
		filename_ = filename;
	}

	void write( std::ostream& out ) const;
	void write( const std::string& filename ) const;

	std::vector<vl::Vec3d> vertexPositions() const { return vertexPositions_; }
	std::vector<vl::Vec3d> textureCoordinates() const { return textureCoords_; }
	std::vector<vl::Vec3d> normals() const { return normals_; }

	MaterialLibraryList materialLibraries() const { return materialLibraries_; }
	void appendMaterialLibraries( MaterialLibraryList libraries );

	size_t faceCount() const;

private:
	void read( const std::string& filename );
	
	typedef std::deque<Group> GroupList;

	std::string uniqueName( std::string name ) const;

	GroupList groups_;
	std::vector< vl::Vec3 > vertexPositions_;
	std::vector< vl::Vec3 > textureCoords_;
	std::vector< vl::Vec3 > normals_;
	std::string filename_;

	MaterialLibraryList materialLibraries_;
};

std::ostream& operator<<( std::ostream &output, const ObjFile::Face<unsigned int>& face);
std::ostream& operator<<( std::ostream &output, const ObjFile::Vertex<unsigned int>& vertex);

#endif
