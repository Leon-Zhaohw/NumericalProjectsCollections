#include "stdafx.h"

#include "twigg/imagelib.h"
#include "twigg/renderUtils.h"

#include "twigg/objfile.h"
#include "twigg/util.h"
#include "twigg/pathUtils.h"

#ifdef __APPLE__
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <boost/spirit/core.hpp>
#include <boost/spirit/attribute.hpp>
#include <boost/spirit/symbols/symbols.hpp>
#include <boost/spirit/utility/chset.hpp>
#include <boost/spirit/utility/escape_char.hpp>
#include <boost/spirit/utility/confix.hpp>
#include <boost/spirit/iterator.hpp>
#include <boost/spirit/error_handling/exceptions.hpp>
#include <boost/thread.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/scoped_array.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/exception.hpp>

#include <typeinfo>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

using std::ifstream;
using std::string;
using namespace vl;
using namespace boost::spirit;

template <typename T>
class SetBool
{
public:
	SetBool( bool& valid, T& value )
		: value_(value), valid_(valid) {}

	void operator()( T value ) const
	{
		value_ = value;
		valid_ = true;
	}

private:
	T& value_;
	bool& valid_;
};

ObjFile::ObjFile()
	: Named( "default" )
{
}

ObjFile::ObjFile( const std::string& filename )
	: Named( filename ), filename_(filename)
{
	read( filename_ );
}

ObjFile::ObjFile( const std::string& filename, const std::string& name )
	: Named(name), filename_(filename)
{
	if( filename.empty() )
		throw CreationException( "Empty filename." );

	read( filename_ );
}

ObjFile::ObjFile( const Group& group,
		std::vector<vl::Vec3> vertexPositions,
		std::vector<vl::Vec3> vertexNormals,
		std::vector<vl::Vec3> textureCoords )
		: Named( group.name() ),
			vertexPositions_(vertexPositions),
			textureCoords_( textureCoords ),
			normals_( vertexNormals )
{
	groups_.push_back(group);
}

void ObjFile::read( const std::string& filename )
{
	std::ostream& detailsOut = std::cerr;

	try
	{
		/*
		boost::filesystem::path fullPath(filename, boost::filesystem::native);
		boost::filesystem::path objPath(fullPath.branch_path());

		boost::filesystem doesn't work worth a crap on Windows
		*/

		const int maxline = 10000;
		std::ifstream ifs(filename.c_str());
		char line[maxline];

		detailsOut << "Parsing .obj file '" << filename << "'." << std::endl;
		if( !ifs )
		{
			std::string message = "Couldn't open .obj file '";
			message.append( filename );
			message.append( "'" );
			throw IOException( message );
		}

		int lineNum = 0;
		size_t faceCount = 0;

		boost::shared_ptr<const ObjMaterial> currentMaterial;

		while( ifs )
		{
			++lineNum;

			ifs.getline( line, maxline );

			std::string command;
			parse_info<char const*> info = parse( line,
				lexeme_d[ +alpha_p[ ::append(command) ] ] );
			if( !info.hit )
				continue;

			const char* curPos = info.stop;

			if( command == "v" ) // vertex
			{
				vl::Vec3 pos;
				parse_info<char const*> info = parse( curPos, 
					real_p[assign(pos[0])] >> real_p[assign(pos[1])] >> real_p[assign(pos[2])],
					space_p );
				if( !info.full )
					throw ParserException("Invalid vertex", filename, lineNum, 1);

				vertexPositions_.push_back( pos );
			}
			else if( command == "vn" )
			{
				vl::Vec3 normal;
				parse_info<char const*> info = parse( curPos, 
					real_p[assign(normal[0])] >> real_p[assign(normal[1])] >> real_p[assign(normal[2])],
					space_p );
				if( !info.full )
					throw ParserException("Invalid normal", filename, lineNum, 1);

				normals_.push_back( normal );
			}
			else if( command == "vt" )
			{
				vl::Vec3 tex = vl_0;	// w defaults to 0
				parse_info<char const*> info = parse( curPos, 
					real_p[assign(tex[0])] >> real_p[assign(tex[1])] >> !(real_p[assign(tex[2])]),
					space_p );
				if( !info.full )
					throw ParserException("Invalid texture coordinate", filename, lineNum, curPos - line);

				textureCoords_.push_back( tex );
			}
			else if( command == "g" )
			{
				std::string name;
				parse_info<char const*> info = parse( curPos, 
					*graph_p[ ::append(name) ], space_p );
				GroupList::iterator group = 
					groups_.insert( groups_.end(), Group(name) );
				group->setMaterial( currentMaterial );
			}
			else if( command == "f" )
			{
				++faceCount;

				if( groups_.empty() )
					groups_.push_back( Group("default") );

				Face<unsigned int> face;

				const char* endLine = line + strlen(line);
				while( curPos < endLine )
				{
					unsigned int pos;
					std::pair< bool, unsigned int > texPos;
					std::pair< bool, unsigned int > normal;
					parse_info<char const*> info = parse( curPos,
						lexeme_d[ 
							uint_p[ assign(pos) ] >> 
								!( ch_p('/') >> !uint_p[SetBool<unsigned int>(texPos.first, texPos.second)] >>
								!( ch_p('/') >> !uint_p[SetBool<unsigned int>(normal.first, normal.second)] ) ) ],
								space_p );

					curPos = info.stop;
					face.addVertex( Vertex<unsigned int>( pos, texPos, normal ) );
				}

				groups_.rbegin()->addFace( face );
			}
			else if( command == "mtllib" )
			{
				std::string name;
				parse_info<char const*> info = parse( curPos, 
					*graph_p[ ::append(name) ], space_p );
				materialLibraries_.push_back( ObjMaterialLibrary( 
					fileInPath( name, pathForFile(filename) ) ) );
				std::sort( materialLibraries_.begin(), materialLibraries_.end(),
					NamedComparator<ObjMaterialLibrary>() );
			}
			else if( command == "usemtl" )
			{
				std::string name;
				parse_info<char const*> info = parse( curPos, 
					*graph_p[ ::append(name) ], space_p );
				for( MaterialLibraryList::const_iterator libItr = materialLibraries_.begin();
					libItr != materialLibraries_.end();
					++libItr )
				{
					boost::shared_ptr<const ObjMaterial> mat = libItr->material(name);
					if( !mat )
						continue;

					currentMaterial = mat;

					if( !groups_.empty() && groups_.rbegin()->faceCount() == 0 )
						groups_.rbegin()->setMaterial(currentMaterial);
				}

				if( !currentMaterial )
				{
					throw CreationException( "Material '" + name + "' not found in any of the material libraries for '"
						+ filename + "'." );
				}
			}
			else if( command == "s" )
			{
				// ignore silently for now
				std::cerr << "Warning: ignoring smoothing groups." << std::endl;
			}
			else if( command == "o" )
			{
				std::cerr << "Warning: ignoring object names." << std::endl;
			}
			else
			{
				std::ostringstream msg;
				msg << "Invalid line in .obj file '" << filename << "': " << line;
				throw ParserException(msg.str(), filename, lineNum, 1);
			}
		}

		std::sort( groups_.begin(), groups_.end(), NamedComparator<Group>() );

		// statistics
		detailsOut << "Parsed obj file '" << filename << "'; statistics:" << std::endl;
		detailsOut << "   " << groups_.size() << " groups," << std::endl;
		detailsOut << "   " << faceCount << " faces," << std::endl;
		detailsOut << "   " << vertexPositions_.size() << " vertices," << std::endl;
		detailsOut << "   " << normals_.size() << " normals, " << std::endl;
		detailsOut << "   " << textureCoords_.size() << " texture coordinates, " << std::endl;
	}
	catch( boost::filesystem::filesystem_error& e )
	{
		throw IOException( e.what() );
	}
}

size_t ObjFile::faceCount() const
{
	size_t result = 0;
	for( GroupList::const_iterator groupItr = groups_.begin();
		groupItr != groups_.end(); ++groupItr )
	{
		result += groupItr->faceCount();
	}

	return result;
}

ObjFile::Group ObjFile::group( std::string name ) const
{
	typedef GroupList::const_iterator GroupListItr;
	typedef std::pair<GroupListItr, GroupListItr> GroupListItrPair;

	GroupListItrPair p = std::equal_range( groups_.begin(), groups_.end(), 
		Group(name), NamedComparator<Group>() );
	if( p.first == p.second )
	{
		std::ostringstream oss;
		oss << "Invalid group name: '" << name << "'.";
		throw CreationException( oss.str() );
	}

	return *p.first;
}

std::string incrementString( std::string str )
{
	std::string oldString = str;
	std::string::reverse_iterator rItr = str.rbegin();
  while( rItr != str.rend() )
  {
    if( isdigit(*rItr) )
    {
      if( *rItr < '9' )
      {
        *rItr = *rItr + 1;
        break;
      }
      else
      {
        *rItr = '0';
      }
    }

		++rItr;
  }

	if( rItr == str.rend() )
		str += "_1";

  return str;
}

std::string ObjFile::uniqueName( std::string name ) const
{
	typedef GroupList::const_iterator GroupListItr;
	typedef std::pair<GroupListItr, GroupListItr> GroupListItrPair;

	GroupListItrPair p = std::equal_range( groups_.begin(), groups_.end(),
		Group(name), NamedComparator<Group>() );
	if( p.first != p.second )
		name += "_1";

	while( p.first != p.second )
	{
		name = incrementString( name );
std::cout << "name: " << name << std::endl;
		p = std::equal_range( groups_.begin(), groups_.end(),
	    Group(name), NamedComparator<Group>() );
	}
	return name;
}

size_t ObjFile::numGroups() const
{
	return groups_.size();
}

ObjFile::Group ObjFile::group( size_t id ) const
{
	return groups_[id];
}

void ObjFile::append( const ObjFile& other )
{
	this->addGroups( other.groups_, 
		other.vertexPositions_,
		other.textureCoords_,
		other.normals_ );

	this->appendMaterialLibraries( other.materialLibraries_ );
}

class NamedEquality
{
public:
	bool operator()( const Named& left, const Named& right )
	{
		return left.name() == right.name();
	}
};

void ObjFile::appendMaterialLibraries( MaterialLibraryList libraries )
{
	std::sort( libraries.begin(), libraries.end(), NamedComparator<ObjMaterialLibrary>() );

	MaterialLibraryList result;
	std::set_union( materialLibraries_.begin(), materialLibraries_.end(), 
		libraries.begin(), libraries.end(), std::back_inserter(result), 
		NamedComparator<ObjMaterialLibrary>() );
	std::swap( result, this->materialLibraries_ );
}

void ObjFile::addGroups( const std::deque<ObjFile::Group>& groups, 
	const std::vector<vl::Vec3>& vertexPositions, 
	const std::vector<vl::Vec3>& textureCoordinates,
	const std::vector<vl::Vec3>& normals )
{
	const size_t positionsOffset = this->vertexPositions_.size();
	const size_t normalsOffset = this->normals_.size();
	const size_t textCoordOffset = this->textureCoords_.size();

	for( std::deque<ObjFile::Group>::const_iterator groupItr = groups.begin();
		groupItr != groups.end(); ++groupItr )
	{
		ObjFile::Group newGroup( uniqueName(groupItr->name()) );

		if( groupItr->material() )
		{
/*
			if( materialLibraries_.empty() )
				materialLibraries_.push_back( ObjMaterialLibrary(this->name() + std::string(".mtl")) );
			materialLibraries_.back().add( groupItr->material() );
*/
			newGroup.setMaterial( groupItr->material() );
		}

		for( unsigned int iFace = 0; iFace < groupItr->faceCount(); ++iFace )
		{
			ObjFile::Face<unsigned int> oldFace = groupItr->face(iFace);
			ObjFile::Face<unsigned int> newFace;
			for( unsigned int iVertex = 0; iVertex < oldFace.vertexCount(); ++iVertex )
			{
				ObjFile::Vertex<unsigned int> oldVertex = oldFace.vertex(iVertex);
	
				std::pair<bool, unsigned int> texCoord(false, 0);
				texCoord.first = oldVertex.hasTextureCoordinate();
				if( texCoord.first )
					texCoord.second = oldVertex.texturePos() + textCoordOffset;
	
				std::pair<bool, unsigned int> normal(false, 0);
				normal.first = oldVertex.hasNormal();
				if( normal.first )
					normal.second = oldVertex.normal() + normalsOffset;
	
				newFace.addVertex(
					ObjFile::Vertex<unsigned int>( oldVertex.position() + positionsOffset, 
						texCoord, normal ) );
			}
	
			newGroup.addFace( newFace );
		}

		GroupList::iterator newGroupItr = std::upper_bound( 
			groups_.begin(), groups_.end(), newGroup, NamedComparator<Group>() );
		groups_.insert( newGroupItr, newGroup );
	}

	std::copy( vertexPositions.begin(), vertexPositions.end(), 
		std::back_inserter(vertexPositions_) );
	std::copy( textureCoordinates.begin(), textureCoordinates.end(), 
		std::back_inserter(textureCoords_) );
	std::copy( normals.begin(), normals.end(), 
		std::back_inserter(normals_) );
}

void writeObjVec( std::ostream& os, const vl::Vec3& vec )
{
	const unsigned int precision = 10;
	os	<< std::setprecision(precision) << vec[0] << " " 
		<< std::setprecision(precision) << vec[1] << " " 
		<< std::setprecision(precision) << vec[2] << "\n";
}

void dumpObjArray( std::ostream& out, const std::vector<vl::Vec3>& array, const std::string& prefix )
{
	for( std::vector<vl::Vec3>::const_iterator itr = array.begin();
		itr != array.end();
		++itr )
	{
		out << prefix << " ";
		writeObjVec( out, *itr );
	}
}

class GroupVertexIndexOrdering
{
public:
	bool operator()( const ObjFile::Group& left, const ObjFile::Group& right ) const
	{
		assert( left.faceCount() != 0 && right.faceCount() != 0 );
		return std::less<unsigned int>()( 
			left.face(0).vertex(0).position(), 
			right.face(0).vertex(0).position() );
	}
};

void ObjFile::write( const std::string& filename ) const
{
	std::ofstream out( filename.c_str() );
	if( !out )
		throw IOException( "Unable to open file '" + filename + "' for writing." );

	write( out );
}

void ObjFile::write( std::ostream& out ) const
{
	for( MaterialLibraryList::const_iterator mtllibIter = materialLibraries_.begin();
		mtllibIter != materialLibraries_.end();
		++mtllibIter )
	{
		out << "mtllib " << mtllibIter->name() << "\n";
	}

	dumpObjArray( out, this->vertexPositions_, "v" );
	dumpObjArray( out, this->textureCoords_, "vt" );
	dumpObjArray( out, this->normals_, "vn" );

	// Only want to export non-empty groups to make things easier on
	//   Rob/Doug
	GroupList groups;
	for( GroupList::const_iterator groupItr = groups_.begin();
		groupItr != groups_.end(); ++groupItr )
	{
		if( groupItr->faceCount() == 0 )
			continue;
		groups.push_back( *groupItr );
	}

	if( groups.empty() )
		throw IOException( "No nonempty groups in .obj file." );

	// Would like to export groups in order lowest vertex->highest vertex
	//   rather than alphabetical
	std::sort( groups.begin(), groups.end(), GroupVertexIndexOrdering() );
	for( GroupList::const_iterator groupItr = groups.begin();
		groupItr != groups.end(); ++groupItr )
	{
		out << "g " << groupItr->name() << "\n";
		if( groupItr->material() )
			out << "usemtl " << groupItr->material()->name() << "\n";

		for( unsigned int iFace = 0; iFace < groupItr->faceCount(); ++iFace )
		{
			out << "f " << groupItr->face( iFace ) << "\n";
		}
	}
}

#ifndef CONDOR
void ObjFile::render( const std::string& groupName ) const
{
	Group group = ObjFile::group( groupName );
	glEnable(GL_LIGHTING);
	boost::scoped_ptr<GLBindTextureHandler> useTexture;

	if( group.material().get() != 0 )
	{
		boost::shared_ptr<const ObjMaterial> mat = group.material();
		vl::Vec3 Ka_v( mat->Ka() );
		vl::Vec3 Kd_v( mat->Kd() );
		vl::Vec3 Ks_v( mat->Ks() );
		
		{
			const GLfloat Kd[] = { Kd_v[0], Kd_v[1], Kd_v[2], 1.0f };
			const GLfloat Ks[] = { Ks_v[0], Ks_v[1], Ks_v[2], 1.0f };
			const GLfloat Ns[] = { mat->Ns() };
			const GLfloat no_mat[] = { 0.0f, 0.0f, 0.0f, 0.0f };

			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, no_mat);
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, Kd);
			glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, Ks);
			glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, Ns);
			glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, no_mat );
		}
		
		if( mat->hasTexture() )
			useTexture.reset( new GLBindTextureHandler( *mat->glTextureName() ) );
	}
	else
	{
		GLfloat color[] = { 0.21f, 0.41f, 0.09f, 1.0f };
		const GLfloat no_mat[] = { 0.0f, 0.0f, 0.0f, 0.0f };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, no_mat);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, no_mat);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, no_mat);
		glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, no_mat );
	}

	GLActionHandler triangles( GL_TRIANGLES );
	for( unsigned int iFace = 0; iFace < group.faceCount(); ++iFace )
	{
		ObjFile::Face<unsigned int> face = group.face( iFace );
		// triangulate if necessary
		for( size_t j = 0; (j+2) < face.vertexCount(); ++j )
		{
			for( unsigned int iVertex = 0; iVertex < 3; ++iVertex )
			{
				Vertex<unsigned int> vertex = face.vertex(
					iVertex == 0 ? iVertex : j + iVertex);

				if( vertex.hasTextureCoordinate() )
					glTexCoord3dv( textureCoordinate( vertex.texturePos() ).Ref() );

				if( vertex.hasNormal() )
					glNormal3dv( normal( vertex.normal() ).Ref() );

				glVertex3dv( vertexPosition( vertex.position() ).Ref() );
			}
		}
	}
}
#endif

ObjMaterialLibrary::ObjMaterialLibrary( const std::string& filename, bool readFile )
	: Named(filename)
{
	std::ostream& detailsOut = std::cerr;

	if( !readFile )
		return;

	std::deque< boost::shared_ptr<ObjMaterial> > materials;

	const int maxline = 1000;
	std::ifstream ifs(filename.c_str());
	std::vector<char> line(maxline);

	detailsOut << "Parsing .mtl file '" << filename << "'." << std::endl;
	if( !ifs )
	{
		std::string message = "Couldn't open .mtl file '";
		message.append( filename );
		message.append( "'" );
		throw CreationException( message );
	}

	int lineNum = 0;
	while( ifs )
	{
		++lineNum;

		ifs.getline( &line[0], maxline );
		const char* curPos = &line[0];
		while( isspace(*curPos) && *curPos != 0 )
			++curPos;

		if( *curPos == 0 || *curPos == '#' )
			continue;

		std::string command;
		while( !isspace(*curPos) && *curPos != 0 )
			command.push_back( *curPos++ );

		if( command == "newmtl" ) // vertex
		{
			std::string name;
			parse_info<char const*> info = parse( curPos, 
				*graph_p[ append(name) ], space_p );
			materials.push_back( 
				boost::shared_ptr<ObjMaterial>( new ObjMaterial(name) ) );
		}
		else
		{
			if( materials.empty() )
			{
				std::ostringstream errMsg;
				errMsg << "Error: encountered material property before 'newmtl' command at line "
					<< lineNum << " of '" << filename << "'.";
				throw CreationException( errMsg.str() );
			}

			if( command == "Ka" || command == "Kd" || command == "Ks" )
			{
				vl::Vec3 k;
				parse_info<char const*> info = parse( curPos, 
					real_p[assign(k[0])] >> real_p[assign(k[1])] >> real_p[assign(k[2])],
					space_p );
				if( !info.hit )
				{
					std::ostringstream errMsg;
					errMsg << "Expected 3-vector for property '" << command << "' at line "
						<< lineNum << " of '" << filename << "'.";
					throw CreationException( errMsg.str() );
				}

				if( command == "Ka" )
					materials.back()->Ka( k );
				else if( command == "Kd" )
					materials.back()->Kd( k );
				else if( command == "Ks" )
					materials.back()->Ks( k );
				else
					assert( false );	// shouldn't ever reach here.
			}
			else if( command == "Tr" || command == "d" || command == "Ns" || command == "Ni" )
			{
				double k;
				parse_info<char const*> info = parse( curPos, 
					real_p[assign(k)],
					space_p );
				if( !info.hit )
				{
					std::ostringstream errMsg;
					errMsg << "Expected double for property '" << command << "' at line "
						<< lineNum << " of '" << filename << "'.";
					throw CreationException( errMsg.str() );
				}

				if( command == "Tr" || command == "d" )
					materials.back()->Tr(k);
				else if( command == "Ns" )
					materials.back()->Ns(k);
				else if( command == "Ni" )
					materials.back()->Ni(k);
				else
					assert( false );	// shouldn't ever reach here.
			}
			else if( command == "illum" )
			{
				unsigned int k;
				parse_info<char const*> info = parse( curPos, 
					uint_p[assign(k)],
					space_p );
				if( !info.hit )
				{
					std::ostringstream errMsg;
					errMsg << "Expected uint for property '" << command << "' at line "
						<< lineNum << " of '" << filename << "'.";
					throw CreationException( errMsg.str() );
				}

				materials.back()->illum(k);
			}
			else if( command == "map_Kd" || command == "map_Ka" )
			{
				enum CurrentState
				{
					None,
					Offset,
					Repeat,
					Unknown
				} currentState = None;

				std::deque<double> offset;
				std::deque<double> repeat;
				std::string texFilename;

				while( *curPos != 0  )
				{
					while( isspace(*curPos) && *curPos != 0 )
						++curPos;

					if( *curPos == 0 )
						break;

					if( *curPos == '-' && isalpha( *(curPos+1) ) )
					{
						++curPos;
						std::string option;
						while( !isspace(*curPos) && *curPos != 0 )
							option.push_back( *curPos++ );
						if( option == "s" )
							currentState = Repeat;
						else if( option == "o" )
							currentState = Offset;
						else
						{
							detailsOut << "Warning: unknown texture option: -" << option << std::endl;
							currentState = Unknown;
						}

						continue;
					}

					if( isdigit(*curPos) || *curPos == '.' || *curPos == '-' )
					{
						// found a number; parse it
						double k;
						parse_info<char const*> info = parse( curPos, 
							real_p[assign(k)],
							space_p );
						if( !info.hit )
						{
							std::ostringstream errMsg;
							errMsg << "Expected: double at line "
								<< lineNum << " of '" << filename << "'.";
							throw CreationException( errMsg.str() );
						}

						switch( currentState )
						{
							case Unknown:
								break;
							case Offset:
								offset.push_back( k );
								break;
							case Repeat:
								repeat.push_back( k );
								break;
							case None:
								{
									std::ostringstream errMsg;
									errMsg << "Unexpected real at line " << 
										lineNum << " of '" << filename << "'.";
									throw CreationException( errMsg.str() );
								}
						}

						curPos = info.stop;
						continue;
					}

					if( !texFilename.empty() )
					{
						std::ostringstream errMsg;
							errMsg << "Unexpected extra texture filename at line "
							<< lineNum << " of '" << filename << "'.";
						throw CreationException( errMsg.str() );
					}

					while( !isspace(*curPos) && *curPos != 0 )
						texFilename.push_back( *curPos++ );
				}

				materials.back()->setTexture( 
					fileInPath( texFilename, pathForFile(filename) ) );

				{
					vl::Vec3d offsetVec( vl::vl_0 );
					for( unsigned int i = 0; i < offset.size(); ++i )
						offsetVec[i] = offset[i];
					materials.back()->textureOffset( offsetVec );
				}

				{
					vl::Vec3d repeatVec( vl::vl_1 );
					for( unsigned int i = 0; i < repeat.size(); ++i )
						repeatVec[i] = repeat[i];
					materials.back()->textureRepeat( repeatVec );
				}
			}
			else if( command == "Tf" )
			{
				detailsOut << "Throwing out mystery material attribute 'Tf'." << std::endl;
			}
			else if( command == "map_opacity" )
			{
				detailsOut << "Throwing out opacity map." << std::endl;
			}
			else if( command == "map_bump" || command == "bump" )
			{
				detailsOut << "Throwing out bump map." << std::endl;
			}
			else if( command == "map_d" )
			{
				detailsOut << "Throwing out unknown map 'map_d'." << std::endl;
			}
			else if( command == "map_kS" )
			{
				detailsOut << "Throwing out specular map." << std::endl;
			}
			else if( command == "map_kA" )
			{
				detailsOut << "Throwing out ambient map." << std::endl;
			}
			else
			{
				std::ostringstream errMsg;
				errMsg << "Error: unknown command '" << command << "' at line " 
					<< lineNum << " of '" << filename << "'.";
				throw CreationException( errMsg.str() );
			}
		}
	}

	std::sort( materials.begin(), materials.end(), 
		NamedPtrComparator() );
	std::copy( materials.begin(), materials.end(),
		std::back_inserter(materials_) );
}

boost::shared_ptr<const ObjMaterial> ObjMaterialLibrary::material( const std::string& name ) const
{
	typedef ObjMaterialPtrList::const_iterator MaterialListItr;
	typedef std::pair<MaterialListItr, MaterialListItr> MaterialListItrPair;
	typedef boost::shared_ptr<const Named> NamedPtr;

	MaterialListItrPair p = std::equal_range( materials_.begin(), materials_.end(), 
		NamedPtr(new Named(name)), NamedPtrComparator() );
	if( p.first == p.second )
		return boost::shared_ptr<const ObjMaterial>();

	/*
	{
		std::ostringstream oss;
		oss << "Material '" + name + "' not found in library '" + ObjMaterialLibrary::name() + "'.";
		throw CreationException( oss.str() );
	}
	*/

	return *p.first;
}


void ObjMaterialLibrary::add( boost::shared_ptr<const ObjMaterial> material )
{
	ObjMaterialPtrList::iterator materialItr = std::lower_bound(
		materials_.begin(), materials_.end(), material, NamedPtrComparator() );
	materials_.insert( materialItr, material );
}

void ObjMaterialLibrary::write( const std::string& filename ) const
{
	std::ofstream out( filename.c_str() );
	if( !out )
		throw IOException( "Unable to open file '" + filename + "' for writing." );

	write( out );
}

void ObjMaterialLibrary::write( std::ostream& out ) const
{
	std::ostream& detailsOut = std::cerr;
	for( ObjMaterialPtrList::const_iterator materialItr = materials_.begin();
		materialItr != materials_.end(); ++materialItr )
	{
		out << "newmtl " << (*materialItr)->name() << "\n";
		out << "Kd "; writeObjVec( out, (*materialItr)->Kd() );
		out << "Ks "; writeObjVec( out, (*materialItr)->Ks() );
		out << "Ka "; writeObjVec( out, (*materialItr)->Ka() );
		out << "Ns " << (*materialItr)->Ns() << "\n";
		out << "Tr " << (*materialItr)->Ns() << "\n";
		out << "Ni " << (*materialItr)->Ns() << "\n";
		out << "illum " << (*materialItr)->illum() << "\n";

		if( (*materialItr)->hasTexture() )
		{
			detailsOut << "Warning: no support currently for exporting textures." << std::endl;
		}
	}
}

ObjMaterial::ObjMaterial( std::string name )
	:	Named( name ),
		Ka_(vl::vl_0),			// default to no ambient,
		Kd_(1.0, 1.0, 1.0),		// diffuse: white,
		Ks_( vl::vl_0 ),		// no specular,
		Tr_( 1.0 ),				// no transparency (alpha = 1),
		Ns_( 100.0 ),			// High specular exponent,
		Ni_( 1.0 ),				// refractive index doesn't matter.
		illuminationModel_(1),	// diffuse-only illumination?
		textureOffset_(vl::vl_0),
		textureRepeat_(vl::vl_1),
		textureAvg_(0.5, 0.5, 0.5)
{
}

void ObjMaterial::textureOffset( const vl::Vec3d& offset )
{
	textureOffset_ = offset;
}

void ObjMaterial::textureRepeat( const vl::Vec3d& repeat )
{
	textureRepeat_ = repeat;
}

void ObjMaterial::setTexture( const std::string& filename )
{
	textureFilename_ = filename;
}

bool ObjMaterial::hasTexture() const
{
	return !textureFilename_.empty();
}

#ifndef CONDOR
boost::shared_ptr<GLTexture> ObjMaterial::glTextureName() const
{
	assert( hasTexture() );

	if( !glTexture_ )
	{
		glTexture_.reset( new GLTexture() );
		GLBindTextureHandler bindTexture( *glTexture_ );

		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);


		try
		{
			Image image( this->textureFilename_ );

			gluBuild2DMipmaps( GL_TEXTURE_2D, 3, image.width(), image.height(), GL_RGB,
				GL_UNSIGNED_BYTE, image.data() );
			TinyVec<unsigned char, 3> avgColor = image.averageColor();
			this->textureAvg_ = vl::Vec3d( avgColor[0], avgColor[1], avgColor[2] );
			this->textureAvg_ /= 255.0;
		}
		catch( CreationException& e )
		{
			std::cerr << e.message();

			// need to create a standin texture:
			std::vector<unsigned char> image( 3, 128 );
			gluBuild2DMipmaps( GL_TEXTURE_2D, 3, 1, 1, GL_RGB,
				GL_UNSIGNED_BYTE, &image[0] );
		}
	}

	return glTexture_;
}
#endif

ObjMaterial::~ObjMaterial()
{
}

std::ostream& operator<<( std::ostream &out, const ObjFile::Face<unsigned int>& face)
{
	for( unsigned int iVertex = 0; iVertex < face.vertexCount(); ++iVertex )
	{
		if( iVertex != 0 )
			out << " ";

		ObjFile::Vertex<unsigned int> vertex = face.vertex(iVertex);
		out << face.vertex(iVertex);
	}

	return out;
}

std::ostream& operator<<( std::ostream &out, const ObjFile::Vertex<unsigned int>& vertex)
{
	out << vertex.position();
	if( vertex.hasNormal() || vertex.hasTextureCoordinate() )
	{
		out << "/";
		if( vertex.hasTextureCoordinate() )
			out << vertex.texturePos();

		if( vertex.hasNormal() )
			out << "/" << vertex.normal();
	}

	return out;
}



