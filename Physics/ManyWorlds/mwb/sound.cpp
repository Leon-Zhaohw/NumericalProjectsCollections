#include "stdafx.h"

#include "sound.h"
#include "physicsObject.h"
#include "voxels.h"
#include "barbicUtils.h"
#include "mechModel.h"

#include "twigg/objFile.h"
#include "twigg/linalg.h"

#include <boost/functional/hash/hash.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <boost/spirit/core.hpp>
#include <boost/spirit/attribute.hpp>
#include <boost/spirit/symbols/symbols.hpp>
#include <boost/spirit/utility/chset.hpp>
#include <boost/spirit/utility/escape_char.hpp>
#include <boost/spirit/utility/confix.hpp>
#include <boost/spirit/iterator.hpp>
#include <boost/spirit/error_handling/exceptions.hpp>

#include <mkl_cblas.h>

#ifdef SOUND
namespace planning {

template <typename T>
struct hash_triple
	: stdext::hash_compare< boost::tuple<T, T, T> >
{
	typedef boost::tuple<T, T, T> Key;
	size_t operator()( const Key& tup ) const
	{
		std::size_t seed = 0;
		boost::hash_combine(seed, tup.get<0>());
		boost::hash_combine(seed, tup.get<1>());
		boost::hash_combine(seed, tup.get<2>());
		return seed;
	}

	bool operator()( const Key& left, const Key& right )
	{
		std::less<Key> comparator;
		return comparator(left, right);
	}
};

struct QuadFace
{
public:
	QuadFace( size_t i, size_t j, size_t k, size_t l )
	{
		vertices[0] = i;
		vertices[1] = j;
		vertices[2] = k;
		vertices[3] = l;
	}

	boost::array<size_t, 4> vertices;
};

struct HashQuadFace
	: stdext::hash_compare< HashQuadFace >
{
public:
	size_t operator()( QuadFace f ) const
	{
		std::size_t seed = 0;

		std::sort( f.vertices.begin(), f.vertices.end() );
		for( size_t i = 0; i < 4; ++i )
            boost::hash_combine(seed, f.vertices[i]);
		return seed;
	}

	bool operator()( QuadFace left, QuadFace right )
	{
		std::sort( left.vertices.begin(), left.vertices.end() );
		std::sort( right.vertices.begin(), right.vertices.end() );
		return std::lexicographical_compare( 
			left.vertices.begin(), left.vertices.end(), 
			right.vertices.begin(), right.vertices.end() );
	}
};

typedef STLEXT hash_set<QuadFace, HashQuadFace> QuadFaceSet;

size_t findVertex( size_t i, size_t j, size_t k, 
				const std::deque<SmallIndexTriple>& vertices )
{
	std::deque<SmallIndexTriple>::const_iterator itr = 
		std::lower_bound( vertices.begin(), vertices.end(), 
			SmallIndexTriple(i, j, k) );
	assert( itr != vertices.end() );

	size_t iVert = std::distance( vertices.begin(), itr );
	return iVert;
}

// adds external face; if it's already in the set, face is declared
//   to be internal and thrown out
void addExternalFace( QuadFaceSet& externalFaces, 
	size_t vert0, size_t vert1, size_t vert2, size_t vert3 )
{
	const QuadFace newFace( vert0, vert1, vert2, vert3 );
	bool inserted;
	QuadFaceSet::iterator itr;
	boost::tie(itr, inserted) = externalFaces.insert( newFace );
	if( !inserted )
		externalFaces.erase( itr );
}

boost::shared_ptr<SquashingCubesModes> computeModes( PhysicsObjectPtr object, Event& event )
{
	if( !object->hasAudio() )
		return boost::shared_ptr<SquashingCubesModes>();

	{
		// first step is to gather up all the materials
		// and determine the bounding box
		std::vector<ConstMaterialPtr> materials;
		BoundingBox3f bounds;
		vl::Mat4f objectTransform = object->transform();
		vl::Mat4f invObjectTransform = vl::inv( objectTransform );

		vl::Mat4f rigidTransform = object->rigidTransform();

		std::deque<ConstPhysicsObjectPtr> allObjects;
		if( boost::dynamic_pointer_cast<const CombinedObject>( object ) )
		{
			std::deque<ConstSceneGraphElementPtr> toHandle( 1, object );
			while( !toHandle.empty() )
			{
				ConstSceneGraphElementPtr current = toHandle.front();
				toHandle.pop_front();

				ConstPhysicsObjectPtr physicsObject =
					boost::dynamic_pointer_cast<const PhysicsObject>( current );
				if( physicsObject && physicsObject->hasAudio() )
				{
					materials.push_back( physicsObject->material() );

					vl::Mat4f localTransform = rigidInverse(rigidTransform) * physicsObject->transform();
					bounds.expand( physicsObject->bounds(localTransform) );

					allObjects.push_back(physicsObject);
				}

				std::copy( current->children_begin(), current->children_end(),
					std::front_inserter(toHandle) );
			}

			std::sort( materials.begin(), materials.end(), MaterialSoundPropertiesOrder() );
			materials.erase( 
				std::unique(materials.begin(), materials.end(), MaterialSoundPropertiesEqual()), 
				materials.end() );
		}
		else
		{
			// want the bounds in the local frame
			// so we factor out the rigid portion of the transform when computing
			// the bounding box
			bounds = object->bounds( rigidInverse(rigidTransform) * objectTransform );
			materials.resize(1, object->material());
			allObjects.resize(1, object);
		}

		// @todo this should be dependent on the bounding box aspect ratio:
		const size_t maxResolution = 5;
		const vl::Vec3d bmin = toVec3d(bounds.minimum());
		const vl::Vec3d diff = toVec3d(bounds.maximum()) - bmin;
		double longestAxis = 0.0;
		for( vl::Int i = 0; i < 3; ++i )
			longestAxis = std::max( longestAxis, diff[i] );

		TinyVec<size_t, 3> resolution;
		for( size_t i = 0; i < 3; ++i )
			resolution[i] = std::max<size_t>( boost::numeric_cast<double>( (diff[i]/longestAxis) 
				* static_cast<double>(maxResolution) + 0.5 ), 1 );
		BoundedVoxelGrid voxels( bounds, resolution, materials.size() );

		std::for_each( allObjects.begin(), allObjects.end(),
			boost::bind( &PhysicsObject::voxelize, _1, boost::ref(voxels), 
				boost::ref(materials), rigidInverse(rigidTransform) ) );

		boost::shared_ptr<SquashingCubesModes> modesStructure
			( new SquashingCubesModes(voxels) );

		// Okay, now we have the voxels
		// count how many voxels
		/*
		size_t nVoxels = 0;
		for( size_t i = 0; i < resolution[0]; ++i )
			for( size_t j = 0; j < resolution[1]; ++j )
				for( size_t k = 0; k < resolution[2]; ++k )
				{
					if( voxels(i, j, k) > 0 )
						++nVoxels;
				}

		// We will need to be able to map back from voxel number to
		//   voxel id, we'll use this simple data structure:
		VoxelGrid voxelToId( resolution[0], resolution[1], resolution[2], nVoxels );
		// start with 1 to distinguish from unused voxels. 
		size_t curVoxel = 1;
		for( size_t i = 0; i < resolution[0]; ++i )
			for( size_t j = 0; j < resolution[1]; ++j )
				for( size_t k = 0; k < resolution[2]; ++k )
				{
					if( voxels(i, j, k) > 0 )
						voxels.set(i, j, k, curVoxel++);
				}
				*/

		STLEXT hash_set<SmallIndexTriple, hash_triple<SmallIndex> > allVerticesSet;

		for( size_t i = 0; i < resolution[0]; ++i )
			for( size_t j = 0; j < resolution[1]; ++j )
				for( size_t k = 0; k < resolution[2]; ++k )
				{
					if( voxels(i, j, k) > 0 )
					{
						allVerticesSet.insert( SmallIndexTriple(i,   j,   k  ) );
						allVerticesSet.insert( SmallIndexTriple(i+1, j,   k  ) );
						allVerticesSet.insert( SmallIndexTriple(i,   j+1, k  ) );
						allVerticesSet.insert( SmallIndexTriple(i,   j,   k+1) );
						allVerticesSet.insert( SmallIndexTriple(i+1, j+1, k  ) );
						allVerticesSet.insert( SmallIndexTriple(i,   j+1, k+1) );
						allVerticesSet.insert( SmallIndexTriple(i+1, j,   k+1) );
						allVerticesSet.insert( SmallIndexTriple(i+1, j+1, k+1) );
					}
				}

		if( allVerticesSet.empty() )
			return boost::shared_ptr<SquashingCubesModes>();

		modesStructure->vertices.swap( 
			std::deque<SmallIndexTriple>( allVerticesSet.begin(), allVerticesSet.end() ) );
		std::sort( modesStructure->vertices.begin(), modesStructure->vertices.end() );

		vl::Vec3d increment;
		for( vl::Int i = 0; i < 3; ++i )
			increment[i] = diff[i] / boost::numeric_cast<double>(resolution[i]);

		std::string prefix( "test" );
		{
			std::ofstream ofs( (prefix + ".inp").c_str() );
			ofs.precision(8);

			ofs << "*HEADING\n";
			ofs << "Model: elastic cubes\n";

			// dump out all the positions
			ofs << "*NODE\n";

			std::vector<vl::Vec3d> vertexPositions( modesStructure->vertices.size() );
			for( size_t iVertex = 0; iVertex < modesStructure->vertices.size(); ++iVertex )
			{
				const SmallIndexTriple& vert = modesStructure->vertices[iVertex];
				vl::Vec3d position = bmin + 
					vl::Vec3d( 
						boost::numeric_cast<double>(vert.get<0>())*increment[0], 
						boost::numeric_cast<double>(vert.get<1>())*increment[1],
						boost::numeric_cast<double>(vert.get<2>())*increment[2] );
				vertexPositions[iVertex] = position;

				ofs << "   " << (iVertex+1) << ",   " 
					<< position[0] << ",   " 
					<< position[1] << ",   " 
					<< position[2] << "\n";
			}

			// vertices done, now for the elements
			ofs << "*ELEMENT, TYPE=C3D8, ELSET=B1\n";

			// we will need to keep track of which elements
			// have which material
			typedef std::pair<size_t, size_t> EltRange;
			std::deque<EltRange> eltRanges;

			QuadFaceSet externalFaces;

			// generate mass matrix
			fortran_matrix eltMass = barbic::ReadMatrixFromDisk("voxel.mass");
			const size_t numElementVertices = 8;
			barbic::SparseMatrix<double> massMatrix( 3*vertexPositions.size(), 3*vertexPositions.size() );
			assert( eltMass.nrows() == numElementVertices );
			assert( eltMass.ncols() == numElementVertices );
			double voxelSpacing3 = increment[0] * increment[1] * increment[2];

			// elements are 1-indexed
			size_t iElement = 1;
			for( size_t iMat = 0; iMat < materials.size(); ++iMat )
			{
				const size_t rangeStart = iElement;
				for( size_t i = 0; i < resolution[0]; ++i )
					for( size_t j = 0; j < resolution[1]; ++j )
						for( size_t k = 0; k < resolution[2]; ++k )
						{
							if( voxels(i, j, k) == (iMat+1) )
							{
								ofs << iElement;
								// locate all 8 vertices
								boost::array<size_t, 8> v;
								v[0] = findVertex( i+0, j+0, k+0, modesStructure->vertices );
								v[1] = findVertex( i+1, j+0, k+0, modesStructure->vertices );
								v[2] = findVertex( i+1, j+1, k+0, modesStructure->vertices );
								v[3] = findVertex( i+0, j+1, k+0, modesStructure->vertices );
								v[4] = findVertex( i+0, j+0, k+1, modesStructure->vertices );
								v[5] = findVertex( i+1, j+0, k+1, modesStructure->vertices );
								v[6] = findVertex( i+1, j+1, k+1, modesStructure->vertices );
								v[7] = findVertex( i+0, j+1, k+1, modesStructure->vertices );

								for( size_t i = 0; i < 8; ++i )
									ofs << "," << (v[i]+1);
								ofs << "\n";

								// build .obj file
								addExternalFace(externalFaces,v[0],v[3],v[2],v[1]);
								addExternalFace(externalFaces,v[4],v[5],v[6],v[7]);
								addExternalFace(externalFaces,v[0],v[1],v[5],v[4]);
								addExternalFace(externalFaces,v[3],v[7],v[6],v[2]);
								addExternalFace(externalFaces,v[1],v[2],v[6],v[5]);
								addExternalFace(externalFaces,v[0],v[4],v[7],v[3]);
								
								// build mass matrix
								double density = materials.at(iMat)->density();
								for( size_t i = 0; i < numElementVertices; ++i )
									for( size_t j = 0; j < numElementVertices; ++j )
										for( size_t kComponent = 0; kComponent < 3; ++kComponent )
											massMatrix.add( 3*v[i] + kComponent, 3*v[j] + kComponent, 
												eltMass(i, j)*density*voxelSpacing3 );

								++iElement;
							}
						}

				// if start > end, range is empty
				// this is safe because start is always >= 1
				const size_t rangeEnd = iElement-1;
				eltRanges.push_back( EltRange(rangeStart, rangeEnd) );
			}

			ObjFile::Group objGroup( "squashingCubes" );
			std::vector<vl::Vec3d> objNormals;
			objNormals.reserve( externalFaces.size() );
			for( QuadFaceSet::const_iterator faceItr = externalFaces.begin();
				faceItr != externalFaces.end(); ++faceItr )
			{
				typedef ObjFile::Face<unsigned int> Face;
				typedef ObjFile::Vertex<unsigned int> Vertex;

				const boost::array<size_t, 4>& index = faceItr->vertices;
				vl::Vec3d normal = vl::cross( 
					toVec3d(vertexPositions[ index[1] ]) - toVec3d(vertexPositions[ index[0] ]),
					toVec3d(vertexPositions[ index[2] ]) - toVec3d(vertexPositions[ index[0] ]) );
				objNormals.push_back( normal );

				std::pair< bool, unsigned int > texPos(false, 0);
				std::pair< bool, unsigned int > objNormal(true, objNormals.size());
				Face newFace1( 
					Vertex( index[0]+1, texPos, objNormal ),
					Vertex( index[1]+1, texPos, objNormal ),
					Vertex( index[2]+1, texPos, objNormal ) );
				Face newFace2( 
					Vertex( index[2]+1, texPos, objNormal ),
					Vertex( index[3]+1, texPos, objNormal ),
					Vertex( index[0]+1, texPos, objNormal ) );

				objGroup.addFace( newFace1 );
				objGroup.addFace( newFace2 );
			}

			ObjFile objFile( objGroup, 
				vertexPositions,
				objNormals,
				std::vector<vl::Vec3d>() );
			objFile.write( "test.obj" );

			// dump out 1 element range for each material
			for( size_t iRange = 0;	iRange < eltRanges.size(); ++iRange )
			{
				if( eltRanges[iRange].first > eltRanges[iRange].second )
					continue;

				ofs << "*ELSET,ELSET=E" << iRange << ",GENERATE\n";
				ofs << "   " << eltRanges[iRange].first << "," << eltRanges[iRange].second << "\n";
			}

			// now, dump out all materials and match them with element sets:
			for( size_t iMaterial = 0; iMaterial < materials.size(); ++iMaterial )
			{
				const Material& mat = *materials[iMaterial];
				ofs << "*MATERIAL,NAME=" << mat.name() << "\n";
				ofs << "*ELASTIC\n";
				ofs << mat.youngsModulus() << "," << mat.poissonRatio() << "\n";
				ofs << "*DENSITY\n";
				ofs << mat.density() << "\n";

				// match element set with material:
				ofs << "*SOLID SECTION,MATERIAL=" << mat.name() 
					<< ",ELSET=E" << iMaterial << "\n";
			}

			const size_t nModes = std::min<size_t>( 200, modesStructure->vertices.size() / 2 );
			ofs << "*STEP\n";
			ofs << "*FREQUENCY\n";
			ofs << nModes << ",100.0,20000.0\n";
			ofs << "*NODE FILE\n";
			ofs << "U\n";
			ofs << "*END STEP\n";
			ofs.close();

			runCalculix( prefix, event );

			boost::tie( modesStructure->frequencies, modesStructure->modes ) = 
				parseFRDFile( prefix + ".frd" );

			// need to normalize columns
			assert( modesStructure->modes.nrows() == 3*vertexPositions.size() );
			for( size_t iMode = 0; iMode < modesStructure->modes.ncols(); ++iMode )
			{
				double* mode = modesStructure->modes.data() + iMode*modesStructure->modes.nrows();
				std::vector<double> massScaledMode = 
					massMatrix.multiply( mode );
				double val = cblas_ddot( modesStructure->modes.nrows(), mode, 1, &massScaledMode[0], 1 );

				// rescale the mode:
				cblas_dscal( modesStructure->modes.nrows(), 1.0/sqrt(val), mode, 1 );
			}

			object->setSquashingCubesModes( modesStructure );

			return modesStructure;

			/*
			barbic::WriteMatrixToDisk( "test.Ulin", modesStructure->modes );
			{
				// write to file for Jernej's reader
				std::ofstream ofs( "test.Ulin.ascii" );
				if( !ofs )
					throw IOException( "Unable to open file 'test.Ulin.ascii' for writing." );

				ofs << modes.nrows() << "\n";
				ofs << modes.ncols();
				for( size_t iRow = 0; iRow < modes.nrows(); ++iRow )
				{
					ofs << "\n";
					for( size_t iCol = 0; iCol < modes.ncols(); ++iCol )
					{
						ofs << modes(iRow, iCol) << " ";
					}
				}
			}
			*/
		}
	}
}

void runCalculix( const std::string& inputFile, Event& event )
{
	STARTUPINFO si;
	PROCESS_INFORMATION pi;
	ZeroMemory( &si, sizeof(si) );
	si.cb = sizeof(si);
	ZeroMemory( &pi, sizeof(pi) );


	// todo make this somehow user-configurable
	std::string calculixLoc( "E:\\Program Files\\CalculiX\\ccx_1.5\\ccx_1.5.exe" );
	std::string commandLine( "\"" + calculixLoc + "\"" );
	commandLine += " " + inputFile;

	std::vector<char> commandLineVec( commandLine.begin(), commandLine.end() );
	commandLineVec.push_back( 0 );

	if( !CreateProcess(
		calculixLoc.c_str(),
		&commandLineVec[0],
		NULL,                // Process handle not inheritable. 
		NULL,                // Thread handle not inheritable. 
		FALSE,               // Set handle inheritance to FALSE. 
		IDLE_PRIORITY_CLASS, // Set to very low priority
		NULL,                // Use parent's environment block. 
		NULL,                // Use parent's starting directory. 
		&si,                 // Pointer to STARTUPINFO structure.
		&pi )                // Pointer to PROCESS_INFORMATION structure.
		)
	{
		std::ostringstream oss;
		oss << "Unable to execute process '" << &commandLineVec[0] << "'." << std::endl;
		oss << "CreateProcess returned error: " << GetLastError() << std::endl;
		throw CreationException( oss.str() );
	}

	HANDLE handles[] = { event.handle(), pi.hProcess };

	// Wait until child process exits.
	DWORD res = WaitForMultipleObjects( 2, handles, FALSE, INFINITE );

	if( res == WAIT_OBJECT_0 )
	{
		BOOL result = TerminateProcess( pi.hProcess, 0 );
		if( result == 0 )
		{
			DWORD err = GetLastError();
		}

		CloseHandle( pi.hProcess );
		CloseHandle( pi.hThread );
		throw Exception( "User interrupt" );
	}

	CloseHandle( pi.hProcess );
	CloseHandle( pi.hThread );
}

std::pair< std::vector<double>, fortran_matrix > parseFRDFile( const std::string& filename )
{
	using namespace boost::spirit;

	std::ifstream ifs( filename.c_str() );
	if( !ifs )
		throw IOException( "Unable to open file '" + filename + "' for reading." );

	std::deque< std::vector<double> > modes;
	std::vector<double> frequencies;

	int currentBlock = 0;
	int lineNum = 0;
	while( ifs )
	{
		// maximum line length is 1024:
		char line[1024];
		ifs.getline( line, sizeof(line) );

		size_t len = ifs.gcount();
		if( len < 6 )
			continue;

		if( line[5] == 'C' )
		{
			if( currentBlock == 100 )
			{
				// check that the mode is ok
				assert( (modes.front().size() % 3) == 0 );
				assert( !modes.empty() );
				if( modes.back().empty() )
					throw ParserException( "Empty mode", filename, lineNum-1, 0 );

				if( modes.back().size() != modes.front().size() )
				{
					std::ostringstream oss;
					oss << "Mismatch in vertex counts; expected " << 
						modes.front().size()/3 << " vertices, but found " <<
						modes.back().size()/3 << ".";
					throw ParserException( oss.str(), filename, lineNum-1, 0 );
				}
			}

			parse_info<char const*> info = parse( &line[1], &line[5],
				int_p[assign(currentBlock)], space_p );
			if( !info.hit )
				throw ParserException("Missing type specifier", filename, lineNum, 1);

			if( currentBlock == 100 )
			{
				double frequency;
				info = parse( &line[12], &line[24],
					real_p[assign(frequency)], space_p );
				if( !info.hit )
					throw ParserException("Invalid frequency", filename, lineNum, 12);

				// Calculix frequencies are in cycles/sec, convert to radians/sec
				frequencies.push_back( frequency * (2.0*M_PI) );

				modes.push_back( std::vector<double>() );
				modes.back().reserve( modes.front().size() );
			}
		}
		else if( isspace(line[1]) && isspace(line[2]) )
		{
		}
		else if( currentBlock == 100 )
		{
			int type;
			parse_info<char const*> info = parse( &line[1], &line[3],
				int_p[assign(type)], space_p );
			if( !info.hit )
				throw ParserException("Missing type specifier", filename, lineNum, 1);

			if( type == -1 )
			{
				for( size_t i = 0; i < 3; ++i )
				{
					size_t start = 13 + 12*i;
					size_t end = start + 12;

					double val;
					info = parse( &line[start], &line[end],
						real_p[assign(val)], space_p );
					if( !info.hit )
						throw ParserException("Expected: floating point", filename, lineNum, start);

					modes.back().push_back( val );
				}
			}
		}

		++lineNum;
	}

	if( modes.size() == 0 )
		throw ParserException( "Didn't find any modes", filename, lineNum, 0 );

	// tossing out the rigid body modes
	/*
	modes.erase( modes.begin(), modes.begin() + 6 );
	frequencies.erase( frequencies.begin(), frequencies.begin() + 6 );
	*/

	fortran_matrix result( modes.front().size(), modes.size() );
	for( size_t iMode = 0; iMode < modes.size(); ++iMode )
	{
		assert( modes[iMode].size() == result.nrows() );
		for( size_t jVert = 0; jVert < modes[iMode].size(); ++jVert )
			result(jVert, iMode) = (modes[iMode])[jVert];
	}


	assert( result.ncols() == frequencies.size() );

	return std::make_pair( frequencies, result );
}

AudioGenerator::AudioGenerator( size_t rate, 
	const std::vector<double>& frequencies, 
	double alpha, double beta )
	: rate_( rate ), dt_( 1.0/boost::numeric_cast<double>(rate) )
{
	oscillators_.reserve( frequencies.size() );
	for( size_t iFrequency = 0; iFrequency < frequencies.size(); ++iFrequency )
	{
		double omega = frequencies[iFrequency];

		// make sure we don't retain modes greater than the Nyquist rate
		double cyclesPerSecond = omega / (2.0*M_PI);
		if( cyclesPerSecond > rate/2 )
			break;

		double damping = 0.5*(alpha/omega + beta*omega);
		if( damping >= 1.0 )
			break;

		oscillators_.push_back( 
			barbic::HarmonicOscillator_IIR(
				1.0,         // mass
				damping,       // damping parameter (todo: change it)
				omega*omega, // frequency is sqrt(stiffness/mass), so stiffness is frequency squared
				dt_ ) );
	}
}

AudioGenerator::~AudioGenerator()
{
}

void AudioGenerator::addImpulses( const std::vector<double>& impulses, double width )
{
	assert( impulses.size() >= oscillators_.size() );
    
	impulses_.push_back( GaussianImpulse(impulses, width/3.0) );
}

std::vector<double> AudioGenerator::integrate( size_t numSteps )
{
	std::vector<double> result( numSteps );
	for( size_t i = 0; i < numSteps; ++i )
	{
		// clear force accumulators
		std::for_each( oscillators_.begin(), oscillators_.end(),
			boost::bind( &barbic::HarmonicOscillator_IIR::SetExternalForce, _1, 0.0 ) );

		for( ImpulseList::iterator impulseItr = impulses_.begin(); 
			impulseItr != impulses_.end(); )
		{
			const double x = impulseItr->currentPosition;
			const double sigma = impulseItr->stdDev;

			double gaussianVal = (1.0 / (sigma*2.0*M_PI)) *
				exp( -(x*x) / (2.0*sigma*sigma) );

			for( size_t iOscillator = 0; iOscillator < oscillators_.size(); ++iOscillator )
				oscillators_[iOscillator].AddExternalForce( gaussianVal * impulseItr->impulses[iOscillator] );

			impulseItr->currentPosition += this->dt_;
			if( impulseItr->currentPosition >= 3.0*impulseItr->stdDev )
				impulseItr = impulses_.erase(impulseItr);
			else
				++impulseItr;
		}

		double val = 0.0;
		for( std::vector<barbic::HarmonicOscillator_IIR>::iterator oscItr = oscillators_.begin();
			oscItr != oscillators_.end(); ++oscItr )
		{
			val += oscItr->Getq();
			oscItr->DoTimestep();
		}

		result[i] = val;
	}

	std::for_each( oscillators_.begin(), oscillators_.end(), 
		boost::bind( &barbic::HarmonicOscillator_IIR::ChopZeroSignal, _1, 1e-10 ) );

	return result;
}

SquashingCubesAudioGenerator::SquashingCubesAudioGenerator( 
	size_t rate,
	boost::shared_ptr<SquashingCubesModes> modes,
	double alpha,
	double beta )
	:	AudioGenerator(rate, modes->frequencies, alpha, beta),
		modes_(modes)
{
	assert( modes->frequencies.size() == modes->modes.ncols() );
}

void SquashingCubesAudioGenerator::addImpulse( const vl::Vec3f& position, const vl::Vec3f& force, float dt )
{
	// first step is to find a nearby voxel
	TinyVec<size_t, 3> voxel = modes_->voxels.voxel( position );
	// todo: do something better in this case, like snap to the
	//   nearest actual vertex
	if( modes_->voxels(voxel) == 0 )
		return;

	BoundingBox3f box = modes_->voxels.box( voxel );
	const vl::Vec3f offset = position - box.minimum();
	const vl::Vec3f weights = offset / (box.maximum() - box.minimum());

	double dtRatio = dt / this->dt();

	vl::Vec3f impulse = force;
	std::vector<double> modeImpulses( modes_->modes.ncols(), 0.0 );
	for( size_t i = 0; i < 8; ++i )
	{
		// pick one of the 8 cube vertices
		TinyVec<size_t, 3> vertexOffsets(
			((i/1)%2 > 0) ? 1 : 0,
			((i/2)%2 > 0) ? 1 : 0,
			((i/4)%2 > 0) ? 1 : 0 );

		size_t iVertex = findVertex( voxel[0] + vertexOffsets[0], 
			voxel[1] + vertexOffsets[1], 
			voxel[2] + vertexOffsets[2], 
			modes_->vertices );
		double weight = 1.0;
		for( size_t i = 0; i < 3; ++i )
			weight *= (vertexOffsets[i] > 0) ? (weights[i]) : (1.0 - weights[i]);

		for( size_t iComponent = 0; iComponent < 3; ++iComponent )
			for( size_t iMode = 0; iMode < modeImpulses.size(); ++iMode )
				modeImpulses[iMode] += fabs( weight * dtRatio
					* modes_->modes( 3*iVertex + iComponent, iMode ) 
					* impulse[iComponent] );
	}

	this->addImpulses( modeImpulses, 5.0*this->dt() );
}

struct FormatChunk
{
  FormatChunk()
  {
	  const char fmtText[] = "fmt ";
	  std::copy( fmtText, fmtText + 4, this->chunkID );
	  this->chunkSize = sizeof( FormatChunk ) - 8;
  }

  char           chunkID[4];
  long           chunkSize;

  short          wFormatTag;
  unsigned short wChannels;
  unsigned long  dwSamplesPerSec;
  unsigned long  dwAvgBytesPerSec;
  unsigned short wBlockAlign;
  unsigned short wBitsPerSample;

/* Note: there may be additional fields here, depending upon wFormatTag. */

};

// WAV file details from here:
// http://www.borg.com/~jglatt/tech/wave.htm
template <typename IntegerType>
void createWavFile( const std::string& filename, const std::vector<IntegerType>& samples )
{
	try
	{
		std::ofstream ofs( filename.c_str(), std::ios::binary );
		if( !ofs )
			throw IOException( "Unable to open file '" + filename + "' for writing." );
		ofs.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );

		const char riffText[] = "RIFF";
		long fullChunkSize = sizeof(FormatChunk) + (8 + samples.size() * sizeof(IntegerType));
		ofs.write( riffText, 4 );
		ofs.write( reinterpret_cast<const char*>(&fullChunkSize), sizeof(long) );

		const char waveText[] = "WAVE";

		// header:
		ofs.write( waveText, 4 );

		// format chunk
		FormatChunk fmt;
		fmt.wFormatTag = 1;
		fmt.wChannels = 2;
		fmt.dwSamplesPerSec = 44100;
		fmt.wBitsPerSample = sizeof(IntegerType) * 8;
		fmt.wBlockAlign = fmt.wChannels * sizeof(IntegerType);
		fmt.dwAvgBytesPerSec = fmt.dwSamplesPerSec * fmt.wBlockAlign;
		ofs.write( reinterpret_cast<const char*>( &fmt ), sizeof( FormatChunk ) );

		// data chunk
		const char dataText[] = "data";
		ofs.write( dataText, 4 );
		unsigned long chunkSize = samples.size() * sizeof(IntegerType);
		ofs.write( reinterpret_cast<const char*>( &chunkSize ), sizeof(unsigned long) );
		ofs.write( reinterpret_cast<const char*>( &samples[0] ), sizeof(IntegerType) * samples.size() );
	}
	catch( std::ifstream::failure& e )
	{
		std::ostringstream oss;
		oss << "Error writing to '" << filename << "': " << e.what();
		throw IOException( oss.str() );
	}
}

template void createWavFile<unsigned char>( const std::string&, const std::vector<unsigned char>& );
template void createWavFile<short>( const std::string&, const std::vector<short>& );
template void createWavFile<int>( const std::string&, const std::vector<int>& );

bool CompareAudioHash::operator()( const AudioHash& left, const AudioHash& right ) const
{
	if( left.first < right.first )
		return true;
	else if( right.first > left.first )
		return false;

	for( size_t i = 0; i < 6; ++i )
	{
		if( !withinEpsilon( left.second[i], right.second[i] ) )
			return left.second[i] < right.second[i];
	}

	return false;
}

AudioHash hashForAudio( ConstPhysicsObjectPtr physicsObject )
{
	boost::shared_ptr<MechModelDeclaredObject> object = physicsObject->toMechModel( MECH_MODEL_AUDIO_ROOT );

	vl::Mat3f transform = toMat3f( physicsObject->transform() );
	vl::Mat3f shear = polarDecomposition( toMat3f(transform) ).second;
	boost::array<float, 6> values = 
		{{ shear[0][0], shear[1][0], shear[1][1], shear[2][0], shear[2][1], shear[2][2] }};

	/*
	// don't actually want to do this since all relevant information is captured
	//   by the scale parameters
	// now need to add the hierarchy up to the root as well.
	ConstSceneGraphElementPtr current = physicsObject;
	while( current->parent().lock() )
	{
		ConstSceneGraphElementPtr parent = current->parent().lock();
		boost::shared_ptr<MechModelDeclaredObject> localObject = 
			parent->toMechModel( MECH_MODEL_AUDIO );
		localObject->add( object );

		object = localObject;
		current = parent;
	}
	*/

	std::ostringstream oss;
	object->dump( oss, "" );
	return AudioHash(oss.str(), values);
}

} // namespace planning

#endif // SOUND
