#include "stdafx.h"

#include "twigg/linalg.h"
#include "twigg/TRAFile.h"

#include <sstream>
#include <cassert>

TRAFile::TRAFile( const std::string& filename )
	: filename_(filename)
{
	std::ifstream ifs( filename.c_str(), std::ios::binary );
	if( !ifs )
		throw CreationException( "Unable to open file '" + filename +
			"' for reading." );

	ifs.exceptions(
		std::ifstream::eofbit |
		std::ifstream::failbit |
		std::ifstream::badbit );

	try
	{
		{
			int endian;
			int numVertices;
			int numCoeffs;
			int numEffectiveVertices;

			ifs.read( reinterpret_cast<char*>(&endian), sizeof(int) );
			ifs.read( reinterpret_cast<char*>(&numVertices), sizeof(int) );
			ifs.read( reinterpret_cast<char*>(&numCoeffs), sizeof(int) );
			ifs.read( reinterpret_cast<char*>(&numEffectiveVertices), sizeof(int) );

			this->nVertices_ = numVertices;
			this->nEffectiveVertices_ = numEffectiveVertices;
			this->nCoeffs_ = numCoeffs;

			if( endian != 1 )
			{
				std::ostringstream oss;
				oss << "Read " << endian << " instead of 1; "
					<< "wrong endian-ness?";
				throw CreationException( oss.str() );
			}
		}

		for( unsigned int iComponent = 0; iComponent < 4; ++iComponent )
			components_[iComponent].resize( this->nEffectiveVertices_, this->nCoeffs_ );

		for( unsigned int iVertex = 0; iVertex < this->nEffectiveVertices_; ++iVertex )
		{
			for( unsigned int iComponent = 0; iComponent < 4; ++iComponent )
			{
				std::vector<float> coeffs;
				readArray( coeffs, ifs, this->nCoeffs_ );
				for( unsigned int iCoeff = 0; iCoeff < coeffs.size(); ++iCoeff )
					(components_[iComponent])( iVertex, iCoeff ) 
						= coeffs[iCoeff];
			}
		}

		monochrome_.resize( nEffectiveVertices_, nCoeffs_ );
		std::fill( monochrome_.data(), monochrome_.data() + 
			monochrome_.nrows()*monochrome_.ncols(), 0.0 );
		add( scaled(components_[0], 0.299), monochrome_ );
		add( scaled(components_[1], 0.587), monochrome_ );
		add( scaled(components_[2], 0.114), monochrome_ );
	}
	catch( std::ifstream::failure& )
	{
		throw CreationException( 
			std::string("Caught iostream exception while trying to ")
			+ std::string("read '") + filename + "'." );
	}
}

TRAFile::TRAFile( const ComponentArray& components, unsigned int actualVertices )
	: components_( components ), nVertices_( actualVertices ) 
{
	assert( this->nVertices_ > 0 );

	this->nCoeffs_ = components[0].ncols();
	assert( this->nCoeffs_ > 0 );

	this->nEffectiveVertices_ = components[0].nrows();
	assert( this->nEffectiveVertices_ > 0 );

	if( nEffectiveVertices_ < nVertices_ )
	{
		std::ostringstream oss;
		oss << "Effective vertex count must be at least as large as vertex count."
			<< "  Found effective vertex count of " << nEffectiveVertices_ 
			<< " but actual vertex count " << nVertices_;
		throw CreationException( oss.str() );
	}

	// verification:
	for( unsigned int iComponent = 0; iComponent < 4; ++iComponent )
	{
		if( components_[iComponent].nrows() != this->nEffectiveVertices_ )
		{
			std::ostringstream oss;
			oss << "Vertex count of " << components_[iComponent].nrows() << " in component " 
				<< iComponent << " does not match vertex count of " << components[0].nrows() 
				<< " in component 0.";
			throw CreationException( oss.str() );
		}

		if( components_[iComponent].ncols() != this->nCoeffs_ )
		{
			std::ostringstream oss;
			oss << "Coeffs count of " << components_[iComponent].ncols() << " in component "
				<< iComponent << " does not match coeffs count of " << components[0].nrows()
				<< " in component 0.";
			throw CreationException( oss.str() );
		}
	}

	monochrome_.resize( nEffectiveVertices_, nCoeffs_ );
	std::fill( monochrome_.data(), monochrome_.data() + 
		monochrome_.nrows()*monochrome_.ncols(), 0.0 );
	add( scaled(components_[0], 0.299), monochrome_ );
	add( scaled(components_[1], 0.587), monochrome_ );
	add( scaled(components_[2], 0.114), monochrome_ );
}

void TRAFile::write( const std::string& filename ) const
{
	std::ofstream ofs( filename.c_str(), std::ios::binary );
	if( !ofs )
		throw IOException( "Unable to open file '" + filename +
			"' for writing." );

	ofs.exceptions(
		std::ifstream::eofbit |
		std::ifstream::failbit |
		std::ifstream::badbit );

	try
	{
		int endian = 1;
		int numVertices = this->nVertices_;
		int numCoeffs = this->nCoeffs_;
		int numEffectiveVertices = this->nEffectiveVertices_;

		ofs.write( reinterpret_cast<const char*>(&endian), sizeof(int) );
		ofs.write( reinterpret_cast<const char*>(&numVertices), sizeof(int) );
		ofs.write( reinterpret_cast<const char*>(&numCoeffs), sizeof(int) );
		ofs.write( reinterpret_cast<const char*>(&numEffectiveVertices), sizeof(int) );

		for( unsigned int iVertex = 0; iVertex < this->nEffectiveVertices_; ++iVertex )
		{
			for( unsigned int iComponent = 0; iComponent < 4; ++iComponent )
			{
				std::vector<float> coeffs( numCoeffs );
				for( unsigned int iCoeff = 0; iCoeff < coeffs.size(); ++iCoeff )
					coeffs[iCoeff] = (components_[iComponent])( iVertex, iCoeff );

				ofs.write( reinterpret_cast<const char*>(&coeffs[0]), sizeof(float)*coeffs.size() );
			}
		}
	}
	catch( std::ifstream::failure& )
	{
		throw IOException( 
			std::string("Caught iostream exception while trying to ")
			+ std::string("write '") + filename + "'." );
	}
}

