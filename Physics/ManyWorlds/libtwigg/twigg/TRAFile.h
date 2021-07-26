#ifndef __TRAFILE_H__
#define __TRAFILE_H__

#include "twigg/linalg.h"
#include <boost/array.hpp>

class TRAFile
{
public:
	typedef boost::array<fortran_matrix, 4> ComponentArray;

	TRAFile( const std::string& filename );
	TRAFile( const ComponentArray& array, unsigned int actualVertices );

	const fortran_matrix& component( unsigned int iComponent ) const
	{
		return components_[iComponent];
	}

	const fortran_matrix& monochrome() const
	{
		return monochrome_;
	}

	void write( const std::string& filename ) const;

	unsigned int numCoefficients() const
	{
		return nCoeffs_;
	}

	unsigned int numVertices() const
	{
		return nVertices_;
	}

	unsigned int numEffectiveVertices() const
	{
		return nEffectiveVertices_;
	}

private:
	boost::array<fortran_matrix, 4> components_;
	fortran_matrix monochrome_;

	unsigned int nVertices_;
	unsigned int nEffectiveVertices_;
	unsigned int nCoeffs_;

	std::string filename_;
};


#endif


