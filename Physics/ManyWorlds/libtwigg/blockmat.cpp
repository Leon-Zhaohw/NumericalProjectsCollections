#include "stdafx.h"

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <stack>
#include <map>
#include <algorithm>
#include <numeric>
#include <functional>
#include <boost/mem_fn.hpp>
#include "twigg/blockmat.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#endif

using namespace std;
using namespace vl;

//#define DEBUG_FORCES

BlockMatrix::BlockMatrix(int size) 
	: rows_( size )
{
	clearEntries();
}

void BlockMatrix::clear()
{
	for( unsigned int iRow = 0; iRow < rows_.size(); ++iRow )
	{
		for( EntryVec::iterator itr = rows_[iRow].begin();
			itr != rows_[iRow].end();
			++itr )
		{
			itr->second = vl::vl_0;
		}
	}
}

void BlockMatrix::clearEntries()
{
	std::for_each( rows_.begin(), rows_.end(),
		boost::mem_fn( &EntryVec::clear ) );

	// Some of the algorithms depend on every row having at least one
	// entry; we'll stick them along the diagonal.
	for( unsigned int iRow = 0; iRow < rows_.size(); ++iRow )
		rows_[iRow].push_back( Entry( iRow, vl::Mat3(vl::vl_0) ) );
}

// Utility class for binary search
class EntryComparator
{
public:
	bool operator() (const BlockMatrix::Entry& left, const BlockMatrix::Entry& right) const {
		return left.first < right.first;
	}
};

const Mat3 BlockMatrix::operator() (unsigned int i, unsigned int j) const
{
	assert( i < size() );
	assert( j < size() );

	EntryVec::const_iterator itr = 
		std::lower_bound(rows_[i].begin(), rows_[i].end(), std::make_pair(j, vl::Mat3()), EntryComparator());
	if (itr == rows_[i].end() || itr->first != j)
		return vl::Mat3(vl::vl_0);
	else
		return itr->second;
}

Mat3& BlockMatrix::operator() (unsigned int i, unsigned int j)
{
	assert( i < size() );
	assert( j < size() );

	EntryVec::iterator itr = 
		std::lower_bound(rows_[i].begin(), rows_[i].end(), std::make_pair(j, vl::Mat3()), EntryComparator());
	if (itr == rows_[i].end() || itr->first != j)
		// Need to insert an entry
		return (rows_[i].insert(itr, std::make_pair(j, vl::vl_0)))->second;
	else
		return itr->second;
}

void BlockMatrix::add( unsigned int i, unsigned int j, const vl::Mat3& value )
{
	assert( !isBad(value) );
	(*this)(i, j) += value;
}

void BlockMatrix::subtract( unsigned int i, unsigned int j, const vl::Mat3& value )
{
	assert( !isBad(value) );
	(*this)(i, j) -= value;
}

const BlockMatrix& BlockMatrix::operator*=(const double rhs)
{
	for (unsigned int i = 0; i < size(); ++i) {
		for (EntryVec::iterator itr = rows_[i].begin(); 
			itr != rows_[i].end();
			++itr)
			itr->second *= rhs;
	}
	return *this;
}

const BlockMatrix& BlockMatrix::operator+=(const BlockMatrix& rhs)
{
	assert(size() == rhs.size());
	for (unsigned int i = 0; i < size(); ++i) {
		EntryVec::iterator lhsIter = rows_[i].begin();
		EntryVec::const_iterator rhsIter = rhs.rows_[i].begin();
		while (rhsIter != rhs.rows_[i].end()) {
			if (lhsIter == rows_[i].end() || 
				lhsIter->first < rhsIter->first) {
				rows_[i].insert(lhsIter, *rhsIter);
				++rhsIter;
			}
			else if (lhsIter->first == rhsIter->first) {
				lhsIter->second += rhsIter->second;
				++lhsIter;
				++rhsIter;
			}
		}
	}

	return *this;
}

const BlockMatrix& BlockMatrix::operator-=(const BlockMatrix& rhs)
{
	assert(size() == rhs.size());
	for (unsigned int i = 0; i < size(); ++i) {
		EntryVec::iterator lhsIter = rows_[i].begin();
		EntryVec::const_iterator rhsIter = rhs.rows_[i].begin();
		while (rhsIter != rhs.rows_[i].end()) {
			if (lhsIter == rows_[i].end() || 
				lhsIter->first < rhsIter->first) {
				(rows_[i].insert(lhsIter, *rhsIter))->second *= -1;
				++rhsIter;
			}
			else if (lhsIter->first == rhsIter->first) {
				lhsIter->second -= rhsIter->second;
				++lhsIter;
				++rhsIter;
			}
		}
	}

	return *this;
}

void BlockMatrix::printDense(std::ostream& os) const
{
	size_t length = numberLength();
	os << "[ ";
	for (unsigned int i = 0; i < size(); ++i) {
		if (i != 0)
			os << "  ";
		printRow(os, i, 0, length);   os << std::endl << "  ";
		printRow(os, i, 1, length);   os << std::endl << "  ";
		printRow(os, i, 2, length);   os << std::endl;
	}
	os << "]" << std::endl;
}

size_t BlockMatrix::numberLength() const
{
	size_t maxLength = 1;
	for (unsigned int k = 0; k < size(); k++) {
		for (EntryVec::const_iterator itr = rows_[k].begin();
		itr != rows_[k].end();
		++itr) {
			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < 3; ++j) {
					std::ostringstream s;
					s << (itr->second)[j][i] << std::flush;
					if (s.str().length() > maxLength)
						maxLength = s.str().length();
				}
			}
		}
	}

	return maxLength;
}

void BlockMatrix::printRow(std::ostream& os, int rowBlock, const unsigned int rowNum, const unsigned int numberLength) const
{
	const static int entrySeparation = 1;
	for (unsigned int i = 0; i < size(); ++i) {
		const Mat3& matrix = (*this)(rowBlock, i);
		for (unsigned int j = 0; j < 3; ++j) {
			std::ostringstream sstream;
			sstream << matrix[rowNum][j] << std::flush;
			std::string s(sstream.str());
			while (s.length() < numberLength + entrySeparation)
				s += " ";
			os << s;
		}
	}
}

unsigned int BlockMatrix::nnz() const
{
	unsigned int nnz = 0;
	for( size_t iRow = 0; iRow < rows_.size(); ++iRow )
		nnz += rows_[iRow].size();
	return nnz;
}

fortran_matrix BlockMatrix::toMatrix() const
{
	fortran_matrix result( 3*size(), 3*size() );
	for( unsigned int iBlockRow = 0; iBlockRow < size(); ++iBlockRow )
	{
		for( std::vector<Entry>::const_iterator blockColItr = rows_[iBlockRow].begin();
			blockColItr != rows_[iBlockRow].end();
			++blockColItr )
		{
			for( unsigned int iCol = 0; iCol < 3; ++iCol )
			{
				const unsigned int actualCol = 3*blockColItr->first + iCol;
				for( int iRow = 0; iRow < 3; ++iRow )
					result( 3*iBlockRow + iRow, actualCol ) = (blockColItr->second)[iRow][iCol];
			}
		}
	}

	return result;
}

void BlockMatrix::toSymmetricCoordinateStorage( long* irow, long* icol, double* a, unsigned int* rowPtrs ) const
{
	typedef long Integer;
	typedef double Real;

	Real* aStart = a;
	for( unsigned int iRow = 0; iRow < rows_.size(); ++iRow )
	{
		unsigned int length = 0;
		{
			for( std::vector<Entry>::const_iterator colItr = rows_[iRow].begin();
				colItr != rows_[iRow].end() && colItr->first < iRow; 
				++colItr )
				++length;
		}

		length *= 3;

		boost::array<Integer*, 3> irows = { { irow, irow + length + 1, irow + length + length + 3 } };
		boost::array<Integer*, 3> icols = { { icol, icol + length + 1, icol + length + length + 3 } };
		boost::array<Real*, 3> as = { { a, a + length + 1, a + length + length + 3 } };

		for( unsigned int i = 0; i < 3; ++i )
			*rowPtrs++ = as[i] - aStart;

		if( rows_[iRow].empty() )
		{
			std::ostringstream oss;
			oss << "Matrix is singular; block row " << iRow << " is empty.";
			throw CreationException( oss.str() );
		}

		std::vector<Entry>::const_iterator colItr = rows_[iRow].begin();
		// First do all the off-diagonal entries.
		while( colItr != rows_[iRow].end() && colItr->first < iRow )
		{
			const vl::Mat3 symmetricPr = (*this)( colItr->first, iRow );
			// This is because output is not blocked.
			for( unsigned int iSplitRow = 0; iSplitRow < 3; ++iSplitRow )
			{
				for( unsigned int iSplitCol = 0; iSplitCol < 3; ++iSplitCol )
				{
					*as[iSplitRow]++ = 0.5*(colItr->second[iSplitRow][iSplitCol] + symmetricPr[iSplitCol][iSplitRow]);
					*irows[iSplitRow]++ = 3*iRow + iSplitRow + 1;
					*icols[iSplitRow]++ = 3*colItr->first + iSplitCol + 1;
				}
			}

			++colItr;
		}

		// Now, do the diagonal as a special case
		if( colItr != rows_[iRow].end() && colItr->first == iRow )
		{
			for( unsigned int iSplitRow = 0; iSplitRow < 3; ++iSplitRow )
			{
				for( unsigned int iSplitCol = 0; iSplitCol <= iSplitRow; ++iSplitCol )
				{
					*as[iSplitRow]++ = 0.5*(colItr->second[iSplitRow][iSplitCol] + colItr->second[iSplitCol][iSplitRow]);
					*irows[iSplitRow]++ = 3*iRow + iSplitRow + 1;
					*icols[iSplitRow]++ = 3*colItr->first + iSplitCol + 1;
				}
			}
		}

		irow = irows[2];
		icol = icols[2];
		a = as[2];
	}
}

void BlockMatrix::toCoordinateStorage( long* irow, long* icol, double* a ) const
{
	typedef long Integer;
	typedef double Real;

	for( size_t iRow = 0; iRow < rows_.size(); ++iRow )
	{
		size_t length = 3 * rows_[iRow].size();
		boost::array<Integer*, 3> irows = { { irow, irow + length, irow + length + length } };
		boost::array<Integer*, 3> icols = { { icol, icol + length, icol + length + length } };
		boost::array<Real*, 3> as = { { a, a + length, a + length + length } };
		for( std::vector<Entry>::const_iterator colItr = rows_[iRow].begin();
			colItr != rows_[iRow].end();
			++colItr )
		{
			for( unsigned int iSplitRow = 0; iSplitRow < 3; ++iSplitRow )
			{
				for( unsigned int iSplitCol = 0; iSplitCol < 3; ++iSplitCol )
				{
					*as[iSplitRow]++ = colItr->second[iSplitRow][iSplitCol];
					*irows[iSplitRow]++ = 3*iRow + iSplitRow + 1;
					*icols[iSplitRow]++ = 3*colItr->first + iSplitCol + 1;
				}
			}
		}

		irow = irows[2];
		icol = icols[2];
		a = as[2];
	}
}


ForceDebuggingInfo::ForceDebuggingInfo( unsigned int particleCount, const std::string& name )
	:	prevEval( new Evaluation(particleCount) ),
		nextEval( new Evaluation(particleCount) ),
		prevValid( false ),
		name_( name ),
		particleCount_( particleCount )
{
}

void ForceDebuggingInfo::evaluate( const VectorView<vl::Vec3>& x, 
		const VectorView<vl::Vec3>& v )
{
#ifdef DEBUG_FORCES
	for( unsigned int iParticle = 0; iParticle < this->nextEval->x.size(); ++iParticle )
		nextEval->x[iParticle] = x[iParticle];
	for( unsigned int iParticle = 0; iParticle < this->nextEval->v.size(); ++iParticle )
		nextEval->v[iParticle] = v[iParticle];

	if( prevValid )
	{
		// Figure out the deltas:
		std::vector<vl::Vec3> delta_x( particleCount() );
		std::transform( nextEval->x.begin(), nextEval->x.end(), prevEval->x.begin(), delta_x.begin(),
			std::minus<vl::Vec3>() );

		std::vector<vl::Vec3> delta_v( particleCount() );
		std::transform( nextEval->v.begin(), nextEval->v.end(), prevEval->v.begin(), delta_v.begin(),
			std::minus<vl::Vec3>() );

		std::vector<vl::Vec3> delta_f( particleCount() );
		std::transform( nextEval->f.begin(), nextEval->f.end(), prevEval->f.begin(), delta_f.begin(),
			std::minus<vl::Vec3>() );

		// We have our delta_f, delta_x, delta_v, now, we need to evaluate the error
		// [(Df/Dx * delta_x) + (Df/Dv * delta_v)] - delta_f
		std::vector<vl::Vec3> error( particleCount(), vl::Vec3(vl::vl_0) );
		prevEval->dfdx.matrixVectorProduct( delta_x, error );
		prevEval->dfdv.matrixVectorProduct( delta_v, error );

		std::transform( error.begin(), error.end(), delta_f.begin(), error.begin(), std::minus<vl::Vec3>() );

		debuggingInfo_.push_back( ForceDebuggingInfo::DebuggingInfo(delta_x, delta_v, error) );

		// Now, need to occasionally dump out these errors to disk:
		const unsigned int frequency = 20;
		if( (debuggingInfo_.size() % frequency) == 0 )
		{
			std::ostringstream filename;
			filename << name_ << ".mat";
			
			unsigned int nrows = 3*particleCount();
			unsigned int ncols = debuggingInfo_.size();
			fortran_matrix dx( nrows, ncols );
			fortran_matrix dv( nrows, ncols );
			fortran_matrix err( nrows, ncols );

			for( unsigned int i = 0; i < debuggingInfo_.size(); ++i )
			{
				const ForceDebuggingInfo::DebuggingInfo& cur = debuggingInfo_[i];

				const double* dx_ptr = reinterpret_cast<const double*>( &cur.delta_x[0] );
				const double* dv_ptr = reinterpret_cast<const double*>( &cur.delta_v[0] );
				const double* err_ptr = reinterpret_cast<const double*>( &cur.error[0] );
				std::copy( dx_ptr, dx_ptr + 3*particleCount(), dx.data() + i*nrows );
				std::copy( dv_ptr, dv_ptr + 3*particleCount(), dv.data() + i*nrows );
				std::copy( err_ptr, err_ptr + 3*particleCount(), err.data() + i*nrows );
			}

			MATFile matFile( filename.str(), "Force error dump" );
			matFile.add( "dx", dx );
			matFile.add( "dv", dv );
			matFile.add( "err", err );
		}
	}

	prevValid = true;
#endif
}

void ForceDebuggingInfo::swap()
{
#ifdef DEBUG_FORCES
	std::swap( prevEval, nextEval );

	nextEval->dfdx.clear();
	nextEval->dfdv.clear();
	std::fill( nextEval->f.begin(), nextEval->f.end(), vl::Vec3(vl::vl_0) );
#endif
}

boost::shared_ptr< AccumulationVector<vl::Vec3> > ForceDebuggingInfo::accum_f( AccumulationVector<vl::Vec3>& other )
{
#ifdef DEBUG_FORCES
	return boost::shared_ptr< AccumulationVector<vl::Vec3> >(
		new DuplicateAccumulationVector( other, *this->nextEval->f_accum ) );
#else
	return boost::shared_ptr< AccumulationVector<vl::Vec3> >(
		new AccumulationVectorWrapper<vl::Vec3>( other ) );
#endif
}

boost::shared_ptr< AccumulationMatrix<vl::Mat3> > ForceDebuggingInfo::accum_dfdx( AccumulationMatrix<vl::Mat3>& other )
{
#ifdef DEBUG_FORCES
	return boost::shared_ptr< AccumulationMatrix<vl::Mat3> >(
		new DuplicateAccumulationMatrix( other, this->nextEval->dfdx ) );
#else
	return boost::shared_ptr< AccumulationMatrix<vl::Mat3> >(
		new AccumulationMatrixWrapper<vl::Mat3>( other ) );
#endif
}

boost::shared_ptr< AccumulationMatrix<vl::Mat3> > ForceDebuggingInfo::accum_dfdv( AccumulationMatrix<vl::Mat3>& other )
{
#ifdef DEBUG_FORCES
	return boost::shared_ptr< AccumulationMatrix<vl::Mat3> >(
		new DuplicateAccumulationMatrix( other, this->nextEval->dfdv ) );
#else
	return boost::shared_ptr< AccumulationMatrix<vl::Mat3> >(
		new AccumulationMatrixWrapper<vl::Mat3>( other ) );
#endif
}

