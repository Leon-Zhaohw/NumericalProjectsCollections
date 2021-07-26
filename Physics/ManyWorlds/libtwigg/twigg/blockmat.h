#ifndef __BLOCKMAT_H__

#define __BLOCKMAT_H__

#include "twigg/util.h"
#include "twigg/vlutil.h"
#include "twigg/linalg.h"

#include <cassert>
#include <vector>
#include <iosfwd>
#include <numeric>
#include <boost/array.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

// This is an "accumulation matrix" because it only has functionality to
// accumulate things like force derivatives, etc.
template <typename T>
class AccumulationMatrix
{
public:
	// Since this is a virtual class, it's not going to work so well
	// to return wrappers.  We'll use the following methods instead.
	virtual void add( unsigned int i, unsigned int j, const T& value ) = 0;
	virtual void subtract( unsigned int i, unsigned int j, const T& value )
	{
		add( i, j, -value );
	}

	virtual ~AccumulationMatrix() {}
};

// Separate out the read-only view from the writeable one:
template <typename T>
class MatrixView
{
public:
	virtual T operator()( unsigned int i, unsigned int j ) const = 0;
	virtual ~MatrixView() {}
};

// Write-only vector view
template <typename T>
class AccumulationVector
{
public:
	virtual void add( unsigned int i, const T& value ) = 0;
	virtual void subtract( unsigned int i, const T& value )
	{
		add( i, -value );
	}

	virtual ~AccumulationVector() {}
};

// Read-only vector view
template <typename T>
class VectorView
{
public:
	virtual T operator[]( unsigned int i ) const = 0;
	virtual ~VectorView() {}
};

class SimpleVectorView
	: public VectorView< vl::Vec3 >
{
public:
	SimpleVectorView( const vl::Vec3* x )
		: x_(x) {}

	SimpleVectorView( const vl::Real* x )
		: x_( reinterpret_cast<const vl::Vec3*>(x) ) {}

	vl::Vec3 operator[]( unsigned int i ) const { assert(!isBad(x_[i]) ); return x_[i]; }

private:
	const vl::Vec3* x_;
};

class SimpleVectorAccumulator
	: public AccumulationVector<vl::Vec3>
{
public:
	SimpleVectorAccumulator( vl::Vec3* x, unsigned int particleCount )
		: x_(x), particleCount_(particleCount) {}

	SimpleVectorAccumulator( vl::Real* x, unsigned int particleCount )
		: x_( reinterpret_cast<vl::Vec3*>(x) ), particleCount_(particleCount) {}

	SimpleVectorAccumulator( std::vector<vl::Vec3>& x )
		: x_( &x[0] ), particleCount_(x.size()) {}

	void add( unsigned int i, const vl::Vec3& value )
	{
		assert( i < particleCount_ );
		x_[i] += value;
	}

	void subtract( unsigned int i, const vl::Vec3& value )
	{
		assert( i < particleCount_ );
		x_[i] -= value;
	}

private:
	vl::Vec3* x_;
	const size_t particleCount_;
};

class BlockMatrix
	:	public AccumulationMatrix<vl::Mat3>
{
public:
	explicit BlockMatrix(int size);

	const BlockMatrix& operator*=(const double rhs);
	const BlockMatrix& operator+=(const BlockMatrix& rhs);
	const BlockMatrix& operator-=(const BlockMatrix& rhs);

	void printDense(std::ostream& os) const;
	void printSparse(std::ostream& os) const;

	fortran_matrix toMatrix() const;

	void toCoordinateStorage( long* irow, long* icol, double* a ) const;
	void toSymmetricCoordinateStorage( long* irow, long* icol, double* a, unsigned int* rowPtrs ) const;

	inline size_t nrows() const { return 3*rows_.size(); }
	inline size_t size() const
		{ return rows_.size(); }

	const vl::Mat3 operator() (unsigned int i, unsigned int j) const;
	vl::Mat3& operator() (unsigned int i, unsigned int j);

	void add( unsigned int i, unsigned int j, const vl::Mat3& value );
	void subtract( unsigned int i, unsigned int j, const vl::Mat3& value );

	void clear();
	void clearEntries();

	typedef std::pair<unsigned int, vl::Mat3> Entry;

	template <typename BlockVecL, typename BlockVecR>
	void matrixVectorProduct( const BlockVecL& rhs, BlockVecR& output )
	{
		for( unsigned int iRow = 0; iRow < rows_.size(); ++iRow )
		{
			for( std::vector<Entry>::const_iterator colItr = rows_[iRow].begin();
				colItr != rows_[iRow].end();
				++colItr )
			{
				output[ iRow ] += colItr->second * rhs[ colItr->first ];
			}
		}
	}

	void matrixVectorProduct( const vl::Vec3* rhs, vl::Vec3* output )
	{
		for( unsigned int iRow = 0; iRow < rows_.size(); ++iRow )
		{
			for( std::vector<Entry>::const_iterator colItr = rows_[iRow].begin();
				colItr != rows_[iRow].end();
				++colItr )
			{
				output[ iRow ] += colItr->second * rhs[ colItr->first ];
			}
		}
	}

	template <typename Matrix>
	void similarityTransform( const Matrix& U, Matrix& out ) const
	{
		// iterate through all the nonzero entries in the matrix
		for( unsigned int iRow = 0; iRow < rows_.size(); ++iRow )
		{
			for( std::vector<Entry>::const_iterator colItr = rows_[iRow].begin();
				colItr != rows_[iRow].end();
				++colItr )
			{
				for( unsigned int k = 0; k < out.nrows(); ++k )
				{
					for( unsigned int l = 0; l < out.ncols(); ++l )
					{
						unsigned int i = 3*iRow;
						unsigned int j = 3*colItr->first;

						out(k,l) += dot(
							vl::Vec3( U(i+0, k), U(i+1, k), U(i+2, k) ),
							colItr->second * vl::Vec3( U(j+0, l), U(j+1, l), U(j+2, l) ) );
					}
				}
			}
		}
	}

	// precondition: The index set _must_ be sorted, for the binary search
	template <typename Matrix1, typename Matrix2, typename Indexer>
	void similarityTransformSparse( const Matrix1& U, Matrix2& out, const Indexer& indexer ) const
	{
		typedef typename Indexer::const_iterator IndexerItr;
		for( IndexerItr iIndex = indexer.begin();
			iIndex != indexer.end(); 
			++iIndex )
		{
			const unsigned int iRow = (*iIndex) / 3;
			const unsigned int i = 3*std::distance( indexer.begin(), iIndex );

			// The assumption we will make for the moment is that the number of entries 
			// in a given row looks something like O(1) (should be the number of adjacent
			// particles).  Hence, it will be cheaper to iterate through all the entries
			// in the row and try to find them in the index set rather than the other way
			// around.  

			const EntryVec& row = rows_[iRow];
			IndexerItr indexItr = indexer.begin();

			for( EntryVec::const_iterator colItr = row.begin();
				colItr != row.end();
				++colItr )
			{
				typedef std::pair< typename Indexer::const_iterator, typename Indexer::const_iterator > ItrPr;
				ItrPr pr = std::equal_range( indexItr, indexer.end(), 3*colItr->first );

				// If this particular entry is actually in the entry set:
				if( pr.first != pr.second )
				{
					const unsigned int j = 3*std::distance( indexer.begin(), pr.first );
					for( unsigned int k = 0; k < out.nrows(); ++k )
					{
						for( unsigned int l = 0; l < out.ncols(); ++l )
						{
							out(k,l) += dot(
								vl::Vec3( U(i+0, k), U(i+1, k), U(i+2, k) ),
								colItr->second * vl::Vec3( U(j+0, l), U(j+1, l), U(j+2, l) ) );
						}
					}
				}

				// Push the current position in the index set.
				indexItr = pr.second;
			}
		}
	}

	unsigned int nnz() const;

private:

	inline void checkIndices(unsigned int i, unsigned int j) const { 
		assert(0 <= i && i < size() && 0 <= j && j < size()); 
	}

	typedef std::vector<Entry> EntryVec;
	std::vector< EntryVec > rows_;
	
	void printRow(std::ostream& os, int rowBlock, const unsigned int rowNum, const unsigned int numberLength) const;
	size_t numberLength() const;
	size_t numberLength(double number) const;
};

class DuplicateAccumulationMatrix
	: public AccumulationMatrix<vl::Mat3>
{
public:
	DuplicateAccumulationMatrix( AccumulationMatrix<vl::Mat3>& left, AccumulationMatrix<vl::Mat3>& right )
		: left_(left), right_(right) {}
	virtual ~DuplicateAccumulationMatrix() {}

	void add( unsigned int i, unsigned int j, const vl::Mat3& value )
	{
		left_.add( i, j, value );
		right_.add( i, j, value );
	}

private:
	AccumulationMatrix<vl::Mat3>& left_;
	AccumulationMatrix<vl::Mat3>& right_;
};

class DuplicateAccumulationVector
	: public AccumulationVector<vl::Vec3>
{
public:
	DuplicateAccumulationVector( AccumulationVector<vl::Vec3>& left, AccumulationVector<vl::Vec3>& right )
		: left_(left), right_(right) {}
	virtual ~DuplicateAccumulationVector() {}

	void add( unsigned int i, const vl::Vec3& value )
	{
		left_.add(i, value);
		right_.add(i, value);
	}

private:
	AccumulationVector<vl::Vec3>& left_;
	AccumulationVector<vl::Vec3>& right_;
};

template <typename T>
class AccumulationVectorWrapper
	: public AccumulationVector<T>
{
public:
	AccumulationVectorWrapper( AccumulationVector<T>& child )
		: child_(child)
	{
	}

	virtual ~AccumulationVectorWrapper() {}

	void add( unsigned int i, const T& value )
	{
		child_.add( i, value );
	}

private:
	AccumulationVector<T>& child_;
};

template <typename T>
class AccumulationMatrixWrapper
	: public AccumulationMatrix<T>
{
public:
	AccumulationMatrixWrapper( AccumulationMatrix<T>& child )
		: child_(child)
	{
	}

	virtual ~AccumulationMatrixWrapper() {}

	void add( unsigned int i, unsigned int j, const T& value )
	{
		child_.add( i, j, value );
	}

private:
	AccumulationMatrix<T>& child_;
};

class DummyAccumulationVector
	: public AccumulationVector<vl::Vec3>
{
public:
	DummyAccumulationVector() {}
	virtual ~DummyAccumulationVector() {}

	void add( unsigned int i, const vl::Vec3& value )
	{
	}
};

class DummyAccumulationMatrix
	: public AccumulationMatrix<vl::Mat3>
{
public:
	DummyAccumulationMatrix() {}
	virtual ~DummyAccumulationMatrix() {}

	void add( unsigned int i, unsigned int j, const vl::Mat3& value )
	{
	}
};

class ForceDebuggingInfo
	: boost::noncopyable
{
public:
	ForceDebuggingInfo( unsigned int particleCount, const std::string& name );

	// Evaluates \del f, \del x, and \del v based on previous stored data
	// and clears out all the df_d(x,v) matrices
	void evaluate( const VectorView<vl::Vec3>& x, 
		const VectorView<vl::Vec3>& v );
	void swap();

	boost::shared_ptr< AccumulationVector<vl::Vec3> > accum_f( AccumulationVector<vl::Vec3>& other );
	boost::shared_ptr< AccumulationMatrix<vl::Mat3> > accum_dfdx( AccumulationMatrix<vl::Mat3>& other );
	boost::shared_ptr< AccumulationMatrix<vl::Mat3> > accum_dfdv( AccumulationMatrix<vl::Mat3>& other );

	unsigned int particleCount() const { return particleCount_; }

private:
	struct Evaluation
	{
		Evaluation( unsigned int particleCount )
			:	x( particleCount, vl::Vec3(vl::vl_0) ),
				v( particleCount, vl::Vec3(vl::vl_0) ),
				f( particleCount, vl::Vec3(vl::vl_0) ),
				dfdx( particleCount ),
				dfdv( particleCount ),
				f_accum( new SimpleVectorAccumulator(f) ) {}

		std::vector<vl::Vec3> x;
		std::vector<vl::Vec3> v;

		std::vector<vl::Vec3> f;
		BlockMatrix dfdx;
		BlockMatrix dfdv;

		boost::scoped_ptr<SimpleVectorAccumulator> f_accum;
	};

	boost::shared_ptr<Evaluation> prevEval;
	boost::shared_ptr<Evaluation> nextEval;
	bool prevValid;

	struct DebuggingInfo
	{
		DebuggingInfo( const std::vector<vl::Vec3>& dx,
			const std::vector<vl::Vec3>& dv,
			const std::vector<vl::Vec3>& err )
			: delta_x( dx ), delta_v( dv ), error( err ) {}

		std::vector<vl::Vec3> delta_x;
		std::vector<vl::Vec3> delta_v;
		std::vector<vl::Vec3> error;
	};
	std::deque<DebuggingInfo> debuggingInfo_;

	std::string name_;
	unsigned int particleCount_;
};

#endif

