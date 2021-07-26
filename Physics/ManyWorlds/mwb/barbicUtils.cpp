#include "stdafx.h"

#include "barbicUtils.h"

#include <cmath>

namespace barbic {

fortran_matrix ReadMatrixFromDisk(const std::string& filename)
{
	try
	{
		std::ifstream ifs( filename.c_str(), std::ios::binary );
		if( !ifs )
			throw IOException( "Unable to open file '" + filename + "' for reading." );
		ifs.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );

		int m, n;
		ifs.read( reinterpret_cast<char*>( &m ), sizeof(int) );
		ifs.read( reinterpret_cast<char*>( &n ), sizeof(int) );

		fortran_matrix result( m, n );
		ifs.read( reinterpret_cast<char*>( result.data() ), m*n*sizeof(double) );

		return result;
	}
	catch( std::ifstream::failure& e )
	{
		std::ostringstream oss;
		oss << "Error reading '" << filename << "': " << e.what();
		throw IOException( oss.str() );
	}
}

void WriteMatrixToDisk(const std::string& filename, const fortran_matrix& mat)
{
	try
	{
		std::ofstream ofs( filename.c_str(), std::ios::binary );
		if( !ofs )
			throw IOException( "Unable to open file '" + filename + "' for reading." );
		ofs.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );

		int m = mat.nrows();
		int n = mat.ncols();
		ofs.write( reinterpret_cast<const char*>( &m ), sizeof(int) );
		ofs.write( reinterpret_cast<const char*>( &n ), sizeof(int) );

		ofs.write( reinterpret_cast<const char*>( mat.data() ), m*n*sizeof(double) );
	}
	catch( std::ifstream::failure& e )
	{
		std::ostringstream oss;
		oss << "Error writing to '" << filename << "': " << e.what();
		throw IOException( oss.str() );
	}
}

template <typename RealType>
void SparseMatrix<RealType>::add( unsigned int i, unsigned int j, const RealType& value )
{
	(*this)(i, j) += value;
}

// Utility class for binary search
template <typename PairType>
struct EntryComparator
{
	bool operator()(const PairType& left, const PairType& right) const
	{
		return left.first < right.first;
	}
};

// this has to be distinct because in this case we're not allowed to
//   modify the matrix
template <typename RealType>
RealType SparseMatrix<RealType>::operator()( size_t i, size_t j ) const
{
	checkIndices(i, j);

	Entry entry(j, RealType(0));
	const EntryVec& row = rows_[i];

	typedef EntryVec::const_iterator Itr;
	typedef std::pair<Itr, Itr> ItrPr;
	ItrPr pr = std::equal_range(row.begin(), row.end(), entry, EntryComparator<Entry>());
	if (pr.first == pr.second)
		return RealType(0);
	else
		return pr.first->second;
}

template <typename RealType>
RealType& SparseMatrix<RealType>::operator()( size_t i, size_t j )
{
	checkIndices(i, j);

	Entry entry(j, RealType(0));
	EntryVec& row = rows_[i];

	typedef EntryVec::iterator Itr;
	typedef std::pair<Itr, Itr> ItrPr;

	ItrPr pr = std::equal_range(row.begin(), row.end(), entry, EntryComparator<Entry>());
	if (pr.first == pr.second)
	{
		Itr newIt = row.insert(pr.first, entry);
		return newIt->second;
	}
	else
		return pr.first->second;
}

template <typename RealType>
std::vector<RealType> SparseMatrix<RealType>::multiply( const RealType* rhs )
{
	std::vector<RealType> result( this->nrows(), RealType(0) );

	for( size_t iRow = 0; iRow < rows_.size(); ++iRow )
	{
		EntryVec& row = rows_.at( iRow );
		for( EntryVec::const_iterator rowItr = row.begin();
			rowItr != row.end(); ++rowItr )
		{
			result[ iRow ] += rowItr->second * rhs[ rowItr->first ];
		}
	}

	return result;
}

template SparseMatrix<float>;
template SparseMatrix<double>;


HarmonicOscillator_IIR::HarmonicOscillator_IIR(double mass, double xi, double stiffness, double dt) 
{
  externalForce = 0;
  q = 0;
  qvel = 0;
  F1 = 0;

  assert( xi > 0.0 && xi < 1.0 );

  const double omega = sqrt(stiffness/mass);
  const double ksi = xi;

  this->omegaDamped = omega * sqrt(1-ksi*ksi);

  // cache -ksi * omega
  this->minKsiOmega = -ksi * omega;

  std::complex<double> alpha(-ksi * omega,omegaDamped);

  std::complex<double> A = exp(dt*alpha);
  std::complex<double> B = 1.0 / (alpha * mass * omegaDamped) * (exp(dt*alpha) - std::complex<double>(1,0));
  
  a = A.real();
  b = A.imag();
  c = B.real();
  d = B.imag();

  assert( !isBad(a) );
  assert( !isBad(b) );
  assert( !isBad(c) );
  assert( !isBad(d) );
}

void HarmonicOscillator_IIR::DoTimestep()
{
	double F1Prev = F1;
	assert( _finite(externalForce) );

	F1 = F1Prev * a - q * b + externalForce * c;
	q = F1Prev * b + q * a + externalForce * d;
	qvel = minKsiOmega * q + omegaDamped * F1;

	assert( _finite(q) );
	assert( _finite(qvel) );
	assert( _finite(F1) );
}

void HarmonicOscillator_IIR::ChopZeroSignal(double threshold)
{
  if ((fabs(qvel) < threshold))
	this->qvel = 0.0;
  if((fabs(q) < threshold))
    this->q = 0.0;
  if((fabs(F1) < threshold))
    this->F1 = 0.0;
}

} // namespace barbic

