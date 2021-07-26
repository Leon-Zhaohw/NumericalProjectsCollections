#ifndef __BARBICUTILS_H__
#define __BARBICUTILS_H__

#include "twigg/linalg.h"

namespace barbic {

// utilities from Jernej Barbic's code for use with the 
//   sound generation stuff: matrix classes, etc.
fortran_matrix ReadMatrixFromDisk(const std::string& filename);
void WriteMatrixToDisk(const std::string& filename, const fortran_matrix& mat);

// Jernej Barbic, CMU, 2005
//
// Integrates a single 1D harmonic oscillator using IIR filters.
// The system must be underdamped.
//
// M * q'' + C * q' + K * q = f,
// 
//   where M,C,K,f are scalars
// 
// The notation (such as for ksi) is taken from:
//   James, D. L., and Pai, D. K., 
//   DyRT: Dynamic Response Textures for Real
//   Time Deformation Simulation with Graphics Hardware. 
//   In Proc. of ACM SIGGRAPH 2002.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

class HarmonicOscillator_IIR
{
public:
  // changed to use dimensionless damping parameter xi
  HarmonicOscillator_IIR(double mass, double xi, double stiffness, double dt);

  inline double Getq() { return q; }
  inline double Getqvel() { return qvel; }
  inline double Setq();
  inline double Setqvel();

  inline void ResetToRestPosition() { q = qvel = F1 = 0.0; }

  // external forces remain in force until explicity changed
  void SetExternalForce(double f) { externalForce = f; }
  void AddExternalForce(double f) { externalForce += f; }

  // performs one step of simulation
  void DoTimestep();

  // reset signal to zero if q below the given threshold
  void ChopZeroSignal(double threshold);

protected:
  double q; // current deformation amplitudes
  double qvel; // current velocities of deformation amplitudes
  double externalForce; 
  double F1;
  double a,b,c,d;
  double omegaDamped;
  double minKsiOmega;
};

// this is actually my class -cdtwigg
template <typename RealType>
class SparseMatrix
{
	typedef std::pair<size_t, RealType> Entry;
	typedef std::vector<Entry> EntryVec;

public:
	explicit SparseMatrix(size_t nrows, size_t ncols)
		: rows_(nrows), nCols_(ncols) {}

	inline size_t nrows() const { return rows_.size(); }
	inline size_t ncols() const { return this->nCols_; }

	void add( size_t i, size_t j, const RealType& value );

	RealType operator()( size_t i, size_t j ) const;
	RealType& operator()( size_t i, size_t j );

	template <typename MatrixType>
	void asDense( MatrixType& mat )
	{
		for( size_t iRow = 0; iRow < rows_.size(); ++iRow )
		{
			EntryVec& row = rows_.at( iRow );
			for( EntryVec::const_iterator rowItr = row.begin();
				rowItr != row.end(); ++rowItr )
			{
				mat( iRow, rowItr->first ) = rowItr->second;
			}
		}
	}

	std::vector<RealType> multiply( const RealType* rhs );

private:
	inline void checkIndices(size_t i, size_t j) const
	{ 
		assert(i < nrows() && j < ncols()); 
	}

	std::vector< EntryVec > rows_;

	size_t nCols_;
};

} // namespace barbic

#endif
