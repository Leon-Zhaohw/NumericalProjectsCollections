/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Conjugate Gradient solver classes, used in Pressure solve
 *
 *****************************************************************************/

#ifndef DDF_CONJUGATEGRADIENT_H
#define DDF_CONJUGATEGRADIENT_H

#include "vectorbase.h"
#include "grid.h"
#include "operators.h"


namespace DDF { 

static const bool CG_DEBUG = false;

// helper functions from conjugategrad.cpp
void InitPrecondition(int method, FlagGrid *Flags,
				Grid<Real> *A0, Grid<Real> *Ai, Grid<Real> *Aj, Grid<Real> *Ak,
				Grid<Real> *orgA0, Grid<Real> *orgAi, Grid<Real> *orgAj, Grid<Real> *orgAk);
void ApplyPrecondition(int method, Grid<Real> *Dst, Grid<Real> *Var1, FlagGrid *Flags,
				Grid<Real> *A0, Grid<Real> *Ai, Grid<Real> *Aj, Grid<Real> *Ak,
				Grid<Real> *orgA0, Grid<Real> *orgAi, Grid<Real> *orgAj, Grid<Real> *orgAk);


// hard coded laplacian operator without stored matix
class goApplyLacplace : public GridOpBase{
	public:
		goApplyLacplace(Grid<Real> *Dst, Grid<Real> *Src, FlagGrid *Flags) : GridOpBase() {
			mpDst = Dst; mpSrc = Src;
			mpFlags = Flags;
			applyOperatorToGrids(this);
		};
		~goApplyLacplace() { };
		void buildCallList() {
			gaDst.gridAccInit(mpDst, AM_WRITE, gaCalls); 
			gaSrc.gridAccInit(mpSrc, AM_READ, gaCalls); 
			setFlags(mpFlags);
		}; 
		inline void operator() (int i, int j, int k) {
			if(!fgIsFluid(getFlagAcc()(i,j,k)) ) {
				gaDst.write(i,j,k) = gaSrc(i,j,k);
				return;
			}

			gaDst.write(i,j,k) = LAPLACEFAC * gaSrc(i,j,k) 
					- gaSrc(i-1,j,k) - gaSrc(i+1,j,k)
					- gaSrc(i,j-1,k) - gaSrc(i,j+1,k)
#if DDF_DIMENSION==3
					- gaSrc(i,j,k-1) - gaSrc(i,j,k+1)
#endif // DDF_DIMENSION==3
				; 
		} 
		void getRidOfWarning() { };
		void reduce(goApplyLacplace &op) { getRidOfWarning(); } 
	protected:
		Grid<Real> *mpDst, *mpSrc;
		GridAccessor<Real,0> gaDst;
		GridAccessor<Real,1> gaSrc;
}; // goApplyLacplace



//*****************************************************************************

// compute residual (init) and add to sigma
class goInitSigma : public GridOpBase{
	public:
		goInitSigma(Grid<Real> *Dst, Grid<Real> *Var1, Grid<Real> *Var2, FlagGrid *Flags, Real &set_sigma) : 
				GridOpBase() {
			mpDst = Dst; // residual
			mpVar1 = Var1; // rhs
			mpVar2 = Var2; // temp
			mpFlags = Flags;
			applyOperatorToGrids( this );
			set_sigma = mSigma;
		}
		~goInitSigma() {};
		virtual void resetVariables() {
			mSigma = 0.;
		}
		void buildCallList() {
			gaResid.gridAccInit( mpDst, AM_WRITE, gaCalls); 
			gaRhs.gridAccInit(   mpVar1, AM_READ, gaCalls); 
			gaTemp.gridAccInit(  mpVar2, AM_READ, gaCalls); 
			setFlags(mpFlags);
		};

		inline void operator() (int i, int j, int k) {
			const Real res = gaRhs(i,j,k) - gaTemp(i,j,k); 
			gaResid.write(i,j,k) = res;

			// only compute residual in fluid region
			if(!fgIsFluid(getFlagAcc()(i,j,k)) ) return; 
			mSigma += res*res;
		} 
		void reduce(goInitSigma &op) { 
			//debMsg("goInitSigma","reduce mSigma"<<mSigma<<" += "<<op.mSigma);
			mSigma += op.mSigma; 
		} 
	protected:
		Real mSigma;
		Grid<Real> *mpDst, *mpVar1, *mpVar2;
		GridAccessor<Real,0> gaResid, gaRhs, gaTemp;
}; // goInitSigma */


class goUpdateSolution : public GridOpBase{
	public:
		goUpdateSolution(Grid<Real> *Dst, Grid<Real> *Var1, FlagGrid *Flags, Real factor) : GridOpBase() {
			mpDst = Dst; // pressure or residual
			mpSrc = Var1; // update vector
			mpFlags = Flags;
			mFactor = factor;
			applyOperatorToGrids( this );
		}
		~goUpdateSolution() {};
		void buildCallList() {
			gaDst.gridAccInit( mpDst, AM_WRITE, gaCalls); 
			gaSrc.gridAccInit( mpSrc, AM_READ, gaCalls); 
			setFlags(mpFlags);
		}; 
		inline void operator() (int i, int j, int k) {
			gaDst.write(i,j,k) += mFactor*gaSrc(i,j,k); 
		} 
		void reduce(goUpdateSolution &op) { };
	protected:
		Real mFactor;
		Grid<Real> *mpDst;
		Grid<Real> *mpSrc;
		GridAccessor<Real,0> gaDst;
		GridAccessor<Real,0> gaSrc;
}; // goUpdateSolution */

class goUpdateSearchVec : public GridOpBase{
	public:
		goUpdateSearchVec(Grid<Real> *Dst, Grid<Real> *Var1, FlagGrid *Flags, Real factor) : GridOpBase() {
			mpDst = Dst; // pressure or residual
			mpSrc = Var1; // update vector
			mpFlags = Flags;
			mFactor = factor;
			applyOperatorToGrids( this );
		}
		~goUpdateSearchVec() {};
		void buildCallList() {
			gaDst.gridAccInit( mpDst, AM_WRITE, gaCalls); 
			gaSrc.gridAccInit( mpSrc, AM_READ, gaCalls); 
			setFlags(mpFlags);
		}; 
		inline void operator() (int i, int j, int k) {
			gaDst.write(i,j,k) = gaSrc(i,j,k) + mFactor*gaDst(i,j,k); 
		} 
		void reduce(goUpdateSearchVec &op) { };
	protected:
		Real mFactor;
		Grid<Real> *mpDst, *mpSrc;
		GridAccessor<Real,0> gaDst, gaSrc;
}; // goUpdateSearchVec */



// currently only copies from source to dest...
class goPreconditionID : public GridOpBase{
	public:
		goPreconditionID(Grid<Real> *Dst, Grid<Real> *Var1, FlagGrid *Flags) : GridOpBase() {
			mpDst = Dst;
			mpVar1 = Var1;
			mpFlags = Flags;
			applyOperatorToGrids( this );
		}
		~goPreconditionID() {};
		void buildCallList() {
			gaDst.gridAccInit(mpDst, AM_WRITE, gaCalls); 
			gaVar1.gridAccInit(mpVar1, AM_READ, gaCalls); 
			setFlags(mpFlags);
		};

		inline void operator() (int i, int j, int k) {
			gaDst.write(i,j,k) = gaVar1(i,j,k);
		} 
		void reduce(goPreconditionID &op) { } 
	protected:
		Grid<Real> *mpDst, *mpVar1;
		GridAccessor<Real,0> gaDst, gaVar1;
}; // goPreconditionID */



//*****************************************************************************
// compute the laplacian of src with right hand side rhs into the tmp grid
// for reference to CG solver
class goJacobiSolver : public GridOpBase{
	public:
		goJacobiSolver(Grid<Real> *Dst, Grid<Real> *Src, Grid<Real> *Rhs, FlagGrid *Flags) : GridOpBase() {
			mpDst = Dst;
			mpSrc = Src;
			mpRhs = Rhs;
			mpFlags = Flags;
			applyOperatorToGrids(this);
		};
		~goJacobiSolver() {};
		void buildCallList() {
			gaDst.gridAccInit(mpDst, AM_WRITE, gaCalls); 
			gaRhs.gridAccInit(mpRhs, AM_READ, gaCalls); 
			gaSrc.gridAccInit(mpSrc, AM_READ, gaCalls); 
			setFlags(mpFlags);
		};

		inline void operator() (int i, int j, int k) {
			//debMsg("OPL","at "<<PRINT_IJK<<" dst="<<gaDst.write(i,j,k)<<" rhs="<<gaRhs(i,j,k)<<", src: "<< gaSrc(i-1,j,k) <<" "<< gaSrc(i+1,j,k)<<" "<< gaSrc(i,j-1,k) <<" "<< gaSrc(i,j+1,k) ); 
			if(fgIsObstacle(getFlagAcc()(i,j,k)) ) return;
			gaDst.write(i,j,k) = ( gaRhs(i,j,k) + 
					+ gaSrc(i-1,j,k) + gaSrc(i+1,j,k)
					+ gaSrc(i,j-1,k) + gaSrc(i,j+1,k) 
#if DDF_DIMENSION==3
					+ gaSrc(i,j,k-1) + gaSrc(i,j,k+1) 
#endif // DDF_DIMENSION==3
				) * (1./LAPLACEFAC); 
		}

		void reduce(goJacobiSolver &op) { }
		
	protected:
		Grid<Real> *mpDst, *mpRhs, *mpSrc;
		GridAccessor<Real,0> gaDst;
		GridAccessor<Real,0> gaRhs;
		GridAccessor<Real,1> gaSrc;
}; // goJacobiSolver



//*****************************************************************************
// compute max values in grid  (no position)
template<class Value>
class goCompMaxSimple : public GridOpBase{
	public:
		goCompMaxSimple(Grid<Value> *Var, FlagGrid *Flags) : GridOpBase() {
			mpVar = Var;
			mpFlags = Flags;
			applyOperatorToGrids( this );
		};
		~goCompMaxSimple() {};
		virtual void resetVariables() {
			mMaxAbs = -1e10;
			mMax = 0.;
		}
		void buildCallList() {
			gaVar.gridAccInit(mpVar, AM_READ, gaCalls); 
			setFlags(mpFlags);
		};

		inline void operator() (int i, int j, int k) {
			const Value s = gaVar(i,j,k); 
			const Real sabs = normHelper(s);

			if(sabs>mMaxAbs) {
				mMaxAbs = sabs;
				mMax = s;
				//mMaxPos = nVec3i(i,j,k);
			}
		}

		void reduce(goCompMaxSimple &op) {
			if(op.mMaxAbs>mMaxAbs) {
				mMaxAbs = op.mMaxAbs;
				mMax = op.mMax;
				//mMaxPos = op.mMaxPos;
			}
		}

		Value getMax() { return mMax; }
		Real getMaxAbs() { return mMaxAbs; }
		//Real getMaxPosValue() { return mMaxPos; }
		string toString() {
			std::ostringstream out;
			out << " max="<<mMax<<" ("<<mMaxAbs<<") ";
			return out.str();
		}
	protected:
		//int outcnt;
		Real mMaxAbs;
		Value mMax;
		Grid<Value> *mpVar;
		GridAccessor<Real,0> gaVar;
}; // goCompMaxSimple

//*****************************************************************************

// basic interface 
class GridCgInterface {
	public:
		GridCgInterface() : mUseResNorm(true) {};
		virtual ~GridCgInterface() {};

		// solving functions
		virtual bool iterate() = 0;
		virtual void solve(int maxIter) = 0;

		// precond
		virtual void setPreconditioner(int method, Grid<Real> *A0, Grid<Real> *Ai, Grid<Real> *Aj, Grid<Real> *Ak) = 0;

		// access
		virtual Real getSigma() const = 0;
		virtual Real getIterations() const = 0;
		virtual Real getResNorm() const = 0;
		virtual void setAccuracy(Real set) = 0;
		virtual Real getAccuracy() const = 0;

		void setUseResNorm(bool set) { mUseResNorm = set; }

	protected:

		// use norm of residual, or max value for threshold?
		bool mUseResNorm; 
};

//! run single iteration of the cg solver
// the template argument determines the type of matrix multiplication,
// typically a goApplyMatrix object (above), another one is needed e.g. for the
// mesh-based wave equation solver
template<class APPLYMAT>
class GridCg : public GridCgInterface {
	public:
		// constructor
		GridCg(Grid<Real>* dst,Grid<Real>* rhs,Grid<Real>* residual,Grid<Real>* search, FlagGrid* flags, Grid<Real>* tmp, 
				Grid<Real> *A0=NULL, Grid<Real> *Ai=NULL, Grid<Real> *Aj=NULL, Grid<Real> *Ak=NULL,
				// optional grids that may be used in the apply-matrix class during matrix myultiplication
				// used eg for surface tension
			  	Grid<Real> *_optionalReal = NULL, Grid<Vec3> *_optionalVec  = NULL  ) :
			GridCgInterface(),
			mInited(false), mIterations(0),
			mpDst(dst), mpRhs(rhs), mpResidual(residual),
			mpSearch(search), mpFlags(flags), mpTmp(tmp),
			mpA0(A0), mpAi(Ai), mpAj(Aj), mpAk(Ak),
			mPcMethod(0), mpPCA0(A0), mpPCAi(Ai), mpPCAj(Aj), mpPCAk(Ak),
			mSigma(0.),
			mAccuracy(gVecEpsilon), mResNorm(1e20)
		{
			if(0) debMsg("PB-CG-Norms","init p"<<sqrt( DDF::GridOpNormNosqrt(mpDst, mpFlags).getValue() ) <<
					" rhs"<<sqrt( DDF::GridOpNormNosqrt(mpRhs, mpFlags).getValue() ) << " search"<<sqrt( DDF::GridOpNormNosqrt(mpSearch, mpFlags).getValue() ) <<
					" res"<<sqrt( DDF::GridOpNormNosqrt(mpResidual, mpFlags).getValue() ) << " tmp"<<sqrt( DDF::GridOpNormNosqrt(mpTmp, mpFlags).getValue() )<<
					" sigma="<<mSigma ); // debug

			GridOpTouchMemory<Real>(NULL, dst, 0.);
			GridOpTouchMemory<Real>(NULL, residual, 0.);
			GridOpTouchMemory<Real>(NULL, search, 0.);
			GridOpTouchMemory<Real>(NULL, tmp, 0.);
			optionalReal = _optionalReal;
			optionalVec = _optionalVec;
		}

		void doInit() {
			mInited = true;

			InitPrecondition(mPcMethod, mpFlags, mpPCA0,mpPCAi,mpPCAj,mpPCAk , mpA0,mpAi,mpAj,mpAk);

			// old APPLYMAT(mpTmp, mpDst, mpA0,mpAi,mpAj,mpAk, mpFlags, optionalReal, optionalVec); 
			// old goInitSigma(mpResidual, mpRhs, mpTmp, mpFlags, mSigma);
			// old ApplyPrecondition(mPcMethod, mpSearch, mpResidual, mpFlags, mpPCA0,mpPCAi,mpPCAj,mpPCAk ); 
			// old GridOpNormNosqrt(mpResidual, mpFlags, &mResNorm);

			goCopyGrid<Real>(mpResidual, mpRhs); // p=0, residual = b

			// tmp == apply preconditioner(residual)
			ApplyPrecondition(mPcMethod, mpTmp, mpResidual, mpFlags, mpPCA0,mpPCAi,mpPCAj,mpPCAk, mpA0,mpAi,mpAj,mpAk );
			
			goCopyGrid<Real>(mpSearch, mpTmp); // search = tmp

			Real dotProdValue = 0.; // sigma = dot(tmp, residual)
			goDotProd(mpTmp, mpResidual, mpFlags, &dotProdValue);
			mSigma = dotProdValue;
			
			if(CG_DEBUG) debMsg("PB-CG-Norms","init p"<<sqrt( DDF::GridOpNormNosqrt(mpDst, mpFlags).getValue() ) <<
					" rhs"<<sqrt( DDF::GridOpNormNosqrt(mpRhs, mpFlags).getValue() ) << " search"<<sqrt( DDF::GridOpNormNosqrt(mpSearch, mpFlags).getValue() ) <<
					" res"<<sqrt( DDF::GridOpNormNosqrt(mpResidual, mpFlags).getValue() ) << " tmp"<<sqrt( DDF::GridOpNormNosqrt(mpTmp, mpFlags).getValue() )<<
					" sigma="<<mSigma ); // debug
		}
		~GridCg() {
		};

		//! perform cg iteration
		bool iterate() {
			if(!mInited) doInit();

			//debMsg("CG","Iterate "<<mIterations<<" "<<mResNorm<<"/"<<mAccuracy);
			//if(mResNorm<mAccuracy) return false;

			mIterations++;

			// create matrix application operator passed as template argument,
			// this could reinterpret the mpA pointers (not so clean right now)
			// tmp = applyMat(search)
			APPLYMAT(mpTmp, mpSearch, mpA0,mpAi,mpAj,mpAk, mpFlags, optionalReal, optionalVec);

			Real dotProdTS;  // alpha = sigma/dot(tmp, search)
			goDotProd(mpTmp, mpSearch, mpFlags, &dotProdTS);
			Real alpha = 0.;
			if(fabs(dotProdTS)>0.) alpha = mSigma / dotProdTS;

			goUpdateSolution(mpDst,      mpSearch, mpFlags,  alpha);  // dst += search * alpha
			goUpdateSolution(mpResidual, mpTmp   , mpFlags, -alpha);  // residual += tmp * -alpha
			ApplyPrecondition(mPcMethod, mpTmp, mpResidual, mpFlags, mpPCA0,mpPCAi,mpPCAj,mpPCAk, mpA0,mpAi,mpAj,mpAk );
			
			// compute norm of the residual?
			if(this->mUseResNorm) {
				GridOpNormNosqrt(mpResidual, mpFlags, &mResNorm);
			} else {
				// TODO test...
				goCompMaxSimple<Real> resMax = goCompMaxSimple<Real>(mpResidual, mpFlags);
				mResNorm = resMax.getMaxAbs();
			}
			if(0 || mIterations % 10 == 9) debMsg("GridCg","Iteration i="<<mIterations<<", resNorm="<<mResNorm<<" accuracy="<<mAccuracy );

			// abort here to safe some work...
			if(mResNorm<mAccuracy) return false;

			Real sigmaNew; // sigma_new = dot(tmp, residual)
			goDotProd(mpTmp, mpResidual, mpFlags, &sigmaNew);
			// abort here?

			Real beta = sigmaNew / mSigma;
			// search =  tmp + beta * srch
			goUpdateSearchVec(mpSearch, mpTmp, mpFlags, beta);

			if(CG_DEBUG) debMsg("PB-Cg","iter2 i="<<mIterations<<" sigma="<<mSigma<<" alpha="<<alpha<<" beta="<<beta<<" ");
			mSigma = sigmaNew;

			if(CG_DEBUG) debMsg("PB-CG-Norms","p"<<sqrt( DDF::GridOpNormNosqrt(mpDst, mpFlags).getValue() ) <<" search"<<sqrt( DDF::GridOpNormNosqrt(mpSearch, mpFlags).getValue() ) 
					<<" res"<<sqrt( DDF::GridOpNormNosqrt(mpResidual, mpFlags).getValue() ) <<" tmp"<<sqrt( DDF::GridOpNormNosqrt(mpTmp, mpFlags).getValue() ) ); // debug
			if(CG_DEBUG) debMsg("PB-CG-Norms","p"<<( DDF::GridOpNormNosqrt(mpDst, mpFlags).getValue() ) <<" search"<<( DDF::GridOpNormNosqrt(mpSearch, mpFlags).getValue() ) 
					<<" res"<<( DDF::GridOpNormNosqrt(mpResidual, mpFlags).getValue() ) <<" tmp"<<( DDF::GridOpNormNosqrt(mpTmp, mpFlags).getValue() ) ); // debug, no sqrt!
			return true;
		}

		void solve(int maxIter) {
			for (int iter=0; iter<maxIter; iter++) {
				if (!iterate()) iter=maxIter;
			} 
			return;
		}

		// init pointers, and copy values from "normal" matrix
		void setPreconditioner(int method, Grid<Real> *A0, Grid<Real> *Ai, Grid<Real> *Aj, Grid<Real> *Ak) {
			mPcMethod = method;
			mpPCA0 = A0;
			mpPCAi = Ai;
			mpPCAj = Aj;
			mpPCAk = Ak;
		}

		// access

		Real getSigma() const { return mSigma; }
		Real getIterations() const { return mIterations; }

		Real getResNorm() const { return mResNorm; }

		void setAccuracy(Real set) { mAccuracy=set; }
		Real getAccuracy() const { return mAccuracy; }

	protected:
		bool mInited;
		int mIterations;
		// grids
		Grid<Real>* mpDst;
		Grid<Real>* mpRhs;
		Grid<Real>* mpResidual;
		Grid<Real>* mpSearch;
		FlagGrid*  mpFlags;
		Grid<Real>* mpTmp;

		Grid<Real> *mpA0, *mpAi, *mpAj, *mpAk;

		// preconditioning method 0=off, 1=incomp cholesky
		int mPcMethod;
		// preconditioning grids
		Grid<Real> *mpPCA0, *mpPCAi, *mpPCAj, *mpPCAk;

		// optional grids that may be used in the apply-matrix class during matrix myultiplication
		// used eg for surface tension
		Grid<Real> *optionalReal;
		Grid<Vec3> *optionalVec;

		// sigma / residual
		Real mSigma;
		// accuracy of solver (max. residuum)
		Real mAccuracy;
		// norm of the residual
		Real mResNorm;
}; // GridCg



} // namespace DDF

#endif // DDF_CONJUGATEGRADIENT_H

