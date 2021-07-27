/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Plugins for calculating pressure correction and diffusion
 *
 *****************************************************************************/

#include "fluidsolver.h"
#include "solverplugin.h"
#include "paramset.h"
#include "conjugategrad.h"

#include <fstream>
#include <sstream>

namespace DDF { 


// *****************************************************************************
// pressure solve
// optional modifications: 
// - global volume correction value (set by compute volume plugin, passed over via mpPlParams)
// - per cell density values (grid), scales RHS values
// - surface tension values (grid), adds
class fsPrMakeRhs : public GridOpBaseFlagbord<1> {
	public:
		fsPrMakeRhs(FlagGrid *flags, Grid<Real> *rhs, Grid<Vec3> *vel, Real corr, 
				Grid<Real> *dens, Grid<Real> *percell) : 
			GridOpBaseFlagbord<1>(), mpRhs(rhs), mpVel(vel), mpPerCellCorrection(percell) , 
			mCorr(corr), mpDensity(dens) {
				mpFlags = flags;
				applyOperatorToGrids(this);
		};
		~fsPrMakeRhs() { }; 
		void resetVariables() { };
		void buildCallList() {
			gaVel.gridAccInit(mpVel, AM_READ, gaCalls); 
			gaRhs.gridAccInit(mpRhs, AM_WRITE, gaCalls); 
			if(mpDensity) {
				gaDens.gridAccInit(mpDensity, AM_READ, gaCalls); 
			}
			if(mpPerCellCorrection) {
				gaPerCellCorr.gridAccInit(mpPerCellCorrection, AM_READ, gaCalls); 
			}
			setFlags(mpFlags);
		};
		inline void operator() (int i, int j, int k) { 
			if (!fgIsFluid(getFlagAcc()(i,j,k)) ) {
				gaRhs.write(i,j,k) = 0.;
				return;
			}

			Real set = 0.;
			set += (gaVel(i,j,k)[0] - gaVel(i+1,j,k)[0]);
			set += (gaVel(i,j,k)[1] - gaVel(i,j+1,k)[1]);
#			if DDF_DIMENSION==3
			set += (gaVel(i,j,k)[2] - gaVel(i,j,k+1)[2]);
#			endif
			gaRhs.write(i,j,k) = set + mCorr;
			
			// multiply RHS by density value, if given...
			if(mpDensity) gaRhs.write(i,j,k) *= gaDens(i,j,k);
			if(mpPerCellCorrection) gaRhs.write(i,j,k) += gaPerCellCorr(i,j,k);
		};

		void reduce(fsPrMakeRhs &op) { };
	protected:
		Grid<Real> *mpRhs;
		Grid<Vec3> *mpVel;
		Grid<Real> *mpPerCellCorrection;
		GridAccessor<Real,0> gaRhs;
		GridAccessor<Vec3,1> gaVel;
		GridAccessor<Real,0> gaPerCellCorr;
		Real mCorr;
		// optional density
		Grid<Real> *mpDensity;
		GridAccessor<Real,0> gaDens;
}; // fsPrMakeRhs */

class fsPrMakeMatrix : public GridOpBaseFlagbord<1> {
	public:
		fsPrMakeMatrix(FlagGrid *flags, Grid<Real> *A0, Grid<Real> *Ai, Grid<Real> *Aj, Grid<Real> *Ak) :
			GridOpBaseFlagbord<1>(), mpA0(A0), mpAi(Ai), mpAj(Aj), mpAk(Ak) { 
			mpFlags = flags;
			applyOperatorToGrids(this);
		};
		~fsPrMakeMatrix() { }; 
		void resetVariables() { };
		void buildCallList() {
			gaA0.gridAccInit(mpA0, AM_WRITE, gaCalls); 
			gaAi.gridAccInit(mpAi, AM_WRITE, gaCalls); 
			gaAj.gridAccInit(mpAj, AM_WRITE, gaCalls); 
			gaAk.gridAccInit(mpAk, AM_WRITE, gaCalls); 
			setFlags(mpFlags);
		};
		inline void operator() (int i, int j, int k) { 
			const int currFlag = getFlagAcc()(i,j,k);
			if (!fgIsFluid(currFlag) ) return;

			if (!fgIsObstacle(getFlagAcc()(i-1,j,k)) ) gaA0.write(i,j,k) += 1.;
			if (!fgIsObstacle(getFlagAcc()(i+1,j,k)) ) gaA0.write(i,j,k) += 1.;
			if (fgIsFluid(getFlagAcc()(i+1,j,k)) )     gaAi.write(i,j,k) = -1.;

			if (!fgIsObstacle(getFlagAcc()(i,j-1,k)) ) gaA0.write(i,j,k) += 1.;
			if (!fgIsObstacle(getFlagAcc()(i,j+1,k)) ) gaA0.write(i,j,k) += 1.;
			if (fgIsFluid(getFlagAcc()(i,j+1,k)) )     gaAj.write(i,j,k) = -1.;

#			if DDF_DIMENSION==3
			if (!fgIsObstacle(getFlagAcc()(i,j,k-1)) ) gaA0.write(i,j,k) += 1.;
			if (!fgIsObstacle(getFlagAcc()(i,j,k+1)) ) gaA0.write(i,j,k) += 1.;
			if (fgIsFluid(getFlagAcc()(i,j,k+1)) )     gaAk.write(i,j,k) = -1.;
#			endif
		};
		void reduce(fsPrMakeMatrix &op) { };
	protected:
		Grid<Real> *mpA0, *mpAi, *mpAj, *mpAk;
		GridAccessor<Real,0> gaA0, gaAi, gaAj, gaAk;
}; // fsPrMakeMatrix */

class fsPrCorrVels : public GridOpBaseFlagbord<1> {
	public:
		fsPrCorrVels(FlagGrid *flags, Grid<Vec3> *vel, Grid<Real> *press, Grid<Real> *dens) : 
			GridOpBaseFlagbord<1>(), mpPress(press), mpVel(vel),
  			mpDensity(dens) { 
			mpFlags = flags;
			applyOperatorToGrids(this);
		};
		~fsPrCorrVels() { }; 
		void resetVariables() { };
		void buildCallList() {
			gaVel.gridAccInit(mpVel, AM_WRITE, gaCalls); 
			gaPress.gridAccInit(mpPress, AM_WRITE, gaCalls); 
			if(mpDensity) {
				gaDens.gridAccInit(mpDensity, AM_READ, gaCalls); 
			}
			setFlags(mpFlags);
		};
		inline void operator() (int i, int j, int k) { 
			const int currFlag = getFlagAcc()(i,j,k);
			if (!fgIsFluid(currFlag) ) return;
			//if (fgIsObstacle(currFlag) ) return;
			Real scale = 1.;
			//if(mpDensity) scale = 1. / gaDens(i,j,k);
			Real densCurr = 1.;
			if(mpDensity) densCurr = gaDens(i,j,k);

			if (fgIsFluid(getFlagAcc()(i-1,j,k)) ) {
				if(mpDensity) {
					scale = 1. / (0.5* (densCurr + gaDens(i-1,j,k)) );
				}
				gaVel.write(i,j,k)[0] -= scale* (gaPress(i,j,k) - gaPress(i-1,j,k) );
			}
			if (fgIsFluid(getFlagAcc()(i,j-1,k)) ) {
				if(mpDensity) {
					scale = 1. / (0.5* (densCurr + gaDens(i,j-1,k)) );
				}
				gaVel.write(i,j,k)[1] -= scale* (gaPress(i,j,k) - gaPress(i,j-1,k) );
			}
			if (fgIsFluid(getFlagAcc()(i,j,k-1)) ) {
				if(mpDensity) {
					scale = 1. / (0.5* (densCurr + gaDens(i,j,k-1)) );
				}
				gaVel.write(i,j,k)[2] -= scale* (gaPress(i,j,k) - gaPress(i,j,k-1) );
			}
		};
		void reduce(fsPrCorrVels &op) { };
	protected:
		Grid<Real> *mpPress;
		Grid<Vec3> *mpVel;
		GridAccessor<Real,1> gaPress;
		GridAccessor<Vec3,0> gaVel;
		// optional density
		Grid<Real> *mpDensity;
		GridAccessor<Real,0> gaDens;
}; // fsPrCorrVels */


// overwrite domain sides with default stencil
// TODO set sides individually?
static void setOpenBound( int bound, Grid<Real>* A0, Grid<Real>* Ai, Grid<Real>* Aj, Grid<Real>* Ak,  Grid<Vec3>* vel ) 
{
	Real fac = 6.;
	if(gDim==2) fac = 4.;

	if(!bound) return;

	// set velocities
	const int iMax = A0->getSizeX()-2;
	const int iMin = 1;
	const int jMax = A0->getSizeY()-2;
	const int jMin = 1;
	const int kMax = A0->getSizeZ()-2;
	const int kMin = 1;

	// TODO on/off per side...
	// TODO smoothen?
	bool dox = false;
	bool doy = true;
	bool doz = false;

	// set velocities
	if(dox)
		for(int j=0;j<A0->getSizeY();j++) 
			for(int k=0;k<A0->getSizeZ();k++) { 
				vel->getGlobal(iMax+1,j,k) = vel->getGlobal(iMax,j,k);
				vel->getGlobal(iMin-1,j,k) = vel->getGlobal(iMin,j,k);
			}

	if(doy)
	for(int i=0;i<A0->getSizeX();i++) 
			for(int k=0;k<A0->getSizeZ();k++) { 
				vel->getGlobal(i,jMax+1,k) = vel->getGlobal(i,jMax,k);
				vel->getGlobal(i,jMin-1,k) = vel->getGlobal(i,jMin,k);
			}

	if(doz)
	for(int i=0;i<A0->getSizeX();i++) 
		for(int j=0;j<A0->getSizeY();j++) { 
				vel->getGlobal(i,j,kMax+1) = vel->getGlobal(i,j,kMax);
				vel->getGlobal(i,j,kMin-1) = vel->getGlobal(i,j,kMin);
			}

	// set matrix stencils at boundary
	if(dox) {
	for(int i=0;i<=       1;i++) 
		for(int j=0;j<A0->getSizeY();j++) 
			for(int k=0;k<A0->getSizeZ();k++) { 
				A0->getGlobal(i,j,k) = fac;
				Ai->getGlobal(i,j,k) = -1.;
				Aj->getGlobal(i,j,k) = -1.;
				if(gDim==3) Ak->getGlobal(i,j,k) = -1.;
			}
	for(int i=A0->getSizeX()-2;i<=A0->getSizeX()-1;i++) 
		for(int j=0;j<A0->getSizeY();j++) 
			for(int k=0;k<A0->getSizeZ();k++) { 
				A0->getGlobal(i,j,k) = fac;
				Ai->getGlobal(i,j,k) = -1.;
				Aj->getGlobal(i,j,k) = -1.;
				if(gDim==3) Ak->getGlobal(i,j,k) = -1.;
			}
	}

	if(doy) {
	for(int i=0;i<A0->getSizeX();i++) 
		for(int j=0;j<=        1;j++) 
			for(int k=0;k<A0->getSizeZ();k++) { 
				A0->getGlobal(i,j,k) = fac;
				Ai->getGlobal(i,j,k) = -1.;
				Aj->getGlobal(i,j,k) = -1.;
				if(gDim==3) Ak->getGlobal(i,j,k) = -1.;
			}
	for(int i=0;i<A0->getSizeX();i++) 
		for(int j=A0->getSizeY()-2;j<=A0->getSizeY()-1;j++) 
			for(int k=0;k<A0->getSizeZ();k++) { 
				A0->getGlobal(i,j,k) = fac;
				Ai->getGlobal(i,j,k) = -1.;
				Aj->getGlobal(i,j,k) = -1.;
				if(gDim==3) Ak->getGlobal(i,j,k) = -1.;
			}
	}

	if(doz) 
		if(gDim==3) {
			for(int i=0;i<A0->getSizeX();i++) 
				for(int j=0;j<A0->getSizeY();j++) 
					for(int k=0;k<=        1;k++) { 
						A0->getGlobal(i,j,k) = fac;
						Ai->getGlobal(i,j,k) = -1.;
						Aj->getGlobal(i,j,k) = -1.;
						if(gDim==3) Ak->getGlobal(i,j,k) = -1.;
					}

			for(int i=0;i<A0->getSizeX();i++) 
				for(int j=0;j<A0->getSizeY();j++) 
					for(int k=A0->getSizeZ()-2;k<=A0->getSizeZ()-1;k++) {
						A0->getGlobal(i,j,k) = fac;
						Ai->getGlobal(i,j,k) = 1.;
						Aj->getGlobal(i,j,k) = 1.;
						if(gDim==3) Ak->getGlobal(i,j,k) = 1.;
					}
		}

	debMsg("FluidSolver::solvePressure","Open boundaries set");
}

static void setBottomOutflow( Grid<Real>* rhs, int height, std::string dir) { // char dir ) {
	const int iMax = rhs->getSizeX()-2;
	const int iMin = 1;
	const int jMax = rhs->getSizeY()-2;
	const int jMin = 1;
	const int kMax = rhs->getSizeZ()-2;
	const int kMin = 1;
	debMsg("setBottomOutflow","Active, using dir='"<<dir<<"' ");

	// remove rhs/divergence at bottom
	for(int trs=0; trs<(int)dir.length(); trs++) {

	switch(dir[trs]) {
		case 'x': // -x
			debMsg("setBottomOutflow","Active, -X");
			for(int i=0; i< height;i++) 
				for(int j=0;j<rhs->getSizeY();j++) 
					for(int k=0;k<rhs->getSizeZ();k++) { 
						rhs->getGlobal(i,j,k) = 0.;
					}
			break;
		case 'X': // +x
			debMsg("setBottomOutflow","Active, +X");
			for(int i=rhs->getSizeX()-1-height; i< rhs->getSizeX();i++) 
				for(int j=0;j<rhs->getSizeY();j++) 
					for(int k=0;k<rhs->getSizeZ();k++) { 
						rhs->getGlobal(i,j,k) = 0.;
					}
			break;
		case 'y': // -y
			debMsg("setBottomOutflow","Active, -Y");
			for(int i=0;i<rhs->getSizeX();i++) 
				for(int j=0;j<= height;j++) 
					for(int k=0;k<rhs->getSizeZ();k++) { 
						rhs->getGlobal(i,j,k) = 0.;
					}
			break;
		case 'Y': // +y
			debMsg("setBottomOutflow","Active, +Y");
			for(int i=0;i<rhs->getSizeX();i++) 
				for(int j=rhs->getSizeY()-1-height; j< rhs->getSizeY();j++) 
					for(int k=0;k<rhs->getSizeZ();k++) { 
						rhs->getGlobal(i,j,k) = 0.;
					}
			break;
		case 'z': // -z
			debMsg("setBottomOutflow","Active, -Z");
			for(int i=0;i<rhs->getSizeX();i++) 
				for(int j=0;j<= rhs->getSizeY();j++) 
					for(int k=0;k<height;k++) { 
						rhs->getGlobal(i,j,k) = 0.;
					}
			break;
		case 'Z': // +z
			debMsg("setBottomOutflow","Active, +Z");
			for(int i=0;i<rhs->getSizeX();i++) 
				for(int j=0;j<rhs->getSizeY();j++) 
					for(int k=rhs->getSizeZ()-1-height; k< rhs->getSizeZ();k++) { 
						rhs->getGlobal(i,j,k) = 0.;
					}
			break;

		default:
			errFatal("setBottomOutflow","Invalid dir "<<dir<<" ",SIMWORLD_GENERICERROR);
	} // dir

	} // string
}

// apply stored symmetric matix 
// simple version - no inclusion of e.g., surface tension forces
class goApplyMatrix : public GridOpBase{
	public:
		goApplyMatrix(Grid<Real> *Dst, Grid<Real> *Src, 
				Grid<Real> *A0, Grid<Real> *Ai, Grid<Real> *Aj, Grid<Real> *Ak, 
				FlagGrid *Flags,
			  	Grid<Real> *unusedReal,
			  	Grid<Vec3> *unusedVec  ) : GridOpBase(),
   			mpA0(A0), mpAi(Ai), mpAj(Aj), mpAk(Ak) {
			mpDst = Dst; mpSrc = Src;
			//mpMat = Mat;
			mpFlags = Flags;
			applyOperatorToGrids(this);
		};
		~goApplyMatrix() { };
		void buildCallList() {
			gaDst.gridAccInit(mpDst, AM_WRITE, gaCalls); 
			gaSrc.gridAccInit(mpSrc, AM_READ, gaCalls); 
			gaA0.gridAccInit(mpA0, AM_READ, gaCalls); 
			gaAi.gridAccInit(mpAi, AM_READ, gaCalls); 
			gaAj.gridAccInit(mpAj, AM_READ, gaCalls); 
			gaAk.gridAccInit(mpAk, AM_READ, gaCalls); 
			setFlags(mpFlags);
		}; 
		inline void operator() (int i, int j, int k) {
			if(!fgIsFluid(getFlagAcc()(i,j,k)) ) {
				//gaDst.write(i,j,k) = 0.; // org
				gaDst.write(i,j,k) = gaSrc(i,j,k);
				//gaDst.write(i,j,k) = gaSrc(i,j,k) * gaA0(i ,j ,k);
				return; 
			}

			gaDst.write(i,j,k) = 
				    gaSrc(i ,j ,k) * gaA0(i ,j ,k)
					+ gaSrc(i-1,j,k) * gaAi(i-1,j,k) 
					+ gaSrc(i+1,j,k) * gaAi(i  ,j,k)
					+ gaSrc(i,j-1,k) * gaAj(i,j-1,k) 
					+ gaSrc(i,j+1,k) * gaAj(i,j  ,k)
					+ gaSrc(i,j,k-1) * gaAk(i,j,k-1) 
					+ gaSrc(i,j,k+1) * gaAk(i,j,k  ); 
		} 

		void getRidOfWarning() { };
		void reduce(goApplyMatrix &op) { getRidOfWarning(); } 
	protected:
		Grid<Real> *mpDst, *mpSrc;
		GridAccessor<Real,0> gaDst;
		GridAccessor<Real,1> gaSrc;

		Grid<Real> *mpA0, *mpAi, *mpAj, *mpAk;
		GridAccessor<Real,1> gaA0, gaAi, gaAj, gaAk;
}; // goApplyMatrix

//*****************************************************************************

// add forces to vels
class spluginSolvePressure : public SolverPlugin {
	public:
		spluginSolvePressure() : SolverPlugin(),
   				mVels("-unn1-"), mFlags("-unn1-"),
   				mPress("-unn1-"), mRhs("-unn1-"),
   				mRes("-unn1-"), mSrch("-unn1-"),
   				mTmp("-unn1-"), mDensity("-unn1-"), mOpenBound(0), 
					mBottomOutflow(0), mBottomOutflowDir("X"),
					mStPerCellCor("-unn2-"), mUseResNorm(true), mPcMethod(0),
   				mNamePcA0("PCA0"), mNamePcAi("PCAi"), mNamePcAj("PCAj"), mNamePcAk("PCAk")
			{ };
		~spluginSolvePressure() { };

		virtual bool parseParams(const ParamSet& params) {
			mVels = params.FindOneString("vel", mVels );
			mPress = params.FindOneString("pressure", mPress );
			mDensity = params.FindOneString("density", mDensity );
			mStPerCellCor = params.FindOneString("st-per-cell-corr", mStPerCellCor );

			// optional grids for preconditining, not all PCs need all four grids!
			mNamePcA0 = params.FindOneString("grid-pca0", mNamePcA0 );
			mNamePcAi = params.FindOneString("grid-pcai", mNamePcAi );
			mNamePcAj = params.FindOneString("grid-pcaj", mNamePcAj );
			mNamePcAk = params.FindOneString("grid-pcak", mNamePcAk );

			mOpenBound = params.FindOneInt("openbound", mOpenBound );
			mPcMethod = params.FindOneInt("precondition", mPcMethod );
			mUseResNorm = params.FindOneInt("use-res-norm", mUseResNorm );
			
			mBottomOutflow = params.FindOneInt("bottom-outflow", mBottomOutflow ); 
			//std::string boStr = std::string("y");
			mBottomOutflowDir = params.FindOneString("bottom-outflow-dir", mBottomOutflowDir );
			//if(boStr.length()>=1) { mBottomOutflowDir = boStr[0]; }

			return true;
		};
		virtual bool initPlugin() {
			debMsg("spluginSolvePressure","init");
			return true;
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("spluginSolvePressure","step "<<dt<<" velsgrid:"<<mVels<<", precond:"<<mPcMethod<<", use-res-norm:"<<mUseResNorm); 

			// TODO fixme, hard coded names for now
			mVels = std::string("vel-curr");
			mFlags = std::string("flags");
			mPress = std::string("pressure");
			mRhs = std::string("rhs");
			mRes = std::string("residual");
			mSrch = std::string("search");
			mTmp = std::string("tmp");

			std::ostringstream pdebug;
			const bool debugSP= false;
			Grid<Vec3> *currVel   = mpPlParams->getGridVec3( mVels );
			FlagGrid*   pFlags     = mpPlParams->getGridInt( mFlags );
			Grid<Real>* pPressure = mpPlParams->getGridReal( mPress );
			Grid<Real>* pRhs      = mpPlParams->getGridReal( mRhs );
			Grid<Real>* pResidual = mpPlParams->getGridReal( mRes );
			Grid<Real>* pSearch   = mpPlParams->getGridReal( mSrch );
			Grid<Real>* pTmp      = mpPlParams->getGridReal( mTmp );

			Grid<Real>* pA0 = mpPlParams->getGridReal( "A0" );
			Grid<Real>* pAi = mpPlParams->getGridReal( "Ai" );
			Grid<Real>* pAj = mpPlParams->getGridReal( "Aj" );
			Grid<Real>* pAk = mpPlParams->getGridReal( "Ak" );

			Grid<Real>* pDens = NULL; 
			if(mpPlParams->haveGridReal( mDensity )) {
				pDens = mpPlParams->getGridReal( mDensity );
				debMsg("spluginSolvePressure","Pressure correction RHS using density grid '"<<mDensity<<"' ");
			}

			Grid<Real>* pStPerCellCorr = NULL; 
			if(mpPlParams->haveGridReal( mStPerCellCor )) {
				pStPerCellCorr = mpPlParams->getGridReal( mStPerCellCor );
				debMsg("spluginSolvePressure","Pressure correction RHS using per cell ST-correction grid '"<<mStPerCellCor<<"' ");
			}

			if (debugSP) {
				pdebug<<"\n PRE \n\n";
				pdebug<<"P \n"<< pPressure->toString() <<" \n";
				pdebug<<"CURRVEL \n"<< currVel->toString() <<" \n";
			}
			
			// reset necessary
			GridOpTouchMemory<Real>(pFlags, pA0, 0.);
			GridOpTouchMemory<Real>(pFlags, pAi, 0.);
			GridOpTouchMemory<Real>(pFlags, pAj, 0.);
			GridOpTouchMemory<Real>(pFlags, pAk, 0.);
			fsPrMakeMatrix(pFlags, pA0,pAi,pAj,pAk);
			setOpenBound(mOpenBound, pA0,pAi,pAj,pAk, currVel);

			// compute divergence and init right hand side
			fsPrMakeRhs(pFlags, pRhs, currVel, 
					mpPlParams->getDivergenceCorrection()  // usually is zero, set by compute-volume-correction plugin
					, pDens , pStPerCellCorr );
			if(mBottomOutflow>0) {
				setBottomOutflow(pRhs, mBottomOutflow, mBottomOutflowDir);
			}

			const int maxIter = (int)(mpPlParams->mCgMaxIterFac *VMAX(pFlags->getSize()) );
			if (1) {
				GridCgInterface *gcg	= NULL;

				GridCg<goApplyMatrix> *mgcg	= new GridCg<goApplyMatrix>(
						pPressure, pRhs, pResidual, pSearch, pFlags, pTmp, pA0,pAi,pAj,pAk );
				gcg = mgcg;
				gcg->setAccuracy( mpPlParams->mCgAccuracy ); 
				gcg->setUseResNorm( mUseResNorm );

				// optional preconditioning, with optional grids
				Grid<Real> *pPCA0 = NULL, *pPCAi=NULL, *pPCAj=NULL, *pPCAk=NULL;
				if(mpPlParams->haveGridReal( mNamePcA0 )) 
					pPCA0 = mpPlParams->getGridReal( mNamePcA0 );
				if(mpPlParams->haveGridReal( mNamePcAi )) 
					pPCAi = mpPlParams->getGridReal( mNamePcAi );
				if(mpPlParams->haveGridReal( mNamePcAj )) 
					pPCAj = mpPlParams->getGridReal( mNamePcAj );
				if(mpPlParams->haveGridReal( mNamePcAk )) 
					pPCAk = mpPlParams->getGridReal( mNamePcAk );
				gcg->setPreconditioner( mPcMethod, pPCA0, pPCAi, pPCAj, pPCAk ); 

				for (int iter=0; iter<maxIter; iter++) {
					if (!gcg->iterate()) iter=maxIter;
				} 
				debMsg("FluidSolver::solvePressure","iterations:"<<gcg->getIterations()<<", res:"<<gcg->getSigma() );
				delete gcg;
			} else {
				// debug tests
				for (int iter=0; iter<4*maxIter; iter++) {
					goJacobiSolver(pTmp, pPressure, pRhs, pFlags);
					goCopyGrid<Real>(pPressure, pTmp);
				}
				debMsg("FluidSolver::solvePressure","jacobi iterations:"<<4*maxIter );
			}
			fsPrCorrVels(pFlags, currVel, pPressure, pDens);

			if (debugSP) {
				goCompMinMax<Real> mm = goCompMinMax<Real>(pPressure, pFlags);
				debMsg("FluidSolver::solvePressure","minmax "<< mm.toString() );

				pdebug<<"\n POST \n\n";
				pdebug<<"RHS \n"<< pRhs->toString() <<" \n";
				pdebug<<"P \n"<< pPressure->toString() <<" \n";
				pdebug<<"CURRVEL \n"<< currVel->toString() <<" \n";

				pdebug<<"A0 \n"<< pA0->toString() <<" \n";
				pdebug<<"Ai \n"<< pAi->toString() <<" \n";
				pdebug<<"Aj \n"<< pAj->toString() <<" \n";
				pdebug<<"Ak \n"<< pAk->toString() <<" \n";
				std::fstream df; // DEBUG
				df.open("debugOut_pressure.txt", std::ios::out);
				df << pdebug.str(); 
				df.close(); // DEBUG
			}

			return true;
		};

	protected:
		// grid names to swap
		std::string mVels, mFlags, mPress, mRhs, mRes, mSrch, mTmp, mDensity;
		// open boundaries?
		int mOpenBound;
		int mBottomOutflow; 
		std::string mBottomOutflowDir;
		// grid name to for optional curvatuere/surface tension grid
		std::string mStPerCellCor;

		// use residual norm or max entry?
		bool mUseResNorm;
		// preconditiong
		int mPcMethod;
		std::string mNamePcA0, mNamePcAi, mNamePcAj, mNamePcAk;
};

//*****************************************************************************

template < class Scalar >
class GridOpDiffuse : public GridOpBase {
	public:
		GridOpDiffuse(Grid<Scalar>* pSrc,Grid<Scalar>* pDst, FlagGrid* flags, Real str)
			: GridOpBase(), mpSrc(pSrc), mpDst(pDst), mStrength(str) {
			mpFlags = flags;
			applyOperatorToGrids( this );
		} 
		~GridOpDiffuse() {};
		void resetVariables() { } 
		void reduce(GridOpDiffuse &op) { }
		
		void buildCallList() {
			gaSrc.gridAccInit(mpSrc, DDF::AM_WRITE, gaCalls); 
			gaDst.gridAccInit(mpDst, DDF::AM_WRITE, gaCalls); 
			setFlags(mpFlags);
		}; 

		inline void operator() (int i, int j, int k) {
			if(!fgIsFluid(getFlagAcc()(i,j,k)) ) {
				gaDst.write(i,j,k) = gaSrc(i,j,k);
			  	return;
			}

			const Scalar curr = gaSrc(i,j,k);
			Scalar avg = 0.;

			if(1 && fgIsFluid(getFlagAcc()(i+1,j,k)) ) {
				avg +=       gaSrc(i+1,j,k);
			} else { avg += curr; }

			if(1 && fgIsFluid(getFlagAcc()(i-1,j,k)) ) {
				avg +=       gaSrc(i-1,j,k);
			} else { avg += curr; }

			if(1 && fgIsFluid(getFlagAcc()(i,j+1,k)) ) {
				avg +=       gaSrc(i,j+1,k);
			} else { avg += curr; }

			if(1 && fgIsFluid(getFlagAcc()(i,j-1,k)) ) {
				avg +=       gaSrc(i,j-1,k);
			} else { avg += curr; }

			if(gDim==3) {
				if(1 && fgIsFluid(getFlagAcc()(i,j,k+1)) ) {
					avg +=       gaSrc(i,j,k+1);
				} else { avg += curr; }

				if(1 && fgIsFluid(getFlagAcc()(i,j,k-1)) ) {
					avg +=       gaSrc(i,j,k-1);
				} else { avg += curr; }
			}

			avg *= (gDim==3 ? 1./6. : 1./4.);
			//gaDst.write(i,j,k) = (1.-mStrength) * curr + mStrength * avg; 
			gaDst.write(i,j,k) = avg*mStrength -curr*mStrength + curr;
		}

	protected:
		Grid<Scalar> *mpSrc, *mpDst;
		GridAccessor<Scalar,1> gaSrc;
		GridAccessor<Scalar,0> gaDst;
		Real mStrength;
}; // GridOpDiffuse */

class spDiffuseGrid : public SolverPlugin {
	public:
		spDiffuseGrid() : 
			SolverPlugin(),
			mDest("real-temp") , mDestVec("vec-temp"), mDiffusivity(1.)
			{ }; 
		~spDiffuseGrid() { };

		virtual bool parseParams(const ParamSet& params) {
			mDiffusivity = params.FindOneFloat("diff", mDiffusivity );
			mSrc = params.FindOneString("src-real", mSrc );
			mSrcVec = params.FindOneString("src-vec3", mSrcVec );
			mDest = params.FindOneString("dst-real", mDest );
			mDestVec = params.FindOneString("dst-vec3", mDestVec );
			return true;
		};
		virtual bool initPlugin() { return true; };


		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			FlagGrid*  pFlags  = mpPlParams->getGridInt( "flags" ); 
			Real diff = mDiffusivity * dt;
			//if(diff<0. || diff>1.) debMsg("spDiffuseGrid","WARNING - diffusion param invalid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
			//CLAMP(diff, (Real)0.,(Real)1.);

			int cnt=0;
			while(diff > 0.) {
				debMsg("spDiffuseGrid","str"<<mDiffusivity<<",dt="<<dt<<", diff="<<diff<<". from "<<mSrc<<","<<mSrcVec<<" to "<<mDest<<","<<mDestVec); 

				if(mSrc.length()>0) {
					Grid<Real>* pSrc = mpPlParams->getGridReal(mSrc); 
					Grid<Real>* pDest = mpPlParams->getGridReal(mDest); 
					GridOpDiffuse<Real> godiv(pSrc, pDest, pFlags, (diff>1.? 1.:diff) );
					cnt++;

					swapGrids(mpPlParams, mSrc, mDest);
				}
				if(mSrcVec.length()>0) {
					Grid<Vec3>* pSrcv = mpPlParams->getGridVec3(mSrcVec); 
					Grid<Vec3>* pDestv = mpPlParams->getGridVec3(mDestVec); 
					GridOpDiffuse<Vec3> godiv(pSrcv, pDestv, pFlags, (diff>1.? 1.:diff) );
					cnt++;

					swapGrids(mpPlParams, mSrcVec, mDestVec);
				}

				diff -= 1.;
			}

			if(cnt==0 && mDiffusivity>0.) {
				errFatal("spDiffuse","No grid for diffusion found or given!", SIMWORLD_PLUGINERROR);
				return false;
			}
			//debMsg("spDiffuseGrid","abs avg="<<godiv.getAbsAvg()<<",abs max="<<godiv.getAbsMax());
			return true;
		};

	protected:
		std::string mSrc, mSrcVec;
		std::string mDest, mDestVec;
		Real mDiffusivity;
};


class fsDiffMakeMatrix : public GridOpBaseFlagbord<1> {
	public:
		fsDiffMakeMatrix(FlagGrid *flags, Grid<Real> *A0, Grid<Real> *Ai, Grid<Real> *Aj, Grid<Real> *Ak, Real coeff) :
			GridOpBaseFlagbord<1>(), mpA0(A0), mpAi(Ai), mpAj(Aj), mpAk(Ak), mCoeff(coeff) { 
			mpFlags = flags;
			applyOperatorToGrids(this);
		};
		~fsDiffMakeMatrix() { }; 
		void resetVariables() { };
		void buildCallList() {
			gaA0.gridAccInit(mpA0, AM_WRITE, gaCalls); 
			gaAi.gridAccInit(mpAi, AM_WRITE, gaCalls); 
			gaAj.gridAccInit(mpAj, AM_WRITE, gaCalls); 
			gaAk.gridAccInit(mpAk, AM_WRITE, gaCalls); 
			setFlags(mpFlags);
		};
		inline void operator() (int i, int j, int k) { 
			const int currFlag = getFlagAcc()(i,j,k);

			//gaA0.write(i,j,k) = 5.; gaAi.write(i,j,k) = -1.; gaAj.write(i,j,k) = -1.; gaAk.write(i,j,k) = 0.; return;

			Real& a0 = gaA0.write(i,j,k);
			a0 = 1.;

			// handle diffusion across free surface boundaries, e.g. for velocity
			if (fgIsEmpty(getFlagAcc()(i,j,k)) ) {
				if (fgIsFluid(getFlagAcc()(i+1,j,k)) ) {
					gaAi.write(i,j,k) -= mCoeff;
				}
				if (fgIsFluid(getFlagAcc()(i,j+1,k)) ) {
					gaAj.write(i,j,k) -= mCoeff;
				}
				if (gDim==3 && fgIsFluid(getFlagAcc()(i,j,k+1)) ) {
					gaAk.write(i,j,k) -= mCoeff;
				}
			}

			if (!fgIsFluid(currFlag) ) return;

			if (!fgIsObstacle(getFlagAcc()(i+1,j,k)) ) {
				a0 += mCoeff;
			}
			if (!fgIsObstacle(getFlagAcc()(i-1,j,k)) ) {
				a0 += mCoeff;
			}
			if(1) {
				gaAi.write(i,j,k) -= mCoeff;
			}

			if (!fgIsObstacle(getFlagAcc()(i,j+1,k)) ) {
				a0 += mCoeff;                     
			}                                   
			if (!fgIsObstacle(getFlagAcc()(i,j-1,k)) ) {
				a0 += mCoeff;
			}
			//if (fgIsFluid(getFlagAcc()(i,j+1,k)) ) {
			//if ( (fgIsFluid(getFlagAcc()(i,j-1,k)) ) &&  (fgIsFluid(getFlagAcc()(i,j+1,k)) )   ) {
			if(1) {
			//if (fgIsFluid(getFlagAcc()(i,j+1,k)) ) {
				gaAj.write(i,j,k) -= mCoeff;
			}

#			if DDF_DIMENSION==3
			if (!fgIsObstacle(getFlagAcc()(i,j,k+1)) ) {
				a0 += mCoeff;                       
			}                                     
			if (!fgIsObstacle(getFlagAcc()(i,j,k-1)) ) {
				a0 += mCoeff;
			}
			//if (fgIsFluid(getFlagAcc()(i,j,k+1)) ) {
			if(1) {
				gaAk.write(i,j,k) -= mCoeff;
			}
#			endif

			//debMsg("D "," at "<<PRINT_IJK<<"a = "<< gaA0.write(i,j,k)<<", "<< gaAi.write(i,j,k)<<", "<< gaAj.write(i,j,k)<<", "<< gaAk.write(i,j,k)<<", "); 
		};
		void reduce(fsDiffMakeMatrix &op) { };
	protected:
		Grid<Real> *mpA0, *mpAi, *mpAj, *mpAk;
		GridAccessor<Real,0> gaA0, gaAi, gaAj, gaAk;
		Real mCoeff;
}; // fsDiffMakeMatrix */

// add forces to vels
class spluginSolveDiffusion : public SolverPlugin {
	public:
		spluginSolveDiffusion() : SolverPlugin(),
   				mDst("-unnamed1-"), mRhs("-unnamed1-"),
   				mRes("-unnamed1-"), mSrch("-unnamed1-"),
   				mTmp("-unnamed1-"), mOpenBound(0), 
					mDiffusion(0.), mAccuracy(1.), mMaxIter(1.) { };
		~spluginSolveDiffusion() { };

		virtual bool parseParams(const ParamSet& params) {
			mDst = params.FindOneString("dst-real", mDst );
			// source = rhs
			mRhs = params.FindOneString("src-real", mRhs );

			mVec3Dst = params.FindOneString("dst-vec3", mVec3Dst );
			mVec3Rhs = params.FindOneString("src-vec3", mVec3Rhs );

			mOpenBound = params.FindOneInt("openbound", mOpenBound );
			mDiffusion = params.FindOneFloat("diffusion", mDiffusion );

			// CG parameters, modify basic fluid solver settings, default = 1
			mAccuracy = params.FindOneFloat("accuracy", mAccuracy );
			mMaxIter = params.FindOneFloat("max-iter", mMaxIter);
			return true;
		};
		virtual bool initPlugin() {
			debMsg("spluginSolveDiffusion","init");
			return true;
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {

			const Real dx = mpPlParams->getDeltaX();
			const Real coeff = dt * mDiffusion / (dx*dx);

			debMsg("spluginSolveDiffusion","step "<<dt<<" diffusionCoeff:"<<mDiffusion<<" -> coeff="<<coeff); 

			// TODO fixme, hard coded names for now
			mRes = std::string("residual");
			mSrch = std::string("search");
			mTmp = std::string("tmp");

			std::ostringstream pdebug;
			const bool debugSP= false;
			FlagGrid*   pFlags     = mpPlParams->getFluidSolver()->getGridFlags(); // getGridInt( mFlags );
			Grid<Real>* pResidual = mpPlParams->getGridReal( mRes );
			Grid<Real>* pSearch   = mpPlParams->getGridReal( mSrch );
			Grid<Real>* pTmp      = mpPlParams->getGridReal( mTmp );

			Grid<Real>* pA0 = mpPlParams->getGridReal( "A0" );
			Grid<Real>* pAi = mpPlParams->getGridReal( "Ai" );
			Grid<Real>* pAj = mpPlParams->getGridReal( "Aj" );
			Grid<Real>* pAk = mpPlParams->getGridReal( "Ak" );

			GridOpTouchMemory<Real>(pFlags, pA0, 0.);
			GridOpTouchMemory<Real>(pFlags, pAi, 0.);
			GridOpTouchMemory<Real>(pFlags, pAj, 0.);
			GridOpTouchMemory<Real>(pFlags, pAk, 0.);
			
			if( (mpPlParams->haveGridVec3( mVec3Rhs )) && (mpPlParams->haveGridVec3( mVec3Dst )) ) {
				// vec3 version
				Grid<Vec3>* pVec3Rhs	= mpPlParams->getGridVec3( mVec3Rhs );
				Grid<Vec3>* pVec3Dst	= mpPlParams->getGridVec3( mVec3Dst );

				// temp real versions, make sure these grids are not used multiply!
				Grid<Real>* pRhs      = mpPlParams->getGridReal( mRhs );
				Grid<Real>* pDst 		 = mpPlParams->getGridReal( mDst );

				// setup diffusion matrix
				fsDiffMakeMatrix(pFlags, pA0,pAi,pAj,pAk, coeff);

				int maxIter = (int)(mMaxIter *mpPlParams->mCgMaxIterFac *VMAX(pFlags->getSize()) ); 
				if(maxIter <= 0) maxIter = 1;


				if(1) {
					goCopyVec3ToScalar<Real,0>(pRhs, pVec3Rhs); 
					//goCopyVec3ToScalar<Real,0>(pDst, pVec3Dst); 
					goCopyVec3ToScalar<Real,0>(pDst, pVec3Rhs); 

					GridCg<goApplyMatrix> gcg = GridCg<goApplyMatrix>(pDst, pRhs, pResidual, pSearch, pFlags, pTmp, pA0,pAi,pAj,pAk );
					gcg.setAccuracy( mAccuracy *mpPlParams->mCgAccuracy ); 
					gcg.solve( maxIter );
					debMsg("FluidSolver::solveDiffusion Vec3","iterations:"<<gcg.getIterations()<<", res:"<<gcg.getSigma() );

					goCopyScalarToVec3<Real,0>(pVec3Dst,pDst); 
				} else {
					goCopyVec3ToScalar<Real,0>(pRhs,pVec3Rhs); 
					goCopyScalarToVec3<Real,0>(pVec3Dst,pRhs); 
				} 


				if(1) {
					goCopyVec3ToScalar<Real,1>(pRhs,pVec3Rhs); 
					//goCopyVec3ToScalar<Real,1>(pDst,pVec3Dst); 
					goCopyVec3ToScalar<Real,1>(pDst,pVec3Rhs); 

					GridCg<goApplyMatrix> gcg	= GridCg<goApplyMatrix>(pDst, pRhs, pResidual, pSearch, pFlags, pTmp, pA0,pAi,pAj,pAk );
					gcg.setAccuracy( mAccuracy *mpPlParams->mCgAccuracy ); 
					gcg.solve( maxIter );
					debMsg("FluidSolver::solveDiffusion Vec3","iterations:"<<gcg.getIterations()<<", res:"<<gcg.getSigma() );

					goCopyScalarToVec3<Real,1>(pVec3Dst,pDst); 
				} else {
					goCopyVec3ToScalar<Real,1>(pRhs,pVec3Rhs); 
					goCopyScalarToVec3<Real,1>(pVec3Dst,pRhs); 
				}


				if(1) {
					goCopyVec3ToScalar<Real,2>(pRhs,pVec3Rhs); 
					goCopyVec3ToScalar<Real,2>(pDst,pVec3Dst); 

					GridCg<goApplyMatrix> gcg	= GridCg<goApplyMatrix>(pDst, pRhs, pResidual, pSearch, pFlags, pTmp, pA0,pAi,pAj,pAk );
					gcg.setAccuracy( mAccuracy *mpPlParams->mCgAccuracy ); 
					gcg.solve( maxIter );
					debMsg("FluidSolver::solveDiffusion Vec3","iterations:"<<gcg.getIterations()<<", res:"<<gcg.getSigma() );

					goCopyScalarToVec3<Real,2>(pVec3Dst,pDst); 
				} else {
					goCopyVec3ToScalar<Real,2>(pRhs,pVec3Rhs); 
					goCopyScalarToVec3<Real,2>(pVec3Dst,pRhs); 
				} 

				swapGrids(mpPlParams, mVec3Dst, mVec3Rhs);

				const bool debugDIFF = false;
				if (debugDIFF) {
					goCompMinMax<Real> mm = goCompMinMax<Real>(pDst, pFlags);
					debMsg("FluidSolver::solvePressure","minmax "<< mm.toString() );

					pdebug<<"\n POST \n\n";
					pdebug<<"RHS \n"<< pRhs->toString() <<" \n";
					pdebug<<"Dst \n"<< pDst->toString() <<" \n";

					pdebug<<"A0 \n"<< pA0->toString() <<" \n";
					pdebug<<"Ai \n"<< pAi->toString() <<" \n";
					pdebug<<"Aj \n"<< pAj->toString() <<" \n";
					pdebug<<"Ak \n"<< pAk->toString() <<" \n";
					std::fstream df; // DEBUG
					df.open("debugOut_diffusionVec.txt", std::ios::out);
					df << pdebug.str(); 
					df.close(); // DEBUG
				}

			} else if( (mpPlParams->haveGridReal( mRhs )) && (mpPlParams->haveGridReal( mDst )) ) {
				// real version
				Grid<Real>* pRhs      = mpPlParams->getGridReal( mRhs );
				Grid<Real>* pDst 		 = mpPlParams->getGridReal( mDst );

				// setup diffusion matrix
				fsDiffMakeMatrix(pFlags, pA0,pAi,pAj,pAk, coeff);

				int maxIter = (int)(mMaxIter *mpPlParams->mCgMaxIterFac *VMAX(pFlags->getSize()) ); 
				if(maxIter <= 0) maxIter = 1;

				GridCg<goApplyMatrix> *gcg	= new GridCg<goApplyMatrix>(pDst, pRhs, pResidual, pSearch, pFlags, pTmp, pA0,pAi,pAj,pAk );
				gcg->setAccuracy( mAccuracy *mpPlParams->mCgAccuracy ); 

				for (int iter=0; iter<maxIter; iter++) {
					if (!gcg->iterate()) iter=maxIter;
				} 
				debMsg("FluidSolver::solveDiffusion Real","iterations:"<<gcg->getIterations()<<", res:"<<gcg->getSigma() );
				delete gcg;

				const bool debugDIFF = false;
				if (debugDIFF) {
					goCompMinMax<Real> mm = goCompMinMax<Real>(pDst, pFlags);
					debMsg("FluidSolver::solvePressure","minmax "<< mm.toString() );

					pdebug<<"\n POST \n\n";
					pdebug<<"RHS \n"<< pRhs->toString() <<" \n";
					pdebug<<"Dst \n"<< pDst->toString() <<" \n";

					pdebug<<"A0 \n"<< pA0->toString() <<" \n";
					pdebug<<"Ai \n"<< pAi->toString() <<" \n";
					pdebug<<"Aj \n"<< pAj->toString() <<" \n";
					pdebug<<"Ak \n"<< pAk->toString() <<" \n";
					std::fstream df; // DEBUG
					df.open("debugOut_diffusionReal.txt", std::ios::out);
					df << pdebug.str(); 
					df.close(); // DEBUG
				}

				swapGrids(mpPlParams, mDst, mRhs);
			}

			return true;
		};

	protected:
		// grid names to swap
		std::string mDst, mRhs, mRes, mSrch, mTmp;
		std::string mVec3Dst, mVec3Rhs;
		// open boundaries?
		int mOpenBound;
		// diffusivity
		Real mDiffusion, mAccuracy, mMaxIter;
};

// create one of the standard hard coded plugins
SolverPlugin* MakePoissionSolverPlugins(std::string name) {

	if(name.compare( "solve-pressure")==0) {
		return new spluginSolvePressure;

	} else if(name.compare( string("diffuse-grid") )==0) {
		return new spDiffuseGrid;

	} else if(name.compare( string("solve-diffusion") )==0) {
		return new spluginSolveDiffusion;

	}
	
	return NULL;
}

} // end namespace DDF 

