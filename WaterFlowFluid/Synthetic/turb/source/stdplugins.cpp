/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Basic plugins (copy, set_bounds, etc.)
 *
 *****************************************************************************/

#include "fluidsolver.h"
#include "solverplugin.h"
#include "paramset.h"

#include "matrixbase.h"

// safety boundary for semi lagrange advection step
#define SLADVBOUND 2

namespace DDF { 

static inline bool doSetVelocity(int nbflag, int myflag) {
	if (fgIsObstacle(nbflag) && !fgIsInflow(myflag) ) {
		return false;
	}
	return true;
}

//*****************************************************************************

class fsSetNoslipBcs : public GridOpBaseFlagbord<1> {
	public:
		fsSetNoslipBcs(FlagGrid *flags, Grid<Vec3> *src) :
			GridOpBaseFlagbord<1>(), mpVecgrid(src) {
			mpFlags = flags;
			applyOperatorToGrids(this);
		};
		~fsSetNoslipBcs() { }; 
		void resetVariables() { };
		void buildCallList() {
			gaVecgrid.gridAccInit(mpVecgrid, AM_WRITE, gaCalls); 
			setFlags(mpFlags);
		};

		// set velocity at NB cell, given current flag myflag?
		inline bool __doSetVelocity(int nbflag, int myflag) {
			if (fgIsObstacle(nbflag) && !fgIsInflow(myflag) ) {
				return false;
			}
			return true;
		}

		// add forces and update empty cells free surface boundaries
		inline void operator() (int i, int j, int k) { 
			int myFlag = getFlagAcc()(i,j,k);
			Vec3& vel = gaVecgrid.write(i,j,k);

			if ( !doSetVelocity(myFlag, myFlag) ) {
				vel = Vec3(0.);
				return;
			}
			
			if (fgIsEmpty(myFlag) ) 
				return;

			if (fgIsFluid(myFlag) ) { 
				// only add for faces that are between fluid/fluid fl/empty cells
				if (  !doSetVelocity(getFlagAcc()(i-1,j,k), myFlag) ) { vel[0] = 0.; }

				if (  !doSetVelocity(getFlagAcc()(i,j-1,k), myFlag) ) { vel[1] = 0.; }

				if ( (!doSetVelocity(getFlagAcc()(i,j,k-1), myFlag) ) && (gDim==3) ) { vel[2] = 0.; }

				return;
			} 
		};
		void reduce(fsSetNoslipBcs &op) { };

	protected:
		Grid<Vec3> *mpVecgrid;
		GridAccessor<Vec3,1> gaVecgrid;
		Vec3 mForce;
		int mModFunc;
		Real mThresh;

}; // fsSetNoslipBcs */

// add forces to vels
class spSetNoslipBcs : public SolverPlugin {
	public:
		spSetNoslipBcs() : SolverPlugin(), mVel("-unnamed1-") { };
		~spSetNoslipBcs() { };
		virtual bool initPlugin() { return true; }; 
		virtual bool parseParams(const ParamSet& params) {
			mVel = params.FindOneString("grid", mVel );
			return true;
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("spSetNoslipBcs","step "<<dt<<" velsgrid:"<<mVel<<" to "<<mVel); 
			Grid<Vec3>* velsrc = mpPlParams->getGridVec3(mVel);
			fsSetNoslipBcs(mpPlParams->getFluidSolver()->getGridFlags(), velsrc );
			return true;
		};
	protected:
		std::string mVel;
};

// copy source grid to destination
class spluginCopyGrid : public SolverPlugin {
	public:
		spluginCopyGrid() : SolverPlugin(),
   				mSrc("-unnamed1-"),
   				mDest("-unnamed2-")
			{ };
		~spluginCopyGrid() { };

		virtual bool parseParams(const ParamSet& params) {
			mDest = params.FindOneString("dest", mDest );
			mSrc = params.FindOneString("src", mSrc );
			return true;
		};
		virtual bool initPlugin() {
			debMsg("spluginCopyGrid","init");
			return true;
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("spluginCopyGrid","step "<<dt);
			if(mpPlParams->haveGridInt(mSrc) && 1) {
				Grid<int>* dest = mpPlParams->getGridInt(mDest);
				Grid<int>* src = mpPlParams->getGridInt(mSrc);
				goCopyGrid<int>(dest, src);
			} else if(mpPlParams->haveGridReal(mSrc) && 1) {
				Grid<Real>* dest = mpPlParams->getGridReal(mDest);
				Grid<Real>* src = mpPlParams->getGridReal(mSrc);
				goCopyGrid<Real>(dest, src);
			} else if(mpPlParams->haveGridVec3(mSrc) && 1) {
				Grid<Vec3>* dest = mpPlParams->getGridVec3(mDest);
				Grid<Vec3>* src = mpPlParams->getGridVec3(mSrc);
				goCopyGrid<Vec3>(dest, src);
			}
			return true;
		};

	protected:
		std::string mSrc, mDest;
};

// copy one element of a vector to a scalar grid
template<class Scalar,int NUM>
class goGetVec3NormToScalar : public GridOpBase {
	public:
		goGetVec3NormToScalar(Grid<Scalar> *gDst, Grid<Vec3> *gSrc) : GridOpBase() { 
			mpDst = gDst;
			mpSrc = gSrc;
			applyOperatorToGridsWithoutFlags(this, mpDst);
		};
		~goGetVec3NormToScalar() { }; 
		void resetVariables() { };
		void buildCallList() {
			gaDst.gridAccInit(mpDst, AM_WRITE, gaCalls); 
			gaSrc.gridAccInit(mpSrc, AM_READ, gaCalls); 
			setNullFlagPointer();
		};
		inline void operator() (int i, int j, int k) { 
			gaDst.write(i,j,k) = norm(gaSrc(i,j,k));
		};
		void reduce(goGetVec3NormToScalar &op) { };
	protected:
		Grid<Scalar> *mpDst;
		Grid<Vec3>   *mpSrc;
		GridAccessor<Scalar,0> gaDst;
		GridAccessor<Vec3,0> gaSrc;
}; // goGetVec3NormToScalar */

class spluginGetNorm : public SolverPlugin {
	public:
		spluginGetNorm() : SolverPlugin(),
   				mSrc("-unnamed-"), mDest("-unnamed2-"),mode(0) { };
		~spluginGetNorm() { };

		virtual bool parseParams(const ParamSet& params) {
			mSrc = params.FindOneString("src", mSrc );
			mDest = params.FindOneString("dest", mDest );
			return true;
		};
		virtual bool initPlugin() { return true; };

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("spluginGetNorm","step "<<dt);
			Grid<int>* flags = mpPlParams->getGridInt("flags");

			Grid<Vec3>* src = mpPlParams->getGridVec3(mSrc);
			Grid<Real>* dest = mpPlParams->getGridReal(mDest);
			goGetVec3NormToScalar<Real,0>(dest,src);
			return true;
		};

	protected:
		// grid names to swap
		std::string mSrc,mDest;
		int mode;
};

// Obtain time-averaged field
class spluginAverage : public SolverPlugin {
	public:
		spluginAverage() : SolverPlugin(), mSummed(0),mStrideCount(0) { };
		~spluginAverage() { };

		virtual bool parseParams(const ParamSet& params) {
			mGrid = params.FindOneString("gridname", "vel-curr");
			mSumGrid = params.FindOneString("sumgrid", "");
			mFrom = params.FindOneInt("from", 0);
			mStride = params.FindOneInt("stride",1);
			mPostQuit = params.FindOneInt("post-quit", 0) != 0;
			mFrames = params.FindOneInt("frames", 100);
			mStrideCount = mStride-1;
			return true;
		};
		virtual bool initPlugin() { return true; };

		virtual void finish() {
			if (mSummed == 0) return;
			Grid<int>* flags = mpPlParams->getGridInt("flags");
			if (mpPlParams->haveGridVec3(mGrid))
				goScale<Vec3>(mpPlParams->getGridVec3(mSumGrid), flags, 1./(Real)mSummed);
			else if (mpPlParams->haveGridReal(mGrid))
				goScale<Real>(mpPlParams->getGridReal(mSumGrid), flags, 1./(Real)mSummed);
			else
				errMsg("spluginAverage","can't find grid " << mGrid << " or wrong type");
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			if (mFrom > 0) {
				debMsg("spluginAverage","will start in " << mFrom << " frames");
				--mFrom;
			} else {
				if (mSummed < mFrames) {
					++mStrideCount;
					if (mStrideCount == mStride) {
						mStrideCount = 0;
						debMsg("spluginAverage", "taking snapshot -- " << mSummed+1 << " frames recorded so far" ); 
				
						Grid<int>* flags = mpPlParams->getGridInt("flags");
						if (mpPlParams->haveGridVec3(mGrid))
							goSum<Vec3>(mpPlParams->getGridVec3(mSumGrid), mpPlParams->getGridVec3(mGrid), flags);
						else if (mpPlParams->haveGridReal(mGrid))
							goSum<Real>(mpPlParams->getGridReal(mSumGrid), mpPlParams->getGridReal(mGrid), flags);
						else
							errMsg("spluginAverage","can't find grid " << mGrid << " or wrong type");
						mSummed++;			
						if (mSummed == mFrames && mPostQuit)
							mpPlParams->setQuit(true);
					}
				}
			}
			return true;
		};

	protected:
		// grid names to swap
		std::string mGrid, mSumGrid;
		int mSummed, mFrom, mFrames,mStride, mStrideCount;
		bool mPostQuit;
};

// Dummy plugin
class spluginDummy : public SolverPlugin {
	public:
		spluginDummy() : SolverPlugin() { };
		~spluginDummy() { };

		virtual bool parseParams(const ParamSet& params) { return true; }
		virtual bool initPlugin() { return true; };

		virtual bool performStep(Real dt) {	return true; };
};


//*****************************************************************************

// solver plugin implementation
SolverPlugin::SolverPlugin() :
	mpPlParams(NULL), mName("-"), mTimeStart(-1.), mTimeStop(1.0e10),
	mWasActive(false), mSingleTime(-1.), mInterval(-1.), mLastActive(0.)
{
};

SolverPlugin::~SolverPlugin() {
};

void SolverPlugin::finish() {
};

// get parameters, e.g. grid names
bool SolverPlugin::parseParams(const ParamSet& params) {
	mTimeStart = params.FindOneFloat("plugin-start", mTimeStart );
	mTimeStop = params.FindOneFloat("plugin-stop", mTimeStop );
	mSingleTime = params.FindOneFloat("plugin-single-time", mSingleTime );
	mInterval = params.FindOneFloat("plugin-interval", mInterval );
	if(mTimeStop<mTimeStart) {
		mTimeStop = mTimeStart = -1.;
	}

	if (mSingleTime>=0.) {
		debMsg("SolverPlugin::parseParams","Single time active "<<mSingleTime);
	} else {
		if(mInterval>=0.) {
			debMsg("SolverPlugin::parseParams","Interval active "<<mInterval); 
			mLastActive = -mInterval;
			if(mTimeStart!=-1.) {
				mLastActive = mTimeStart-mInterval;
			}
		}
	}

	return true;
} 

// check for start & end time of plugin
bool SolverPlugin::isCurrentlyActive(double t) {
	// all timings off, return true by default
	//if( mTimeStart<= -1. && mTimeStop>=1.0e10 && mSingleTime <= -1. ) { return true; }

	// else: single time?
	if(mSingleTime >= 0.) {
		if(t >= mSingleTime && !mWasActive) {
			mWasActive = true;
			return true;
		} else {
			return false;
		}
	}
	
	// else: times given & valid, check:
	if(t < mTimeStart || t > mTimeStop) {
		return false;
	}

	if(mInterval >= 0.) {
		if(t > mLastActive+mInterval) {
			mLastActive += mInterval;
			return true;
		} else {
			return false;
		}
	}

	return true;
}

//*****************************************************************************

bool swapGrids(SolverParams* params, std::string mGrid1, std::string mGrid2)
{
	if(0) debMsg("swapGrids","DEBUG: "<<params->printGridMaps() );

	// distinguish data types
	// replace pointers in maps, and then change
	// names in grid classes too!
	GridBase *g1=NULL, *g2=NULL;
	if(params->haveGridInt(mGrid1) && 1) {
		Grid<int>* tmp = params->getGridInt(mGrid1);
		Grid<int>* grid2 = params->getGridInt(mGrid2);
		params->setGridInt(mGrid1, grid2);
		params->setGridInt(mGrid2, tmp);
		g1 = tmp; g2 = grid2;

		//grid2->setName(mGrid1);
		//tmp->setName(mGrid2);
	} else if(params->haveGridReal(mGrid1) && 1) {
		Grid<Real>* tmp = params->getGridReal(mGrid1);
		Grid<Real>* grid2 = params->getGridReal(mGrid2);
		params->setGridReal(mGrid1, grid2);
		params->setGridReal(mGrid2, tmp);
		g1 = tmp; g2 = grid2;

		//grid2->setName(mGrid1);
		//tmp->setName(mGrid2);
	} else if(params->haveGridVec3(mGrid1) && 1) {
		Grid<Vec3>* tmp = params->getGridVec3(mGrid1);
		Grid<Vec3>* grid2 = params->getGridVec3(mGrid2);
		params->setGridVec3(mGrid1, grid2);
		params->setGridVec3(mGrid2, tmp);
		g1 = tmp; g2 = grid2;
	} else { 
		return false;
	}

	int tmp1 = g1->getDisplayFlags();
	g1->setDisplayFlags( g2->getDisplayFlags() );
	g2->setDisplayFlags( tmp1 );

	g2->setName(mGrid1);
	g1->setName(mGrid2);

	if(0) debMsg("swapGrids","DEBUG: "<<params->printGridMaps() );
	return true;
}



// create one of the standard hard coded plugins
SolverPlugin* MakeStandardPlugin(std::string name) {
	if(name.compare( "copy-grid")==0) {
		return new spluginCopyGrid;
	} else if(name.compare( "dummy")==0) {
		return new spluginDummy;
	} else if(name.compare( "set-noslip-bcs")==0) {
		return new spSetNoslipBcs;
	} else if(name.compare( "average" )==0) {
		return new spluginAverage;
	} else if(name.compare( string("get-vec3-norm") )==0) {
		return new spluginGetNorm;
	}	
	return NULL;
}


} // end namespace DDF 

