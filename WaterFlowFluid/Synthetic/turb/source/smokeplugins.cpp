/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Plugins for seeding smoke
 *
 *****************************************************************************/

#include "fluidsolver.h"
#include "solverplugin.h"
#include "paramset.h"

#include "waveletnoise.h"
using WAVELETNOISE::WaveletNoiseField;


namespace DDF { 


//*****************************************************************************
// helper plugin

// optionally, write out inflow values for file export
class fsInitDensityInflow : public GridOpBase {
	public:
		fsInitDensityInflow(FlagGrid *flags, Grid<Real> *dst , Grid<Real> *writeGrid, Real target, int flag,
			WaveletNoiseField* noise = NULL) :
				GridOpBase(), mpDst(dst), mpWriteGrid(writeGrid),
  				mTargetValue(target), mFlag(flag), mpNoise(noise) { 
			mpFlags = flags;
			applyOperatorToGrids(this);
		};
		~fsInitDensityInflow() { }; 
		void resetVariables() { };
		void reduce(fsInitDensityInflow &op) { };
		void buildCallList() {
			gaDst.gridAccInit(mpDst, AM_WRITE, gaCalls); 
			if(mpWriteGrid) gaWriteGrid.gridAccInit(mpWriteGrid, AM_WRITE, gaCalls); 
			setFlags(mpFlags);
		};

		// add forces and update empty cells free surface boundaries
		inline void operator() (int i, int j, int k) { 
			if(mpWriteGrid) gaWriteGrid.write(i,j,k) = 0.;

			// safety - always write
			if( (getFlagAcc()(i,j,k) & mFlag) == 0) return;
			Real& v = gaDst.write(i,j,k);

			if(!mpNoise) {
				// constant
				if(v<mTargetValue) v = mTargetValue;
			} else {
				const Real targ = mpNoise->evaluate( Vec3(i, j, k) ) * mTargetValue;
				if(v<targ) v = targ;
			}

			//debMsg("at"," "<<PRINT_IJK<<" flag ="<< getFlagAcc()(i,j,k)<<" smoke "<< gaDst.write(i,j,k));
			if(mpWriteGrid) gaWriteGrid.write(i,j,k) = v;
		};

	protected:
		Grid<Real> *mpDst, *mpWriteGrid;
		GridAccessor<Real,0> gaDst;
		GridAccessor<Real,0> gaWriteGrid;
		Real mTargetValue;
		int mFlag;
		// optional noise eval
		WaveletNoiseField* mpNoise;
}; // fsInitDensityInflow */

class fsInflowVecApply : public GridOpBase {
	public:
		fsInflowVecApply(FlagGrid *flags, Grid<Vec3> *dst , Vec3 target, 
			WaveletNoiseField* noise = NULL, bool set=false, Real noiseThresh=0.) :
				GridOpBase(), mpDst(dst), 
  				mTargetValue(target), mpNoise(noise), mSet(set), mNoiseThreshold(noiseThresh) { 
			debMsg("fsInflowVecApply","noise thresh "<<mNoiseThreshold);
			if(mNoiseThreshold>=0. && !mpNoise) {
				errFatal("fsInflowVecApply","Noise treshhold needs noise field!",SIMWORLD_PLUGINERROR);
			}
			mpFlags = flags;
			applyOperatorToGrids(this);
		};
		~fsInflowVecApply() { }; 
		void resetVariables() { };
		void buildCallList() {
			gaDst.gridAccInit(mpDst, AM_WRITE, gaCalls); 
			setFlags(mpFlags);
		};

		// add forces and update empty cells free surface boundaries
		inline void operator() (int i, int j, int k) { 
			if(! fgIsInflow(getFlagAcc()(i,j,k)) ) return;
			Vec3& v = gaDst.write(i,j,k);

			if(!mpNoise) {
				// constant
				if(!mSet) v += mTargetValue;
				else v = mTargetValue;
			} else if(mNoiseThreshold<=0.) {
				// use noise...

				if(!mSet) { 
					// add 
					v[0] += mpNoise->evaluate( Vec3(i, j, k) ) * mTargetValue[0];
					v[1] += mpNoise->evaluate( Vec3(i+100., j, k) ) * mTargetValue[0];
					v[2] += mpNoise->evaluate( Vec3(i, j+200., k) ) * mTargetValue[0]; 
				} else { 
					// overwrite old
					v[0] = mpNoise->evaluate( Vec3(i, j, k) ) * mTargetValue[0];
					v[1] = mpNoise->evaluate( Vec3(i+100., j, k) ) * mTargetValue[0];
					v[2] = mpNoise->evaluate( Vec3(i, j+200., k) ) * mTargetValue[0]; 
				}
			} else {
				// use noise threshold...

				if(!mSet) { 
					// add 
					errFatal("fsInflowVecApply","Noise treshhold only works with set=1",SIMWORLD_PLUGINERROR);
				} else { 
					// overwrite old
					if( mpNoise->evaluate( Vec3(i, j, k) ) > mNoiseThreshold ) {
						v = mTargetValue;
						//debMsg("vec"," "<<PRINT_IJK);
					}
				}
			}
		};
		void reduce(fsInflowVecApply &op) { };

	protected:
		Grid<Vec3> *mpDst;
		GridAccessor<Vec3,0> gaDst;
		Vec3 mTargetValue;
		WaveletNoiseField* mpNoise;
		bool mSet;
		// dont set noise values, but set constant value if noise is smaller than threshold
		Real mNoiseThreshold;
}; // fsInflowVecApply */

class spluginInitDensityInflow : public SolverPlugin {
	public:
		spluginInitDensityInflow() : SolverPlugin(),
   			mDensName("-unnamed-"),  mWriteName(""),
				mNoiseName(""),
				mAnimOutCounter(0), mTargetValue(1.), mTargetVec(0.), mSet(false),
  				mNoiseThresh(0.)	{
			debMsg("spluginInitDensityInflow","cons");
		};
		~spluginInitDensityInflow() {
			debMsg("spluginInitDensityInflow","des");
		};

		virtual bool parseParams(const ParamSet& params) {
			debMsg("spluginInitDensityInflow","parse");
			mDensName = params.FindOneString("density", mDensName );
			mFlag = params.FindOneInt("flag", FINFLOW);
			mVecName = params.FindOneString("vecgridname", mVecName );
			mWriteName = params.FindOneString("write-grid", mWriteName );
			mTargetValue = params.FindOneFloat("target-value", mTargetValue );

			mSet = 0 < params.FindOneInt("set", mSet );
			
			// preinit vec in case only real value is given
			mTargetVec = Vec3(mTargetValue);
			mTargetVec = params.FindOneVector("target-vec", mTargetVec );

			mNoiseName = params.FindOneString("noise", mNoiseName );
			mNoiseThresh = params.FindOneFloat("noise-threshold", mNoiseThresh );
			return true;
		};
		virtual bool initPlugin() {
			debMsg("spluginInitDensityInflow","init");
			return true;
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("spluginInitDensityInflow"," dt="<<dt<<" dest:"<<mDensName<<"/"<<mVecName );
			//Grid<Real>* dens = mpPlParams->getGridReal(mDensName);
			Grid<Real>* writeGrid = NULL;
			if(mWriteName.length()>0) {
				writeGrid = mpPlParams->getGridReal(mWriteName);
			}

			WaveletNoiseField* noise = NULL;
			if(mNoiseName.length()>0) 
				noise = mpPlParams->getNoiseField(mNoiseName);

			Grid<Real>* dens = NULL;
			if(mpPlParams->haveGridReal(mDensName))    dens    = mpPlParams->getGridReal(mDensName);
			Grid<Vec3>* srcVecGrid = NULL;
			if(mpPlParams->haveGridVec3(mVecName)) srcVecGrid = mpPlParams->getGridVec3(mVecName);

			GridBase* baseGrid = NULL;
			if(dens) {
				baseGrid = dens;
			} else if(srcVecGrid) {
				baseGrid = srcVecGrid;
			} else {
				errFatal("spluginInitDensityInflow","Provide either 'density'='"<<mDensName<<"', or 'vecgridname'='"<<mVecName<<"' ",SIMWORLD_PLUGINERROR);
			}

			if(dens) {
				fsInitDensityInflow(mpPlParams->getFluidSolver()->getGridFlags(), dens, 
							writeGrid, mTargetValue, mFlag,noise);				

				/* writing out grids should now be done with the "dump-grip-scalar" plugin */
			} 
			if(srcVecGrid) {
				// sanity checks
				if(dens) {
					errFatal("spluginInitDensityInflow","Error, can't do real-grid and vec-grid handling at same time!", SIMWORLD_PLUGINERROR);
				}
				
				// apply...
				fsInflowVecApply(mpPlParams->getFluidSolver()->getGridFlags(), srcVecGrid, 
						mTargetVec, noise, mSet, mNoiseThresh);
			}

			return true;
		};

	protected:
		// grid names to swap
		std::string mDensName, mWriteName, mVecName;
		std::string mNoiseName;
		int mAnimOutCounter, mFlag;
		Real mTargetValue;
		Vec3 mTargetVec;
		bool mSet;
		int mDoLevelset;
		Real mNoiseThresh;
};


//*****************************************************************************

SolverPlugin* MakeSmokePlugin(std::string name) {

	if(name.compare( string("init-density-inflow") )==0) {
		return new spluginInitDensityInflow;
	}
	return NULL;
}


} // end namespace DDF 


