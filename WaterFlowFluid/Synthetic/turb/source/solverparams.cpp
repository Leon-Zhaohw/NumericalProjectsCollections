/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Class to store solver parameters
 *
 *****************************************************************************/

#include "solverparams.h"
#include "paramset.h"
#include "waveletnoise.h"
namespace DDF {

//*****************************************************************************
// default init
SolverParams::SolverParams(const std::string& name) :
			mName(name),
			mMaxTimestep(0.5), mCflFactor(1.),
			//mSurfaceTracker(0), // by default, no surface
			mTimestepLevelsetFactor(1.0),
			mCgMaxIterFac(1.), mCgAccuracy(5.0 * gVecEpsilon),
			mTimestepAnim(0.5), mAnimOutputFile("animout_"), mOutputProgram(0),			
			mStepFactor(1), 
			mGridsReal(), mGridsVec3(), mGridsInt(),
			mInited(false), mQuit(false), mSimTime(-1.), 
			mTimestep(0.5), 
			mGridSize(-1), mDeltaX(-1.), 
			mDimMax(-1), // updated in setSize
			mAdaptInit(false), mAdaptDtScale(1.),
			mAdapGridMultiplier(-1.),
			mNoiseFields(), 
			mDivergenceCorrection(0.),
			mpFluidSolver(NULL)
{
	mU0 = Vec3(0.);
}; 	

// adapt init from another solver (eg set timestep)
bool SolverParams::adaptInit(const SolverParams* const otherParams)
{ 
	// use mAdapGridMultiplier ?
	if(mAdapGridMultiplier > 0.) {
		// round size
		nVec3i s = otherParams->getGridSize();
		Vec3 ds = Vec3(s[0],s[1],s[2]);
	  	ds *= mAdapGridMultiplier;
		setGridSize( nVec3i( (int)ds[0], (int)ds[1], (int)ds[2] ));
	}

	mAdaptDtScale = 1. / (double)otherParams->mDimMax * (double)mDimMax;
	mTimestep = otherParams->mTimestep * mAdaptDtScale;
	debMsg("SolverParams","adaptInit, gridsize="<<mGridSize<<", dt="<<mTimestep );

	mAdaptInit = true;

	return true;
}
		 
//*****************************************************************************
// free grids etc.
SolverParams::~SolverParams() {
		for(std::map<string, Grid<Real>* >::iterator iter=mGridsReal.begin(); 
				iter!=mGridsReal.end(); iter++) {
			if( (*iter).second ) 
				delete (*iter).second;
		}
		mGridsReal.clear();
		for(std::map<string, Grid<int>* >::iterator iter=mGridsInt.begin(); 
				iter!=mGridsInt.end(); iter++) {
			if( (*iter).second ) 
				delete (*iter).second;
		}
		mGridsInt.clear();
		for(std::map<string, Grid<Vec3>* >::iterator iter=mGridsVec3.begin(); 
				iter!=mGridsVec3.end(); iter++) {
			if( (*iter).second ) 
				delete (*iter).second;
		}
		mGridsVec3.clear();
}

bool SolverParams::verifyInit(string caller) {
	if(!mInited) { 
		errMsg("SolverParams::verifyInit","Not inited at all! "<<caller);
		return false;
	}
	if(mGridSize[0]<=0 || 
	   mGridSize[1]<=0 || 
	   mGridSize[2]<=0 || 
	   mGridSize[0]>10000 || 
	   mGridSize[1]>10000 || 
	   mGridSize[2]>10000 || 
	   false ) { 
		errMsg("SolverParams::verifyInit","Invalid gridsize "<<mGridSize);
		return false;
	}
	return true;
}


// initialize from parsed paramset
bool SolverParams::initFromParamSet(const ParamSet &params) {
	// reset values, might be copied from previous settings
	mInited =false;
	mGridsReal.clear();
	mGridsInt.clear();
	mGridsVec3.clear();
	mpFluidSolver = NULL;

	// either init from vector or int
	DDF::Vec3 sizeVec = params.FindOneVector("gridsize", vec2R(mGridSize) );
	int sizeInt = params.FindOneInt("gridsize", -1);
	if(sizeInt>0) {
		setGridSize( DDF::nVec3i(sizeInt) );
	} else {
		setGridSize( vec2I(sizeVec) );
	}
	// use adapt init multiplier? -1 for off
	mAdapGridMultiplier = params.FindOneFloat("adap-grid-multiplier", mAdapGridMultiplier);

	mU0 = params.FindOneVector("u0", mU0);
	mHostVorticitySystem = params.FindOneInt("host-vorticity-system", 0) != 0;

	mTimestep = params.FindOneFloat("timestep", mTimestep );
	mMaxTimestep = mTimestep;
	mCflFactor = params.FindOneFloat("cflfactor", mCflFactor );
	mTimestepLevelsetFactor = params.FindOneFloat("timestep_levelset_factor", mTimestepLevelsetFactor );

	mCgMaxIterFac = params.FindOneFloat("cg-max-iter-fac", mCgMaxIterFac );
	mCgAccuracy = params.FindOneFloat("cg-accuracy", mCgAccuracy );

	mOutputProgram = params.FindOneInt("outputprogram", mOutputProgram);
	if(mOutputProgram == 1) debMsg("Params","Output program set to blender...");

	mTimestepAnim = params.FindOneFloat("timestep_anim", mTimestepAnim );
	if (mAdapGridMultiplier>1) mTimestepAnim*=mAdapGridMultiplier;

	mAnimOutputFile = params.FindOneString("animfile", mAnimOutputFile);
	mStepFactor = params.FindOneInt("step-factor", mStepFactor);
	
	setInited(true);
	// check sanity... (requires mInited = true)
	if(!verifyInit("ddfParseGlobals")) {
		setInited(false);
	}

	return isInited();
}

//*****************************************************************************
// get grids

Grid<int>*  SolverParams::getGridInt(std::string name) {
	if(mGridsInt.find(name) != mGridsInt.end() && mGridsInt[name]!=NULL) {
		return mGridsInt[name];
	} else {
		errFatal("SolverParams::getGridInt","'"<< name <<"' int grid required! "<<this->toString(), SIMWORLD_GRIDERROR);
	}
	return NULL;
}

Grid<Real>* SolverParams::getGridReal(std::string name) {
	if(mGridsReal.find(name) != mGridsReal.end() && mGridsReal[name]!=NULL) {
		return mGridsReal[name];
	} else {
		errFatal("SolverParams::getGridReal","'"<< name <<"' real grid required! "<<this->toString(), SIMWORLD_GRIDERROR);
	}
	return NULL;
}

Grid<Vec3>* SolverParams::getGridVec3(std::string name) {
	if(mGridsVec3.find(name) != mGridsVec3.end() && mGridsVec3[name]!=NULL) {
		return mGridsVec3[name];
	} else {
		errFatal("SolverParams::getGridVec3","'"<< name <<"' vec3 grid required! "<<this->toString(), SIMWORLD_GRIDERROR);
	}
	return NULL;
}


void SolverParams::setGridInt(std::string name, Grid<int>* set) {
	mGridsInt[name] = set;
}
void SolverParams::setGridReal(std::string name, Grid<Real>* set) {
	mGridsReal[name] = set;
} 
void SolverParams::setGridVec3(std::string name, Grid<Vec3>* set) {
	mGridsVec3[name] = set;
} 

bool SolverParams::haveGridInt(std::string name) {
	if(mGridsInt.find(name) != mGridsInt.end() && mGridsInt[name]!=NULL) {
		return true;
	}
	return false;
	//return ( mGridsInt[name] != NULL);
}
bool SolverParams::haveGridReal(std::string name) {
	if(mGridsReal.find(name) != mGridsReal.end() && mGridsReal[name]!=NULL) {
		return true;
	}
	return false;
	//return ( mGridsReal[name] != NULL);
} 
bool SolverParams::haveGridVec3(std::string name) {
	if(mGridsVec3.find(name) != mGridsVec3.end() && mGridsVec3[name]!=NULL) {
		return true;
	}
	return false;
	//return ( mGridsVec3[name] != NULL);
} 
bool SolverParams::haveGridAny(std::string name) {
	return ( haveGridInt(name) ||
			 haveGridReal(name) ||
			 haveGridVec3(name) );
} 

//*****************************************************************************
// noise fields

void SolverParams::addNoiseField(WAVELETNOISE::WaveletNoiseField* f, std::string name) 
{ 
	mNoiseFields[name] = f; 
	Real dt = mNoiseFields[name]->getTimeAnim();
	// rescaling not necessary anymore
	//debMsg("ADAP"," dt="<<dt<<" scale="<<mAdaptDtScale );
	mNoiseFields[name]->setTimeAnim( dt ); 
}

WAVELETNOISE::WaveletNoiseField* SolverParams::getNoiseField(std::string name) {
	if(mNoiseFields.find(name) != mNoiseFields.end() && mNoiseFields[name]!=NULL) {
		return mNoiseFields[name];
	} else {
		errFatal("SolverParams::getNoiseField","'"<< name <<"' noise field required! "<<this->toString(), SIMWORLD_GRIDERROR);
	}
	return NULL;
}


// increase noise times
void SolverParams::addToNoiseTime(Real dt) {

	// scale by size because time anim is essentially
	// a diagonal movement...
	const Real dtFac = 1./(Real)mDimMax;

	for(std::map<std::string, WAVELETNOISE::WaveletNoiseField*>::iterator 
			iter=mNoiseFields.begin(); iter != mNoiseFields.end(); iter++) {
		WAVELETNOISE::WaveletNoiseField* pntval = (WAVELETNOISE::WaveletNoiseField*)( (*iter).second );

		//debMsg("SolverParams","addToNoiseTime, increased by "<<dt); 
		pntval->increaseTime(dt * dtFac);
		//pntval->increaseTime(1.);
	}
}

//*****************************************************************************
// helper

//void helperPrintMap(std::map<std::string, T*>& pmap, std::string caller) {
std::string SolverParams::printGridMaps() {
	std::ostringstream out;
	out << helperPrintMap(mGridsInt,"SolverParams::printGridMaps - int");
	out << helperPrintMap(mGridsReal,"SolverParams::printGridMaps - real");
	out << helperPrintMap(mGridsVec3,"SolverParams::printGridMaps - vec3");
	return out.str();
}

std::string SolverParams::toString() {
	std::ostringstream out;
	out <<	"SolverParams: gridsize="<<mGridSize<<" dt="<<mTimestep<<" cflfac="<<mCflFactor ;
	return out.str();
}

} // DDF

