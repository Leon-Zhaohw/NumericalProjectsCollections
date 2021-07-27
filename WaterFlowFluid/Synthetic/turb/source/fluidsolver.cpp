/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Class encapsuling the solver
 *
 *****************************************************************************/

// lsdebug
#include "globals.h"
#include "vectorbase.h"

#include "solverparams.h"
#include "solverplugin.h"
#include "paramset.h"

#include "vortexpart.h"
//#include "boundbox.h"
#include "randomstream.h"
#include <fstream>
#include <sstream>

// define global solver list
map<string, DDF::FluidSolver*> gFluidSolvers;
DDF::FluidSolver *gpFsolver = NULL;

namespace DDF { 
// search for fluid solve with given name
FluidSolver* ddfWorldFindSolver(const std::string solverName) {
	if(gFluidSolvers.find(solverName) != gFluidSolvers.end() ) {
		return gFluidSolvers[solverName];
	}

	errMsg("ddfWorldFindSolver","Solver '"<<solverName<<"' not found!");
	std::cout<<helperPrintMap(gFluidSolvers,"ddfWorldFindSolver")<<"\n";
	return NULL;
}

// ctor
FluidSolver::FluidSolver(const std::string& name) : 
	mpParams(NULL),
	mDt(0.), mDim(0), mGravity(0.),
	mpFlags(NULL), 
	mStepCounter(0),
	mName(name), mSolverInited(false)	
{
	mpVorticitySys = NULL;
	// debug test inits
	mGravity = Vec3( 0., 0., 0.);
};

// dtor
FluidSolver::~FluidSolver() {
	// free plugins
	for (int i=0; i< (int)mPlugins.size(); i++)
		delete mPlugins[i];
	for (int i=0; i< (int)mInitPlugins.size(); i++)
		delete mInitPlugins[i];
	for (int i=0; i< (int)mEndPlugins.size(); i++)
		delete mEndPlugins[i];

	mPlugins.clear();
	mEndPlugins.clear();
	mInitPlugins.clear();
	mPluginTimes.clear();
	delete mpParams;
	if (mpVorticitySys != NULL) delete mpVorticitySys;

	//note, handle free'ing other pointers...
};

// *****************************************************************************
//! init solver
bool FluidSolver::initFluid() {
	if(!mpParams) {
		errMsg("FluidSolver::initFluid","Fatal - no params given...");
		return false;
	}
	mpParams->setFluidSolver( this );

	// restrict z-axis in 2D
	if (gDim==2) {
		nVec3i nn = mpParams->getGridSize();
		nn[2] = gPatchSize;
		mpParams->setGridSize(nn);
	}

	const nVec3i n = mpParams->getGridSize();
	// set time to zero
	mpParams->addToSimTime( -mpParams->getSimTime() );
	mpFlags = mpParams->getGridInt("flags");

	// init from unittests
	const int nMax = mpParams->getDimMax();

	mTwodKs = 0;
	mTwodKe = n[2];
	if (gDim==2) {
		mTwodKs = gPatchSize/2;
		mTwodKe = mTwodKs+1;
	}
	const int ks = mTwodKs;
	const int ke = mTwodKe;

	// setup boundaries
	for (int k=ks; k< ke; k++) 
		for (int j=0; j< n[1]-0; j++) 
			for (int i=0; i< n[0]-0; i++) {
				int set = FFLUID;
				if( k<=0 || j<=0 || i<=0 ||
					 (k>=n[2]-1) || (j>=n[1]-1) || (i>=n[0]-1)) { set = FOBSTACLE; }
				mpFlags->getGlobal(i,j,k) = set;
			}
	mDim = n;

	if (mpParams->mHostVorticitySystem)
		mpVorticitySys = new VorticitySystem(this);
	
	return true;
};

//! run init plugins once, similar to simulateFluid
void FluidSolver::runInitPlugins() {
	mDt = mpParams->getTimestep();
	DDF::myTime_t tStart = DDF::getTime();

	// run plugins
	for (int i=0; i< (int)mInitPlugins.size(); i++) {

		// execute
		DDF::myTime_t plStart = DDF::getTime();
		if (! mInitPlugins[i]->performStep(mDt) ) {
			errFatal("FluidSolver::runInitPlugins","Plugin step "<<i<<" failed!", SIMWORLD_PLUGINERROR);
		}
		DDF::myTime_t plEnd = DDF::getTime();
	} 

	// timing...
	DDF::myTime_t tEnd = DDF::getTime();
	debMsg("FluidSolver::runInitPlugins",mpParams->mName<<" t="<<mpParams->getSimTime() <<", dt="<<mDt<<", took "<< DDF::getTimeString(tEnd-tStart)<<" ");
};

//! run end plugins once, similar to simulateFluid
void FluidSolver::finalize() {
	mDt = mpParams->getTimestep();
	DDF::myTime_t tStart = DDF::getTime();

	// finalize main loop plugins
	for (int i=0; i< (int)mPlugins.size(); i++)
		mPlugins[i]->finish();

	// run end plugins
	for (int i=0; i< (int)mEndPlugins.size(); i++) {

		// execute
		DDF::myTime_t plStart = DDF::getTime();
		if (! mEndPlugins[i]->performStep(mDt) ) {
			errFatal("FluidSolver::runEndPlugins","Plugin step "<<i<<" failed!", SIMWORLD_PLUGINERROR);
		}
		DDF::myTime_t plEnd = DDF::getTime();
		//mPluginTimes[i] += (double)(plEnd-plStart);
	} 

	// timing...
	DDF::myTime_t tEnd = DDF::getTime();
	debMsg("FluidSolver::finalize",mpParams->mName<<" t="<<mpParams->getSimTime() <<", dt="<<mDt<<", took "<< DDF::getTimeString(tEnd-tStart)<<" ");
};

// setup plugin stack
void FluidSolver::addPlugins(vector<SolverPlugin*>& plugins) {
	debMsg("addPlugins","adding "<<plugins.size() ); 
	mPluginTimes.clear();

	// store plugin stack, and perform first init
	for (int i=0; i< (int)plugins.size(); i++) {
		mPlugins.push_back( plugins[i] );
		mPluginTimes.push_back( 0. );

		if(! mPlugins[i]->initPlugin() ) {
			errFatal("FluidSolver::addPlugins","Init "<<i<<" failed!", SIMWORLD_PLUGINERROR);
		}
	}
}

// setup init plugin stack
void FluidSolver::addInitPlugins(vector<SolverPlugin*>& plugins) {
	debMsg("addInitPlugins","adding "<<plugins.size() ); 
	mPluginTimes.clear();

	// store plugin stack, and perform first init
	for (int i=0; i< (int)plugins.size(); i++) {
		mInitPlugins.push_back( plugins[i] );

		if(! mInitPlugins[i]->initPlugin() ) {
			errFatal("FluidSolver::addInitPlugins","Init "<<i<<" failed!", SIMWORLD_PLUGINERROR);
		}
	}
}

// setup end plugin stack
void FluidSolver::addEndPlugins(vector<SolverPlugin*>& plugins) {
	debMsg("addInitPlugins","adding "<<plugins.size() ); 
	
	// store plugin stack, and perform first init
	for (int i=0; i< (int)plugins.size(); i++) {
		mEndPlugins.push_back( plugins[i] );

		if(! mEndPlugins[i]->initPlugin() ) {
			errFatal("FluidSolver::addEndPlugins","Init "<<i<<" failed!", SIMWORLD_PLUGINERROR);
		}
	}
}


//! advance by time t
bool FluidSolver::simulateFluid() {
	mDt = mpParams->getTimestep();
	DDF::myTime_t tStart = DDF::getTime();

	debMsg("\nsimulateFluid"," for '"<<this->getName()<<"' ");
	if(!isSolverInited()) {
		errFatal("simulateFluid","Not inited!? 'SolverEnd' statement missing...?", SIMWORLD_GENERICERROR);
	}

	if (mPlugins.size() == 0) {
		mpParams->setQuit(true);
		return false;
	}

	// run plugins
	for (int i=0; i< (int)mPlugins.size(); i++) {
		// plugin active? start&stop time
		if(! mPlugins[i]->isCurrentlyActive( mpParams->getSimTime() ) ) {
			continue;
		}

		// execute
		DDF::myTime_t plStart = DDF::getTime();
		if (! mPlugins[i]->performStep(mDt) ) {
			errFatal("FluidSolver::simulateFluid","Plugin step "<<i<<" failed!", SIMWORLD_PLUGINERROR);
		}
		DDF::myTime_t plEnd = DDF::getTime();
		mPluginTimes[i] += (double)(plEnd-plStart);
	}

	mpParams->addToSimTime( mDt * mpParams->getDeltaX() );
	mpParams->addToNoiseTime(mDt);
	mStepCounter++;
	
	// timing...
	DDF::myTime_t tEnd = DDF::getTime();
	//debMsg("SG-Fl","simtime="<<sgSimTime<<", Time = "<< DDF::getTimeString(tEnd-tStart) );
	debMsg("FluidSolver::simulateFluid",mpParams->mName<<" t="<<mpParams->getSimTime() <<", dt="<<mDt<<", took "<< DDF::getTimeString(tEnd-tStart)<<" #steps:"<<mStepCounter); 

	return true;
};

// velocity interpolation helpers
#		define PINIT_INTERPOL(offa,offb,offc) \
				Real srcpi = pos[0]+(offa); \
				Real srcpj = pos[1]+(offb); \
				Real srcpk = pos[2]+(offc); \
				int srci = (int)srcpi; \
			 	int	srcj = (int)srcpj; \
			 	int	srck = (int)srcpk; \
				CLAMP_TO_GRID(srci,srcj,srck,vels); \
				const float s1 = srcpi-(float)srci, s0 = 1.-s1; \
				const float t1 = srcpj-(float)srcj, t0 = 1.-t1; \
				const float f1 = srcpk-(float)srck, f0 = 1.-f1; 
		/* end init */

#		define VELACC vels->getGlobal
#		if DDF_DIMENSION==3
#		define PADV_INTERPOLATE(index) \
				vel[index] = f0*( \
					s0*(t0*VELACC(srci  ,srcj,srck)[index]  + t1*VELACC(srci  ,srcj+1,srck)[index] ) + \
					s1*(t0*VELACC(srci+1,srcj,srck)[index]  + t1*VELACC(srci+1,srcj+1,srck)[index] ) ) \
					+ f1 * ( \
					s0*(t0*VELACC(srci  ,srcj,srck+1)[index]  + t1*VELACC(srci  ,srcj+1,srck+1)[index] ) + \
					s1*(t0*VELACC(srci+1,srcj,srck+1)[index]  + t1*VELACC(srci+1,srcj+1,srck+1)[index] ) ) ;
		/* end interp */
#		else
		// NOTE! srck not used for 2d!
#		define PADV_INTERPOLATE(index) \
				vel[index] = f0*( \
					s0*(t0*VELACC(srci  ,srcj,srck)[index]  + t1*VELACC(srci  ,srcj+1,srck)[index] ) + \
					s1*(t0*VELACC(srci+1,srcj,srck)[index]  + t1*VELACC(srci+1,srcj+1,srck)[index] ) ) ;
		/* end interp */
#		endif // DDF_DIMENSION==3

Vec3 FluidSolver::interpolateVelocity(Vec3 pos, Grid<Vec3>* vels, const Real dx) {
	Vec3 vel;
	pos /=  dx; // mpParams->getDeltaX();
	if (gDim==2) pos[2] = gPatchSize/2;

	Real koff= -0.5;
	if (gDim==2) koff=0.;

	{ PINIT_INTERPOL(0., -0.5, koff);
	  PADV_INTERPOLATE(0); }

	{ PINIT_INTERPOL(-0.5, 0., koff);
	  PADV_INTERPOLATE(1); }

#	if DDF_DIMENSION==3
	{ PINIT_INTERPOL(-0.5, -0.5, 0.);
	  PADV_INTERPOLATE(2); }
#	else
	vel[2] = 0.;
#	endif

	//debMsg("interpolateVelocity","At "<<pos<<" = "<<vel ); // DEBUG
	return vel;
} // interpolateVelocity


Vec3 FluidSolver::interpolateVpVelocity(Vec3 pos, Grid<Vec3>* vels) {
	Vec3 vel;
	if (gDim==2) pos[2] = gPatchSize/2;

	{ PINIT_INTERPOL(0.5, 0., 0.);
	  PADV_INTERPOLATE(0); }

	{ PINIT_INTERPOL(0., 0.5, 0.);
	  PADV_INTERPOLATE(1); }

#	if DDF_DIMENSION==3
	{ PINIT_INTERPOL(0., 0., 0.5);
	  PADV_INTERPOLATE(2); }
#	else
	vel[2] = 0.;
#	endif

	return vel;
}

#undef PIOFF 
#undef PINIT_INTERPOL 
#undef PADV_INTERPOLATE
#undef VELACC

// set velocity, check flags
// similar to force add (but slower!)
void FluidSolver::setGlobVelNoobs(int i,int j, int k, Vec3 &set, Grid<Vec3>* velGrid) {
	if ( fgIsObstacle(mpFlags->getGlobal(i,j,k)) ) return;

	// new test, only set empty ones
	if (!fgIsEmpty(mpFlags->getGlobal(i,j,k)) ) return;
	Vec3& velE = velGrid->getGlobal(i,j,k);
	if (fgIsEmpty(mpFlags->getGlobal(i-1,j,k)) ) {
		velE[0] = set[0]; 
	} 
	if (fgIsEmpty(mpFlags->getGlobal(i,j-1,k)) ) {
		velE[1] = set[1]; 
	} 
	if (fgIsEmpty(mpFlags->getGlobal(i,j,k-1)) ) {
		velE[2] = set[2]; 
	}
	return;


	//Vec3& vel = mpCurrVel->getGlobal(i,j,k);
	Vec3& vel = velGrid->getGlobal(i,j,k);
	if (!fgIsObstacle(mpFlags->getGlobal(i-1,j,k)) ) {
		vel[0] = set[0]; 
	}

	if (!fgIsObstacle(mpFlags->getGlobal(i,j-1,k)) ) {
		vel[1] = set[1]; 
	}

	if (!fgIsObstacle(mpFlags->getGlobal(i,j,k-1)) ) {
		vel[2] = set[2]; 
	}
}

// return timing statistics string
std::string FluidSolver::getTimingStats(int sort) {
	std::ostringstream ret;
	ret << "FluidSolver::getTimingStats, #steps="<< (mStepCounter+1) <<"\n";
	const double istep = 1  / (double)(mStepCounter+1.);

	double total = 0.;
	double max = -1.;
	for(int i=0; i<(int)mPluginTimes.size(); i++) {
		total += mPluginTimes[i]*istep;
		if(mPluginTimes[i]*istep > max) max = mPluginTimes[i]*istep;
	}

	if(sort==0) {	
		for(int i=0; i<(int)mPluginTimes.size(); i++) {
			ret << i <<", "<< mPlugins[i]->getName() << " = "<< getTimeString( (myTime_t)(mPluginTimes[i]*istep) ) 
				<<"; "  << (mPluginTimes[i]*istep / total *100.) <<"%"	 // per cent
				<<"\n";
		}
	} else {
		std::vector<bool> done;
		int printed = 0;
		done.resize( mPluginTimes.size() );
		for(int i=0; i<(int)mPluginTimes.size(); i++) {
			done[i] = false;
		}
		// sort the slow way...
		while(printed< (int)mPluginTimes.size()) {
			double nextMax = -1.;
			for(int i=0; i<(int)mPluginTimes.size(); i++) {
				if(mPluginTimes[i]*istep==max) {
					ret << i <<", "<< mPlugins[i]->getName() << " = "<< getTimeString( (myTime_t)(mPluginTimes[i]*istep) ) 
						<<"; "  << (mPluginTimes[i]*istep / total *100.) <<"%"	 // per cent
						<<"\n";
					done[i] = true;
				}
				if(!done[i] && mPluginTimes[i]*istep > nextMax) nextMax = mPluginTimes[i]*istep;
			}
			max = nextMax;
			printed++;
		}
	}
	ret <<"Total time = "  << getTimeString( (myTime_t)total ) <<"\n";

	return ret.str();
}

}// end namespace DDF 
