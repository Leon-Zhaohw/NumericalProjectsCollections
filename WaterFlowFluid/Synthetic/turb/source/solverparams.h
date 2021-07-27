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

#ifndef DDF_SOLVERPARAMS_H
#define DDF_SOLVERPARAMS_H

#include "grid.h"
class ParamSet;

namespace WAVELETNOISE { class WaveletNoiseField; }

namespace DDF { 

class GeomFile; 
class FluidSolver; 


class SolverParams {
	public:
		// cons
		SolverParams(const std::string& name="");
		// des
		~SolverParams();
		// check that globals have been inited, and are ok
		bool verifyInit(string caller);
		// initializer from parsed paramset
		bool initFromParamSet(const ParamSet &params);
		// adapt init from another solver (eg set timestep)
		bool adaptInit(const SolverParams* const otherParams);

		// solver name (from config)
		std::string mName;

		// time step (maximal value, might be smaller due to cflfac)
		Real mMaxTimestep;
		// cfl multiplier to change time step according to max velocitites
		Real mCflFactor;

		// size of timestep for level set tracker (mTimestep * factor)
		Real mTimestepLevelsetFactor;

		// cg pressure solver control: max iteration factor
		Real mCgMaxIterFac;
		// desired accuracy
		Real mCgAccuracy;

		//! time step for animation output (<0 for no output)
		Real mTimestepAnim;
		//! output name prefix
		std::string mAnimOutputFile;
		//! output for blender or other program?
		int mOutputProgram;

		//! step factor for fixed dt solver comparison
		int mStepFactor;
		
		// VP stuff
		Vec3 mU0;
		bool mHostVorticitySystem;

		// grids
		std::map<string, Grid<Real>* > mGridsReal;
		std::map<string, Grid<Vec3>* > mGridsVec3;
		std::map<string, Grid<int>* >  mGridsInt;

		//! access functions
		void setInited(bool set) { mInited=set; }
		bool isInited() const { return mInited; } 
		void setQuit(bool set) { mQuit=set; }
		bool getQuit() const { return mQuit; }
		void setFluidSolver(FluidSolver* set) { mpFluidSolver=set; }
		FluidSolver* getFluidSolver() const { return mpFluidSolver; }

		Real getDt() const       { return mTimestep; }
		Real getTimestep() const { return mTimestep; }
		void setTimestep(Real set) { mTimestep = set; }

		//! get/set grids in a secure(er) way
		Grid<int>*  getGridInt(std::string name); 
		Grid<Real>* getGridReal(std::string name);
		Grid<Vec3>* getGridVec3(std::string name);
		void setGridInt(std::string name, Grid<int>* set); 
		void setGridReal(std::string name, Grid<Real>* set);
		void setGridVec3(std::string name, Grid<Vec3>* set);
		void getGridType(std::string name, Grid<int>* &set)  { set = getGridInt(name); };
		void getGridType(std::string name, Grid<Real>* &set) { set = getGridReal(name); };
		void getGridType(std::string name, Grid<Vec3>* &set) { set = getGridVec3(name); };
		
		//! get noise field
		WAVELETNOISE::WaveletNoiseField* getNoiseField(std::string name);
		//! add noise field
		void addNoiseField(WAVELETNOISE::WaveletNoiseField* f, std::string name);
		void addToNoiseTime(Real dt);

		//! debug output
		std::string printGridMaps();
		std::string toString();
		//! check if gridname exists
		bool haveGridInt(std::string name);
		bool haveGridReal(std::string name);
		bool haveGridVec3(std::string name);
		bool haveGridAny(std::string name);

		void setGridSize(nVec3i set) { mGridSize=set; mDimMax=VMAX(set); mDeltaX=1./(Real)mDimMax; }
		nVec3i getGridSize() const { return mGridSize; }
		// compute spatial discretization, currently given by dimMax,
		// might differ in future
		inline Real getDeltaX() const  { return mDeltaX; }
		// max dim.
		inline int  getDimMax() const  { return mDimMax; }

		// simulation time
		inline Real getSimTime() const { return mSimTime; }
		inline void addToSimTime(Real add) { mSimTime += add; }

		// div. corr.
		inline Real getDivergenceCorrection() const { return mDivergenceCorrection; }
		inline void setDivergenceCorrection(Real set) { mDivergenceCorrection = set; }

		inline Real getMultiplier() { return mAdapGridMultiplier; }
		inline void setMultiplier(Real v) { mAdapGridMultiplier = v; }

	protected: 
		// global settings
		bool mInited;
		// quit?
		bool mQuit;
		// current time
		Real mSimTime;

		// time step (maximal value, might be smaller due to cflfac)
		Real mTimestep;

		// size of simulation grid
		nVec3i mGridSize;
		// size of a cell, spatial discr. dx, set by gridsize
		Real mDeltaX;
		// max. dimension
		int mDimMax;

		// adapt init 
		bool mAdaptInit;
		double mAdaptDtScale;
		// gridsize multiplier for adapted init
		Real mAdapGridMultiplier;

		//! noise field storage
		std::map<string, WAVELETNOISE::WaveletNoiseField* > mNoiseFields;

		//! pass divergence correction on to poisson solver, e.g. for volume correction
		Real mDivergenceCorrection;

		//! simplify access to "parent" fluid solver
		FluidSolver* mpFluidSolver;
};

}; // DDF

#endif // DDF_SOLVERPARAMS_H

