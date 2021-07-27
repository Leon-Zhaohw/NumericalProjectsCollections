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

#ifndef DDF_FLUIDSOLVER_H
#define DDF_FLUIDSOLVER_H

#include "globals.h"
#include "grid.h"
#include "solverparams.h"

namespace DDF { 
class GeomFile;
class SolverPlugin;
class VorticitySystem;

extern FluidSolver* ddfWorldFindSolver(const std::string solverName);

//! patch grid based cg solver
class FluidSolver {
	public:
		// constructor
		FluidSolver(const std::string& name="fluidSolverUnnamed");
		~FluidSolver();

		//! setup plugin stack
		void addPlugins(vector<SolverPlugin*>& plugins);
		//! setup init plugin stack
		void addInitPlugins(vector<SolverPlugin*>& plugins);
		//! setup end plugin stack
		void addEndPlugins(vector<SolverPlugin*>& plugins);

		//! init solver
		bool initFluid();
		//! init level set
		void initLevelSet(int eiksolver);
		//! advance by time t
		bool simulateFluid();
		bool advanceParticles(Real t, Grid<Vec3>* vels);
		//! run init plugins
		void runInitPlugins();
		
		//! run final plugins
		void finalize();

		// access functions
		inline Real getDt() const        { return mDt; }
		inline Real getSimTime() const   { return mpParams->getSimTime(); }
		inline nVec3i getDim() const     { return mDim; }
		inline VorticitySystem* getVorticitySys() { return mpVorticitySys; }
		
		inline Vec3 getGravity() const   { return mGravity; }
		inline void setGravity(Vec3 set) { mGravity=set; }
		inline Real getCellSize() const   { return mCellSize; }

		inline int getStepCounter() const { return mStepCounter; }
		inline FlagGrid*   getGridFlags()      { return getParams()->getGridInt("flags"); }
		inline Grid<Vec3>* getGridCurrVel()    { return getParams()->getGridVec3("vel-curr"); }
		inline Grid<Vec3>* getGridOtherVel()   { return getParams()->getGridVec3("vel-old"); }
		inline Grid<Real>* getGridPressure()   { return getParams()->getGridReal("pressure"); }

		inline int get2dKstart() const        { return mTwodKs; }
		inline int get2dKend()   const        { return mTwodKe; }
		
		inline void setName(std::string set) { mName = set; }
		inline std::string getName() { return mName; }

		void setParams(SolverParams *params) { mpParams = params; }
		SolverParams* getParams() { return mpParams; }

		bool isSolverInited() { return mSolverInited; }
		void setSolverInited(bool set) { mSolverInited=set; }

		int getNumInitSteps() { return  mInitPlugins.size(); }
		int getNumPluginSteps() { return  mPlugins.size(); }

		// return timing statistics string
		std::string getTimingStats(int sort); 

		// helper info functions
		int getFrameNum() { 
			return (int)(getSimTime() / (getParams()->getDt() * getParams()->getDeltaX())); }
		int getAniFrameNum() { 
			return (int)(getSimTime() / (getParams()->mTimestepAnim ) ); }

	protected:
		typedef struct {
			std::string name;
			vector<Real> data;
		} DataVector;

		// parameter set
		SolverParams* mpParams;

		// time step
		Real mDt;
		// dimensions
		nVec3i mDim;
		// size of a single cell
		Real mCellSize;
		// gravity force
		Vec3 mGravity;

		// flag grid
		FlagGrid*  mpFlags;

		VorticitySystem* mpVorticitySys;
		
		// some global info vars from init
		// 2d k region start & end
		int mTwodKs, mTwodKe;
		// for debugging steps
		int mStepCounter;

		// actual steps to be performed
		vector<SolverPlugin*> mPlugins;
		// steps for init
		vector<SolverPlugin*> mInitPlugins;
		// final steps
		vector<SolverPlugin*> mEndPlugins;
		// profiling of steps
		vector<double> mPluginTimes;

		// debug, store name from init
		std::string mName;
		bool mSolverInited;

	public:

		// interpolate velocity at a given position (note! velocity is face-centered, 
		// so "normal" grid->getInerpolated(..) isnt correct 
		static Vec3 interpolateVelocity(Vec3 pos, Grid<Vec3>* vels, const Real dx); 
		Vec3 interpolateVpVelocity(Vec3 pos, Grid<Vec3>* vels); 

		// set velocity using global accesses, and
		// respecting obstacle flags (only !fgIsObstacle cells are modified)
		void setGlobVelNoobs(int i,int j, int k, Vec3 &set, Grid<Vec3>* vel);
		
		// get velocity component interpolated at grid center
		inline Vec3 getCenteredVel(Grid<Vec3> *vel, int x, int y, int z) { 
			Vec3 v = vel->getGlobal(x,y,z);
			v[0] +=  vel->getGlobal(x+1,y  ,z  )[0];
			v[1] +=  vel->getGlobal(x  ,y+1,z  )[1];
			v[2] +=  vel->getGlobal(x  ,y  ,z+1)[2];
			v *= 0.5;
			return v;
		};

}; // FluidSolver


} // namespace DDF 

#endif // DDF_FLUIDSOLVER_H
