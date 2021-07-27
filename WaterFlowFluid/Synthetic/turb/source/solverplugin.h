/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Base class for Plugins
 *
 *****************************************************************************/

#ifndef DDF_SOLVERPLUGIN_H
#define DDF_SOLVERPLUGIN_H

#include "globals.h"
class ParamSet;

namespace DDF { 
class SolverParams;

// base functions for a plugin
// plugin handling implementation in solver.cpp, 
// base methods implemented in stdplugins.cpp
class SolverPlugin {
	public:
		// cons
		SolverPlugin();
		// des
		virtual ~SolverPlugin();

		// get parameters, e.g. grid names
		virtual bool parseParams(const ParamSet& params);
		// init plugin, return failure
		virtual bool initPlugin() = 0;
		// optional method, called after main loop is finished
		virtual void finish();

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) = 0;

		inline void setPluginParams(SolverParams* set) { mpPlParams=set; }

		inline void setName(std::string set) { mName=set; }
		inline std::string getName() { return mName; }

		// check for start & end time of plugin
		bool isCurrentlyActive(double t);

		static SolverPlugin* makePlugin(const std::string& name);

	protected: 
		// parameter storage
		SolverParams* mpPlParams;

		// store plugin name for debugging
		std::string mName;

		// debugging - enable plugin only in an interval
		Real mTimeStart, mTimeStop;

		bool mWasActive;
		Real mSingleTime;

		Real mInterval, mLastActive;

}; // SolverPlugin


// create one of the plugins, e.g., advection or init steps

// stdplugins.cpp
SolverPlugin* MakeStandardPlugin(std::string name);

// stdplugins.cpp
SolverPlugin* MakeAdvectionPlugin(std::string name);

// initplugins.cpp
SolverPlugin* MakeInitPlugin(std::string name);

// vortexplugins.cpp
SolverPlugin* MakeVortexPlugin(std::string name);

// fileio.cpp
SolverPlugin* MakeIoPlugin(std::string name);

// poissonsolvers.cpp
SolverPlugin* MakePoissionSolverPlugins(std::string name);

// animplugins.cpp
SolverPlugin* MakeAnimPlugin(std::string name);

// smokeplugins.cpp
SolverPlugin* MakeSmokePlugin(std::string name);

// helper function for plugins moving data from one grid to another
bool swapGrids(SolverParams* params, std::string mGrid1, std::string mGrid2);

//************************************
// invoke plugins
inline SolverPlugin* SolverPlugin::makePlugin(const std::string& name)
{
	SolverPlugin *plugin = MakeStandardPlugin(name);
	if(!plugin) plugin = MakeInitPlugin(name);
	if(!plugin) plugin = MakeAdvectionPlugin(name);
	if(!plugin) plugin = MakeIoPlugin(name);
	if(!plugin) plugin = MakePoissionSolverPlugins(name);
	if(!plugin) plugin = MakeVortexPlugin(name);
	if(!plugin) plugin = MakeAnimPlugin(name);
	if(!plugin) plugin = MakeSmokePlugin(name);

	if(plugin) plugin->setName(name);
	return plugin;
}
}; // DDF


#endif // DDF_SOLVERPLUGIN_H


