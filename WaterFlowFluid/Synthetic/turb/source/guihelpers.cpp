/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Helper function for the GUI
 *
 *****************************************************************************/

#include "globals.h"
#include "vectorbase.h"
#include "solverparams.h"
#include "paramset.h"

#include "fluidsolver.h"

extern std::string gPhiRefName;
extern std::string gShowVelgridName;
extern std::string gShowRealgridName;

// solver object from main
extern DDF::FluidSolver *gpGlutFsolver;

// from parser.cpp
extern std::map<std::string, DDF::FluidSolver*> gFluidSolvers;

// values from gui
extern std::map<std::string,float> gScaleRealDisplay;

namespace DDF { 

extern Real gRealgridMin, gRealgridmMax;
extern Vec3 gVecgridMin, gVecgridMax;

#if DDF_GLUTGUI==1

/******************************************************************************/
// switch solvers & grids


template<class Scalar>
void glutgSelectNextGrid(std::map<std::string, Grid<Scalar>*> &gridMap, std::string &gname) {
	bool debSwitch = false;
	SolverParams* prms = gpGlutFsolver->getParams();
	Grid<Scalar> *retgrid = NULL;

	// clean up list, might contain entries with NULL pointers...
	if(1) {
		typename std::map<std::string, Grid<Scalar>*>::iterator iter = gridMap.begin(); 
		for ( ; iter!= gridMap.end(); iter++) {
			if(!(*iter).second) {
				gridMap.erase(iter);
				iter--;
				continue;
			}
			if (debSwitch) debMsg("deb disp","'"<<(*iter).second->getName()<<"' "); 
		}

		// sanity check, swap messes up names?
		iter = gridMap.begin(); 
		for ( ; iter!= gridMap.end(); iter++) {
			string name = (*iter).second->getName();
			prms->getGridType(name, retgrid); 
			//debMsg("glutGui","Scalar-grid '"<<retgrid->getName()<<"', '"<<name<<"' selected for display");

			if(retgrid->getName().length() != name.length() ) {
				errFatal("glutgSelectNextGrid","Invalid grid names!? "<<retgrid->getName()<<" vs. "<<name, SIMWORLD_GENERICERROR );
			}
			retgrid = NULL;
		}
	}

	typename std::map<std::string, Grid<Scalar>*>::iterator iter = gridMap.begin();
	if (!(gname.length()>0)) {
		// select first
		if (debSwitch) debMsg("glutgSelectNextGrid","Not inited, getting first");

		while((*iter).second->getDisplayFlags() & 1) {
			if (debSwitch) debMsg("glutgSelectNextGrid","1 Skipping "<<" "<<(*iter).second->getName() );
			iter++;
			if (iter == gridMap.end()) { iter = gridMap.begin(); }
		}

	} else { 
		// iterate through map
		//typename std::map<std::string, Grid<Scalar>*>::iterator iter = gridMap.begin();
		int cnt=0;
		for (; iter!= gridMap.end(); iter++) {
			if (debSwitch) debMsg("glutgSelectNextGrid","Searching "<<(*iter).second->getName() );
			if(gname == (*iter).second->getName()) break;
			cnt++;
		}
		if (iter==gridMap.end()) {
			if (debSwitch) debMsg("glutgSelectNextGrid","Not found, getting first "<<cnt);
			iter = gridMap.begin();
		} else {
			if (debSwitch) debMsg("glutgSelectNextGrid","Found, getting next "<<cnt<<" "<<(*iter).second->getName() );
			iter++;

			if (iter==gridMap.end()) {
				if (debSwitch) debMsg("glutgSelectNextGrid","Last, getting first "<<cnt );
				iter = gridMap.begin();
			}
		}

		// skip hidden grids... (getDisplayFlags&1)
		while((*iter).second->getDisplayFlags() & 1) {
			if (debSwitch) debMsg("glutgSelectNextGrid","Skipping "<<" "<<(*iter).second->getName() );

			iter++;
			if (iter == gridMap.end()) { iter = gridMap.begin(); }

			if (debSwitch) debMsg("glutgSelectNextGrid","Now "<<" "<<(*iter).second->getName() );
		}
	}

	gname = (*iter).second->getName();
	prms->getGridType(gname, retgrid);

	debMsg("glutGui","Scalar-grid '"<<retgrid->getName()<<"', '"<<gname<<"' selected for display");
}

void glutgSelectNextFloatArray() {
	SolverParams* prms = gpGlutFsolver->getParams();
	glutgSelectNextGrid<Real>( prms->mGridsReal, gShowRealgridName );
}

// init display, called eg after glutgSelectNextFloatArray
void initFloatArrayDisplay() {
	// new grid? set standard scale
	if( gScaleRealDisplay[gShowRealgridName] <= 0. ) gScaleRealDisplay[gShowRealgridName] = 1.;

	// show min/max values for selected grid
	if(1) {
		Grid<Real> *realGrid = NULL; 
		realGrid = gpGlutFsolver->getParams()->getGridReal(gShowRealgridName); // DEBrealGrid
		goFindMinMax<Real> gomm = goFindMinMax<Real>( gpGlutFsolver->getGridFlags(), realGrid );

		gRealgridMin = gomm.mMinVal;
		gRealgridmMax = gomm.mMaxVal;
		debMsg("MINMAX-REAL","min="<<gRealgridMin<<" at "<<gomm.mMinPos<<" max="<<gRealgridmMax<<" at "<<gomm.mMaxPos);
	} 
}

void glutgSelectNextVec3Array() {
	SolverParams* prms = gpGlutFsolver->getParams();
	glutgSelectNextGrid<Vec3>( prms->mGridsVec3, gShowVelgridName );
}

void glutgSelectNextFloatArrayO() {
	//if (gFluidSolvers.find(solverName) != gFluidSolvers.end() ) { gpGlutFsolver = gFluidSolvers[solverName]; }
	bool debSwitch = false;
	SolverParams* prms = gpGlutFsolver->getParams();
	Grid<Real> *grid = NULL;

	// clean up list, might contain entries with NULL pointers...
	for ( std::map<std::string, Grid<Real>*>::iterator iter = prms->mGridsReal.begin(); iter!= prms->mGridsReal.end(); iter++) {
		if(!(*iter).second) {
			prms->mGridsReal.erase(iter);
			iter--;
			continue;
		}
		//if (debSwitch) debMsg("deb disp","'"<<(*iter).second->getName()<<"' "); 
	}

	if (!(gShowRealgridName.length()>0)) {
		// select first
		if (debSwitch) debMsg("glutgSelectNextFloatArray","Not inited, getting first");
		gShowRealgridName = (*prms->mGridsReal.begin()).second->getName();
		grid = prms->getGridReal(gShowRealgridName);

	} else { 
		// iterate through map
		std::map<std::string, Grid<Real>*>::iterator iter = prms->mGridsReal.begin();
		int cnt=0;
		for (; iter!= prms->mGridsReal.end(); iter++) {
			if (debSwitch) debMsg("glutgSelectNextFloatArray","Searching "<<(*iter).second->getName() );
			if(gShowRealgridName == (*iter).second->getName()) break;
			cnt++;
		}
		if (iter==prms->mGridsReal.end()) {
			if (debSwitch) debMsg("glutgSelectNextFloatArray","Not found, getting first "<<cnt);
			iter = prms->mGridsReal.begin();
		} else {
			if (debSwitch) debMsg("glutgSelectNextFloatArray","Found, getting next "<<cnt<<" "<<(*iter).second->getName() );
			iter++;

			if (iter==prms->mGridsReal.end()) {
				if (debSwitch) debMsg("glutgSelectNextFloatArray","Last, getting first "<<cnt );
				iter = prms->mGridsReal.begin();
				//return;
			}
		}

		gShowRealgridName = (*iter).second->getName();
		grid = prms->getGridReal(gShowRealgridName);
	}

	debMsg("glutGui","Real-grid '"<<grid->getName()<<"' selected for display");
}

void glutgSelectNextSolver() {
	//if (gFluidSolvers.find(solverName) != gFluidSolvers.end() ) { gpGlutFsolver = gFluidSolvers[solverName]; }
	bool debSwitch = false;

	if (!gpGlutFsolver) {
		if (debSwitch) debMsg("glutgSelectNextSolver","Not inited, getting first");
		gpGlutFsolver = (*gFluidSolvers.begin()).second;
		return;
	} else {

		std::map<std::string, DDF::FluidSolver*>::iterator iter = gFluidSolvers.begin();
		int cnt=0;
		for (; iter!= gFluidSolvers.end(); iter++) {
			//debMsg("---"," = "<<(long int)(gpGlutFsolver)<<" c "<<cnt<<","<< (long int)((*iter).second) );
			if (gpGlutFsolver == (*iter).second) break;
			cnt++;
		}
		if (iter==gFluidSolvers.end()) {
			if (debSwitch) debMsg("glutgSelectNextSolver","Not found, getting first "<<cnt);
			gpGlutFsolver = (*gFluidSolvers.begin()).second;
		} else {
			iter++;
			if (iter==gFluidSolvers.end()) {
				if (debSwitch) debMsg("glutgSelectNextSolver","Last, getting first "<<cnt);
				gpGlutFsolver = (*gFluidSolvers.begin()).second;
			} else {
				if (debSwitch) debMsg("glutgSelectNextSolver","Found, getting next "<<cnt);
				gpGlutFsolver = (*iter).second; 
			}
		}
	}
	debMsg("glutGui","Solver '"<<gpGlutFsolver->getName()<<"' selected for display");

	SolverParams* prms = gpGlutFsolver->getParams();
	if(!prms->haveGridReal(gShowRealgridName)) {
		glutgSelectNextFloatArray();
	}
	if(!prms->haveGridVec3(gShowVelgridName)) {
		glutgSelectNextVec3Array();
	}
}


#endif // DDF_GLUTGUI==1
}; // namespace DDF

