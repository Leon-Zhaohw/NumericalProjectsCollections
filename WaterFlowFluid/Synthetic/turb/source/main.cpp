/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Main function, solver description
 *
 *****************************************************************************/

#include "globals.h"
#include "fluidsolver.h"

#include "solverparams.h"
#include "paramset.h"
#include "solverplugin.h"
#include "solverinit.h"
using namespace std;

namespace DDF
{
	extern void closeDatabaseWriter();
	extern std::vector<Vec3> initDatabaseWriter(const string&, int, int);

	extern int showGlutGui(const char* text); // from glutgui.cpp
	SolverObject *solver;
	vector<SolverObject*> solvers;

	// functions taken over from standard main
	bool advanceAllSolvers()
	{
		for (unsigned i=0;i<solvers.size();i++) {
			if (solvers[i]->performStep()) return true;
		}
		return false;
	}

	void finalizeAllSolvers()
	{
		for (unsigned i=0;i<solvers.size();i++)
			solvers[i]->finalize();
	}

	void run(bool allowGui=true)
	{
		solver = solvers[0];
		// run solver
		for (unsigned i=0;i<solvers.size();i++)
			solvers[i]->forceInit();
				
		if ( DDF_GLUTGUI==1 && allowGui )
			showGlutGui("DDF noparser");
		else
			while ( !advanceAllSolvers() ) {};
		
		finalizeAllSolvers();
		for (unsigned i=0;i<solvers.size();i++)
			delete solvers[i];

		solvers.clear();
	}

};
using namespace DDF;

void makeScenes()
{
	{ // static scene
		SolverObject* solver = new SolverObject( "makescene", nVec3i ( 200,80,150 ) );
		solver->createVec3Grid ( "normal" );
		solver->createRealGrid ( "dist" );
		solver->addInitPlugin ( "init-box-domain",  IntArg ( "flag-inside",FFLUID ) + IntArg ( "flag-border",FINFLOW ) + IntArg("flag-floor", FOBSTACLE) );
		// obstacles
		solver->addInitPlugin ( "init-sphere", VecArg("center",Vec3(0.25,0.6,0.35)) + RealArg("radius",0.1) + StringArg("norm","normal")+StringArg("dist","dist"));	
		solver->addInitPlugin ( "init-box", VecArg("pos1",Vec3(0.2,0.0,0.3)) + VecArg("pos2",Vec3(0.3,0.6,0.4)) + StringArg("norm","normal")+StringArg("dist","dist"));	
		solver->addInitPlugin ( "init-box", VecArg("pos1",Vec3(0.3,0.0,0.65)) + VecArg("pos2",Vec3(0.5,0.3,0.9)) + StringArg("norm","normal")+StringArg("dist","dist"));	
		solver->addInitPlugin ( "init-box", VecArg("pos1",Vec3(0.6,0.0,0.3)) + VecArg("pos2",Vec3(0.65,0.4,0.5)) + StringArg("norm","normal")+StringArg("dist","dist"));	
		// smoke seed
		solver->addInitPlugin ( "init-sphere", VecArg("center",Vec3(0.3,0.65,0.35)) + RealArg("radius",0.1) + IntArg("type", FDENSITYSOURCE));	
		solver->addInitPlugin ( "init-box", VecArg("pos1",Vec3(0.25,0.0,0.3)) + VecArg("pos2",Vec3(0.35,0.65,0.4)) + IntArg("type", FDENSITYSOURCE));	
		solver->addInitPlugin ( "init-box", VecArg("pos1",Vec3(0.35,0.0,0.65)) + VecArg("pos2",Vec3(0.55,0.35,0.9)) + IntArg("type", FDENSITYSOURCE));	
		solver->addInitPlugin ( "init-box", VecArg("pos1",Vec3(0.65,0.0,0.3)) + VecArg("pos2",Vec3(0.7,0.45,0.5)) + IntArg("type", FDENSITYSOURCE));	

		solver->addInitPlugin ( "dump-universal", StringArg ( "grid", "flags" ) + StringArg ( "override-name","scene/static-flags" )  + IntArg ( "single-dump", 1 ) );	
		solver->addInitPlugin ( "dump-universal", StringArg ( "grid", "dist" ) + StringArg ( "override-name","scene/static-dist" )  + IntArg ( "single-dump", 1 ) );	
		solver->addInitPlugin ( "dump-universal", StringArg ( "grid", "normal" ) + StringArg ( "override-name","scene/static-normal" )  + IntArg ( "single-dump", 1 ) );	
		solvers.push_back(solver);
		run(false);
	}
	{ // moving objects
		SolverObject* solver = new SolverObject( "makescene", nVec3i ( 80,80,80 ) );
		solver->createVec3Grid ( "normal" );
		solver->createRealGrid ( "dist" );
		solver->addInitPlugin ( "init-box-domain",  IntArg ( "flag-inside",FFLUID ) + IntArg ( "flag-border",FINFLOW ) );
		solver->addInitPlugin ( "init-box", VecArg("pos1",Vec3(0.2,0.4,0.4)) + VecArg("pos2",Vec3(0.8,0.5,0.6)) + StringArg("norm","normal")+StringArg("dist","dist"));	
		solver->addInitPlugin ( "init-box", VecArg("pos1",Vec3(0.4,0.5,0.4)) + VecArg("pos2",Vec3(0.6,0.6,0.6)) + StringArg("norm","normal")+StringArg("dist","dist"));	
		solver->addInitPlugin ( "init-sphere", VecArg("center",Vec3(0.4,0.5,0.5)) + RealArg("radius",0.1) + StringArg ( "norm","normal" ) + StringArg ( "dist","dist" ));
		solver->addInitPlugin ( "dump-universal", StringArg ( "grid", "flags" ) + StringArg ( "override-name","scene/dyn-flags" )  + IntArg ( "single-dump", 1 ) );	
		solver->addInitPlugin ( "dump-universal", StringArg ( "grid", "dist" ) + StringArg ( "override-name","scene/dyn-dist" )  + IntArg ( "single-dump", 1 ) );	
		solver->addInitPlugin ( "dump-universal", StringArg ( "grid", "normal" ) + StringArg ( "override-name","scene/dyn-normal" )  + IntArg ( "single-dump", 1 ) );	
		solvers.push_back(solver);
		run(false);
	}
}

void doPrecompute(const std::string& name, const Vec3& inflow, int frames, bool dynamic, const Vec3& rotAxis = Vec3(0.), Real rotSpeed = 0.)
{
	bool rotate = (rotSpeed != 0.);
	SolverObject* solver = new SolverObject( "precompute", "scene/" + name + "-flags.gz" );

	// create grids
	solver->createVec3Grid ( "normal" );
	solver->createRealGrid ( "dist" );
	solver->createVec3Grid ( "mean-vel" );
	solver->createVec3Grid ( "abl" );
	solver->addStandardSolverGrids();
	
	// additional grids for rot. precomputation
	if (rotate) {
		solver->createIntGrid ("obstacle-flags");
		solver->createIntGrid ("empty-flags");
		solver->addInitPlugin ( "load-universal", StringArg("grid","obstacle-flags") + StringArg("file","scene/" + name + "-flags.gz"));
		solver->addInitPlugin ( "init-box-domain",  StringArg("gridname","empty-flags") + IntArg ( "flag-inside",FFLUID ) + IntArg ( "flag-border",FINFLOW ) );		
	}

	// load grids, initialize fluid velocities
	solver->addInitPlugin ( "load-universal", StringArg("grid","dist") + StringArg("file","scene/" + name + "-dist.gz"));
	solver->addInitPlugin ( "load-universal", StringArg("grid","normal") + StringArg("file","scene/" + name + "-normal.gz"));
	solver->addInitPlugin ( "set-conditional", StringArg ( "gridname","vel-curr" ) + VecArg ( "target-vec",inflow ) + IntArg ( "flag", FFLUID ) );

	// program solver main loop
	solver->addPlugin ( "set-conditional", StringArg ( "gridname","vel-curr" ) + VecArg ( "target-vec",inflow ) + IntArg ( "flag", FINFLOW ) );
	solver->addPlugin ( "maccormack-advect-vec3", StringArg ( "vel-src", "vel-curr" ) );
	solver->addPlugin ( "set-noslip-bcs", StringArg ( "grid","vel-curr" ) );
	solver->addPlugin ( "diffuse-grid", StringArg ( "src-vec3", "vel-curr" ) + RealArg ( "diff", 0.3 ) );
	if (rotate)
		solver->addPlugin ("set-moving-obs-bcs", StringArg("obstacle","obstacle-flags") + StringArg("flags-src","empty-flags") +
							VecArg("obs-rot-axis", rotAxis) + RealArg("obs-rot-vel", rotSpeed) + VecArg("obs-center", Vec3(0.5,0.5,0.5)));
	solver->addPlugin ( "solve-pressure", IntArg ( "openbound",0 ) );
	if (rotate)
		solver->addPlugin ( "average", StringArg ( "gridname","vel-curr" ) + StringArg ( "sumgrid","mean-vel" ) + IntArg ( "from", frames ) + IntArg ( "frames", 3 ) + IntArg ( "post-quit",1 ) + IntArg("stride", frames) );
	else
		solver->addPlugin ( "average", StringArg ( "gridname","vel-curr" ) + StringArg ( "sumgrid","mean-vel" ) + IntArg ( "from", frames ) + IntArg ( "frames", frames ) + IntArg ( "post-quit",1 ) );
	
	// program final steps
	solver->addEndPlugin ( "calc-abl", StringArg ("mean-vel","mean-vel") + StringArg ("dist","dist") + StringArg("normal","normal") + StringArg("abl","abl") + RealArg("d", 1.7));
	
	if (dynamic)
		solver->addEndPlugin ( "add-database", StringArg("grid","abl") + StringArg("normal","normal") + VecArg("u0",rotate ? (rotSpeed*rotAxis):inflow));
	else {
		solver->addEndPlugin ( "dump-universal", StringArg ( "grid","abl" ) + StringArg ( "override-name","scene/" + name + "-abl" )  + IntArg ( "single-dump", 1 ) );	
		solver->addEndPlugin ( "dump-universal", StringArg ( "grid","mean-vel" ) + StringArg ( "override-name","scene/" + name + "-mean" ) + IntArg ( "single-dump", 1 ) );
	}
		
	solvers.push_back(solver);		
	run();
}

void runStatic(const Vec3& inflow)
{
	SolverObject* solver = new SolverObject( "run_static", "scene/static-flags.gz" );
	solver->getParams().mU0 = inflow;
	solver->getParams().mTimestepAnim = 0.005;
		
	// create grids
	solver->createVec3Grid ( "mean-flow" );
	solver->createRealGrid ( "dist" );
	solver->createVec3Grid ( "vorticity" );
	solver->createVec3Grid ( "ABL" );
	solver->createVec3Grid ( "pre-ABL" );
	solver->createVec3Grid ( "vort" );
	solver->createRealGrid ( "pdf" );
	solver->createRealGrid ( "density" );
	solver->addStandardSolverGrids();
	solver->createNoiseField("noise", Vec3(0.), Vec3(50,50,50), -0.4, 2.0, 0.002);

	// program solver initialization process
	solver->addInitPlugin ( "load-universal", StringArg("grid","dist") + StringArg("file","scene/static-dist.gz"));
	solver->addInitPlugin ( "load-universal", StringArg("grid","mean-flow") + StringArg("file","scene/static-mean.gz"));
	solver->addInitPlugin ( "load-universal", StringArg("grid","pre-ABL") + StringArg("file","scene/static-abl.gz"));
	
	// program solver main loop
	solver->addPlugin ( "copy-grid", StringArg ( "src","mean-flow" ) + StringArg ( "dest","vel-curr") );
	solver->addPlugin ("init-density-inflow", StringArg("density","density") + RealArg("target-value",0.7) + IntArg("flag", FDENSITYSOURCE) + StringArg("noise","noise")); 
	solver->addPlugin ( "gen-vpart", StringArg ("source","pre-ABL") + StringArg ("flow","ABL") + StringArg("dist","dist") + StringArg("pdf","pdf") +
						RealArg("thres-vort", 2e-2) + RealArg("thres-pdf",5e-5) + RealArg("mult-pdf",1) + RealArg("scale-flow",0.94) + RealArg("max-bl",0.15) +
						RealArg("min-dist", 3) + RealArg("min-rad", 3) + RealArg("max-rad", 7) + RealArg("vortex-gain", 2.5) + RealArg("fade-in", 0));
 
	solver->addPlugin ( "semi-lagr-advect-vec3", StringArg ( "vel-src","ABL" ) + IntArg ( "mac", 0) );
	solver->addPlugin ( "apply-vpart", StringArg ( "vorticity", "vorticity" ) );
	solver->addPlugin ( "advect-vpart");
	solver->addPlugin ( "merge-vpart", StringArg("ndist","dist") + RealArg("init-time", 30) + RealArg("decay-time",450) + RealArg("merge-dist",0.8) +
   						RealArg("dissipate-radius",1) + RealArg("radius-cascade",1.5) );
	solver->addPlugin ("compute-vorticity", StringArg("vorticity","vort"));	
	solver->addPlugin ("maccormack-advect-real", StringArg("real-src","density"));
	solver->addPlugin( "dump-df3", IntArg("max-frames",200) + StringArg("gridname","density") + StringArg("prefix","sta") + RealArg("start-time",0) + IntArg("pbrt",1) + StringArg ( "override-name", "render/static" ));

	solvers.push_back(solver);
	run();
}


void runDynamic()
{
	SolverObject* sMain = new SolverObject( "main", nVec3i(100,40,100));
	SolverObject* sUpscale = new SolverObject( "upscale", *sMain, 2);
	SolverObject* sAnim = new SolverObject( "anim", "scene/dyn-flags.gz");
	sMain->getParams().mU0 = Vec3(.1,0,0);
	sUpscale->getParams().mHostVorticitySystem = false;
	sUpscale->getParams().mTimestepAnim = 0.01;
	sAnim->getParams().mHostVorticitySystem = false;
		
	// define grids
	sMain->createIntGrid ( "static-flags" );
	sMain->createVec3Grid ( "abl" );
	sMain->createVec3Grid ( "abl-pre" );
	sMain->createVec3Grid ( "abl-pre-static" );
	sMain->createRealGrid ( "ndist" );
	sMain->createRealGrid ( "ndist-static" );
	sMain->createRealGrid ( "pdf" );
	sMain->addStandardSolverGrids();
	sUpscale->createIntGrid ( "static-flags" );
	sUpscale->createVec3Grid ( "vort" );
	sUpscale->createRealGrid ( "density" );
	sUpscale->addStandardSolverGrids();
	sUpscale->createNoiseField("noise", Vec3(0.), Vec3(50,50,50), -0.3, 2.0, 0.01);
	sAnim->createRealGrid ( "ndist-static" );
	
	// main solver init
	sMain->addInitPlugin ("init-box-domain", IntArg("flag-inside",FFLUID ) + IntArg("flag-border",FINFLOW) + IntArg("flag-floor", FOBSTACLE) +
							 StringArg("gridname","static-flags"));
	sMain->addInitPlugin ("init-box-domain", IntArg("flag-inside",FFLUID ) + IntArg("flag-border",FINFLOW) + IntArg("flag-floor", FOBSTACLE));
	
	// main solver loop
	sMain->addPlugin ("interpolate-grid-from", StringArg("name", "s_upscale") + StringArg("src","vel-curr") + StringArg("dst","vel-curr"));
	sMain->addPlugin ("maccormack-advect-vec3", StringArg("vel-src", "vel-curr"));
	sMain->addPlugin ("animate", StringArg("file", "scene/anim.txt") + StringArg("params","s_anim") + StringArg("flags-src","static-flags") + StringArg("flags-dst","flags") +
					StringArg("ndist-src","ndist-static") + StringArg("ndist-dst","ndist") + StringArg("vort-src","abl-pre-static") + StringArg("vort-dst","abl-pre") + 
					StringArg("vort-file","scene/dyn-database.gz"));

	sMain->addPlugin ("gen-vpart", StringArg("source","abl-pre") + StringArg("flow","abl") + StringArg("dist","ndist") + StringArg("pdf","pdf") +
					RealArg("thres-vort",5e-2) + RealArg("thres-pdf",4e-5) + RealArg("mult-pdf",0.5) + RealArg("scale-flow",0.99) +
					RealArg("min-dist",1) + RealArg("min-rad",2) + RealArg("max-rad",5) + RealArg("vortex-gain",0.4));

	sMain->addPlugin ("semi-lagr-advect-vec3", StringArg("vel-src","abl") + IntArg("mac",0));
	sMain->addPlugin ("merge-vpart", StringArg("ndist","ndist") + RealArg("init-time",40.0) + RealArg("decay-time",200.0) + RealArg("merge-dist", 0.5) +
					RealArg("dissipate-radius",0.5) + RealArg("radius-cascade",1.3));

	sMain->addPlugin ( "set-noslip-bcs", StringArg ( "grid","vel-curr" ) );
	sMain->addPlugin ("set-conditional", StringArg("gridname","vel-curr") + VecArg("target-vec", Vec3(0.)) + IntArg("flag",FINFLOW));
	sMain->addPlugin ("solve-pressure", IntArg("openbound",0));

	// upscale solver init
	sUpscale->addInitPlugin ("init-box-domain", IntArg("flag-inside",FFLUID ) + IntArg("flag-border",FINFLOW) + IntArg("flag-floor", FOBSTACLE) +
							 StringArg("gridname","static-flags"));
	sUpscale->addInitPlugin ("init-box-domain", IntArg("flag-inside",FFLUID ) + IntArg("flag-border",FINFLOW) + IntArg("flag-floor", FOBSTACLE));
	
	// upscale solver loop
	sUpscale->addPlugin("interpolate-grid-from", StringArg("name","s_main") + StringArg("src","vel-curr") + StringArg("dst","vel-curr"));
	sUpscale->addPlugin("animate", StringArg("file", "scene/anim.txt") + StringArg("params","s_anim") + StringArg("flags-src","static-flags") +
					StringArg("flags-dst","flags") + IntArg("do-forces",0));
		
	sUpscale->addPlugin ("init-density-inflow", StringArg("density","density") + RealArg("target-value",0.7) + IntArg("flag", FDENSITYSOURCE) + StringArg("noise","noise")); 

	sUpscale->addPlugin ("compute-vorticity", StringArg("vorticity","vort"));
	sUpscale->addPlugin ("apply-vpart", StringArg("vorticity","vort") + StringArg("from-solver", "s_main"));

	sUpscale->addPlugin ("maccormack-advect-real", StringArg("real-src","density"));
	sUpscale->addPlugin ("advect-vpart", StringArg("from-solver", "s_main"));
	sUpscale->addPlugin( "dump-df3", IntArg("max-frames",200) + StringArg("gridname","density") + StringArg("prefix","dyn") + RealArg("start-time",0) + IntArg("pbrt",1)  + StringArg ( "override-name", "render/dyn" ));

	// anim "solver"
	sAnim->addInitPlugin ("load-universal", StringArg("grid","ndist-static") + StringArg("file","scene/dyn-dist.gz"));
	sAnim->addPlugin ("dummy");
	
	solvers.push_back(sMain);
	solvers.push_back(sUpscale);
	solvers.push_back(sAnim);
	run();
}


void createDatabase(const string& name, int nTheta, int nPhi)
{	
	const int nFrames = 50;
	const int rotSpeed = 5;

	// obtain frames to simulate from database plugin
	vector<Vec3> axes = initDatabaseWriter("scene/" + name + "-database", nTheta, nPhi);
	
	for(unsigned i=0;i<axes.size();++i)
	{
		// translational component		
		Vec3 inflow = axes[i] * .3;
		cout << endl << endl << " PRECOMPUTING TRANSLATION ELEMENT " << axes[i] << " (" << i+1 << " of "<< axes.size() << ")" << endl << endl;
		doPrecompute (name, inflow, nFrames, true);
		
		// rotational component
		cout << endl << endl << " PRECOMPUTING ROTATION ELEMENT " << axes[i] << " (" << i+1 << " of "<< axes.size() << ")" << endl << endl;
		doPrecompute (name, Vec3(0.), 360/rotSpeed, true, axes[i], rotSpeed);		
	}

	// closes database file
	closeDatabaseWriter();
}

//*****************************************************************************
// main...
int main ( int argc,const char* argv[] )
{
	if (argc==2 && !strcmp(argv[1],"scene")) {
		makeScenes();
	}
	else if (argc==2 && !strcmp(argv[1],"pre-static")) {
		// PRECOMPUATION STATIC SCENE
		Vec3 inflow (0.4, 0, 0);
		const int nFrames = 50;
		doPrecompute("static", inflow, nFrames, false);		
	}
	else if (argc==2 && !strcmp(argv[1],"pre-dynamic")) {
		// PRECOMPUATION DYNAMIC SCENE
		const int nTheta = 5, nPhi = 8;
		createDatabase("dyn", nTheta, nPhi);	
	}  
	else if (argc==2 && !strcmp(argv[1],"run-static")) {
		// STATIC SOLVER		
		Vec3 inflow (0.4, 0, 0);
		runStatic(inflow);
	}
	else if (argc==2 && !strcmp(argv[1],"run-dynamic")) {
		runDynamic();
	}
	else {
		cout << "usage : turbCmd/turbGui mode" << endl <<
		        "  where mode is " << endl << 
			"    scene       : generate scene files. has to be run first." << endl <<
			"    pre-static  : precompute static scene" << endl <<
			"    pre-dynamic : precompute dynamic scene" << endl <<
			"    run-static  : run static scene" << endl <<
			"    run-dynamic : run dynamic scene" << endl;	  
	}
	
	return 0;
}


