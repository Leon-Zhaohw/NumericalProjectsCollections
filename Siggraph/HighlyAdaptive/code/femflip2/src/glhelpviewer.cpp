/*
 *	glhelpviewer.cpp
 *	
 *	Created by Ryoichi Ando on 11/7/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "opengl.h"
#include "glhelpviewer.h"
#include "glviewer.h"
#include "flip2.h"
#include "testcase.h"
#include <string>

#define MAX_RES		16384
#define MIN_RES		8

using namespace std;
using namespace testcase;

static const char visibilityKeys[] = {
	't',
	'y',
	'g',
	'l',
	'u',
	'p',
	'n',
	'z',
	'k',
	'j',
};

static const char controllerKeys[] = {
	'h',
	'r',
	'+',
	'-',
	'd',
	'b',
	'v',
	's',
	' ',
    'a',
	'm',
	'i',
	'x',
};

enum stateName {
	S_HELP,
	S_RES,
	S_SURF,
	S_VARIATION,
	S_CORRECTION,
	S_SCENE,
    S_ADAPTIVE,
	S_SOLVER,
	S_REMESH,
	S_MESHER,
};

enum keyName {
	K_HELP,
	K_RESET,
	K_INCR,
	K_DECR,
	K_RESTORE,
	K_SURF,
	K_VARIATION,
	K_CORRECTION,
	K_SCENE,
    K_ADAPTIVE,
	K_SOLVER,
	K_REMESH,
	K_MESHER
};

static const char *controllerString[] = {
	"Push \"%c\" to toggle help. %s",
	"Push \"%c\" to reset the simulation. %s",
	"Push \"%c\" to increase the resolution. %s",
	"Push \"%c\" to decrease the resolution. %s",
	"Push \"%c\" to restore defaults. %s",
	"Push \"%c\" to change the accuracy of boundary condition. %s",
	"Push \"%c\" to toggle variational pressure solver. %s",
	"Push \"%c\" to toggle spring correction. %s",
	"Push \"%c\" to change the test case. %s",
    "Push \"%c\" to toggle adaptive sampling. %s",
	"Push \"%c\" to change fluid solver. %s",
	"Push \"%c\" to toggle grid remeshing. %s",
	"Push \"%c\" to toggle external surface reconstructor. %s",
};

static const char *visibilityString[] = {
	"Push \"%c\" to toggle sim time visibility. %s",
	"Push \"%c\" to toggle grid line visibility. %s",
	"Push \"%c\" to toggle grid velocity visibility. %s",
	"Push \"%c\" to toggle liquid levelset visibility. %s",
	"Push \"%c\" to toggle particle velocity visibility. %s",
	"Push \"%c\" to toggle pressure visibility. %s",
	"Push \"%c\" to toggle neighborhood visibility. %s",
	"Push \"%c\" to toggle particle radius visibility. %s",
	"Push \"%c\" to toggle matrix connection visibility. %s",
	"Push \"%c\" to toggle external mesh visibility. %s",
};

static const char *visibilityStatusStr[2] = {
	"Current: shown", "Current: hidden",
};

static void (glviewer::*funcPtr[])(bool) = {
	&glviewer::setSimTimeVisibility,
	&glviewer::setGridlineVisibility,
	&glviewer::setGridVelocityVisibility,
	&glviewer::setLevelsetVisibility,
	&glviewer::setParticleVelocityVisibility,
	&glviewer::setPressureVisibility,
	&glviewer::setNeighborConnectionVisibility,
	&glviewer::setParticleAnisotropyVisibility,
	&glviewer::setMatrixConnectionVisibility,
	&glviewer::setMeshGeneratorVisibility,
};

static size_t numVisibilityKeys = sizeof(visibilityString)/sizeof(const char *);
static size_t numControllerKeys = sizeof(controllerKeys)/sizeof(char);

glhelpviewer::glhelpviewer( glviewer& viewer, flip2& sim ) : viewer(viewer), sim(sim) {
	visibilityStates.push_back(false);	// show sim time
	visibilityStates.push_back(true);	// grid line
	visibilityStates.push_back(false);	// grid velocity
	visibilityStates.push_back(true);	// liquid levelset
	visibilityStates.push_back(false);	// particle velocity
	visibilityStates.push_back(false);	// pressure
	visibilityStates.push_back(false);	// neighborhood
	visibilityStates.push_back(false);	// radius view
	visibilityStates.push_back(false);	// matrix connection
	visibilityStates.push_back(true);	// external mesh view
	default_visibilityStates = visibilityStates;

	controllerStates.push_back(0);		// show help
	controllerStates.push_back(64);		// resolution
	controllerStates.push_back(2);		// free surface boundary accuracy
	controllerStates.push_back(1);		// variational pressure solver
	controllerStates.push_back(1);		// spring correction
	controllerStates.push_back(0);		// test case
    controllerStates.push_back(1);		// adaptive sampling
	controllerStates.push_back(1);		// fluid solver
	controllerStates.push_back(1);		// grid remesh
	controllerStates.push_back(0);		// external surface reconstruction method
	default_controllerStates = controllerStates;
	controllerStatusStr.resize(numControllerKeys);
	setControllerString();
}

void glhelpviewer::setControllerString() {
	char str[64];
	
	controllerStatusStr[K_SURF] = controllerStates[S_SURF] == 2 ? 
	"Current: 2nd order" : "Current: 1st order";
	
	controllerStatusStr[K_VARIATION] = controllerStates[S_VARIATION] ? 
	"Current: Variational solver" : "Current: Voxelized solver";
	
	int state = controllerStates[S_CORRECTION];
	if( state == 0 ) controllerStatusStr[K_CORRECTION] = "Current: Disabled";
	else if( state == 1 ) controllerStatusStr[K_CORRECTION] = "Current: Enabled";
	
	state = controllerStates[S_SOLVER];
	if( state == 0 ) controllerStatusStr[K_SOLVER] = "Current: MAC solver";
	else if( state == 1 ) controllerStatusStr[K_SOLVER] = "Current: FEM solver";
	else if( state == 2 ) controllerStatusStr[K_SOLVER] = "Current: FVM solver";
	else if( state == 3 ) controllerStatusStr[K_SOLVER] = "Current: Octree solver";
	
	controllerStatusStr[K_REMESH] = controllerStates[S_REMESH] ? 
	"Current: Enabled" : "Current: Disabled";
	
    controllerStatusStr[K_ADAPTIVE] = controllerStates[S_ADAPTIVE] ? 
	"Current: Enabled" : "Current: Disabled";
	
	int num = controllerStates[S_SCENE]+1;
	string suffix = "th";
	if( num == 1 ) suffix = "st";
	else if( num == 2 ) suffix = "nd";
	else if( num == 3 ) suffix = "rd";
	sprintf( str, "Current %d%s case", num, suffix.c_str() );
	controllerStatusStr[K_SCENE] = str;
	
	uint gn = sim.gn;
	sprintf( str, "To %d^2", imin(MAX_RES,gn*2) );
	controllerStatusStr[K_INCR] = str;
	
	sprintf( str, "To %d^2", imax(MIN_RES,gn/2) );
	controllerStatusStr[K_DECR] = str;
	
	state = controllerStates[S_MESHER];
	if( state == 0 ) controllerStatusStr[K_MESHER] = "Current: Ours";
	else if( state == 1 ) controllerStatusStr[K_MESHER] = "Current: Adams07";
	else if( state == 2 ) controllerStatusStr[K_MESHER] = "Current: Zhu & Bridson 2005";
	else if( state == 3 ) controllerStatusStr[K_MESHER] = "Current: Yu and Turk 2010";
}

void glhelpviewer::applyStates( bool forceReset ) {
	for( size_t i=0; i<numVisibilityKeys; i++ ) {
		if( funcPtr[i] ) {
			(viewer.*funcPtr[i])(visibilityStates[i]);
		}
	}
	
	int order = controllerStates[S_SURF];
	sim.setBoundaryAccuracy(order);
	sim.setVariation(controllerStates[S_VARIATION]);
	int state = controllerStates[S_CORRECTION];
	sim.setCorrection(state);
	sim.setAdaptiveSampling(controllerStates[S_ADAPTIVE]);
	sim.setFluidSolver(controllerStates[S_SOLVER]);
	sim.setRemesh(controllerStates[S_REMESH]);
	viewer.setSurfaceExtractorMethod(controllerStates[S_MESHER]);
	setControllerString();
	
	bool reset = forceReset;
	if( sim.gn != controllerStates[S_RES] ) reset = true;
	if( reset ) resetSim();
}

void glhelpviewer::resetSim() {
	int sceneNum = controllerStates[S_SCENE];
	sim.init(controllerStates[S_RES],scenes[sceneNum][0],scenes[sceneNum][1],1.0);
}

void glhelpviewer::restoreDefaults() {
	bool force = controllerStates[S_SCENE] != 0;
	visibilityStates = default_visibilityStates;
	controllerStates = default_controllerStates;
	applyStates(force);
}

bool glhelpviewer::keyDown( char key ) {
	bool keyHandled = false;
	// Detect vibility key change
	for( int i=0; i<numVisibilityKeys; i++ ) {
		if( key == visibilityKeys[i] ) {
			visibilityStates[i] = ! visibilityStates[i];
			keyHandled = true;
		}
	}
	bool reset = false;
	
	if( key == controllerKeys[K_HELP] ) {
		controllerStates[K_HELP] = ! controllerStates[K_HELP];
		keyHandled = true;
	} else if( key == controllerKeys[K_RESTORE] ) {
		restoreDefaults();
		keyHandled = true;
	} else if( key == controllerKeys[K_RESET] ) {
		reset = true;
		keyHandled = true;
	} else if( key == controllerKeys[K_INCR] ) {
		uint gn = controllerStates[S_RES];
		if( gn < MAX_RES ) {
			controllerStates[S_RES] = 2*gn;
		}
		keyHandled = true;
	} else if( key == controllerKeys[K_DECR] ) {
		uint gn = controllerStates[S_RES];
		if( gn > MIN_RES ) {
			controllerStates[S_RES] = gn/2;
		}
		keyHandled = true;
	} else if( key == controllerKeys[K_SURF] ) {
		int order = controllerStates[S_SURF];
		order = order == 1 ? 2 : 1;
		controllerStates[S_SURF] = order;
		keyHandled = true;
	} else if( key == controllerKeys[K_VARIATION] ) {
		int enabled = controllerStates[S_VARIATION];
		enabled = enabled == 1 ? 0 : 1;
		controllerStates[S_VARIATION] = enabled;
		keyHandled = true;
	} else if( key == controllerKeys[K_CORRECTION] ) {
		int state = controllerStates[S_CORRECTION];
		state = (state+1)%2;
		controllerStates[S_CORRECTION] = state;
		keyHandled = true;
	} else if( key == controllerKeys[K_SCENE] ) {
		int num = controllerStates[S_SCENE];
		num = (num+1) % NUMSCENE;
		controllerStates[S_SCENE] = num;
		reset = true;
		keyHandled = true;
	} else if( key == controllerKeys[K_ADAPTIVE] ) {
		int enabled = controllerStates[S_ADAPTIVE];
		enabled = enabled == 1 ? 0 : 1;
		controllerStates[S_ADAPTIVE] = enabled;
		keyHandled = true;
	} else if( key == controllerKeys[K_SOLVER] ) {
		int name = controllerStates[S_SOLVER];
		name = (name+1)%4;
		controllerStates[S_SOLVER] = name;
		keyHandled = true;
	} else if( key == controllerKeys[K_REMESH] ) {
		int enabled = controllerStates[S_REMESH];
		enabled = enabled == 1 ? 0 : 1;
		controllerStates[S_REMESH] = enabled;
		keyHandled = true;
	} else if( key == controllerKeys[K_MESHER] ) {
		int num = controllerStates[S_MESHER];
		num = (num+1) % 4;
		controllerStates[S_MESHER] = num;
		keyHandled = true;
	}
	if( reset ) resetSim();
	
	// Apply changes
	applyStates();
	return keyHandled;
}

void glhelpviewer::drawGL() {
#if ! RECORD
	// Setup display spacing info
	int winsize[2] = { glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) };
	double dw = 1.0/(double)winsize[0];
	double dh = 1.0/(double)winsize[1];
	int peny = 15;
	char text[256];
	char statusText[256];
	
	// Draw half shadow
	if( controllerStates[S_HELP] ) {
		glColor4f(0.0,0.0,0.0,0.5);
		glRectf(0.0,0.0,1.0,1.0);
	}
	
	// Now show texts...
	glColor4f(1.0,1.0,1.0,1.0);
	size_t i=0;
	for( i=0; i<numControllerKeys; i++ ) {
		if( i==0 || controllerStates[S_HELP] ) {
			glRasterPos2d(30*dw, 1.0-(i+1)*(peny)*dh-20*dh);
			const char *statusStr = controllerStatusStr[i].c_str();
			if( statusStr[0] ) {
				sprintf( statusText, "( %s )", statusStr );
			} else {
				statusText[0] = 0;
			}
			
			sprintf( text, controllerString[i], controllerKeys[i], statusText );
			glviewer::drawBitmapString(text,GLUT_BITMAP_HELVETICA_10);
		}
	}
	
	for( size_t j=0; j<numVisibilityKeys; j++ ) {
		if( i==0 || controllerStates[S_HELP] ) {
			glRasterPos2d(30*dw, 1.0-(j+i+1.5)*(peny)*dh-20*dh);
			const char *statusStr = visibilityStatusStr[1-visibilityStates[j]];
			if( statusStr[0] ) {
				sprintf( statusText, "( %s )", statusStr );
			} else {
				statusText[0] = 0;
			}
			
			sprintf( text, visibilityString[j], visibilityKeys[j], statusText );
			glviewer::drawBitmapString(text,GLUT_BITMAP_HELVETICA_10);
		}
	}
#endif
}