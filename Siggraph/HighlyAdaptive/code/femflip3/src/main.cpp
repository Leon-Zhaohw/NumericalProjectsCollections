/*
 *	main.cpp
 *	
 *	Created by Ryoichi Ando on 11/20/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "util3.h"
#include "email.h"
#include "opengl.h"
#include "glviewer.h"
#include "exporter3.h"
#include "flip3.h"
#include "object3.h"
#include "testcase.h"
#include "octhash3.h"
#include "extsurf3.h"
#include <vector>
using namespace std;

#if defined(__APPLE__) || defined(MACOSX) // For Mac OS X
#define USE_OPENGL		1
#define RESOLUTION		32			// 128
#define FRAMERATE		0			// (1.0/120.0)
#define DURATION		2			// seconds
#define MINSTEP			3
#else // Mostly for Linux
#define USE_OPENGL		0
#define RESOLUTION		256
#define FRAMERATE		(1.0/240.0) //(1.0/60.0)
#define DURATION		2			// seconds
#define MINSTEP			3
#endif

// Window size
unsigned int width = 512;
unsigned int height = 512;

// Useful variables for mouse control...
static int prev_x, prev_y;
static int state = 1;
static int paused = 0;
static FLOAT64 elapsedTime = 0.0;
static FLOAT64 lastTime = 0.0;
static FLOAT64 video_msec = 0.0;

// Sim instances
const uint gsize = RESOLUTION;									// Simulator resolution
const FLOAT64 timePerFrame = FRAMERATE;							// Recoding intervals
static flip3 sim;												// Simulator
static glviewer viewer(sim);									// Viewer
static exporter3 exporter;										// Exporter

static void display(void) {
	// Clear display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	// Visualize simulation
	viewer.drawGL(width,height);
	
	// Swap buffer
	glutSwapBuffers();
}

static void init() {
	// Activate email alert if on recording mode
	email::setTitle("Flip3D");
	if( FRAMERATE ) {
		// Email alert is permanently disabled...
		//email::activate();
	}
	
	// Set background color
	glClearColor(0.0, 0.0, 0.0, 1.0);
	
	// Turn on blending
	glEnable(GL_BLEND);
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	
	// Turn on depth
	glEnable(GL_DEPTH_TEST);
}

static void reshape(int w, int h) {
	// Record new window size
	width = w;
	height = h;
	
	// Set 2% margin
	double margin = 0.02;
	glViewport(0, 0, w, h);
	glLoadIdentity();
	glOrtho(-margin,1.0+margin,-margin,1.0+margin,-1.0,1.0);
}

static void writeFrame( FLOAT64 lastTime, FLOAT64 curTime, FLOAT64 timeInterval, FLOAT64 maxtime ) {
	if( timeInterval ) {
		int prevFrame = lastTime / timeInterval;
		int curFrame = curTime / timeInterval;
		
		// If frame has changed, output an image
		if( curFrame != prevFrame ) {
			const mesher3 *mesher = sim.fluidSolver->getMesh();
			const octree3 *octree = sim.fluidSolver->getOctree();
			exporter.writeMesh(curFrame,sim.surf,sim.particles,sim.dpx,mesher,octree,sim.camera,&sim.objects,true);
#if USE_OPENGL
			viewer.writeImage(curFrame);
#endif
			sim.setVideoFrameNumber(curFrame);
			if( maxtime && maxtime < curTime ) {
				email::print( "Maximum time reached.\n" );
				email::send();
				exit(0);
			}
			FLOAT64 msec = util3::getMilliseconds()-video_msec;
			video_msec = util3::getMilliseconds();
			writeNumber("video_simtime", curTime);
			writeNumber("video_time", msec);
			writeNumber("video_step", curFrame);
			dump( "Took %s to compute %s video frame.\n", tstr(msec), nth(curFrame) );
		}
	}
}

static void idle ( void ) {
	if( ! paused || FRAMERATE ) {
		tick();
		
		// Simulate a step
		sim.setElapsedTime(elapsedTime);
		FLOAT64 maxdt = FRAMERATE ? FRAMERATE-fmod(lastTime,FRAMERATE) : 1.0/30.0;
		if( MINSTEP == 1 ) maxdt = FRAMERATE;
		FLOAT64 curTime = sim.simStep(FRAMERATE,maxdt,FRAMERATE/(FLOAT64)MINSTEP);
		if( testcase::flycamera ) {
			testcase::flycamera->setTime(curTime);
		}
		writeFrame(lastTime, curTime, FRAMERATE, DURATION);
		lastTime = curTime;		
		elapsedTime += tock();
		
		// Send display command
		glutPostRedisplay ();
	}
}

static void keyboard( unsigned char key, int x, int y ) {
	// Hidden keys...
	switch( key ) {
		case '\e':
			exit(0);
			break;
		case '/':
			paused = ! paused;
			break;
		case '.': {
			sim.setElapsedTime(elapsedTime);
			tick();
			FLOAT64 curTime = sim.simStep(1.0,1.0,0.0);
			elapsedTime += tock();
			if( testcase::flycamera ) testcase::flycamera->setTime(curTime);
			glutPostRedisplay();
			break;
		}
		case 'e': {
			if( ! FRAMERATE ) {
				const mesher3 *mesher = sim.fluidSolver->getMesh();
				const octree3 *octree = sim.fluidSolver->getOctree();
				exporter.writeMesh(sim.step,sim.surf,sim.particles,sim.dpx,mesher,octree,sim.camera,&sim.objects,true);
				break;
			}
		}
	}
}

static void passive( int x, int y ) {
	glutPostRedisplay();
}

static void mouse ( int button, int _state, int x, int y ) {
	prev_x = x;
	prev_y = y;
	state = _state;
}

static void motion ( int x, int y ) {
	if( state == 0 ) {
	}
	prev_x = x;
	prev_y = y;
	passive( x, y );
}

int main (int argc, char * argv[]) {
	// Starting recording elapsed time
	tick();
	
	// Build testcase
	vector<object3 *> objs = testcase::buildTestcase();
	
	// See if we want to generate particle mesh
	if( argc > 1 ) {
		// Set root path
		setRootPath("external");
		
		uint from = 0;
		sscanf(argv[1],"%d",&from);
		uint to = from;
		if( argc > 2 ) sscanf(argv[2],"%d",&to);
		uint interval = 1;
		if( argc > 3 ) sscanf(argv[3],"%d",&interval);
		if( interval <= 0 ) {
			dump( "Interval must be > 0");
		}
		uint cnum = from;
		if( argc > 4 ) sscanf(argv[4],"%d",&cnum);
		setConsoleNumber(cnum);
		
		levelset3 *solid = new objLevelset3(object3::SOLID,objs);
		for( uint frame=from; frame<=to; ) {
			setTimestepNumber(frame);
			const char *path = format_str("render/particles/%d_particles.dat.gz",frame);
			if( is_exist(path)) {
				tick(); dump( ">>> Generating %s mesh...\n", nth(frame));
				extsurf3 extsurf;
				extsurf.loadFile(path);
				extsurf.writeMesh(frame,solid);
				surf3 surf;
				surf.vertices = extsurf.vertices;
				surf.normals = extsurf.normals;
				surf.faces = extsurf.faces;	
				exporter.writeMesh(frame,surf,sim.particles,extsurf.dpx,&extsurf.mesher,&extsurf.octree,sim.camera,&sim.objects,false);
				dump( "<<< Done. Took %s.\n", stock("extmesher_time"));
				frame += interval;
			} else {
				// Wait 10 seconds until the file to come out.
				sleep(10);
			}
		}
		exit(0);
	} else {
		// Set root path
		setRootPath("render");
	}
	
	// Delete render folder and recreate it
	if( is_exist(root_path)) {
		run( "rm -rf %s", root_path);
		run( "mkdir %s", root_path );
	}
	
	video_msec = util3::getMilliseconds();
	
	// Init simulator
	if( testcase::flycamera ) {
		testcase::flycamera->setEnabled(true);
		testcase::flycamera->setTime(0.0);
	}
	sim.init(gsize,objs,testcase::gravity,testcase::flycamera);
	
	// Write out first mesh
	if( FRAMERATE ) {
		const mesher3 *mesher = sim.fluidSolver->getMesh();
		const octree3 *octree = sim.fluidSolver->getOctree();
		exporter.writeMesh(0,sim.surf,sim.particles,sim.dpx,mesher,octree,sim.camera,&sim.objects,true);
	}
	elapsedTime += tock();
	sim.setElapsedTime(elapsedTime);
	sim.setVideoFrameNumber(0);
	
#if USE_OPENGL
	// Setup OpenGL stuff..
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GL_DOUBLE | GLUT_DEPTH);
	glutInitWindowPosition ( 100, 100 );
	glutInitWindowSize ( width, height );
	glutCreateWindow(argv[0]);
	glutIdleFunc(idle);
	if( ! FRAMERATE ) {
		glutKeyboardFunc(keyboard);
		glutMouseFunc(mouse);
		glutMotionFunc (motion);
		glutPassiveMotionFunc (passive);
	}
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	init();
	glutMainLoop();
#else
	// Simulate a step
	while(true) {
		// Simulate a step
		sim.setElapsedTime(elapsedTime);
		tick();
		FLOAT64 maxdt = FRAMERATE ? FRAMERATE-fmod(lastTime,FRAMERATE) : 1.0/60.0;
		if( MINSTEP == 1 ) maxdt = FRAMERATE;
		FLOAT64 curTime = sim.simStep(FRAMERATE,maxdt,FRAMERATE/(FLOAT64)MINSTEP);
		if( testcase::flycamera ) {
			testcase::flycamera->setTime(curTime);
		}
		writeFrame(lastTime, curTime, FRAMERATE, DURATION);
		lastTime = curTime;
		elapsedTime += tock();
	}
#endif
	return 0;
}
