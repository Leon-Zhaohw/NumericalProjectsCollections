/*
 *  main.cpp
 *	
 *	Created by Ryoichi Ando on 2012/05/17
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "email.h"
#include "opengl.h"
#include "testcase.h"
#include "flip2.h"
#include "glviewer.h"
#include "glhelpviewer.h"

#define FRAMERATE		0	//(1.0/120.0)  // framerate of the video
#define DURATION		5	// seconds
#define MINSTEP			2
#define PATH			"/Users/mscp/Desktop/img"		// Change here to change where to export images.

// Window size
static uint win_x = 600;
static uint win_y = 600;

// Maginification variables
FLOAT64 scale = 1.0;
vec2d center(0.5,0.5);
FLOAT64 trans = 0.0;
FLOAT64 goal_scale = 1.0;
vec2d goal_center = center;

// Useful variables for mouse control...
static int prev_x, prev_y;
static int cur_x, cur_y;
static int state = 1;
static int paused = 0;

// Debug level defined in util2.h
int debugLevel = 1;

// Simulator variables
static flip2 sim;							// Simulator
static glviewer viewer(sim);				// Viewer
static glhelpviewer helpviewer(viewer,sim);	// Help viewer

// Variables to record time
FLOAT64 lastTime = 0.0;

static void display(void) {
	// Clear display
	glClear(GL_COLOR_BUFFER_BIT);
	
	// Visualize simulation
	viewer.drawGL();
	
	// Show help
	if( ! FRAMERATE ) helpviewer.drawGL();

	// Swap buffer
	glutSwapBuffers();
}

static void transition ( FLOAT64 r=0.2 ) {
	center = (goal_center-center)*r+center;
	scale = (goal_scale-scale)*r+scale;
}

static void writeFrame( FLOAT64 lastTime, FLOAT64 curTime, FLOAT64 timeInterval, FLOAT64 maxtime, const char *path ) {
	if( timeInterval ) {
		int prevFrame = lastTime / timeInterval;
		int curFrame = curTime / timeInterval;
		
		// If frame has changed, output an image
		if( curFrame != prevFrame ) {
			if( curFrame == 1 ) {
				run( "rm -rf %s/*", path );
			}
			viewer.writeImage(format_str("%s/frame_%d.bmp", path, curFrame));
			//viewer.writeSVG(format_str("%s/frame_%d.svg", path, curFrame), PATH);
			if( maxtime && maxtime < curTime ) {
				email::print( "Maximum time reached.\n" );
				email::send();
				exit(0);
			}
		}
	}
}

static void init( const char *path ) {
	if( path ) {
		printf( "Path = %s\n", path );
		exit(0);
	}
	if( FRAMERATE ) paused = 0;
	
	// Activate email alert if on recording mode
	email::setTitle("Flip2D");
	if( FRAMERATE ) {
		//email::activate();
	}
	
	// Set background color
	glClearColor(0.0, 0.0, 0.0, 1.0);
	
	// Turn on blending
	glEnable(GL_BLEND);
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

	// Setup helpview
	helpviewer.applyStates(true);
}

static void reshape(int w, int h) {
	// Reshape window to default size
	if( FRAMERATE && (w!=win_x || h!=win_y) ) {
		glutReshapeWindow(win_x,win_y);
	} else {
		// FRAMERATE new window size
		win_x = w;
		win_y = h;
		
		// Set 2% margin
		double margin = 0.02;
		glViewport(0, 0, w, h);
		glLoadIdentity();
		glOrtho(-margin+center[0]-0.5/scale,center[0]+0.5/scale+margin,
				-margin+center[1]-0.5/scale,center[1]+0.5/scale+margin,-1.0,1.0);
	}
}

static void idle ( void ) {
	if( ! paused ) {
		// Simulate a step
		FLOAT64 maxdt = FRAMERATE ? 1.001*FRAMERATE-fmod(lastTime,FRAMERATE) : 1.0/60.0;
		FLOAT64 curTime = sim.simStep(FRAMERATE,maxdt,FRAMERATE ? FRAMERATE/(FLOAT64)MINSTEP : 0.0015);
		writeFrame(lastTime, curTime, FRAMERATE, DURATION, PATH);
		lastTime = curTime;
	} else {
		if( fabs(scale-goal_scale) < sqr(1e-2)) glutIdleFunc(NULL);
	}
	// Send display command
	reshape(win_x,win_y);
	transition();
	glutPostRedisplay ();
}

static void keyboard( unsigned char key, int x, int y ) {
	cur_x = x;
	cur_y = y;
	// Hidden keys...
	vec2d p( x/(GLdouble)win_x, 1.0 - y/(GLdouble)win_y );
	switch( key ) {
		case '\e':
			exit(0);
			break;
		case '/':
			paused = ! paused;
			if( ! paused ) glutIdleFunc(idle);
			break;
		case '.':
			sim.simStep(1.0,0.01,0.0001);
			glutPostRedisplay();
			break;
		case ',':
			goal_scale = goal_scale != 1.0 ? 1.0 : 10.0;
			goal_center = goal_scale != 1.0 ? p : vec2d(0.5,0.5);
			if( paused ) glutIdleFunc(idle);
			break;
		// Export for SVG file
		case '\\':
			viewer.writeSVG( "/Users/mscp/Desktop/particles.svg", "/Users/mscp/Desktop/" );
			break;
		// Controls for debug
		case ']':
			debugLevel = ! debugLevel;
			break;
		default:
			helpviewer.keyDown(key);
			break;
	}
	viewer.setMousePosition(p);
	if( paused ) glutPostRedisplay ();
}

static void mouse ( int button, int _state, int x, int y ) {
	prev_x = x;
	prev_y = y;
	cur_x = x;
	cur_y = y;
	state = _state;
	viewer.setMousePressed(!state);
	if( paused ) glutPostRedisplay ();
}

static void passive( int x, int y ) {
	cur_x = x;
	cur_y = y;
	vec2d p( x/(GLdouble)win_x, 1.0 - y/(GLdouble)win_y );
	viewer.setMousePosition(p);
	if( paused ) glutPostRedisplay ();
}

static void motion ( int x, int y ) {
	cur_x = x;
	cur_y = y;
	if( state == 0 ) {
		// Compute normalized position
		vec2d p( x/(GLdouble)win_x, 1.0 - y/(GLdouble)win_y );
		vec2d f( (x-prev_x)/(GLdouble)win_x, -(y-prev_y)/(GLdouble)win_y );
		FLOAT64 k = 5.0e+2;
		// Add a mouse driven external force
		sim.addExternalForce(p,k*f);
	}
	prev_x = x;
	prev_y = y;
	passive( x, y );
}

int main (int argc, char * argv[]) {
	// Setup OpenGL stuff..
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GL_DOUBLE);
	glutInitWindowPosition ( 100, 50 );
	glutInitWindowSize ( win_x, win_y );
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
	init( argc>1 ? argv[1] : NULL );
	glutMainLoop();
	return 0;
}
