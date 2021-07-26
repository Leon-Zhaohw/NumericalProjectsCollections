#include <GLUT/glut.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "fastlevelset.h"

#ifndef M_PI
#define M_PI 3.1415926535897
#endif

#define sqr(x) ((x)*(x))

const char *writefile = "sculpt.fls";

void handleDisplay (void);
void handleReshape (int w, int h);
void handleKeypress (unsigned char key, int x, int y);
void handleClick (int button, int state, int x, int y);
void handleDrag (int x, int y);
void handleIdle (void);
void initLighting();
void initMaterial();
void setLighting();
void drawAxes();
void drawTarget();
void generateRay (int x, int y, float *ox, float *oy, float *oz, float *dx, float *dy, float *dz);
float sphere (float x, float y, float z);
float weirdphi (float x, float y, float z);


FastLevelSet fls(30,sphere);
float tx = .5, ty = .5, tz = .5;
int width = 512, height = 512;
int oldmousex = 0, oldmousey = 0, mousebutton = 0;
enum {NONE, MOUSETOOL, MOUSEROTATE, MOUSEZOOM, MOUSEMOVEXY, MOUSEMOVEXZ} mousemode = NONE;
bool mouseswitch=false;
float heading = 0, pitch = 0, camdist=2;
float lookx=0, looky=0, lookz=-1;
float nearclip = .0625, farclip = 4;
int draw_axes = 1, draw_target = 0;
float tool[3];
bool tool_on = false;
float toolradius = .5*fls.maxToolSize();
float toolamount = .0005;
enum {TOOLADD, TOOLSMOOTH} tooltype = TOOLADD;



int main (int argc, char **argv)
{
	glutInit (&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize (width,height);
	glutCreateWindow ("sculpt");
	glutDisplayFunc (handleDisplay);
	glutReshapeFunc (handleReshape);
	glutKeyboardFunc (handleKeypress);
	glutMouseFunc (handleClick);
	glutMotionFunc (handleDrag);
	glutIdleFunc (handleIdle);

	if (argc > 1) fls.read (argv[1]);

	glEnable (GL_DEPTH_TEST);
	glCullFace (GL_BACK);
	glClearColor (0, 0, 0, 0);
	glClearDepth (1);
	initLighting();
	initMaterial();

	glutMainLoop ();
}


void initLighting()
{
	glEnable (GL_LIGHTING);
	{
		GLfloat ambient[4] = {.3, .3, .3, 1};
		glLightModelfv (GL_LIGHT_MODEL_AMBIENT,ambient);
	}
	{
		GLfloat color[4] = {.8, .8, .8, 1};
		glLightfv (GL_LIGHT0, GL_DIFFUSE, color);
		glLightfv (GL_LIGHT0, GL_SPECULAR, color);
		glEnable (GL_LIGHT0);
	}
	{
		GLfloat color[4] = {.4, .4, .4, 1};
		glLightfv (GL_LIGHT1, GL_DIFFUSE, color);
		glLightfv (GL_LIGHT1, GL_SPECULAR, color);
		glEnable (GL_LIGHT1);
	}
	{
		GLfloat color[4] = {.2, .2, .2, 1};
		glLightfv (GL_LIGHT2, GL_DIFFUSE, color);
		glLightfv (GL_LIGHT2, GL_SPECULAR, color);
		glEnable (GL_LIGHT2);
	}
}


void setLighting()
{
	{
		GLfloat position[4] = {.3, 1, 1, 0};
		glLightfv (GL_LIGHT0, GL_POSITION, position);
	}
	{
		GLfloat position[4] = {-1, .2, .5, 0};
		glLightfv (GL_LIGHT1, GL_POSITION, position);
	}
	{
		GLfloat position[4] = {1, -1, .05, 0};
		glLightfv (GL_LIGHT2, GL_POSITION, position);
	}
}


void handleDisplay()
{
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity();
	{
		float aspect = width/(float)height;
		glFrustum (-.3*aspect*nearclip, .3*aspect*nearclip, -.3*nearclip, .3*nearclip, nearclip, farclip);
	}

	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity();
	setLighting();
	gluLookAt (tx-camdist*lookx, ty-camdist*looky, tz-camdist*lookz, tx, ty, tz, 0, 1, 0);

	if (draw_axes) drawAxes();
	if (draw_target) drawTarget();
	fls.glRender();

	if (tool_on) {
		glMatrixMode (GL_MODELVIEW);
		glPushMatrix();
		glTranslatef (tool[0], tool[1], tool[2]);
		glDisable (GL_LIGHTING);
		if (tooltype == TOOLADD) glColor3f (1, 0, 0); else glColor3f (0, 1, 0);
		glutWireSphere (toolradius, 8, 5);
		glEnable (GL_LIGHTING); 
		glPopMatrix();
	}
	glutSwapBuffers();
}


void generateRay (int x, int y, float orig[3], float dir[3])
{
	float ch = cos(heading), sh = sin(heading), cp = cos(pitch), sp = sin(pitch);
	float mx = .6*nearclip*width/(float)height*((x/(float)width)-.5), my = .6*nearclip*(.5-y/(float)height), mz = -nearclip;
	float overmag;

	orig[0] = camdist*sh*cp + tx;
	orig[1] = camdist*sp + ty;
	orig[2] = camdist*ch*cp + tz;

	dir[0] = (mz*cp-my*sp)*sh + mx*ch;
	dir[1] = mz*sp + my*cp;
	dir[2] = (mz*cp-my*sp)*ch - mx*sh;
	overmag = 1./sqrt(sqr(dir[0])+sqr(dir[1])+sqr(dir[2]));
	dir[0] *= overmag; dir[1] *= overmag; dir[2] *= overmag;
}


void handleReshape (int w, int h)
{
	width = w;
	height = h;
	glViewport (0, 0, (GLsizei)width, (GLsizei)height);
	glutPostRedisplay();
}


void handleKeypress (unsigned char key, int x, int y)
{
	switch (key) {
		case 'q': std::exit(0);
		case 'r': fls.refine(); toolradius *= .5; toolamount *= .5; glutPostRedisplay(); break;
		case 's': tooltype = TOOLSMOOTH; break;
		case 'a': tooltype = TOOLADD; break;
		case 'w': fls.write (writefile); break;
		case '[': toolamount *= .8; break;
		case ']': toolamount *= 1.25; break;
		case '-': toolradius *= .8; break;
		case '=': toolradius *= 1.25; if (toolradius > fls.maxToolSize()) toolradius = fls.maxToolSize(); break;
		case 'd': fls.redistance(); glutPostRedisplay(); break;
		case 'e': fls.enforceSharedValues(); glutPostRedisplay(); break;
		case 'S': fls.smoothRefined(); glutPostRedisplay(); break;
		case 'm': mouseswitch=!mouseswitch; break;
	}
}


void handleClick (int button, int state, int x, int y)
{
	int mods = glutGetModifiers();

	if (state == GLUT_UP)				 {mousemode = NONE;if (tool_on) {tool_on = false;glutPostRedisplay();}}
	else if (mods & GLUT_ACTIVE_SHIFT) {
		if (mouseswitch) {
			if (button == GLUT_LEFT_BUTTON)	 mousemode = MOUSEMOVEXY;
			else							 mousemode = MOUSEMOVEXZ;
		} else {
			if (button == GLUT_LEFT_BUTTON)	 mousemode = MOUSEROTATE;
			else							 mousemode = MOUSEZOOM;
		}
	} else {
		mousemode = MOUSETOOL;
		float orig[3], dir[3];
		generateRay (x, y, orig, dir);
		if (fls.traceRay (orig, dir, tool)) {
			tool_on = true;
			if (tooltype == TOOLADD) fls.localAdd (tool, toolradius, toolamount*(button==GLUT_LEFT_BUTTON ? 1 : -1));
			else fls.localSmooth (tool, toolradius, toolamount*(button==GLUT_LEFT_BUTTON ? 1 : -1));
		} else
			tool_on = false;
		glutPostRedisplay();
	}

	oldmousex = x;
	oldmousey = y;
	mousebutton = button;
}


void handleDrag (int x, int y)
{
	switch (mousemode) {
		case MOUSEROTATE:
		{
			float ch, sh, cp, sp;

			heading += 0.007*(oldmousex-x);
			if (heading < -M_PI) heading += 2*M_PI;
			else if (heading > M_PI) heading -= 2*M_PI;

			pitch += 0.007*(y-oldmousey);
			if (pitch < -.45*M_PI) pitch = -.45*M_PI;
			else if (pitch > .45*M_PI) pitch = .45*M_PI;

			ch = cos(heading); sh = sin(heading);
			cp = cos(pitch);   sp = sin(pitch);
			lookx = -sh*cp; looky = -sp; lookz = -ch*cp;
		} break;
		case MOUSEZOOM:
		{
			double factor = std::pow(1.01, y-oldmousey);
			camdist *= factor;
			nearclip = .0625*camdist;
			farclip = 4*camdist;
		} break;
		case MOUSEMOVEXY:
		{
			double dx = 0.002*(oldmousex-x);
			double dy = 0.002*(y-oldmousey);
			double ch = cos(heading), sh = sin(heading);
			tx += camdist*ch*dx;
			ty += camdist*dy;
			tz -= camdist*sh*dx;
		} break;
		case MOUSEMOVEXZ:
		{
			double dx = 0.002*(oldmousex-x);
			double dz = 0.002*(y-oldmousey);
			double ch = cos(heading), sh = sin(heading);
			tx += camdist*(ch*dx-sh*dz);
			tz += camdist*(-ch*dz-sh*dx);
		} break;
		case MOUSETOOL:
		{
			float orig[3], dir[3];
			generateRay (x, y, orig, dir);
			if (fls.traceRay (orig, dir, tool)) {
				tool_on = true;
				if (tooltype == TOOLADD) fls.localAdd (tool, toolradius, toolamount*(mousebutton==GLUT_LEFT_BUTTON ? 1 : -1));
				else fls.localSmooth (tool, toolradius, toolamount*(mousebutton==GLUT_LEFT_BUTTON ? 1 : -1));
			} else
				tool_on = false;
		} break;
		default: ;
	} 
	oldmousex = x;
	oldmousey = y;
	glutPostRedisplay();
}


void handleIdle ()
{
	if (tool_on) {
		float orig[3], dir[3];
		generateRay (oldmousex, oldmousey, orig, dir);
		if (fls.traceRay (orig, dir, tool)) {
			if (tooltype == TOOLADD) fls.localAdd (tool, toolradius, toolamount*(mousebutton==GLUT_LEFT_BUTTON ? 1 : -1));
			else fls.localSmooth (tool, toolradius, toolamount*(mousebutton==GLUT_LEFT_BUTTON ? 1 : -1));
		} else
			tool_on = false;
		glutPostRedisplay();
	}
}

void drawTarget ()
{
	glMatrixMode (GL_MODELVIEW);
	glPushMatrix();
	glTranslatef (tx, ty, tz);
	glDisable (GL_LIGHTING);
	glColor3f (1, 1, 1);
	glutWireCube (.1);
	glEnable (GL_LIGHTING); 
	glPopMatrix();
}


void drawAxes ()
{
	int i;

	glDisable (GL_LIGHTING);

	glBegin (GL_LINES);
	glColor3f (1, .25, .25);
	for (i = 0; i < 10; ++i) {
		glVertex3f (0, 0, .1*i);
		glVertex3f (1, 0, .1*i);
	}
	glColor3f (.25, 1, .25);
	glVertex3f (0, 0, 0);
	glVertex3f (0, 1, 0);
	glColor3f (.25, .25, 1);
	for (i = 0; i < 10; ++i) {
		glVertex3f (.1*i, 0, 0);
		glVertex3f (.1*i, 0, 1);
	}
	glEnd();

	glEnable (GL_LIGHTING); 
}


void initMaterial()
{
	GLfloat ambient[4] = {.7, .7, 0};
	GLfloat diffuse[4] = {.6, .8, 1};
	GLfloat specular[4] = {0, .2, .3};

	glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
	glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
	glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, specular);
	glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, 50);
}


float sphere (float x, float y, float z)
{
	return sqrt(sqr(x-.501234)+sqr(y-.4999923)+sqr(z-.500234))-.4013214;
}


float weirdphi (float x, float y, float z)
{
	float sphere1 = sqrt(sqr(x-.5)+sqr(y-.5)+sqr(z-.5))-.332;
	float sphere2 = sqrt(sqr(x-.25)+sqr(y-.4)+sqr(z-.4))-.21;
	float deathstar = (sphere1 > -sphere2) ? sphere1 : -sphere2;
	float cylinder = sqrt(sqr(x-.46)+sqr(z-.35))-.1;
	return ((deathstar > -cylinder) ? deathstar : -cylinder);
}


