#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <eigen/Dense>
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include <cassert>

#ifndef M_PI
#define	M_PI		3.14159265358979323846	/* pi */
#endif

double size = 0.5;

class Light {
public:
  GLfloat ambient[4]; // ambient color, RGBA
  GLfloat diffuse[4]; // diffuse color, RGBA
  GLfloat specular[4]; // specular color, RGBA
  GLfloat pos[4]; // light position, XYZW
  GLenum id; // light identifier
  Light() {}; // constructor
  void apply() {
	glLightfv(id, GL_AMBIENT, ambient);
	glLightfv(id, GL_DIFFUSE, diffuse);
	glLightfv(id, GL_SPECULAR, specular);
	glLightfv(id, GL_POSITION, pos);
	glEnable(id);
  }
};

// Global Variables
std::vector<Light> lights;
int vWidth, vHeight;
int currentFrame;
bool singleStep;
std::vector<std::vector<Eigen::Vector3d> > frames;
char *framestring;

void readFrame(const char *fname, std::vector<Eigen::Vector3d> &particles) {
  char ch;
  Eigen::Vector3d p;
  int i;
  int nparts;
  double foo;
  
  std::ifstream in(fname, std::ios::in);
  std::cout<<"reading "<<fname<<std::endl;
  if (!in.good()) return;
  std::cout<<"good"<<std::endl;
  in>>nparts;
  while(!in.eof()) {
	in>>p[0]>>p[1]>>p[2]>>foo>>foo>>foo;
	if (!in.eof()) particles.push_back(p);
  }
  in.close();
  std::cout<<"inputfile read"<<std::endl;
}

///////////////////////////////////////////////////
// Begin Class Function Definitions
///////////////////////////////////////////////////
void draw() {
  std::cout<<currentFrame<<" / "<<frames.size()<<std::endl;

  if (currentFrame >= frames.size()) {
	if (currentFrame > frames.size()) currentFrame = frames.size();
	char fname[80];
	sprintf(fname, framestring, currentFrame);
	std::vector<Eigen::Vector3d> particles;
	readFrame(fname, particles);
	if (particles.size() > 0) frames.push_back(particles);
	else currentFrame = 0;
  }

  glPushMatrix();
  GLfloat c[4] = {0.0, 1.0, 0.0, 1.0};
  glMaterialfv(GL_FRONT, GL_AMBIENT, c);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, c);
  glMaterialfv(GL_FRONT, GL_SPECULAR, c);

  //glBegin(GL_POINTS);
  for (unsigned int i=0; i<frames[currentFrame].size(); i++) {
	glPushMatrix();
	glTranslatef(frames[currentFrame][i][0], frames[currentFrame][i][1], frames[currentFrame][i][2]);
	glutSolidSphere(0.0025, 10, 10);
	glPopMatrix();
	//glVertex3f(frames[currentFrame][i][0], frames[currentFrame][i][1], frames[currentFrame][i][2]);
  }

  //glEnd();
  
  glPopMatrix();
}

///////////////////////////////////////////////////
// End Class Definitions
///////////////////////////////////////////////////

void myReshape(int w, int h) {
  glViewport (0,0,w,h);
  vWidth = w;
  vHeight = h;
}

void myDisplay() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST); /* enable depth testing; required for z-buffer */
  glEnable(GL_CULL_FACE); /* enable polygon face culling */ 
  glCullFace(GL_BACK); /* tell opengl to cull the back face*/ 
  glEnable(GL_NORMALIZE);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60, ((float)vWidth)/vHeight, 0.01, 20);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0.0, 0.5, 0.0,
	  0.0, 0.0, 0.0,
	  0.0, 0.0, 1.0);
  for (int i=0; i<lights.size(); i++) lights[i].apply();
  draw();

  glutSwapBuffers();
}

void mySpecial(int key, int x, int y) {
  switch(key) {
  case GLUT_KEY_LEFT:
	if (singleStep) {
	  currentFrame--;
	  if (currentFrame < 0) currentFrame = 0;
	  std::cout<<"frame: "<< currentFrame<<std::endl;
	}
	break;
  case GLUT_KEY_RIGHT:
	if (singleStep) {
	  currentFrame++;
	  std::cout<<"frame: "<< currentFrame<<std::endl;
	}
	break;
  default:
	std::cerr<<"That key is not recognized"<<std::endl;
	break;
  }
  glutPostRedisplay();
}

void myKeyboard(unsigned char key, int x, int y) {
  int i;
  switch(key) {
  case 'q':
  case 'Q':
	exit(0);
  case 'p':
  case 'P':
	singleStep = false;
	break;
  case 's':
  case 'S':
	singleStep = true;
	break;
  default:
	std::cerr<<"key "<<key<<" not supported"<<std::endl;
  }
  glutPostRedisplay();
}

void myTimerFunc(int id) {
  if (!singleStep) {
	currentFrame++;
	//if (currentFrame > frames.size()-1) currentFrame = 0;
  }
  glutTimerFunc(33, myTimerFunc, 0);
  glutPostRedisplay();
}

int main(int argc, char *argv[]) {
  glutInit(&argc, argv);
  framestring = argv[1];
  
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(1440, 810);
  glutInitWindowPosition(0,0);
  glutCreateWindow(argv[0]);
  
  glutDisplayFunc(myDisplay);
  glutReshapeFunc(myReshape);
  glutKeyboardFunc(myKeyboard);
  glutSpecialFunc(mySpecial);
  glutTimerFunc(1000, myTimerFunc, 0);
  
  for (int i=0; i<1; i++) {
	Light l;
	l.ambient[0] = 0.2; l.ambient[1] = 0.2; l.ambient[2] = 0.2; l.ambient[3] = 0.2;
	l.diffuse[0] = 0.7; l.diffuse[1] = 0.7; l.diffuse[2] = 0.7; l.diffuse[3] = 0.7;
	l.specular[0] = 0.4; l.specular[1] = 0.4; l.specular[2] = 0.4; l.specular[3] = 0.4;
	switch(i) {
	case 0:
	  l.pos[0] = 0.0; l.pos[1] = 0.0; l.pos[2] = 1.0; l.pos[3] = 1.0;
	  break;
	case 1:
	  l.pos[0] = 1.0; l.pos[1] = 0.0; l.pos[2] = 0.0; l.pos[3] = 1.0;
	  break;
	case 2:
	  l.pos[0] = -1.0; l.pos[1] = 0.0; l.pos[2] = 0.0; l.pos[3] = 1.0;
	  break;
	case 3:
	  l.pos[0] = 0.0; l.pos[1] = 0.0; l.pos[2] = -1.0; l.pos[3] = 1.0;
	  break;
	case 4:
	  l.pos[0] = 0.0; l.pos[1] = 0.0; l.pos[2] = 1.0; l.pos[3] = 1.0;
	  break;
	}
	l.id = GL_LIGHT0+i;
	lights.push_back(l);
  }
  
  singleStep = false;
  currentFrame = 0;

  glutMainLoop();
  
  return 0;
}

