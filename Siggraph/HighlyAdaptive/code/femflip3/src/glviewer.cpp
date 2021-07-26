/*
 *	glviewer.cpp
 *	
 *	Created by Ryoichi Ando on 1/9/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "glviewer.h"
#include "opengl.h"
#include "flip3.h"
#include "particle3.h"
#include "write_bmp.h"

glviewer::glviewer(const flip3 &sim) : sim(sim) {
}

static void drawWiredBox() {
	glPushMatrix();
	glTranslatef(0.5,0.5,0.5);
	glutWireCube(1.0);
	glPopMatrix();
}

static void setCamera( int w, int h, int dim ) {
	if( dim == 3 ) {
		// Projection Matrix
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective( 35, w/(double)h, 1, 1000 );
		
		// Model Matrix
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	
		// Set Viewport
		glViewport( 0, 0, w, h );
		
		// Set Rest of Matrix
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glScaled(-1.0, 1.0, 1.0);		
		gluLookAt( -1.0, 1.3, -1.5, 0.5, 0.35, 0.5, 0, 1, 0 );
	} else if( dim == 2 ) {
		// Set 2D coordinate
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0.0,1.0,0.0,1.0,-1.0,1.0);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}
}

static void setLighting( int enabled ) {
	if( enabled ) {
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);
		
		GLfloat light0pos[] = { 0.0, 3.0, 5.0, 1.0 };
		GLfloat light1pos[] = { 5.0, 3.0, 0.0, 1.0 };
		
		glLightfv(GL_LIGHT0, GL_POSITION, light0pos);
		glLightfv(GL_LIGHT1, GL_POSITION, light1pos);
	} else {
		glDisable(GL_LIGHTING);
		glDisable(GL_LIGHT0);
		glDisable(GL_LIGHT1);
	}
}

static void setMaterial( int enabled ) {
	if( enabled ) {
		GLfloat water[] = { 0.4, 0.5, 0.8, 1.0 };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, water);
	}
}

static void drawBitmapString( const char *string, void *font ) {
	if( ! font ) font = GLUT_BITMAP_HELVETICA_12;
	while (*string) glutBitmapCharacter(font, *string++);
}

void glviewer::drawGL( int width, int height ) {
	// Set camera 3D
	setCamera(width,height,3);
	
	// Draw wired box
	glColor4f(1.0,1.0,1.0,1.0);
	drawWiredBox();
	
	// Draw mesh
	glColor4d(1.0,1.0,1.0,0.3);
	sim.surf.drawGL(sim.solid,0);
	sim.solid->draw();
	
	// Draw particles
#if 1
	glPointSize(2.0);
	glBegin(GL_POINTS);	
	FOR_EACH_PARTICLE(sim.particles) {
		if( p.isolated ) {
			glColor4d(1.0,0.5,0.4,0.6);
		} else if( p.levelset > -p.r*sim.dx ) {
			glColor4d(0.4,0.7,1.0,0.6);
		} else {
			glColor4d(0.4,1.0,0.4,0.6);
		}
		if( ! sim.camera->isStatic() && sim.camera->isEnabled()) {
			bool visible = true;
			for( uint n=0; n<sim.objects.size(); n++ ) {
				const object3 &object = *sim.objects[n];
				if( object.type == object3::SOLID ) {
					vec3d camera_pos = sim.camera->origin;
					if( object.intersect(p.p, camera_pos)) {
						visible = false;
						break;
					}
				}
			}
			vec3d cam = sim.camera->getProjectedCoord(p.p);
			if( visible && fabs(cam[0]) < 1.0 && fabs(cam[1]) < 1.0 && cam[2] > 0.0 ) {
				glColor4d(1.0,0.5,0.4,1.0);
			}
		}
		glVertex3dv(p.p.v);
	} END_FOR
	glEnd();
	glPointSize(1.0);
#endif
	
	// Draw mesh
	sim.fluidSolver->render(fluid3::MESH);
	
#if 0
	// Draw simulation camera
	if( sim.camera->isEnabled()) {
		sim.camera->draw();
	}
#endif
	
	// Set camera 2D
	setCamera(width,height,2);
	
	// Draw name
	sim.fluidSolver->render(fluid3::NAME);
	
	// Back to 3D
	setCamera(width,height,3);
}

void glviewer::writeImage( uint step ) {
	int width = glutGet(GLUT_WINDOW_WIDTH);
	int height = glutGet(GLUT_WINDOW_HEIGHT);
	unsigned char *buffer = new unsigned char[width*height*4];
	
	glFinish();
	glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
	glReadPixels( 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, buffer );
	
	if( ! is_exist("%s/screenshots/", root_path)) {
		system(format_str("mkdir %s/screenshots/", root_path));
	}
	
    // ffmpeg -qscale 4 -r 60 -b 9600 -i frame_%d.bmp out.mp4
	write_bmp( format_str("%s/screenshots/frame_%d.bmp", root_path, step), buffer, width, height, true );
	
	delete [] buffer;
}
