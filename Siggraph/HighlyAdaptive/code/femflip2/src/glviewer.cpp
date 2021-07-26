/*
 *	glviewer.cpp
 *	
 *	Created by Ryoichi Ando on 11/5/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "flip2.h"
#include "opengl.h"
#include "glviewer.h"
#include "levelset2.h"
#include "particle2.h"
#include "util2.h"
#include "macfluid2.h"
#include "write_bmp.h"
#include "octhash2.h"
#include "femfluid2.h"
#include "extsurf2.h"
#include <math.h>
#include <vector>
using namespace std;

glviewer::glviewer(const flip2 &sim) : sim(sim) {
	// Feel free to manually set these value to visualize
	showGridline = true;
	showGridVelocity = false;
	showVolumeFraction = false;
	showLevelset = true;
	showParticles = true;
	showNeighbors = false;
	showPressure = false;
	showSimTime = true;
	showParticleVelocity = false;
	showAnisotropy = false;
	showMatrixConnection = false;
	showExternalMeshGen = false;
	showDivergence = false;
	pressed = false;
	lastFrame = -1;
}

void glviewer::setGridlineVisibility( bool visible ) {
	showGridline = visible;
}

void glviewer::setGridVelocityVisibility( bool visible ) {
	showGridVelocity = visible;
}

void glviewer::setLevelsetVisibility( bool visible ) {
	showLevelset = visible;
}

void glviewer::setVolumeFractionVisibility( bool visible ) {
	showVolumeFraction = visible;
}

void glviewer::setParticleVisibility( bool visible ) {
	showParticles = visible;
}

void glviewer::setNeighborConnectionVisibility( bool visible ) {
	showNeighbors = visible;
}

void glviewer::setPressureVisibility( bool visible ) {
	showPressure = visible;
}

void glviewer::setDivergenceVisibility( bool visible ) {
	showDivergence = visible;
}

void glviewer::setSimTimeVisibility( bool visible ) {
	showSimTime = visible;
}

void glviewer::setParticleVelocityVisibility( bool visible ) {
	showParticleVelocity = visible;
}

void glviewer::setParticleAnisotropyVisibility( bool visible ) {
	showAnisotropy = visible;
}

void glviewer::setMatrixConnectionVisibility( bool visible ) {
	showMatrixConnection = visible;
}

void glviewer::setMeshGeneratorVisibility( bool visible ) {
	showExternalMeshGen = visible;
}

void glviewer::setMousePosition( vec2d p ) {
	mouse = p+vec2d(0.01,0.01);
}

void glviewer::drawBitmapString( const char *string, void *font ) {
	if( ! font ) font = GLUT_BITMAP_HELVETICA_12;
	while (*string) glutBitmapCharacter(font, *string++);
}

void glviewer::setMousePressed( bool pressed ) {
	this->pressed = pressed;
}

void glviewer::setSurfaceExtractorMethod( int method ) {
	extsurf.setMethod(method);
}

void glviewer::drawGL() {	
	// Gather simulator info
	int winsize[2] = { glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) };
	uint gn = sim.gn;
	FLOAT64 dx = 1.0/gn;
	
	const vector<particle2 *> &particles = sim.particles;
	const fluid2 *fluidSolver = sim.fluidSolver;
	
	if( showLevelset ) {
		// Draw solid levelset
		fluidSolver->render(fluid2::SOLID,mouse);
		sim.surf.drawGL(sim.solid,0);
	}
	
	if( showPressure ) {
		// Draw pressure
		fluidSolver->render(fluid2::PRESSURE,mouse);
	} else if( showDivergence ) {
		// Draw divergence
		fluidSolver->render(fluid2::DIVERGENCE,mouse);
	}
	
	if( showGridline ) {// && ! showExternalMeshGen ) {
		fluidSolver->render(fluid2::MESH,mouse);
	}
	
	if( showExternalMeshGen ) {
		const octree2 *octree = fluidSolver->getOctree();
		extsurf.loadParticles(sim.particles,sim.dpx,octree);
		extsurf.generateMesh(sim.solid);
		glColor4d(1.0,1.0,1.0,0.2);
		extsurf.draw();
	}
	
	if( showGridVelocity ) {
		// Draw flow
		fluidSolver->render(fluid2::VELOCITY,mouse);
	}
	
	if( showMatrixConnection ) {
		fluidSolver->render(fluid2::MATRIX_CONNECTION,mouse);
	}

	if( showParticles ) {
		glPointSize(2.0);
		glBegin(GL_POINTS);
		FOR_EACH_PARTICLE(particles) {
			if( p.isolated ) glColor4d(1.0,0.2,0.2,0.9);
			else if( p.levelset < -p.r*sim.dpx ) glColor4d(0.4,1.0,0.4,0.9);
			else glColor4d(0.4,0.7,1.0,1.0);
			glVertex2d(p.p.v[0],p.p.v[1]);
		} END_FOR
		glEnd();
	}
	
	if( showParticles && showAnisotropy ) {
		FOR_EACH_PARTICLE(particles) {
			if( p.isolated ) glColor4d(1.0,0.2,0.2,0.9);
			else if( p.levelset < -p.r*sim.dpx ) glColor4d(0.4,1.0,0.4,0.9);
			else glColor4d(0.4,0.7,1.0,1.0);			
			uint dth = 20;
			FLOAT64 r = 0.5*sim.dpx*p.r;
			glBegin(GL_LINE_LOOP);
			for( int t=0; t<dth; t++ ) {
				FLOAT64 theta = 2.0*t*PI/dth;
				vec2d rpos( r*cos(theta),r*sin(theta) );
				glVertex2f(rpos[0]+p.p[0],rpos[1]+p.p[1]);
			}
			glEnd();
		} END_FOR
	}
	
	if( showLevelset && ! showExternalMeshGen ) {
		sim.surf.drawGL(sim.solid,1);
	}

	if( showParticleVelocity ) {
		glColor4d(1.0,1.0,1.0,0.6);
		FOR_EACH_PARTICLE(particles) {
			if( p.n > 1 ) {
				glBegin(GL_LINES);
				glVertex2f(p.p[0],p.p[1]);
				glVertex2f(p.p[0]+dx*p.u[0],p.p[1]+dx*p.u[1]);
				glEnd();
			}
		} END_FOR
	}
	
	if( showNeighbors ) {
		// Draw Neighbors
		glColor4d(1.0,1.0,1.0,0.6);
		glBegin(GL_LINES);
		FLOAT64 r = sim.dpx*powf(2,sim.discreteDepthLevelset->evalLevelset(mouse)-1);
		std::vector<particle2 *> neighbors = sim.sorter->getNeighbors(mouse,r);
		for( std::vector<particle2 *>::iterator it=neighbors.begin(); it!=neighbors.end(); it++ ) {
			glVertex2f(mouse[0],mouse[1]);
			glVertex2f((*it)->p[0],(*it)->p[1]);
		}
		glEnd();
		
		uint dth = 20;
		glBegin(GL_LINE_LOOP);
		for( int t=0; t<dth; t++ ) {
			FLOAT64 theta = 2.0*t*PI/dth;
			glVertex2f(r*cos(theta)+mouse[0],r*sin(theta)+mouse[1]);
		}
		glEnd();
	}

	if( showSimTime ) {
		glColor4d(1.0,1.0,1.0,1.0);
		glRasterPos2d((winsize[0]-170)/(double)winsize[0],(winsize[1]-40)/(double)winsize[1]);
		drawBitmapString(format_str("Time: %.2f msec", sim.msec),GLUT_BITMAP_HELVETICA_18);
	}
	
#if 0
	// Draw camera
	sim.camera.draw();
#endif
}

void glviewer::writeSVG( const char *path, const char *dir ) {
	FILE *fp;
	if( ! (fp = fopen( path, "w" ))) exit(1);
	
	// Write SVG header
	fprintf( fp, "<?xml version=\"1.0\" standalone=\"no\"?>\n" );
	fprintf( fp, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n" );
	fprintf( fp, "<svg width=\"30cm\" height=\"30cm\" viewBox=\"0 0 1 1\"\n");
	fprintf( fp, "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n" );
	const char *colors[] = { "mediumslateblue", "orange", "orangered", "mediumseagreen", "aqua" };
	
	// Generate a surface
	const octree2 *octree = sim.fluidSolver->getOctree();
	extsurf.loadParticles(sim.particles,sim.dpx,octree);
	extsurf.generateMesh(sim.solid);
	
	// Fetch geometry
	const std::vector<vec2d> &vertices = extsurf.vertices;
	const std::vector<vec2i> &faces = extsurf.faces;
	
	// Compute connections
	std::vector<std::vector<uint> > connections;
	connections.resize(vertices.size());
	for( uint n=0; n<faces.size(); n++ ) {
		connections[faces[n][0]].push_back(faces[n][1]);
		connections[faces[n][1]].push_back(faces[n][0]);
	}
	std::vector<bool> picked(vertices.size());
	for( uint n=0; n<vertices.size(); n++ ) {
		picked[n] = false;
	}
	while( true ) {
		int unpicked = -1;
		for( uint n=0; n<vertices.size(); n++ ) {
			if( ! picked[n] ) {
				unpicked = n;
				break;
			}
		}
		if( unpicked == -1 ) break;
		std::vector<uint> points;
		points.push_back(unpicked);
		uint previous = unpicked;
		uint current = connections[unpicked][0];
		picked[previous] = true;
		picked[current] = true;
		while(true) {
			bool found = false;
			if( connections[current].size() == 2 ) {
				for( uint k=0; k<2; k++ ) {
					uint next = connections[current][k];
					if( next != previous ) {
						if( next == unpicked ) break;
						previous = current;
						points.push_back(next);
						picked[next] = true;
						current = next;
						found = true;
						break;
					}
				}
			}
			if( ! found ) {
				break;
			}
		}
		// Write this polygon
		fprintf(fp, "<polygon fill=\"rgb(190,210,255)\" stroke=\"rgb(40,40,240)\" stroke-width=\"%f\"\n", 0.1*sim.dpx);
		fprintf(fp, "points= \"" );
		for( uint n=1; n<points.size(); n++ ) {
			fprintf( fp, "%lf,%lf ", vertices[points[n]][0], 1.0-vertices[points[n]][1] );
		}
		fprintf(fp, "\"" );
		fprintf(fp,"/>\n");
	}
	
	// Write mesh
	if( sim.fluidSolver->getMesh()) {
		// Draw a mesh
		const mesher2 &g = *sim.fluidSolver->getMesh();
		for( uint n=0; n<g.node2node.size(); n++) {
			uint idx0 = n;
			for( uint m=0; m<g.node2node[n].size(); m++) {
				uint idx1 =g.node2node[n][m];
				if( idx0 < idx1 ) {
					vec2d p0 = g.nodes[idx0];
					vec2d p1 = g.nodes[idx1];
					fprintf( fp,"<line x1=\"%lf\" y1=\"%lf\" x2=\"%lf\" y2=\"%lf\" style=\"stroke:rgb(60,60,60);stroke-width:%lf\"/>\n",
							p0[0], 1.0-p0[1], p1[0], 1.0-p1[1], 0.05*sim.dpx );
				}
			}
		}
	} else if( sim.fluidSolver->getOctree()) {
		const octree2 &octree = *sim.fluidSolver->getOctree();
		for( uint n=0; n<octree.terminals.size(); n++ ) {
			octree2::leaf2 *leaf = octree.terminals[n];
			FLOAT64 dx = leaf->dx/(FLOAT64)octree.resolution;
			vec2d center = vec2d(leaf->center)/(FLOAT64)octree.resolution;
			fprintf(fp, "<polygon fill=\"none\" stroke=\"rgb(60,60,60)\" stroke-width=\"%f\"\n", 0.05*sim.dpx);
			fprintf(fp, "points= \"" );
			fprintf( fp, "%lf,%lf ", center[0]-0.5*dx,1.0-(center[1]-0.5*dx) );
			fprintf( fp, "%lf,%lf ", center[0]-0.5*dx,1.0-(center[1]+0.5*dx) );
			fprintf( fp, "%lf,%lf ", center[0]+0.5*dx,1.0-(center[1]+0.5*dx) );
			fprintf( fp, "%lf,%lf ", center[0]+0.5*dx,1.0-(center[1]-0.5*dx) );
			fprintf(fp, "\"" );
			fprintf(fp,"/>\n");
		}
	} else {
		// Draw a uniform grid
		FLOAT64 dx = 1.0/(sim.gn);
		for( uint i=0; i<sim.gn+1; i++ ) {
			fprintf( fp,"<line x1=\"%lf\" y1=\"%lf\" x2=\"%lf\" y2=\"%lf\" style=\"stroke:rgb(60,60,60);stroke-width:%lf\"/>\n",
					i*dx, 0.0, i*dx, 1.0, 0.05*sim.dpx );
		}
		for( uint j=0; j<sim.gn+1; j++ ) {
			fprintf( fp,"<line x1=\"%lf\" y1=\"%lf\" x2=\"%lf\" y2=\"%lf\" style=\"stroke:rgb(60,60,60);stroke-width:%lf\"/>\n",
					0.0, j*dx, 1.0, j*dx, 0.05*sim.dpx );
		}
	}
	
	// Write particles
	for( uint n=0; n<sim.particles.size(); n++ ) {
		particle2 &p = *sim.particles[n];
		fprintf( fp, "<circle cx=\"%lf\" cy=\"%lf\" fill=\"%s\" stroke=\"none\" r=\"%lf\"/>\n",
				p.p[0], 1.0-p.p[1], (p.n>1 ? colors[3] : colors[0]), 0.4*sim.dpx*p.r );
	}	
	// Write SVG footer
	fprintf( fp, "</svg>\n" );
	fclose(fp);
	if( dir ) {
		run( "qlmanage -t -s 1280 -o %s %s", dir, path);
		run( "convert %s.png -crop 1280x720+0+560 %s.png", path, path );
	}
}

void glviewer::writeImage( const char *path ) {
	int width = glutGet(GLUT_WINDOW_WIDTH);
	int height = glutGet(GLUT_WINDOW_HEIGHT);
	unsigned char *buffer = new unsigned char[width*height*4];
	
	glFlush();
	glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
	glReadPixels( 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, buffer );
	write_bmp( path, buffer, width, height, true );
	
	delete [] buffer;
}
