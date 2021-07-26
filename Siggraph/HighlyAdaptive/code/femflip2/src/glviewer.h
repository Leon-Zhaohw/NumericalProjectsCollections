/*
 *	glviewer.h
 *	
 *	Created by Ryoichi Ando on 11/5/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "vec2.h"
#include "array2.h"
#include "extsurf2.h"

#ifndef _GLVIEWER_H
#define _GLVIEWER_H

class flip2;
class levelset2;
class glviewer {
public:
	glviewer(const flip2 &sim);
	void drawGL();	// Drawn in { (0,0) - (1,1) } coordinate
	void writeImage( const char *path );
	void writeSVG( const char *path, const char *dir );
	
	// View options...
	void setGridlineVisibility( bool visible );
	void setGridVelocityVisibility( bool visible );
	void setLevelsetVisibility( bool visible );
	void setVolumeFractionVisibility( bool visible );
	void setParticleVisibility( bool visible );
	void setNeighborConnectionVisibility( bool visible );
	void setPressureVisibility( bool visible );
	void setSimTimeVisibility( bool visible );
	void setParticleVelocityVisibility( bool visible );
	void setParticleAnisotropyVisibility( bool visible );
	void setMatrixConnectionVisibility( bool visible );
	void setMeshGeneratorVisibility( bool visible );
	void setDivergenceVisibility( bool visible );
	
	void setMousePosition( vec2d p );
	void setMousePressed( bool pressed );
	void setSurfaceExtractorMethod( int method );
	
	// Public utility functions
	static void drawBitmapString( const char *string, void *font=NULL );
protected:
	const flip2& sim;
	extsurf2 extsurf;
	vec2d mouse;
	bool pressed;
	bool showGridline;
	bool showGridVelocity;
	bool showLevelset;
	bool showVolumeFraction;
	bool showParticles;
	bool showNeighbors;
	bool showPressure;
	bool showDivergence;
	bool showSimTime;
	bool showParticleVelocity;
	bool showAnisotropy;
	bool showMatrixConnection;
	bool showExternalMeshGen;
	int lastFrame;
};

#endif