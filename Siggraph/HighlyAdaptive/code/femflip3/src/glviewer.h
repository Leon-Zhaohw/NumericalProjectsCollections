/*
 *	glviewer.h
 *	
 *	Created by Ryoichi Ando on 1/9/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"

#ifndef _GLVIEWER_H_
#define _GLVIEWER_H_

class flip3;
class glviewer {
public:
	glviewer(const flip3 &sim);
	void drawGL( int width, int height );
	void writeImage( uint step );
protected:
	const flip3& sim;
};

#endif
