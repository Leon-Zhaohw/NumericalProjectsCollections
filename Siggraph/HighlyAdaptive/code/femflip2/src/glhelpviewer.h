/*
 *	glhelpviewer.h
 *	
 *	Created by Ryoichi Ando on 11/7/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <vector>
#include <string>

#ifndef _GLHELPVIEWER_H
#define _GLHELPVIEWER_H

#define MAX_SLOT	64
class glviewer;
class flip2;
class glhelpviewer {
public:
	glhelpviewer( glviewer& viewer, flip2& sim );
	void applyStates( bool forceReset=false );
	void restoreDefaults();
	bool keyDown( char key );
	void drawGL();
protected:
	glviewer& viewer;
	flip2& sim;
	
	std::vector<bool> visibilityStates;
	std::vector<bool> default_visibilityStates;
	
	std::vector<int> controllerStates;
	std::vector<int> default_controllerStates;
	
	std::vector<std::string> controllerStatusStr;
	void setControllerString();
private:
	void resetSim();
};

#endif