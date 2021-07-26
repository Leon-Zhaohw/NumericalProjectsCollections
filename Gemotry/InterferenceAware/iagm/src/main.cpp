// main.cpp
//

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <iostream>

#include "visualizer.h"
#include "state.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <AntTweakBar.h>



#define WindowW 1000
#define WindowH 750

Visualizer *visualizer;
int tw_my_var = 0;
bool movierunning;
int moviecount;

void display ()
{
    visualizer->paintGL ();
    TwDraw ();
    glutSwapBuffers ();
}

void reshape (int w, int h)
{
    visualizer->resizeGL (w, h);
    TwWindowSize (w, h);
    glutPostRedisplay ();
}
 
void mouse_click (int button, int state, int x, int y)
{
    if (TwEventMouseButtonGLUT (button, state, x, y))
    {
		glutPostRedisplay ();
        return;
    }

	visualizer->mousePress(button, state, x, y);
    display();
	glutPostRedisplay ();
}

void mouse_move (int x, int y)
{
    if (TwEventMouseMotionGLUT (x, y))
    {
		//glutPostRedisplay ();
        display();
        return;
    }
	
    visualizer->mouseMove (x, y);
    display();
}

void keyboard (unsigned char k, int x, int y)
{
    TwEventKeyboardGLUT (k, x, y);
    visualizer->keyPress(k);
    display();
}

void special (int k, int x, int y)
{
    TwEventSpecialGLUT (k, x, y);
	glutPostRedisplay();
}

int main(int argc, char *argv[])
{
    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize (WindowW, WindowH);
    glutCreateWindow ("visualizer");
	
    TwInit (TW_OPENGL, NULL);
    TwWindowSize (WindowW, WindowH);
    
    visualizer = new Visualizer (WindowW, WindowH);

	State* st = State::Instance();
    
    visualizer->initializeGL();
    visualizer->initMesh();

    glutDisplayFunc (display);
    glutReshapeFunc (reshape);
    glutMouseFunc (mouse_click);
    glutMotionFunc (mouse_move);
    glutPassiveMotionFunc (mouse_move);
    glutKeyboardFunc (keyboard);
    glutSpecialFunc (special);
    TwGLUTModifiersFunc (glutGetModifiers);
	
    TwBar *cBar;
    cBar = TwNewBar ("Spline2D");
    TwDefine("Spline2D position='750 100'");
    TwDefine("Spline2D size='200 420'");
    TwDefine("Spline2D valueswidth=75");
    //TwDefine("Visualizer color='192 255 192' text=dark ");
	
    TwAddVarRW(cBar, "Show Control Line", TW_TYPE_BOOLCPP, &st->m_showControlPoints, "");
    TwAddVarRW(cBar, "Show Final   Spline", TW_TYPE_BOOLCPP, &st->m_showFinalPoints, "");
    TwAddSeparator(cBar, "", "");
    TwAddVarRW(cBar, "Show Control Points", TW_TYPE_BOOLCPP, &st->m_paintP, "");
    TwAddVarRW(cBar, "Show Final   Points", TW_TYPE_BOOLCPP, &st->m_paintQ, "");
    TwAddVarRW(cBar, "Move all / selected / unselected vertices", TW_TYPE_CHAR, &st->m_moveVertices, " min=48 max=50 ");
    
    TwAddSeparator(cBar, "", "");
    
    TwEnumVal UIStateEV[] = { 
        {State::UI_SELECT, "Select"},
        {State::UI_NEWSPLINE, "New Spline"},
        {State::UI_DELETESPLINE, "Delete Spline"},
        {State::UI_ADDCONTROLPOINT, "Insert Control Point"},
        {State::UI_MOVE, "Move Control Points"}, 
        {State::UI_ROTATE, "Rotate Control Points"},
        {State::UI_SCALE, "Scale Control Points"}
    };
    TwType UIStateType = TwDefineEnum("UIStateType", UIStateEV, 7);
    TwAddVarRW(cBar, "Tool", UIStateType, &st->m_uistate, NULL);
    
    glutMainLoop ();
    return 1;
}
