// visualizer.cpp
//

#include "visualizer.h"

#include <iomanip>
#include <AntTweakBar.h>
#include "state.h"

using namespace Eigen;



Visualizer::Visualizer(int w, int h)
{
    width = w;
    height = h;
    isSelecting = false;
    initialPos[0] = 0;
    initialPos[1] = 0;
    isDown = false;
}

void Visualizer::initMesh ()
{
}

void Visualizer::updateGL ()
{
    glutPostRedisplay ();
}

void Visualizer::initializeGL ()
{
	glDisable(GL_CULL_FACE);
    glShadeModel(GL_SMOOTH);
}

void Visualizer::resizeGL (int w, int h)
{
    width = w;
    height = h;
    glViewport (0, 0, (GLsizei) w, (GLsizei) h);
    initializeGL ();
}

void Visualizer::paintGL ()
{
    State* st = State::Instance();

    glClearColor(223.0/255.0,223.0/255.0,223.0/255.0,255.0);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glDisable(GL_DEPTH_TEST);
	
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
	
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity ();
    
    st->paint();
    
    // paint the selection rectangle
    if (isSelecting)
    {
        glColor3b(0,0,0);
        Vector2d m(std::min(initialPos[0],mousePos[0]),
                   std::min(initialPos[1],mousePos[1]));
        Vector2d M(std::max(initialPos[0],mousePos[0]),
                   std::max(initialPos[1],mousePos[1]));
        
        glBegin(GL_LINE_LOOP);
        glVertex2f(m[0], m[1]);
        glVertex2f(M[0], m[1]);
        glVertex2f(M[0], M[1]);
        glVertex2f(m[0], M[1]);
        glEnd();
	}
	
}

void Visualizer::mousePress (int button, int state, int x, int y)
{
    State* st = State::Instance();

	Point3d center = unproject(x, height-y, 0);
    int km = glutGetModifiers();

    if (state == GLUT_DOWN) // Button pressed
    {
        isDown = true;
        initialPos[0] = center[0];
        initialPos[1] = center[1];

        switch(st->m_uistate)
        {
            case State::UI_NEWSPLINE:
            {
                Spline s;
                s.addP(Vector2d(center[0],center[1]));
                st->m_splines.push_back(s);
                st->m_selected = pair<int,int>(st->m_splines.size()-1,0);
                st->m_uistate = State::UI_ADDCONTROLPOINT;
                st->updateM();
            }
                break;
            case State::UI_ADDCONTROLPOINT:
            {
                if (st->isSelected())
                {
                    pair<int,int> p = st->m_selected;
                    st->m_splines[p.first].addP(Vector2d(center[0],center[1]),p.second);
                    st->m_selected = pair<int,int>(p.first,p.second+1);
                    st->updateM();
                }
            }
                break;
            case State::UI_SELECT:
            {
                isSelecting = true;
            }
                break;

			default:
				break;
        }
        
//        if (km & GLUT_ACTIVE_CTRL)
//        {
//            mouseMove(x,y);
//        }
//        else if (km & GLUT_ACTIVE_ALT)
//        {
//            mouseMove(x,y);
//        }
//        else
//        {
//            if (button == GLUT_RIGHT_BUTTON)
//            {
//                
//            }
//        }
    }
    else // Button released
    {
        isDown = false;

        switch(st->m_uistate)
        {
            case State::UI_SELECT:
            {
                isSelecting = false;
                if (!(km & GLUT_ACTIVE_SHIFT))
                    st->clearAllSelected();
                st->addSelected(std::min(initialPos[0],center[0]),
                                std::min(initialPos[1],center[1]),
                                std::max(initialPos[0],center[0]),
                                std::max(initialPos[1],center[1])
                                );
                st->setSelected();

            }
                break;

            default:
				break;
        }
    }

}

void Visualizer::mouseMove (int x, int y)
{
    Point3d center3d = unproject(x, height-y, 0);
    Point2d center(center3d[0],center3d[1]);
    
    State* st = State::Instance();
    mousePos = center;
    
    switch(st->m_uistate)
    {
        case State::UI_MOVE:
        {
            if (isDown)
            {
                Point2d diff = mousePos - initialPos;
                initialPos = mousePos;
                
                vector<pair<int,int> > cps = st->getSelected();
                for (size_t i=0; i<cps.size(); ++i)
                {
                    int a = cps[i].first;
                    int b = cps[i].second;
                    
                    st->m_splines[a].m_p[b] += diff;
                }
                st->updateQ();
            }
        }
            break;

        case State::UI_SCALE:
        {
            if (isDown)
            {
                Point2d diff = mousePos - initialPos;
                double scale = 1.0 + diff[0];
                initialPos = mousePos;
                
                vector<pair<int,int> > cps = st->getSelected();
                
                // compute barycenter
                Point2d bary(0,0);
                for (size_t i=0; i<cps.size(); ++i)
                {
                    int a = cps[i].first; int b = cps[i].second;
                    bary += st->m_splines[a].m_p[b];
                }
                bary /= cps.size();
                
                for (size_t i=0; i<cps.size(); ++i)
                {
                    int a = cps[i].first;
                    int b = cps[i].second;
                    
                    st->m_splines[a].m_p[b] -= bary;
                    st->m_splines[a].m_p[b] *= scale;
                    st->m_splines[a].m_p[b] += bary;
                    
                }
                st->updateQ();
            }
        }
            break;
        case State::UI_ROTATE:
        {
            if (isDown)
            {
                Point2d diff = mousePos - initialPos;
                double t = diff[0];
                initialPos = mousePos;
                
                vector<pair<int,int> > cps = st->getSelected();
                
                // compute barycenter
                Point2d bary(0,0);
                for (size_t i=0; i<cps.size(); ++i)
                {
                    int a = cps[i].first; int b = cps[i].second;
                    bary += st->m_splines[a].m_p[b];
                }
                bary /= cps.size();

                Matrix2d rot;
                rot << cos(t), -sin(t), sin(t), cos(t);
                cerr << rot << endl;

                for (size_t i=0; i<cps.size(); ++i)
                {
                    int a = cps[i].first;
                    int b = cps[i].second;
                    
                    st->m_splines[a].m_p[b] -= bary;
                    st->m_splines[a].m_p[b] = rot * st->m_splines[a].m_p[b];
                    st->m_splines[a].m_p[b] += bary;
                    
                }
                st->updateQ();
            }
        }
            break;

		default:
			break;
    }

}

void Visualizer::keyPress(unsigned char k)
{
    State* st = State::Instance();

    switch(k)
    {
        case 'q':
            st->m_uistate = State::UI_NEWSPLINE;
            break;
        case 'w':
            st->m_uistate = State::UI_SELECT;
            break;
        case 'e':
            st->m_uistate = State::UI_ADDCONTROLPOINT;
            break;
        case 'a':
            st->m_uistate = State::UI_MOVE;
            break;
        case 's':
            st->m_uistate = State::UI_SCALE;
            break;
        case 'd':
            st->m_uistate = State::UI_ROTATE;
            break;
        case 'z':
        {
            st->deleteSelected();
            st->clearSelected();
            st->clearAllSelected();
            st->updateM();
            st->m_uistate = State::UI_SELECT;
        }
            break;
            
    }
}

void Visualizer::renderText(Point3d p, string str, bool _small) 
{  
	const char *c = str.c_str();
	glRasterPos3f(p[0], p[1], p[2]);
	for (; *c != '\0'; c++) {
		glutBitmapCharacter(_small ? GLUT_BITMAP_TIMES_ROMAN_10 : GLUT_BITMAP_TIMES_ROMAN_24, *c);
	}
}

string Visualizer::convertInt(int number)
{
	stringstream ss;//create a stringstream
	ss << number;//add number to the stream
	return ss.str();//return a string with the contents of the stream
}

string Visualizer::convertDouble(double number)
{
	stringstream ss;//create a stringstream
	ss << setprecision(4) << number;//add number to the stream
	return ss.str();//return a string with the contents of the stream
}

Point3d Visualizer::project(const Point3d& p)
{
	double model_view_matrix[16];
	double projection_matrix[16];
	int    view_port[4];
	glGetDoublev(GL_MODELVIEW_MATRIX,  model_view_matrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
	glGetIntegerv(GL_VIEWPORT, view_port);
	double x, y, z;
	gluProject(p[0], p[1], p[2], model_view_matrix, projection_matrix,
			   view_port, &x, &y, &z);
	return Point3d(x, y, z);
}

Point3d Visualizer::unproject(double sx, double sy, double sz)
{
	double model_view_matrix[16];
	double projection_matrix[16];
	int    view_port[4];
	glGetDoublev(GL_MODELVIEW_MATRIX,  model_view_matrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
	glGetIntegerv(GL_VIEWPORT, view_port);
	double x, y, z;
	gluUnProject(sx, sy, sz, model_view_matrix, projection_matrix,
				 view_port, &x, &y, &z);
	return Point3d(x, y, z);
}

