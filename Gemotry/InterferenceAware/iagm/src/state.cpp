// state.cpp
//

#include "timer.h"
#include "state.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#undef min
#undef max
#include <limits>
#include <cmath>
#include <fstream>
#include <iostream>

extern void display ();



void Spline::updateM(int lvl)
{
    if (m_p.size() < 4)
    {
        m_q.clear();
        return;
    }
    
    // Every segment is subdivided lvl times
	//
    int f = pow(2.0,lvl)+1;
    VectorXd ranges(f);
    
    // Ranges from 0 to 1
	//
    for (int i=0; i<f; ++i)
        ranges(i) = (double)i/(f-1.0);
    
    // Compute the subdivision matrix for a single segment
	//
    MatrixXd sM = MatrixXd::Zero(f,4);

    MatrixXd A(4,4);
    A << -1,  3, -3, 1,
          3, -6,  3, 0,
         -3,  0,  3, 0,
          1,  4,  1, 0;

    for (int i=0; i<f; ++i)
    {
        double t = ranges(i);
        MatrixXd tm(4,1); tm << t*t*t, t*t, t, 1;
        MatrixXd temp = tm.transpose()*(1.0/6.0)*A;
        for(int j=0;j<4;++j)
            sM(i,j) = temp(j);
    }
    
    int segmentCount = (int)m_p.size()-3;
    int finalPointCount = (segmentCount * f) - (segmentCount-1);
    
    M = MatrixXd::Zero(finalPointCount, m_p.size());
    
    // Produce the set of final points and the global subdivision matrix
	//
    for (size_t m=1; m<m_p.size()-2; ++m)
    {
        int blockN = m-1;
        
        int i = blockN*(f-1);
        int j = blockN;
        M.block(i, j, sM.rows(), sM.cols()) = sM;
    }
    
    m_q.resize(finalPointCount);
    
    updateQ();
}

void Spline::updateQ()
{
	MatrixXd b = MatrixXd(m_p.size(),2);
	for (size_t i=0; i<m_p.size();++i)
		b.row(i) = m_p[i];

	MatrixXd finalPos = M*b;
	for (size_t i=0; i<m_q.size(); ++i)
		m_q[i] = finalPos.row(i);
}

void Spline::addP(Vector2d p)
{
    m_p.push_back(p);
    updateM();
}

void Spline::addP(Vector2d p, int i)
{
    m_p.insert(m_p.begin()+i+1, p);
    updateM();
}

void Spline::removeP(int i)
{
    m_p.erase(m_p.begin()+i);
}

void paint_aux(vector<Vector2d>& v, int r, int g, int b, bool paintPoints)
{
    // Inside
    glLineWidth(1);
    glColor3ub(60, 60, 60);
    glBegin(GL_LINE_STRIP);
    for(size_t i=0;i<v.size();++i)
        glVertex3f(v[i](0),v[i](1),1);
    glEnd();

    if (paintPoints)
    {
        // Outer Border
        glPointSize(10);
        glColor3ub(r-20, g-20, b-20);
        glBegin(GL_POINTS);
        for(size_t i=0;i<v.size();++i)
            glVertex3f(v[i](0),v[i](1),1);
        glEnd();
        
        // Inside
        glPointSize(5);
        glColor3ub(r, g, b);
        glBegin(GL_POINTS);
        for(size_t i=0;i<v.size();++i)
            glVertex3f(v[i](0),v[i](1),1);
        glEnd();
    }
}

void Spline::paint()
{
    State* st = State::Instance();

    if (st->m_showControlPoints)
        paint_aux(m_p,200,50,50,st->m_paintP);

    if (st->m_showFinalPoints)
        paint_aux(m_q,50,200,50,st->m_paintQ);
    
    glColor3ub(50,50,250);
    glPointSize(15);
    glBegin(GL_POINTS);
    for(size_t i=0; i<m_selected.size(); ++i)
    {
        glVertex3f(m_p[m_selected[i]](0),m_p[m_selected[i]](1),1);
    }
    glEnd();
}


// Global static pointer used to ensure a single instance of the class.
State* State::m_pInstance = NULL; 

/** This function is called to create an instance of the class.
 Calling the constructor publicly is not allowed. The constructor
 is private and is only called by this Instance function.
 */

State::State()
 : m_xStart(0), m_xEnd(0), m_xCtrl(0), m_edges(0)
{
    m_uistate = UI_NEWSPLINE;
    m_showControlPoints = true;
    m_showFinalPoints = true;
    clearAllSelected();
    clearSelected();
    m_paintQ = false;
    m_paintP = true;
}


State* State::Instance()
{
    if (!m_pInstance)   // Only allow one instance of class to be generated.
        m_pInstance = new State;
    return m_pInstance;
}

void State::clearSelected()
{
    m_selected = pair<int,int>(-1,-1);
}

void State::clearAllSelected()
{
    for (size_t i=0; i<m_splines.size(); ++i)
        m_splines[i].m_selected.clear();

    m_selected = pair<int,int>(-1,-1);
}

vector<pair<int,int> > State::getSelected()
{
    vector<pair<int,int> > ret;
    for (size_t i=0; i<m_splines.size(); ++i)
    {
        Spline& s = m_splines[i];
        for (size_t j=0; j<s.m_selected.size(); ++j)
            ret.push_back(pair<int,int>(i,s.m_selected[j]));
    }

    return ret;
}

void State::deleteSelected()
{
    for (size_t i=0; i<m_splines.size(); ++i)
    {
        Spline& s = m_splines[i];
        int count = 0;
        for (size_t j=0; j<s.m_selected.size(); ++j)
            s.m_p.erase(s.m_p.begin() + count++ + s.m_selected[j]);
    }
}

void State::paint()
{
    for (size_t i=0; i<m_splines.size(); ++i)
        m_splines[i].paint();
}

void State::load(char* file)
{
    // TODO
}

void State::save(char* file)
{
    // TODO
}

// Set the selected vertex if only one is selected
//
void State::setSelected()
{
    vector<pair<int,int> > s = getSelected();
    if (s.size() == 1)
    {
        m_selected.first  = s[0].first;
        m_selected.second = s[0].second;
    }
    else
    {
        m_selected.first  = -1;
        m_selected.second = -1;
    }
}

// Add new vertices to the current selection
//
void State::addSelected(double xmin, double ymin, double xmax, double ymax)
{
    for (size_t i=0; i<m_splines.size(); ++i)
    {
        Spline& s = m_splines[i];
        for (size_t j=0; j<s.m_p.size(); ++j)
        {
            Vector2d& p = s.m_p[j];
            if (p[0] >= xmin && p[0] < xmax && p[1] >= ymin && p[1] < ymax)
                s.m_selected.push_back(j);
        }
        
        // Clean duplicate selected
		//
        std::sort(s.m_selected.begin(), s.m_selected.end());
        s.m_selected.erase(std::unique(s.m_selected.begin(), s.m_selected.end()), s.m_selected.end());
    }
}

void State::updateM()
{
    m_psize = 0;
    m_qsize = 0;

	// First update each spline individually
	//
    for (size_t i=0; i<m_splines.size(); ++i)
    {
        Spline& s = m_splines[i];

        s.updateM();
        m_psize += s.m_p.size();
        m_qsize += s.m_q.size();
    }

	// Construct the global subdivision matrix
	//
    M = MatrixXd::Zero(m_qsize, m_psize);
    size_t currenti = 0;
    size_t currentj = 0;
    
    for (size_t i=0; i<m_splines.size(); ++i)
    {
        // For every spline copy the subdivision matrix in the diagonal of M
		//
        Spline& s = m_splines[i];
        
        M.block(currenti,currentj,s.M.rows(),s.M.cols()) = s.M;
        currenti += s.M.rows();
        currentj += s.M.cols();
    }

	if (m_qsize)
	{
		// Re-allocate data for iagm
		//
		delete[] m_xStart;
		delete[] m_xEnd;
		delete[] m_xCtrl;
		delete[] m_edges;

		m_xStart = new double[2 * m_qsize];
		m_xEnd   = new double[2 * m_qsize];
		m_xCtrl  = new double[2 * m_psize];

		m_edges  = new unsigned int[2 * (m_qsize - m_splines.size())];

		// Fill out edge data structure
		//
		for (size_t i=0, currPt=0; i<m_splines.size(); ++i)
		{
			size_t qsize = m_splines[i].m_q.size();
			if (qsize == 0)
				continue;

			for (size_t j=0; j<(qsize-1); ++j, ++currPt)
			{
				m_edges[2*currPt  ] = currPt;
				m_edges[2*currPt+1] = currPt + 1;
			}
		}
	}
}

void State::updateQ()
{
	// Copy updated control mesh data
	//
	for (size_t i=0, currPt=0; i<m_splines.size(); ++i)
	{
		Spline& s = m_splines[i];
		for (size_t j=0; j<s.m_p.size(); ++j, ++currPt)
		{
			m_xCtrl[2*currPt  ] = s.m_p[j][0];
			m_xCtrl[2*currPt+1] = s.m_p[j][1];
		}
	}

	// Current positions are the start positions for the trajectory
	//
    for (size_t i=0, currPt=0; i<m_splines.size(); ++i)
    {
        Spline& s = m_splines[i];
        for (size_t j=0; j<s.m_q.size(); ++j, ++currPt)
		{
			m_xStart[2*currPt  ] = s.m_q[j][0];
			m_xStart[2*currPt+1] = s.m_q[j][1];
		}
    }

	updateSlaveMesh();
	
	// Call iagm library for interference detection / response
	//
	findAndRemoveInterference();
}

void State::updateSlaveMesh()
{
	// Update each spline with new control points
	//
    for (size_t i=0; i<m_splines.size(); ++i)
		m_splines[i].updateQ();
    
	// Copy updated positions as candidate end trajectories
	//
	for (size_t i=0, currPt=0; i<m_splines.size(); ++i)
    {
        Spline& s = m_splines[i];
        for (size_t j=0; j<s.m_q.size(); ++j, ++currPt)
		{
			m_xEnd[2*currPt  ] = s.m_q[j][0];
			m_xEnd[2*currPt+1] = s.m_q[j][1];
		}
    }
}

void State::getSubspaceVector(double *in, double *out)
{
//	Map<VectorXd, 0, InnerStride<1> >(out  , m_psize) = M.transpose() * Map<VectorXd, 0, InnerStride<1> >(in  , m_qsize);
//	Map<VectorXd, 0, InnerStride<1> >(out+1, m_psize) = M.transpose() * Map<VectorXd, 0, InnerStride<1> >(in+1, m_qsize);

    MatrixXd q(m_qsize, 2);
	for (size_t i=0; i<m_qsize; ++i)
		q.row(i) = Vector2d(in[2*i], in[2*i+1]);
	
	MatrixXd p = M.transpose() * q;

	for (size_t i=0; i<m_psize; ++i)
	{
		out[2*i+0] = p.row(i)[0];
		out[2*i+1] = p.row(i)[1];
	}
}

void State::getStartConfiguration(double *&x)
{
	x = m_xStart;
}

void State::getEndConfiguration(double *&x)
{
	x = m_xEnd;
}

void State::getControlConfiguration(double *&x)
{
	x = m_xCtrl;
}

void State::copyControlConfiguration()
{
	// Update each spline with new control configuration
	//
    for (size_t i=0, currPt=0; i<m_splines.size(); ++i)
    {
        Spline& s = m_splines[i];
		for (size_t j=0; j<s.m_p.size(); ++j, ++currPt)
		{
			s.m_p[j][0] = m_xCtrl[2*currPt  ];
			s.m_p[j][1] = m_xCtrl[2*currPt+1];
		}
    }
}

void State::getMovableControlVertices(set<size_t> &movable)
{
	movable.clear();

	switch (m_moveVertices)
	{

	case '1': // Move only selected
		for (size_t i=0, currTotal=0; i<m_splines.size(); ++i)
		{
			Spline& s = m_splines[i];
			for (size_t j=0; j<s.m_selected.size(); ++j)
				movable.insert(currTotal + s.m_selected[j]);

			currTotal += s.m_p.size();
		}
		break;
	
	case '2': // Move only unselected
	{
		set<size_t> selected, all;
		for (size_t i=0, currTotal=0; i<m_splines.size(); ++i)
		{
			Spline& s = m_splines[i];
			for (size_t j=0; j<s.m_selected.size(); ++j)
				selected.insert(currTotal + s.m_selected[j]);
			
			for (size_t j=0; j<s.m_p.size(); ++j)
				all.insert(currTotal+j);

			currTotal += s.m_p.size();
		}
	
		set_difference(all.begin(), all.end(), selected.begin(), selected.end(),
					   inserter(movable, movable.end()));

		break;
	}
	default: // Move all vertices
		for (size_t i=0; i<m_psize; ++i)
			movable.insert(i);

		break;
	}
}

void State::getEdgeIndices(size_t &nbr, unsigned int *&e)
{
	e = m_edges;

	nbr = m_qsize - m_splines.size();
}

