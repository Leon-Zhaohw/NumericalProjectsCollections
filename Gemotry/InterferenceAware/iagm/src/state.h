// state.h
//

#include <Eigen/Core>
#include <vector>

#include <Curve.h>

using namespace std;
using namespace Eigen;



class Spline
{
public:
    vector<Vector2d> m_p; // Control Points
    vector<Vector2d> m_q; // Subdivided points
    MatrixXd M;           // Subdivision matrix (q = Mp)

    vector<int> m_selected; // indices of selected control points
    
    // Update the subdivision matrix
	//
    void updateM(int lvl = 3);
    
    // Recompute the position of the subdivided points
	//
    void updateQ();
    
    // Insert a new control point
	//
    void addP(Vector2d p);
    void addP(Vector2d p, int i);
    
    // Delete a control point
	//
    void removeP(int i);
    
    // Paint the spline
	//
    void paint();
};



class State : public IAGM::Curve
{
public:

    typedef enum 
    {
        UI_SELECT,
        UI_NEWSPLINE,
        UI_DELETESPLINE,
        UI_ADDCONTROLPOINT,
        UI_MOVE,
        UI_ROTATE,
        UI_SCALE
    } UIState;
    
    static State* Instance();
    vector<Spline> m_splines;
    
    pair<int,int> m_selected;
    bool isSelected()
	{ return m_selected != pair<int,int>(-1,-1); }
    
    // Global subdivision matrix
	//
    MatrixXd M;
    size_t m_qsize;
    size_t m_psize;
    
    // Paint Control Vertices
	//
    bool m_paintP;
    
    // Paint Final Vertices
	//
    bool m_paintQ;

    // Enable visualization of control points
	//
    bool m_showControlPoints;
    
    // Enable visualization of finer points
	//
    bool m_showFinalPoints;

	// Move all (0), selected (1), or unselected (2) vertices
	//
	char m_moveVertices;
    
    // UI State
	//
    UIState m_uistate;

	double *m_xStart;
	double *m_xEnd;
	double *m_xCtrl;

	unsigned int *m_edges;

    
    // Clear the selected flag for all control points
	//
    void clearAllSelected();
    void clearSelected();
    
    // Set the selected vertex if only one is selected
	//
    void setSelected();

    // Delete the selected vertices
	//
    void deleteSelected();
    
    // Add new vertices to the current selection
	//
    void addSelected(double xmin, double ymin, double xmax, double ymax);
    
    // Return all the selected control points
	//
    vector<pair<int,int> > getSelected();
    
    // Paint all splines
	//
    void paint();

    // Load from disk
	//
    void load(char* file);

    // Save to disk
	//
    void save(char* file);

    // Compute subdivision matrix
	//
    void updateM();
    
    // Update final positions
	//
    void updateQ();


	// Derived methods
	//
	void updateSlaveMesh();

	void getSubspaceVector(double *in, double *out);

	size_t getNbrVertices() { return m_qsize; }
	size_t getNbrCtrlVertices() { return m_psize; }

	void getStartConfiguration(double *&x);
	void getEndConfiguration(double *&x);

	void getControlConfiguration(double *&x);
	void copyControlConfiguration();

	void getMovableControlVertices(std::set<size_t> &movable);

	void getEdgeIndices(size_t &nbr, unsigned int *&e);
   
    
private:
	State();  // Private so that it can not be called

    State(State const&) { };             // copy constructor is private
    void operator=(State const&) { };    // assignment operator is private
    static State* m_pInstance;
};

