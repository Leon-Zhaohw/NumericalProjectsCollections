#ifdef WIN32
#pragma once
#endif

#include "simulationTree.h"

#include <wx/glcanvas.h>
#include <wx/event.h>
#include <wx/scrolbar.h>
#include <wx/notebook.h>
#include <wx/choice.h>

BEGIN_DECLARE_EVENT_TYPES()
	DECLARE_EVENT_TYPE(wxEVT_CHANGE_NODE_COUNT, -1)
END_DECLARE_EVENT_TYPES()

namespace planning {


class BestScoresCanvas;
class SimulationTree;
class MouseCapture;

class ObjectiveComputer
{
public:
	virtual ~ObjectiveComputer();
	virtual float operator()( const Path* path ) const = 0;

	virtual std::vector<size_t> objects() const = 0;
	virtual bool operator==( const ObjectiveComputer& other ) const = 0;

private:
};

class SelectScore
	: public ObjectiveComputer
{
public:
	SelectScore( size_t metric, const std::vector<size_t>& objects, const Objective& objective );
	~SelectScore();

	float operator()( const Path* path ) const;

	size_t sortBy() const					{ return this->which_; }
	std::vector<size_t> objects() const		{ return this->objects_; }
	bool operator==( const ObjectiveComputer& other ) const;

private:
	const Objective& objective_;
	size_t which_;
	std::vector<size_t> objects_;
};

class CountConstraints
	: public ObjectiveComputer
{
public:
	CountConstraints( ConstSimulationTreePtr tree );
	~CountConstraints();

	float operator()( const Path* path ) const;
	std::vector<size_t> objects() const;
	bool operator==( const ObjectiveComputer& other ) const;

private:
	ConstSimulationTreePtr tree_;
};

typedef boost::shared_ptr<const ObjectiveComputer> ObjectiveComputerPtr;

class BestScoresPanel;

class BestScoresFrame
	: public wxFrame
{
public:
	BestScoresFrame(wxFrame *frame, wxGLContext* sharedContext, const wxString& title, const wxPoint& pos, const wxSize& size,
        long style = wxDEFAULT_FRAME_STYLE);
	~BestScoresFrame();

	void OnClose(wxCloseEvent& event);

	void updateShowAllPaths();
	void updateSelected();
	void newPanel(boost::shared_ptr<SimulationTree> tree, const std::string& title);
	void removePanel( boost::shared_ptr<SimulationTree> tree );

	void clearPanels();
	void setPath( SimulationTreePtr tree, const Path* path );

	void preDestroy();

private:
	enum {
		NOTEBOOK,
	};

	wxGLContext* sharedContext_;
	wxNotebook* notebook_;

	DECLARE_EVENT_TABLE()
};

class BestScoresPanel
	: public wxPanel
{
public:
	BestScoresPanel( SimulationTreePtr tree,
		wxWindow *parent, 
		wxGLContext* sharedContext, 
		const wxPoint& pos, 
		const wxSize& size,
        long style = wxDEFAULT_FRAME_STYLE);
	~BestScoresPanel();

	void OnScroll(wxScrollEvent& event);
	void OnSetSortBy(wxCommandEvent& event);
	void OnSetSortOrder(wxCommandEvent& event);
	void OnSetNodeCount(wxCommandEvent& event);
	void setScrollSize( int range, int pageSize );

	int getScrollPosition() const;
	void setScrollPosition(int pos);

	void setPath( const Path* path );
	SimulationTreePtr tree();

	void updateSelected();
	void updateShowAllPaths();


	enum SortOrder
	{
		SORT_ORDER_ASCENDING,
		SORT_ORDER_DESCENDING,
	};

	SortOrder sortOrder() const;
	size_t sortBy() const;

	void preDestroy();

private:
	SortOrder sortOrder_;

	enum {
		SCROLL_BAR,
		CHOICE_SORT_BY,
		CHOICE_SORT_DIRECTION,
	};

	int scrollPosition_;

	BestScoresCanvas* canvas_;
	wxScrollBar* scrollBar_;

	wxChoice *sortByChoice_;
	wxChoice *sortDirectionChoice_;

	DECLARE_EVENT_TABLE()
};

class BestScoresCanvas
	: public wxGLCanvas, public SimulationNodeCallback
{
public:
	BestScoresCanvas(SimulationTreePtr tree, wxWindow *parent, wxGLContext* sharedContext, const wxWindowID id = -1, const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize, long style = 0, const wxString& name = wxT("GraphCanvas"));
	~BestScoresCanvas();

    void OnPaint(wxPaintEvent& event);
	void OnMouse( wxMouseEvent& event );
    void OnSize(wxSizeEvent& event);

	void OnEraseBackground(wxEraseEvent& event);

	ObjectiveComputerPtr selectFunction() const { return this->selectScore_; }
	void setSortBy( size_t value );

	bool updateObjectiveComputer();
	void updateSort();

	void updateSelected();
	void updateShowAllPaths();

	void addNodes( const SimulationTree* tree, const std::deque<const Path*>& nodes, const boost::dynamic_bitset<>& active );
	void setNodes( const SimulationTree* tree, const boost::dynamic_bitset<>& nodes );

	std::deque<const Path*> nodes() const;

	void setPath( const Path* path );
	SimulationTreePtr tree();

	void preDestroy();

	size_t numShown() const;

private:
	SimulationTreePtr tree_;
	boost::scoped_ptr<SimulationNodeCallbackRegistrar> callbackRegistrar_;

	ObjectiveComputerPtr selectScore_;

	std::deque<const Path*> nodes_;
	mutable boost::mutex nodesMutex_;

	size_t subImageWidth_;
	size_t subImageHeight_;

	size_t gridWidth_;
	size_t gridHeight_;

	size_t border_;

	size_t selected_;

	int currentWheelRot_;

	DECLARE_EVENT_TABLE()
};

} // namespace planning

