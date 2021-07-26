#include "stdafx.h"

#include "bestScoresView.h"
#include "planningUI.h"
#include "GLView.h"
#include "constraints.h"
#include "walkTreeHandler.h"

#include "twigg/renderUtils.h"
#include <wx/stattext.h>

DEFINE_EVENT_TYPE(wxEVT_CHANGE_NODE_COUNT)

namespace planning {

BEGIN_EVENT_TABLE(BestScoresFrame, wxFrame)
	EVT_CLOSE(BestScoresFrame::OnClose)
END_EVENT_TABLE()

BestScoresFrame::BestScoresFrame(wxFrame *frame, wxGLContext* sharedContext, const wxString& title, const wxPoint& pos,
    const wxSize& size, long style)
	: wxFrame(frame, -1, title, pos, size, style)
{
	sharedContext_ = sharedContext;
	notebook_ = new wxNotebook(this, NOTEBOOK, wxDefaultPosition, wxDefaultSize, wxNB_LEFT);
}

BestScoresFrame::~BestScoresFrame()
{
}

void BestScoresFrame::newPanel(boost::shared_ptr<SimulationTree> tree, const std::string& name)
{
	BestScoresPanel* newPanel = new BestScoresPanel(tree, notebook_, sharedContext_, wxDefaultPosition, wxDefaultSize);
	this->notebook_->AddPage( newPanel, wxT(name.c_str()), true );
}

void BestScoresFrame::removePanel( boost::shared_ptr<SimulationTree> tree )
{
	for( size_t iPage = 0; iPage < this->notebook_->GetPageCount(); ++iPage )
	{
		wxWindow* page = this->notebook_->GetPage(iPage);
		BestScoresPanel* panel = dynamic_cast<BestScoresPanel*>( page );
		if( panel->tree() == tree )
		{
			panel->preDestroy();
			notebook_->DeletePage(iPage);
			return;
		}		
	}
}

void BestScoresFrame::clearPanels()
{
	for( size_t iPage = 0; iPage < this->notebook_->GetPageCount(); ++iPage )
	{
		wxWindow* page = this->notebook_->GetPage(iPage);
		BestScoresPanel* panel = dynamic_cast<BestScoresPanel*>( page );
		panel->preDestroy();
	}

	this->notebook_->DeleteAllPages();
}

void BestScoresFrame::setPath( SimulationTreePtr tree, const Path* path )
{
	for( size_t i = 0; i < this->notebook_->GetPageCount(); ++i )
	{
		wxWindow* page = this->notebook_->GetPage(i);
		BestScoresPanel* panel = dynamic_cast<BestScoresPanel*>( page );
		if( panel->tree() == tree )
		{
			panel->setPath( path );
			if( i != notebook_->GetSelection() )
				notebook_->SetSelection(i);
			return;
		}
	}
}

void BestScoresFrame::updateSelected()
{
	for( size_t i = 0; i < this->notebook_->GetPageCount(); ++i )
	{
		wxWindow* page = this->notebook_->GetPage(i);
		BestScoresPanel* panel = dynamic_cast<BestScoresPanel*>( page );
		panel->updateSelected();
	}
}

void BestScoresFrame::updateShowAllPaths()
{
	for( size_t i = 0; i < this->notebook_->GetPageCount(); ++i )
	{
		wxWindow* page = this->notebook_->GetPage(i);
		BestScoresPanel* panel = dynamic_cast<BestScoresPanel*>( page );
		panel->updateShowAllPaths();
	}
}

void BestScoresFrame::OnClose(wxCloseEvent& event)
{
	if( event.CanVeto() )
	{
		this->Hide();
		event.Veto();
	}
	else
	{
		this->Destroy();
	}
}

void BestScoresFrame::preDestroy()
{
	for( size_t i = 0; i < this->notebook_->GetPageCount(); ++i )
	{
		wxWindow* page = this->notebook_->GetPage(i);
		BestScoresPanel* panel = dynamic_cast<BestScoresPanel*>( page );
		panel->preDestroy();
	}
}

BEGIN_EVENT_TABLE(BestScoresPanel, wxPanel)
	EVT_COMMAND_SCROLL(BestScoresPanel::SCROLL_BAR, BestScoresPanel::OnScroll)
	EVT_CHOICE(BestScoresPanel::CHOICE_SORT_BY, BestScoresPanel::OnSetSortBy)
	EVT_CHOICE(BestScoresPanel::CHOICE_SORT_DIRECTION, BestScoresPanel::OnSetSortOrder)
	EVT_COMMAND(wxID_ANY, wxEVT_CHANGE_NODE_COUNT, BestScoresPanel::OnSetNodeCount)
END_EVENT_TABLE()

BestScoresPanel::BestScoresPanel(SimulationTreePtr tree, wxWindow *frame, wxGLContext* sharedContext, const wxPoint& pos,
    const wxSize& size, long style)
	: wxPanel(frame, -1, pos, size, style)
{
	wxBoxSizer *topsizer = new wxBoxSizer( wxVERTICAL );

	wxBoxSizer *mainsizer = new wxBoxSizer( wxHORIZONTAL );

	wxBoxSizer* choiceSizer = new wxBoxSizer(wxHORIZONTAL);
	sortByChoice_ = new wxChoice(this, CHOICE_SORT_BY, wxDefaultPosition, wxDefaultSize, 0, 0);

	sortDirectionChoice_ = new wxChoice(this, CHOICE_SORT_DIRECTION, wxDefaultPosition, wxDefaultSize, 0, 0);
	sortDirectionChoice_->Append( wxT("Ascending") );
	sortDirectionChoice_->Append( wxT("Descending") );
	sortDirectionChoice_->SetSelection(1);
	sortOrder_ = SORT_ORDER_DESCENDING;

	choiceSizer->Add( 
		new wxStaticText(this, -1, "Sort by:", wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT) );
	choiceSizer->Add( sortByChoice_ );
	choiceSizer->Add( sortDirectionChoice_ );
	topsizer->Add( choiceSizer, 0, wxEXPAND );

	scrollBar_ = new wxScrollBar( this, SCROLL_BAR, wxDefaultPosition, wxDefaultSize, wxSB_VERTICAL );
	this->scrollPosition_ = scrollBar_->GetThumbPosition();

	canvas_ = new BestScoresCanvas(tree, this, sharedContext, -1, wxPoint(0, 0), wxSize(200, 200), wxSUNKEN_BORDER);
	mainsizer->Add( canvas_,
		1,
		wxEXPAND );
	mainsizer->Add( scrollBar_, 
		0,                      // don't expand horizontally
		wxEXPAND );             // but expand vertically is okay

	topsizer->Add( mainsizer, 1, wxEXPAND );
	SetSizer( topsizer );


	{
		std::vector<ObjectivePtr> objectives = tree->objectives();

		for( std::vector<ObjectivePtr>::const_iterator iter = objectives.begin();
			iter != objectives.end(); ++iter )
		{
			sortByChoice_->Append( wxT((*iter)->name().c_str()) );
		}
		this->sortByChoice_->SetSelection(0);

		sortByChoice_->Append( wxT("Count constraints") );
	}
}

void BestScoresPanel::preDestroy()
{
	this->canvas_->preDestroy();
}

void BestScoresPanel::updateSelected()
{
	this->canvas_->updateSelected();
}

void BestScoresPanel::updateShowAllPaths()
{
	this->canvas_->updateShowAllPaths();
}

BestScoresPanel::~BestScoresPanel()
{
}

void BestScoresPanel::setPath(const Path* path)
{
	this->canvas_->setPath( path );
}

void BestScoresPanel::OnScroll(wxScrollEvent& event)
{
	scrollPosition_ = event.GetInt();
	this->canvas_->Refresh(FALSE);
}

void BestScoresPanel::OnSetNodeCount(wxCommandEvent& event)
{
	int num = event.GetInt();
	long scrollPos = event.GetExtraLong();

	scrollBar_->SetScrollbar(
		scrollPos,
		canvas_->numShown(),
		num,
		canvas_->numShown(),
		true);
	this->scrollPosition_ = scrollPos;

	this->canvas_->Refresh(FALSE);
}

void BestScoresPanel::setScrollSize( int range, int pageSize )
{
	scrollBar_->SetScrollbar(
		this->scrollPosition_,
		pageSize,
		range,
		pageSize,
		true);
}

void BestScoresPanel::OnSetSortBy(wxCommandEvent& event)
{
	this->canvas_->setSortBy( event.GetInt() );
}

size_t BestScoresPanel::sortBy() const
{
	return sortByChoice_->GetSelection();
}

void BestScoresPanel::OnSetSortOrder(wxCommandEvent& event)
{
	if( this->sortDirectionChoice_->GetSelection() == 1 )
		this->sortOrder_ = BestScoresPanel::SORT_ORDER_DESCENDING;
	else
		this->sortOrder_ = BestScoresPanel::SORT_ORDER_ASCENDING;

	this->canvas_->updateSort();
}

BestScoresPanel::SortOrder BestScoresPanel::sortOrder() const
{
	return sortOrder_;
}

int BestScoresPanel::getScrollPosition() const
{
	return scrollPosition_;
}

void BestScoresPanel::setScrollPosition(int pos)
{
	this->scrollBar_->SetThumbPosition(pos);
}

SimulationTreePtr BestScoresPanel::tree()
{
	return this->canvas_->tree();
}


BEGIN_EVENT_TABLE(BestScoresCanvas, wxGLCanvas)
    EVT_PAINT(BestScoresCanvas::OnPaint)
    EVT_MOUSE_EVENTS(BestScoresCanvas::OnMouse)
    EVT_SIZE(BestScoresCanvas::OnSize)
	EVT_ERASE_BACKGROUND(BestScoresCanvas::OnEraseBackground)
END_EVENT_TABLE()

class SelectFunctionComparator
{
public:
	SelectFunctionComparator( ObjectiveComputerPtr select, BestScoresPanel::SortOrder order )
		: select_(select), order_(order) {}

	bool operator()( const Path* left, 
		const Path* right ) const
	{
		if( order_ == BestScoresPanel::SORT_ORDER_ASCENDING )
		{
			std::less<float> comparator;
			return comparator((*select_)( left ), (*select_)( right ) );
		}
		else
		{
			std::greater<float> comparator;
			return comparator((*select_)( left ), (*select_)( right ) );
		}
	}

private:
	ObjectiveComputerPtr select_;
	BestScoresPanel::SortOrder order_;
};

ObjectiveComputer::~ObjectiveComputer()
{
}

SelectScore::SelectScore( size_t metric, const std::vector<size_t>& objects, const Objective& objective )
	: objective_(objective), which_(metric), objects_(objects)
{
	std::sort( objects_.begin(), objects_.end() );
}

SelectScore::~SelectScore()
{
}

float SelectScore::operator()( const Path* path ) const
{
	float result = this->objective_.startValue();
	for( std::vector<size_t>::const_iterator itr = this->objects_.begin();
		itr != this->objects_.end(); ++itr )
	{
		result = this->objective_.combine( 
			result, 
			path->score( *itr, this->which_ ) );
	}
	return result;
}

bool SelectScore::operator==( const ObjectiveComputer& other ) const
{
	const SelectScore* selectScore = dynamic_cast<const SelectScore*>( &other );
	if( selectScore == 0 )
		return false;
	return (which_ == selectScore->which_) && (objects_ == selectScore->objects_);
}

CountConstraints::CountConstraints( ConstSimulationTreePtr tree )
	: tree_(tree)
{
}

CountConstraints::~CountConstraints()
{
}

bool CountConstraints::operator==( const ObjectiveComputer& other ) const
{
	const CountConstraints* countConstraints = dynamic_cast<const CountConstraints*>( &other );
	return (countConstraints != 0) && (countConstraints->tree_ == this->tree_);
}


float CountConstraints::operator()( const Path* path ) const
{
	float result = 0.0f;

	size_t id = tree_->idForPath( path );
	for( size_t iConstraint = 0; iConstraint < tree_->constraintCount(); ++iConstraint )
	{
		if( tree_->satisfiesConstraint( id, iConstraint ) )
			result += 1.0f;
	}

	return result;
}

std::vector<size_t> CountConstraints::objects() const
{
	std::deque<ConstraintPtr> constraints = tree_->constraints();
	std::vector<size_t> result;
	for( std::deque<ConstraintPtr>::const_iterator itr = constraints.begin();
		itr != constraints.end(); ++itr )
	{
		std::copy( (*itr)->objects_begin(), (*itr)->objects_end(), 
			std::back_inserter(result) );
	}

	std::sort( result.begin(), result.end() );
	result.erase( std::unique(result.begin(), result.end()), result.end() );
	return result;
}


BestScoresCanvas::BestScoresCanvas(SimulationTreePtr tree, wxWindow *parent, wxGLContext* sharedContext, wxWindowID id,
    const wxPoint& pos, const wxSize& size, long style, const wxString& name)
	:	wxGLCanvas(parent, sharedContext, id, pos, size, style, name), currentWheelRot_(0), selected_(0)
{
	subImageHeight_ = 130;
	subImageWidth_ = 3*subImageHeight_ / 2;
	border_ = 5;

	this->tree_ = tree;
	this->callbackRegistrar_.reset( 
		new SimulationNodeCallbackRegistrar(tree_, this) );

	this->setSortBy( 0 );
}

void BestScoresCanvas::preDestroy()
{
	this->callbackRegistrar_.reset();
	this->tree_.reset();
}

BestScoresCanvas::~BestScoresCanvas()
{
}

SimulationTreePtr BestScoresCanvas::tree()
{
	return tree_;
}

void InitGL()
{
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glDisable(GL_AUTO_NORMAL);
	glDisable(GL_NORMALIZE);
	glDisable(GL_TEXTURE_2D);
	glShadeModel(GL_SMOOTH);

	{
		// set up lighting
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
		glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);

		// set up lighting
		vl::Vec4f lightAmbient0( 0.2, 0.2, 0.2, 1.0 );
		vl::Vec4f lightPosition0( 1, -2, 1, 0 );
		vl::Vec4f lightDiffuse0( 0.5, 0.5, 0.5, 1.0 );
		vl::Vec4f lightPosition1( -2, -1, 5, 0 );
		vl::Vec4f lightDiffuse1( 0.5, 0.5, 0.5, 1 );
		vl::Vec4f lightPosition2 = vl::Vec4f(1, 1, 2, 0.0);
		vl::Vec4f lightDiffuse2( 0.3, 0.3, 0.3, 1 );

		glLightfv( GL_LIGHT0, GL_POSITION, norm( lightPosition0 ).Ref() );
		glLightfv( GL_LIGHT0, GL_DIFFUSE, lightDiffuse0.Ref() );
		glLightfv( GL_LIGHT0, GL_AMBIENT, lightAmbient0.Ref() );
		glLightfv( GL_LIGHT1, GL_POSITION, norm( lightPosition1 ).Ref() );
		glLightfv( GL_LIGHT1, GL_DIFFUSE, lightDiffuse1.Ref() );
		glLightfv( GL_LIGHT2, GL_POSITION, norm(lightPosition2).Ref() );
		glLightfv( GL_LIGHT2, GL_DIFFUSE, lightDiffuse2.Ref() );

		glEnable( GL_LIGHT0 );
		glEnable( GL_LIGHT1 );
		glEnable( GL_LIGHT2 );
	}
}

void BestScoresCanvas::OnEraseBackground(wxEraseEvent& event)
{
    // Do nothing, to avoid flashing on MSW
}

struct BestScoresRenderObjectHandler
{
	void operator()( ConstSceneGraphElementPtr object, const vl::Mat4f& transform, const UsefulBits& bits, const char* statePtr )
	{
		SceneGraphElement::RenderState renderState;
		renderState.alpha = 1.0;
		renderState.polygonMode = PhysicsObject::RenderState::FILL;
		object->render(renderState, transform, statePtr);
	}
};

void BestScoresCanvas::OnPaint( wxPaintEvent& event )
{
	BestScoresPanel* parent = dynamic_cast<BestScoresPanel*>(this->GetParent());

	// must always be here
    wxPaintDC dc(this);
#ifndef __WXMOTIF__
    if (!GetContext()) return;
#endif
    SetCurrent();
	InitGL();

	MainFrame* mainFrame = dynamic_cast<MainFrame*>(parent->GetParent()->GetParent()->GetParent());
	vl::Vec3f backgroundColor = mainFrame->backgroundColor();
	glClearColor( backgroundColor[0], backgroundColor[1], backgroundColor[2], 0.0 );
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if( !this->tree_ )
		return;

	{
		boost::mutex::scoped_lock lock( this->nodesMutex_ );

		int position = parent->getScrollPosition();
		for( size_t iPos = position; iPos < std::min(position + gridHeight_, nodes_.size()); ++iPos )
		{
			if( selected_ == iPos )
			{
				glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();

				int viewportWidth, viewportHeight;
				GetClientSize(&viewportWidth, &viewportHeight);
				glViewport( 0, 0, viewportWidth, viewportHeight );

				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();
				glOrtho(0.0, viewportWidth, 0.0, viewportHeight, -1, 1);

				int actualPos = iPos - position;
				int left = 0;
				int right = gridWidth_*(this->subImageWidth_ + border_) + border_;
				int top = viewportHeight - actualPos*(this->subImageHeight_ + 2*border_) - border_;
				int bottom = top - this->subImageHeight_ - 2*border_;


				GLDisableHandler depthTest( GL_DEPTH_TEST );
				GLDisableHandler lighting( GL_LIGHTING );
				{
					glColor3ub(255, 251, 149);
					GLActionHandler quads(GL_QUADS);
					glVertex2f(right, bottom);
					glVertex2f(right, top);
					glVertex2f(left, top);
					glVertex2f(left, bottom);
				}
				{
					glColor3ub(255, 246, 0);
					glLineWidth(2.0);
					GLActionHandler lines(GL_LINE_LOOP);
					glVertex2f(left, bottom);
					glVertex2f(left, top);
					glVertex2f(right, top);
					glVertex2f(right, bottom);
				}


				glLineWidth(1.0);

			}

			const Path* path = nodes_.at(iPos);
			float score = (*this->selectScore_)(path);

			std::vector<Simulation::State> states = 
				path->samples( gridWidth_ );
			assert( !states.empty() );
			size_t numStoredStates = states.size();

			for( size_t iTimeSlice = 0; iTimeSlice < gridWidth_; ++iTimeSlice )
			{
				if( iTimeSlice >= states.size() )
					continue;

				// Now we have the correct node to render.

				// need to use the correct view and projection matrices:
				Camera& cam = mainFrame->GetCanvas()->camera();
				cam.setRatio( 1.5 );
				cam.applyProjectiveTransform( 0 );
				cam.applyViewingTransform();

				{
					int viewportWidth, viewportHeight;
					GetClientSize(&viewportWidth, &viewportHeight);

					int actualPos = iPos - position;
					// now for the viewport transform.  
					glViewport( border_ + (gridWidth_ - iTimeSlice - 1)*(this->subImageWidth_ + border_), 
						viewportHeight - (actualPos+1)*(this->subImageHeight_ + 2*border_),
						subImageWidth_,
						subImageHeight_ );

					/*
					int viewportWidth, viewportHeight;
					GetClientSize(&viewportWidth, &viewportHeight);
					glViewport(0, 0, viewportWidth, viewportHeight);
					*/
				}

				PhysicsObject::RenderState renderState;
				renderState.alpha = 1.0;
				renderState.polygonMode = PhysicsObject::RenderState::FILL;

				size_t iState = numStoredStates - (iTimeSlice * numStoredStates) / gridWidth_ - 1;
				BestScoresRenderObjectHandler handler;
				mainFrame->walkTree( handler, states.at(iState) );

				{
					GLDisableHandler depthCheck(GL_DEPTH_TEST);
					glMatrixMode(GL_MODELVIEW);
					glLoadIdentity();

					glMatrixMode(GL_PROJECTION);
					glLoadIdentity();

					glColor4f( 0.1, 0.1, 0.1, 1.0 );
					GLActionHandler lines( GL_LINE_LOOP );
					glVertex3f( -1.0, -1.0, 0.0 );
					glVertex3f( -1.0, 0.999, 0.0 );
					glVertex3f( 0.999, 0.999, 0.0 );
					glVertex3f( 1.0, -1.0, 0.0 );
				}

				if( iTimeSlice == (gridWidth_ - 1) )
				{
					std::string rankString( "Rank: " );
					rankString.append( boost::lexical_cast<std::string>(iPos + 1) );

					std::string scoreString( "Score: " );
					scoreString.append( boost::lexical_cast<std::string>(score) );

					glColor4f( 1.0, 1.0, 1.0, 1.0 );
					BoundingBox3d box = printString(
						rankString.c_str(), 
						vl::Vec3d(-1.0, -1.0, 0.0), 
						ALIGN_HORIZ_LEFT,
						ALIGN_VERT_BASELINE,
						false );
					double height = (box.maximum() - box.minimum())[1];

					printString( 
						rankString.c_str(), 
						vl::Vec3d( -1.0 + 0.3*height, -1.0 + 0.3*height, 0.0 ),
						ALIGN_HORIZ_LEFT,
						ALIGN_VERT_BASELINE,
						true );
					printString( 
						rankString.c_str(), 
						vl::Vec3d( -1.0 + 0.3*height, -1.0 + 0.3*height, 0.0 ),
						ALIGN_HORIZ_LEFT,
						ALIGN_VERT_BASELINE,
						true );
					printString( 
						scoreString.c_str(), 
						vl::Vec3d( -1.0 + 0.3*height, -1.0 + 1.5*height, 0.0 ),
						ALIGN_HORIZ_LEFT,
						ALIGN_VERT_BASELINE,
						true );
				}
			}
		}
	}

	glFlush();
	SwapBuffers();
}

void BestScoresCanvas::OnMouse( wxMouseEvent& event )
{
	BestScoresPanel* parent = dynamic_cast<BestScoresPanel*>(this->GetParent());
	MainFrame* mainFrame = dynamic_cast<MainFrame*>(parent->GetParent()->GetParent()->GetParent());
	assert( mainFrame != 0 );

	int wheelRot = event.GetWheelRotation();
	currentWheelRot_ += wheelRot;
	if( currentWheelRot_ > event.GetWheelDelta() || currentWheelRot_ < -event.GetWheelDelta() )
	{
		int sign = (currentWheelRot_ > 0) ? 1 : -1;
		int numLines = sign*(sign*currentWheelRot_) / 120;
		currentWheelRot_ = sign*((sign*currentWheelRot_) % 120);
		parent->setScrollPosition( parent->getScrollPosition() + numLines );
	}

	if( event.ButtonUp() )
	{
		this->SetFocus();
		long vpos = event.GetY();
		int pos = vpos / (subImageHeight_ + 2*border_);

		size_t newSelected = pos + parent->getScrollPosition();
		if( newSelected < this->nodes_.size() )
		{
			selected_ = newSelected;
			boost::shared_ptr<Action> setPathAction( 
				new SetPathSelectionAction(
					std::make_pair(this->tree_, nodes_[newSelected]),
					std::make_pair(mainFrame->getTree(), mainFrame->getPath()),
					mainFrame ) );
			mainFrame->performAction( setPathAction );
			Refresh( FALSE );
		}
		return;
	}

//	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent()->GetParent());
	event.Skip();
	//Refresh( FALSE );
}

void BestScoresCanvas::OnSize(wxSizeEvent& event)
{
	BestScoresPanel* parent = dynamic_cast<BestScoresPanel*>(this->GetParent());
	// I don't understand why I need this:
	if( !parent )
		return;

	wxSize sz = event.GetSize();
	gridWidth_ = std::max<size_t>( (sz.GetWidth() - 2*border_) / (subImageWidth_ + border_), 1 );
	gridHeight_ = std::max<size_t>( (sz.GetHeight() - 2*border_) / (subImageHeight_ + 2*border_) + 1, 1 );

	{
		boost::mutex::scoped_lock lock( this->nodesMutex_ );
		parent->setScrollSize( nodes_.size(), gridHeight_ - 1);
	}

    wxGLCanvas::OnSize(event);
	Refresh( FALSE );
}

size_t BestScoresCanvas::numShown() const
{
	return this->gridHeight_ - 1;
}

std::deque<const Path*> BestScoresCanvas::nodes() const
{
	boost::mutex::scoped_lock lock( this->nodesMutex_ );
	std::deque<const Path*> result( this->nodes_ );
	return result;
}

void BestScoresCanvas::updateShowAllPaths()
{
	this->setNodes( tree_.get(), tree_->activePaths() );
}

void BestScoresCanvas::setNodes( const SimulationTree* tree, const boost::dynamic_bitset<>& nodes )
{
	BestScoresPanel* parent = dynamic_cast<BestScoresPanel*>(this->GetParent());
	MainFrame* mainFrame = dynamic_cast<MainFrame*>(parent->GetParent()->GetParent()->GetParent());

	bool showAllPaths = mainFrame->showAllPaths();

	std::deque<const Path*> newPaths;
	for( size_t i = 0; i < nodes.size(); ++i )
	{
		if( !nodes[i] && !showAllPaths )
			continue;

		newPaths.push_back( tree->path(i) );
	}

	SelectFunctionComparator comparator(this->selectScore_, parent->sortOrder());
	std::sort( newPaths.begin(), newPaths.end(), comparator );

	const Path* selectedNode = 0;
	{
		boost::mutex::scoped_lock lock( this->nodesMutex_ );
		int scrollPos = parent->getScrollPosition();
		if( this->selected_ < this->nodes_.size() )
			selectedNode = this->nodes_[this->selected_];
	}

	// need to find where our current selection has been moved to:
	size_t newSelected = selected_;
	for( size_t i = 0; i < newPaths.size(); ++i )
	{
		if( newPaths[i] == selectedNode )
		{
			newSelected = i;
			break;
		}
	}

	{
		boost::mutex::scoped_lock lock( this->nodesMutex_ );
		this->nodes_.swap( newPaths );
		this->selected_ = newSelected;
	}

	{
		wxCommandEvent event(wxEVT_CHANGE_NODE_COUNT, this->GetId());
		event.SetEventObject( this );
		event.SetInt( nodes_.size() );
		event.SetExtraLong( newSelected );
		::wxPostEvent( parent, event );
	}
}

void BestScoresCanvas::addNodes( const SimulationTree* tree, 
								const std::deque<const Path*>& inNodes, 
								const boost::dynamic_bitset<>& active )
{
	assert( tree == this->tree_.get() );
	BestScoresPanel* parent = dynamic_cast<BestScoresPanel*>(this->GetParent());
	MainFrame* mainFrame = dynamic_cast<MainFrame*>(parent->GetParent()->GetParent()->GetParent());

	bool showAllPaths = mainFrame->showAllPaths();

	std::deque<const Path*> nodes;
	if( showAllPaths )
	{
		nodes = inNodes;
	}
	else
	{
		for( size_t iPath = 0; iPath < inNodes.size(); ++iPath )
		{
			if( active[iPath] )
				nodes.push_back( inNodes[iPath] );
		}
	}

	SelectFunctionComparator comparator(this->selectScore_, parent->sortOrder());
	std::sort( nodes.begin(), nodes.end(), comparator );

	std::deque<const Path*> origNodes;
	{
		boost::mutex::scoped_lock lock( this->nodesMutex_ );
		origNodes = this->nodes_;
	}

	std::deque<const Path*> newNodes;
	std::merge( nodes.begin(), nodes.end(),
		origNodes.begin(), origNodes.end(),
		std::back_inserter(newNodes),
		comparator );

	// need to find where our current selection has been moved to:
	size_t newSelected = selected_;
	if( selected_ < origNodes.size() )
	{
		int lowPos = std::max<int>( boost::numeric_cast<int>(selected_) - boost::numeric_cast<int>(nodes.size()), 0 );
		int highPos = std::min<int>( boost::numeric_cast<int>(selected_) + boost::numeric_cast<int>(nodes.size()), newNodes.size()-1 );
		for( size_t i = lowPos; i <= highPos; ++i )
		{
			if( newNodes[i] == origNodes[selected_] )
			{
				newSelected = i;
				break;
			}
		}
	}

	int scrollPos = parent->getScrollPosition();
	if( scrollPos < origNodes.size() )
	{
		int lowPos = std::max<int>( scrollPos - boost::numeric_cast<int>(nodes.size()), 0 );
		int highPos = std::min<int>( scrollPos + boost::numeric_cast<int>(nodes.size()), newNodes.size()-1 );
		for( int i = lowPos; i <= highPos; ++i )
		{
			if( newNodes[i] == origNodes[scrollPos] )
			{
				scrollPos = i;
				break;
			}
		}
	}

	{
		boost::mutex::scoped_lock lock( this->nodesMutex_ );
		this->nodes_.swap( newNodes );
		this->selected_ = newSelected;
	}

	{
		wxCommandEvent event(wxEVT_CHANGE_NODE_COUNT, this->GetId());
		event.SetEventObject( this );
		event.SetInt( nodes_.size() );
		event.SetExtraLong( scrollPos );
		::wxPostEvent( parent, event );
	}
	/*
	parent->setScrollSize( newNodes.size(), gridHeight_ );
	if( newSelected != selected_ )
	{
		parent->setScrollPosition( scrollPos );
		selected_ = newSelected;
	}
	*/

//	this->Refresh( FALSE );
}

void BestScoresCanvas::setSortBy( size_t value )
{
	if( !updateObjectiveComputer() )
		return;

	this->updateSort();
}

bool BestScoresCanvas::updateObjectiveComputer()
{
	BestScoresPanel* parent = dynamic_cast<BestScoresPanel*>(this->GetParent());
	size_t sortBy = parent->sortBy();
	boost::shared_ptr<ObjectiveComputer> newSelectScore;
	if( sortBy >= tree_->objectives().size() )
	{
		newSelectScore.reset( new CountConstraints(this->tree_) );
	}
	else
	{
		MainFrame* mainFrame = dynamic_cast<MainFrame*>( 
			this->GetParent()->GetParent()->GetParent()->GetParent() );
		std::vector<size_t> dynamicObjectIds = 
			mainFrame->objectListToDynamicObjectIds( 
			mainFrame->getSelected(), this->tree_->dynamicObjects() );
		std::sort( dynamicObjectIds.begin(), dynamicObjectIds.end() );
		if( dynamicObjectIds.empty() )
		{
			dynamicObjectIds.resize( tree_->dynamicObjects().size() );
			for( size_t i = 0; i < dynamicObjectIds.size(); ++i )
				dynamicObjectIds[i] = i;
		}

		newSelectScore.reset(
			new SelectScore( 
				sortBy,
				dynamicObjectIds, 
				*this->tree_->objectives().at(sortBy) ) );
	}

	if( this->selectScore_ && *newSelectScore == *this->selectScore_ )
		return false;

	this->selectScore_ = newSelectScore;
		return true;
}


void BestScoresCanvas::updateSelected()
{
	if( !updateObjectiveComputer() )
		return;

	this->updateSort();
}

void BestScoresCanvas::updateSort()
{
	BestScoresPanel* parent = dynamic_cast<BestScoresPanel*>(this->GetParent());
	const Path* selected = 0;
	if( this->selected_ < nodes_.size() )
		selected = nodes_[this->selected_];

	{
		boost::mutex::scoped_lock lock( this->nodesMutex_ );
		std::sort( nodes_.begin(), nodes_.end(), SelectFunctionComparator(this->selectScore_, parent->sortOrder()) );

		for( size_t i = 0; i < nodes_.size(); ++i )
		{
			if( nodes_[i] == selected )
				selected_ = i;
		}
	}

	this->Refresh( FALSE );
}

void BestScoresCanvas::setPath( const Path* path )
{
	if( selected_ < nodes_.size() && nodes_[selected_] == path )
		return;

	for( size_t i = 0; i < this->nodes_.size(); ++i )
	{
		if( this->nodes_[i] == path )
		{
			selected_ = i;
			break;
		}
	}

	BestScoresPanel* parent = dynamic_cast<BestScoresPanel*>(this->GetParent());
	parent->setScrollPosition( this->selected_ );	
	this->Refresh(FALSE);
}

} // namespace planning

