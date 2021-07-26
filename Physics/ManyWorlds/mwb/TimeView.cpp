#include "stdafx.h"

#include "TimeView.h"
#include "planningUI.h"

#include <iostream>
#include <iomanip>

namespace planning {

BEGIN_EVENT_TABLE(TimeView, wxPanel)
    EVT_SIZE(TimeView::OnSize)
    EVT_PAINT(TimeView::OnPaint)
    EVT_ERASE_BACKGROUND(TimeView::OnEraseBackground)
    EVT_MOUSE_EVENTS(TimeView::OnMouse)
	EVT_KEY_DOWN(TimeView::OnKey)
END_EVENT_TABLE()


TimeView::TimeView(wxWindow* parent, const wxWindowID id, const wxPoint& pos,
    const wxSize& size, long style, const wxString& name)
	:	wxPanel(parent, id, pos, size, style, name), fps_(30)
{
	increment_ = 1.0/boost::numeric_cast<double>(fps_);

	timeRange_ = std::make_pair( -1.0, 30.0 );
	activeTimeRange_ = std::make_pair( 0.0, 10.0 );
	mouseAction_ = MOUSE_NOTHING;
}

TimeView::~TimeView()
{
}

std::pair<double, double> TimeView::activeTimeRange() const
{
	return this->activeTimeRange_;
}

void TimeView::OnKey(wxKeyEvent& event)
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());

	switch( event.GetKeyCode() )
	{
	case WXK_LEFT:
	case WXK_NUMPAD_LEFT:
		this->setCurrentTime( parent->getTime() - increment_, false );
		break;
	case WXK_RIGHT:
	case WXK_NUMPAD_RIGHT:
		this->setCurrentTime( parent->getTime() + increment_, false );
		break;
	}
}


const int selectedRangeHeight = 10;
const int selectedPolygonWidth = 6;
const int minCenterGrippersSize = 20;

void TimeView::OnPaint(wxPaintEvent& event)
{
	const MainFrame* parent = dynamic_cast<const MainFrame*>(this->GetParent());

	// must always be here
	wxBufferedPaintDC dc(this);

	dc.Clear();
	assert( timeRange_.first < timeRange_.second );

	wxSize size = dc.GetSize();

	double diff = timeRange_.second - timeRange_.first;
	scale_ = boost::numeric_cast<double>(size.GetWidth()) / diff;

	double minBlock = 10.0;
	double minLabelledBlock = 100.0;

	wxFont font( 8, wxFONTFAMILY_MODERN, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_LIGHT, false, wxString(), wxFONTENCODING_SYSTEM );
	dc.SetFont(font);
	dc.SetMapMode( wxMM_TEXT );

	bool labeled = false;

	size_t maxLineHeight = 12;

	// line separating the top half from the bottom
	// color the bottom half gray

	wxPen lightGrayPen( wxColor(212, 208, 200) );
	wxPen darkGrayPen( wxColor(128, 128, 128) );
	wxPen paleGrayPen( wxColor(229, 229, 229) );

	wxBrush lightGrayBrush( wxColor(212, 208, 200) );
	wxBrush darkGrayBrush( wxColor(128, 128, 128) );
	wxBrush paleGrayBrush( wxColor(229, 229, 229) );

	dc.SetBrush( lightGrayBrush );
	dc.DrawRectangle( 0, selectedRangeHeight, size.GetWidth(), size.GetHeight() );

	dc.SetPen( *wxBLACK_PEN );
	dc.DrawLine( 0, selectedRangeHeight, size.GetWidth(), selectedRangeHeight );


	int activeRangeLeft = timeToPosition(activeTimeRange_.first);
	int activeRangeRight = timeToPosition(activeTimeRange_.second);
	if( activeRangeLeft == activeRangeRight )
		activeRangeRight += 1;

	// draw selected region
	{
		int xPos = activeRangeLeft;
		wxPoint points[] = 
			{ 
				wxPoint(xPos, 0), 
				wxPoint(xPos, selectedRangeHeight),
				wxPoint(xPos - selectedPolygonWidth, selectedRangeHeight),
				wxPoint(xPos - selectedPolygonWidth, selectedPolygonWidth)
			};
		dc.SetBrush( lightGrayBrush );
		dc.SetPen( *wxBLACK_PEN );
		dc.DrawPolygon(4, points);

		wxPoint innerPoints[] = 
		{ 
			wxPoint(xPos - 1, 2), 
			wxPoint(xPos - 1, selectedRangeHeight - 1),
			wxPoint(xPos - selectedPolygonWidth + 1, selectedRangeHeight - 1),
			wxPoint(xPos - selectedPolygonWidth + 1, selectedPolygonWidth)
		};

		dc.SetPen( *wxWHITE_PEN );
		dc.DrawLine(innerPoints[2].x, 
			innerPoints[2].y,
			innerPoints[3].x, 
			innerPoints[3].y );
		dc.DrawLine(innerPoints[3].x, 
			innerPoints[3].y,
			innerPoints[0].x, 
			innerPoints[0].y );

		dc.SetPen( darkGrayPen );
		dc.DrawLine(innerPoints[0].x, 
			innerPoints[0].y,
			innerPoints[1].x, 
			innerPoints[1].y );
		dc.DrawLine(innerPoints[1].x, 
			innerPoints[1].y,
			innerPoints[2].x, 
			innerPoints[2].y );
	}

	{
		int xPos = activeRangeRight;
		wxPoint points[] = 
			{ 
				wxPoint(xPos, 0), 
				wxPoint(xPos + selectedPolygonWidth, selectedPolygonWidth),
				wxPoint(xPos + selectedPolygonWidth, selectedRangeHeight),
				wxPoint(xPos, selectedRangeHeight)
			};
		dc.SetBrush( lightGrayBrush );
		dc.SetPen( *wxBLACK_PEN );
		dc.DrawPolygon(4, points);

		wxPoint innerPoints[] = 
		{ 
			wxPoint(xPos + 1, 1), 
			wxPoint(xPos + selectedPolygonWidth - 1, selectedPolygonWidth + 1),
			wxPoint(xPos + selectedPolygonWidth - 1, selectedRangeHeight - 1),
			wxPoint(xPos + 1, selectedRangeHeight - 1)
		};

		dc.SetPen( *wxWHITE_PEN );
		dc.DrawLine(innerPoints[3].x, 
			innerPoints[3].y,
			innerPoints[0].x, 
			innerPoints[0].y );

		dc.SetPen( darkGrayPen );
		dc.DrawLine(innerPoints[1].x, 
			innerPoints[1].y,
			innerPoints[2].x, 
			innerPoints[2].y );
		dc.DrawLine(innerPoints[2].x, 
			innerPoints[2].y,
			innerPoints[3].x, 
			innerPoints[3].y );
	}

	// draw in-between region
	{
		int left = activeRangeLeft;
		int right = activeRangeRight;
		dc.SetPen( *wxTRANSPARENT_PEN );
		dc.SetBrush( paleGrayBrush );
		dc.DrawRectangle( left + 1,
			0,
			right - left - 1,
			selectedRangeHeight );

		if( right - left > minCenterGrippersSize )
		{
			int center = timeToPosition(0.5*(activeTimeRange_.first + activeTimeRange_.second));
			int left = center - 5;
			for( size_t i = 0; i < 4; ++i )
			{
				int startPos = left + i*3;
				dc.SetPen( *wxWHITE_PEN );
				dc.DrawLine( startPos, 2, startPos, selectedRangeHeight - 3 );

				dc.SetPen( darkGrayPen );
				dc.DrawLine( startPos+1, 3, startPos+1, selectedRangeHeight - 2 );
			}
		}
	}

	dc.SetPen( *wxBLACK_PEN );
	for( int exponent = -2; exponent < 2; ++exponent )
	{
		for( int offset = 1; offset <= 5; offset += 4 )
		{
			double unit = boost::numeric_cast<double>(offset)*pow( 10.0, boost::numeric_cast<double>(exponent) );
			if( scale_ * unit < minBlock )
				continue;

			double lineLength = pow( 2.0, boost::numeric_cast<double>(exponent-2) );
			if( offset == 5 )
				lineLength *= 1.5;

			int start = timeRange_.first / unit - 1;
			int end = timeRange_.second / unit + 1;

			{
				for( int i = start; i <= end; ++i )
				{
					double xPos = scale_ * (boost::numeric_cast<double>(i) * unit - timeRange_.first);
					dc.DrawLine( xPos, lineLength*maxLineHeight + selectedRangeHeight, 
						xPos, selectedRangeHeight );
				}
			}
			if( scale_ * unit > minLabelledBlock && !labeled )
			{
				labeled = true;
				for( int i = start; i <= end; ++i )
				{
					double xPos = static_cast<float>(i) * unit;
					std::ostringstream oss;
					oss << std::fixed << std::setprecision( std::max<int>(1.0, -exponent) ) << xPos;
					wxString str( wxT(oss.str().c_str()) );
					wxCoord w, h;
					dc.GetTextExtent( str, &w, &h );
					dc.DrawText( str, scale_*(xPos - timeRange_.first) - w/2, size.GetHeight() - h );
				}
			}
		}
	}

	{
		wxPen redPen( *wxRED, 3 );
		dc.SetPen( redPen );
		
		int xPos = timeToPosition(parent->getTime());
		dc.DrawLine( xPos, selectedRangeHeight, xPos, size.GetHeight() );
	}

    if( this->mouseAction_ == MOUSE_ADDITIVE_CONSTRAINT ||
        this->mouseAction_ == MOUSE_SUBTRACTIVE_CONSTRAINT ||
        this->mouseAction_ == MOUSE_REFINE )
    {
        wxColour color;
        switch( this->mouseAction_ )
        {
        case MOUSE_ADDITIVE_CONSTRAINT:
            color = wxColour(40, 255, 40);
            break;
        case MOUSE_SUBTRACTIVE_CONSTRAINT:
            color = wxColour(255, 40, 40);
            break;
        case MOUSE_REFINE:
            color = wxColour(255, 255, 40);
            break;
        }

        // Need to draw a box
        wxPen pen( color, 3 );
        wxBrush brush( color, wxBDIAGONAL_HATCH );
        dc.SetBrush( brush );
        dc.SetPen( pen );

        double minVal = clickedPosition_[0];
        double maxVal = currentPosition_[0];
        if( minVal > maxVal )
            std::swap( minVal, maxVal );

		dc.DrawRectangle( minVal,
			selectedRangeHeight + 1,
			maxVal - minVal,
			size.GetHeight() - selectedRangeHeight - 1 );
    }
}

int TimeView::timeToPosition( double time )
{
	return scale_ * (time - timeRange_.first);
}


void TimeView::OnSize(wxSizeEvent& event)
{
    // this is also necessary to update the context on some platforms
    wxPanel::OnSize(event);
    
    {
//        SetCurrent();
		this->Refresh(FALSE);
    }

}

void TimeView::OnEraseBackground(wxEraseEvent& event)
{

}

void TimeView::OnMouse( wxMouseEvent& event )
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());
	double currentTime = event.GetX() / scale_ + timeRange_.first;

    currentPosition_ = vl::Vec2d( event.GetX(), event.GetY() );

	double clickedTime = clickedPosition_[0] / scale_ + clickedTimeRange_.first;
	double newTime = event.GetX() / scale_ + clickedTimeRange_.first;

	if( event.ButtonDown() )
	{
		this->SetFocus();

		mouseCapture_.reset( new MouseCapture(this) );

		clickedPosition_ = vl::Vec2d( event.GetX(), event.GetY() );
		clickedTimeRange_ = timeRange_;
		clickedActiveTimeRange_ = activeTimeRange_;
		this->backupTime_ = parent->getTime();

		int activeRangeStart = timeToPosition( this->activeTimeRange_.first );
		int activeRangeEnd = timeToPosition( this->activeTimeRange_.second );
		int activeRangeMid = timeToPosition( 0.5*(this->activeTimeRange_.second+this->activeTimeRange_.first) );

		BoundingBox2d activeRangeStartBox( vl::Vec2d(activeRangeStart - selectedPolygonWidth, 0),
			vl::Vec2d(activeRangeStart, selectedRangeHeight) );
		BoundingBox2d activeRangeEndBox( vl::Vec2d(activeRangeEnd, 0),
			vl::Vec2d(activeRangeEnd + selectedPolygonWidth, selectedRangeHeight) );
		BoundingBox2d centerGrippersBox( vl::Vec2d(activeRangeMid - 10, 0),
			vl::Vec2d(activeRangeMid + 10, selectedRangeHeight) );

		if( contains( activeRangeStartBox, clickedPosition_ ) )
			mouseAction_ = MOUSE_MOVE_LEFT_ACTIVE_ENDPOINT;
		else if( contains( activeRangeEndBox, clickedPosition_ ) )
			mouseAction_ = MOUSE_MOVE_RIGHT_ACTIVE_ENDPOINT;
		else if( activeRangeEnd - activeRangeStart > minCenterGrippersSize 
			&& contains( centerGrippersBox, clickedPosition_ ) )
			mouseAction_ = MOUSE_MOVE_CENTER_GRIPPERS;
		else if( event.AltDown() )
		{
			if( event.LeftIsDown() )
				mouseAction_ = MOUSE_NOTHING;
			else if( event.MiddleIsDown() )
				mouseAction_ = MOUSE_TRANSLATE;
			else if( event.RightIsDown() )
				mouseAction_ = MOUSE_SCALE;
		}
		else if( event.ControlDown() )
		{
		}
        else if( parent->getMode() == MainFrame::MODE_SOLVE && 
            (parent->mouseMode() == MainFrame::MODE_MOUSE_ADDITIVE_CONSTRAINT ||
            parent->mouseMode() == MainFrame::MODE_MOUSE_SUBTRACTIVE_CONSTRAINT ||
            parent->mouseMode() == MainFrame::MODE_MOUSE_REFINE) )
        {
            switch( parent->mouseMode() )
            {
            case MainFrame::MODE_MOUSE_ADDITIVE_CONSTRAINT:
                this->mouseAction_ = MOUSE_ADDITIVE_CONSTRAINT;
                break;
            case MainFrame::MODE_MOUSE_SUBTRACTIVE_CONSTRAINT:
                this->mouseAction_ = MOUSE_SUBTRACTIVE_CONSTRAINT;
                break;
            case MainFrame::MODE_MOUSE_REFINE:
                this->mouseAction_ = MOUSE_REFINE;
                break;
            default:
                break;
            }
        }
        else
		{
			mouseAction_ = MOUSE_CHANGING;
			setCurrentTime( currentTime, event.Button(wxMOUSE_BTN_LEFT) );
		}
	}
	else if( event.Dragging() )
	{
	    switch( mouseAction_ )
	    {
	    case MOUSE_TRANSLATE:
		    {
			    timeRange_ = clickedTimeRange_;
			    double translate = newTime - clickedTime;
			    timeRange_.second -= translate;
			    timeRange_.first -= translate;
			    break;
		    }

	    case MOUSE_SCALE:
		    {
			    double clickedTime = clickedPosition_[0] / scale_ - clickedTimeRange_.first;
			    double newTime = event.GetX() / scale_ - clickedTimeRange_.first;

			    double middle = 0.5*(clickedTimeRange_.first + clickedTimeRange_.second);
			    std::pair<double, double> offsetRange( clickedTimeRange_.first - middle, 
				    clickedTimeRange_.second - middle );
			    double scale = 1.0 + (clickedTime - newTime) / (clickedTimeRange_.second - clickedTimeRange_.first);
			    timeRange_ = std::make_pair( middle + scale*offsetRange.first, 
				    middle + scale*offsetRange.second );
			    break;
		    }

	    case MOUSE_CHANGING:
		    {
			    setCurrentTime( currentTime, event.LeftIsDown() );
			    break;
		    }

	    case MOUSE_MOVE_LEFT_ACTIVE_ENDPOINT:
		    {
			    this->activeTimeRange_.first = std::min( currentTime, this->activeTimeRange_.second );
			    break;
		    }

	    case MOUSE_MOVE_RIGHT_ACTIVE_ENDPOINT:
		    {
			    this->activeTimeRange_.second = std::max( currentTime, this->activeTimeRange_.first );
			    break;
		    }

	    case MOUSE_MOVE_CENTER_GRIPPERS:
		    {
			    double translate = newTime - clickedTime;
			    this->activeTimeRange_ = std::make_pair(
				    clickedActiveTimeRange_.first + translate,
				    clickedActiveTimeRange_.second + translate );
			    break;
		    }

	    case MOUSE_NOTHING:
		    break;
	    }
	}
	else if( event.ButtonUp() )
	{
		mouseCapture_.reset();

		if( mouseAction_ == MOUSE_CHANGING )
		{
			setCurrentTime( currentTime, event.Button(wxMOUSE_BTN_LEFT) );
			ActionPtr action( new SetTimeAction(
				parent, this->backupTime_, parent->getTime() ) );
			parent->performAction( action );
		}
		else if( mouseAction_ == MOUSE_MOVE_LEFT_ACTIVE_ENDPOINT
			|| mouseAction_ == MOUSE_MOVE_RIGHT_ACTIVE_ENDPOINT
			|| mouseAction_ == MOUSE_MOVE_CENTER_GRIPPERS )
		{
			parent->updateActiveTimeRange( this->activeTimeRange_ );
		}
        else if( mouseAction_ == MOUSE_REFINE )
        {
            parent->applyTimeSelection( clickedTime, newTime );
        }

		mouseAction_ = MOUSE_NOTHING;
	}

	Refresh( FALSE );
}

void TimeView::setCurrentTime( double currentTime, bool changing )
{
	/*
	typedef boost::numeric::converter<double, long, 
		boost::numeric::conversion_traits<double, long>, 
		boost::numeric::def_overflow_handler,
		boost::numeric::RoundEven<double> > Converter;

	Converter converter;
	double currentFrame = 
		currentTime / increment_;
	long curTime = converter.convert( currentFrame );
	currentTime_ = std::max( static_cast<double>(curTime) * increment_, 0.0 );
	*/

	double currentFrame = 
		currentTime / increment_;
	double intPart;
	double frac = std::modf( currentFrame, &intPart );
	if( frac > 0.5 )
		intPart += 1.0;

	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());
	parent->setTime( intPart * increment_, changing );
}

} // namespace planning
