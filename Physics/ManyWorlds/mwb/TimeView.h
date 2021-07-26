#ifndef __TIMEVIEW_H__
#define __TIMEVIEW_H__

#include <wx/defs.h>
#include <wx/dcbuffer.h>
#include <wx/dcclient.h>

namespace planning {

class MouseCapture;


class TimeView : public wxPanel
{
public:
	TimeView(wxWindow* parent, const wxWindowID id = -1, const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize, long style = 0, const wxString& name = wxT("TimeView"));
	~TimeView();

    void OnPaint(wxPaintEvent& event);
    void OnSize(wxSizeEvent& event);
    void OnEraseBackground(wxEraseEvent& event);
    void OnMouse( wxMouseEvent& event );

	void OnKey(wxKeyEvent& event);

	std::pair<double, double> activeTimeRange() const;

private:
	void setCurrentTime( double time, bool move );
	int timeToPosition( double time );

	std::pair<double, double> timeRange_;
	double scale_;

	enum MouseAction
	{
		MOUSE_NOTHING,
		MOUSE_TRANSLATE,
		MOUSE_SCALE,
		MOUSE_CHANGING,
		MOUSE_MOVE_LEFT_ACTIVE_ENDPOINT,
		MOUSE_MOVE_RIGHT_ACTIVE_ENDPOINT,
		MOUSE_MOVE_CENTER_GRIPPERS,
        MOUSE_ADDITIVE_CONSTRAINT,
        MOUSE_SUBTRACTIVE_CONSTRAINT,
        MOUSE_REFINE,
	} mouseAction_;

	vl::Vec2d clickedPosition_;
    vl::Vec2d currentPosition_;
	std::pair<double, double> clickedTimeRange_;
	std::pair<double, double> clickedActiveTimeRange_;

	// marks the range that has been active for e.g. rendering
	std::pair<double, double> activeTimeRange_;
	double increment_;
	double backupTime_;
	unsigned int fps_;

	boost::scoped_ptr<MouseCapture> mouseCapture_;

    DECLARE_EVENT_TABLE()
};

} // namespace planning

#endif