
#ifndef __GLVIEW_H__
#define __GLVIEW_H__

#include "novodexWrappers.h"
#include "physicsFwd.h"
#include "action.h"
#include "threads.h"

#include "twigg/keyable.h"
#include "twigg/camera.h"

#include <wx/glcanvas.h>

class TileRendererWrapper;

namespace planning {

class MouseCapture;

class CameraTransformConverter
	: public ScreenSpaceConverter
{
public:
	CameraTransformConverter( const Camera& camera, int viewportWidth, int viewportHeight );
	vl::Vec3d toScreenSpace( const vl::Vec3d& worldSpace ) const;
	vl::Vec3d toWorldSpace( const vl::Vec3d& screenSpace ) const;

private:
	vl::Mat4d projectiveTransform_;
	vl::Mat4d modelViewTransform_;
	vl::Mat4d invModelViewTransform_;
	vl::Mat4d invProjectiveTransform_;
	float viewportWidth_;
	float viewportHeight_;
};

class GLView: public wxGLCanvas
{
public:
 	typedef enum { PerspectiveCamera, OrthographicCamera } CameraType;

	
	GLView(wxWindow *parent, const wxWindowID id = -1, const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize, long style = 0, const wxString& name = wxT("GLView"));
 	GLView(wxWindow *parent, wxGLContext* sharedContext, const wxWindowID id = -1, const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize, long style = 0, const wxString& name = wxT("GLView"));
	~GLView();

    void OnPaint(wxPaintEvent& event);
    void OnSize(wxSizeEvent& event);
    void OnEraseBackground(wxEraseEvent& event);
    void OnMouse( wxMouseEvent& event );
	void OnKey(wxKeyEvent& event);
	void OnKeyUp(wxKeyEvent& event);

	void setupViewport() const;
	void clearCache();

	const Camera& camera() const;
	Camera& camera();
	void setCamera( CameraWrapperPtr camera );
	CameraWrapperPtr cameraWrapper();

	void errMsg(const std::string& message);
	float fps() const;

	void dumpFrame( 
		int width, 
		int height, 
		const std::string& filename );

	void dumpFrame( 
		int width, 
		int height, 
		BYTE* buffer );

	// Useful for the dumpPaths stuff: we want to be able to dump just the
	// dynamic objects at the specified states
	wxImage dumpFrame(
		int width,
		int height,
		const std::vector< Simulation::State >& states );

	void nextCameraKey();
	void prevCameraKey();
	void setCameraKey();
	void deleteCameraKey();

private:
	int viewportWidth_;
	int viewportHeight_;

	CameraWrapperPtr camera_;

	boost::scoped_ptr<TileRendererWrapper> tileRenderer_;

	std::deque<time_t> frameTimes_;
	double frameTime_;

	bool headsUp_;

	boost::scoped_ptr<MouseCapture> mouseCapture_;

	bool selecting_;
	vl::Vec2 selectingAnchor_;
	vl::Vec2 selectingFree_;

    std::vector<vl::Vec2f> sketchPoints_;
    std::vector<LARGE_INTEGER> times_;

	bool spacePressed_;
	time_t keyPressTime_;

	// HACK HACK HACK
	std::vector< Simulation::State > statesToPaint_;

	/*
	// these are for debugging the select code
	vl::Vec3d nearTestPt_;
	vl::Vec3d farTestPt_;
	*/

	DECLARE_EVENT_TABLE()
};

} // namespace planning

#endif

