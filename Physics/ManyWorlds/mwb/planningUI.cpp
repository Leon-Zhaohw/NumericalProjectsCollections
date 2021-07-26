#include "stdafx.h"

#include "physicsObject.h"
#include "scene.h"
#include "planningUI.h"
#include "parser.h"
#include "action.h"
#include "mechModel.h"
#include "constraints.h"
#include "simulationTree.h"
#include "sphere_tree.h"
#include "ribFile.h"
#include "quat.h"
#include "compression.h"
#include "threads.h"
#include "OBB.h"
#include "renderLines.h"
#include "sound.h"

#include "GLView.h"
#include "TimeView.h"
#include "bestScoresView.h"
#include "hierarchyView.h"

#include "ImfRgbaFile.h"
#include "half.h"
#ifdef DRAW_LINES
#include <cairo.h>
#endif

#include <crtdbg.h>
#include <cmath>

#include "twigg/ioutil.h"
#include "twigg/imagelib.h"
#include "twigg/renderUtils.h"
#include "twigg/gslWrappers.h"
#include "twigg/prt.h"
#include "twigg/imagelib.h"
#include "twigg/objfile.h"
#include "twigg/tr.h"
#include "twigg/magicWrappers.h"
#include "twigg/stlext.h"

#include "Wm4Ray3.h"
#include "Wm4DistVector2Line2.h"
#include "Wm4DistVector2Segment2.h"
#include "Wm4DistLine3Line3.h"
#include "Wm4DistLine3Segment3.h"
#include "Wm4DistLine3Ray3.h"
#include "Wm4DistLine3Circle3.h"
#include "Wm4DistVector3Circle3.h"
#include "Wm4DistVector3Line3.h"

/*
#ifdef WIN32
#include "Psapi.h"
#endif
*/

#include <wx/mstream.h>
#include <wx/dcbuffer.h>
#include <wx/progdlg.h>
#include <wx/msgdlg.h>
#include <wx/toolbar.h>
#include <wx/filedlg.h>
#include <wx/textdlg.h>
#include <wx/stattext.h>

#include <boost/logic/tribool.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/bind.hpp>
#include <boost/mem_fn.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/once.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/graph/connected_components.hpp>

#include <GL/glut.h>

#ifdef WINDOWS_MEDIA
#include "Mmreg.h"
#endif


#include <cassert>
#include <iostream>
#include <iomanip>
#include <bitset>


#include "bitmaps/play.xpm"
#include "bitmaps/pause.xpm"
#include "bitmaps/greenRect.xpm"
#include "bitmaps/redRect.xpm"
#include "bitmaps/selectRect.xpm"
#include "bitmaps/refine.xpm"
#include "bitmaps/arc.xpm"

#include "bitmaps/cursor.xpm"
#include "bitmaps/translate.xpm"
#include "bitmaps/rotate.xpm"
#include "bitmaps/scale.xpm"
#include "bitmaps/changePivot.xpm"
#include "bitmaps/points.xpm"

#include "bitmaps/splash.h"

#if defined(__WXMSW__) && !defined(__WXWINCE__)
#pragma comment(linker, "\"/manifestdependency:type='win32' name='Microsoft.Windows.Common-Controls' version='6.0.0.0' processorArchitecture='X86' publicKeyToken='6595b64144ccf1df'\"")
#endif

DEFINE_EVENT_TYPE(wxEVT_FINISHED_COMPUTING_MODES)
DEFINE_EVENT_TYPE(wxEVT_LOG_ERROR)
DEFINE_EVENT_TYPE(wxEVT_REPAINT_ALL_VIEWS)

namespace planning {

// todo fix this hack hack hack
void changeDirectory( const std::string& filename )
{
	boost::filesystem::path fn( filename );
	boost::filesystem::path filePath( fn.branch_path() );
	if( !filePath.empty() )
	{
#ifdef WIN32
		std::string win32String = filePath.native_file_string();
		_chdir( win32String.c_str() );
#endif
	}
}

size_t numberOfProcessors = 1;

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


template <typename Pr>
class Compare1st
{
public:
	bool operator()( const Pr& left, const Pr& right ) const
	{
		return std::less<typename Pr::first_type>()( left.first, right.first );
	}
};

#ifdef LOGGING
Logger::Logger( const std::string& filename )
	: filename_(filename), ofs_(filename.c_str())
{
	if( !ofs_ )
		exit(1);
	startTime_ = GetTickCount();
	this->log( "Starting." );
}

void Logger::log( const std::string& message )
{
	boost::mutex::scoped_lock lock( this->lock_ );
	DWORD currentTime = GetTickCount() - startTime_;
	DWORD ms = currentTime % 1000;
	DWORD s = (currentTime / 1000) % 60;
	DWORD m = currentTime / (60*1000);
	ofs_ << m << ":" << s << "." << ms << " " << message << "\n";
}

Logger::~Logger()
{
	this->log( "Exiting." );
}
#endif

typedef std::basic_string<WCHAR> wchar_string;

// Convert to Unicode string using UTF-8
wchar_string toWideString( const std::string& str )
{
	if( str.empty() )
		return wchar_string();

	const size_t size = str.size();
	int requiredSize = MultiByteToWideChar(
		CP_UTF8,            // not sure about this
		0,
		&str[0],
		size,
		NULL,
		0 );
	assert( requiredSize > 0 );

	wchar_string result( requiredSize, ' ' );
	int error = MultiByteToWideChar(
		CP_UTF8,            // not sure about this
		0,
		&str[0],
		size,
		&result[0],
		requiredSize );
	assert( error > 0 );

	return result;
}

// Convert to Unicode string using UTF-8
std::string toString( const wchar_string& str )
{
	if( str.empty() )
		return std::string();

	const size_t size = str.size();
	int requiredSize = WideCharToMultiByte(
		CP_UTF8,            // not sure about this
		0,
		&str[0],
		size,
		0,
		0,
		NULL,
		NULL );
	assert( requiredSize > 0 );

	std::string result( requiredSize, ' ' );
	int error = WideCharToMultiByte(
		CP_UTF8,            // not sure about this
		0,
		&str[0],
		size,
		&result[0],
		requiredSize,
		NULL,
		NULL );
	assert( error > 0 );

	return result;
}


MouseCapture::MouseCapture( wxWindow* window )
	: window_(window)
{
#ifdef WIN32
	// We need to do this ourselves because in certain cases 
	// the OS may take away mouse capture and wxWidgets blows
	// up if we try to release capture without having it
	HWND handle = (HWND) window_->GetHandle();
	SetCapture( handle );
#else
	window_->CaptureMouse();
#endif
}

MouseCapture::~MouseCapture()
{
#ifdef WIN32
	ReleaseCapture();
#else
	window_->ReleaseMouse();
#endif
}

void errMsg(const std::string& message, wxWindow* parent )
{
//	assert( parent != 0 );
	boost::shared_ptr<wxMessageDialog> dlg( new wxMessageDialog(parent,
		message.c_str(), 
		"Error reading file",
		wxOK | wxICON_ERROR),
		boost::mem_fn(&wxMessageDialog::Destroy) );
	dlg.get()->ShowModal();
}

#ifdef USE_NOVODEX
class wxNxOutputStream
	: public NxUserOutputStream
{
	void reportError (NxErrorCode code, const char *message, const char *file, int line)
	{
		//this should be routed to the application
		//specific error handling. If this gets hit
		//then you are in most cases using the SDK
		//wrong and you need to debug your code!
		//however, code may just be a warning or
		//information.
		if (code < NXE_DB_INFO)
			wxLogError("NX: error in file %s, line %d: %s", file, line, message);
		else
			wxLogWarning("NX: warning in file %s, line %d: %s", file, line, message);
	}

	NxAssertResponse reportAssertViolation (const char *message, const char *file, int line)
	{
		//this should not get hit by
		// a properly debugged SDK!
		assert(0);
		return NX_AR_BREAKPOINT;
	}

	void print( const char *message )
	{
		wxLogVerbose(message);
	}
};
#endif

class ProgressDialogUpdater
	: public ProgressUpdater
{
public:
	ProgressDialogUpdater( const std::string& title, wxWindow* parent )
		: parent_(parent), title_(title)
	{
	}

	void create( const std::string& message, int maxValue )
	{
		progressDialog_.reset(
			new wxProgressDialog(
				wxT(title_.c_str()),
				wxT(message.c_str()),
				maxValue,
				parent_),
			boost::mem_fn(&wxProgressDialog::Destroy) );
	}

	void update( int value )
	{
		progressDialog_->Update( value );
	}

private:
	wxWindow* parent_;
	std::string title_;
	boost::shared_ptr<wxProgressDialog> progressDialog_;
};

std::vector<size_t> MainFrame::objectListToDynamicObjectIds( 
	const std::deque<SceneGraphElementPtr>& objects,
	const std::vector<ConstPhysicsObjectPtr>& dynamicObjects ) const
{
	std::vector<size_t> dynamicObjectIds;
	dynamicObjectIds.reserve( std::min(dynamicObjects.size(), objects.size()) );
	for( std::deque<SceneGraphElementPtr>::const_iterator iter = objects.begin();
		iter != objects.end(); ++iter )
	{
		SceneGraphElementPtr current = *iter;

		// find object with the same name in the dynamic object list
		std::vector<ConstPhysicsObjectPtr>::const_iterator currentItr = 
			std::find_if( dynamicObjects.begin(), dynamicObjects.end(), NameEquals(current->name()) );
		if( currentItr == dynamicObjects.end() )
			continue;

		if( isModifiedFromTree(current) )
			continue;

		dynamicObjectIds.push_back( std::distance(dynamicObjects.begin(), currentItr) );
	}

	std::sort( dynamicObjectIds.begin(), dynamicObjectIds.end() );
	return dynamicObjectIds;
}


BEGIN_EVENT_TABLE(MainFrame, wxFrame)
	EVT_IDLE(MainFrame::OnIdle)
#ifdef SOUND
	EVT_COMMAND(wxID_ANY, wxEVT_FINISHED_COMPUTING_MODES, MainFrame::OnFinishedComputingModes)
#endif
	EVT_COMMAND(wxID_ANY, wxEVT_LOG_ERROR, MainFrame::OnLogError)
	EVT_COMMAND(wxID_ANY, wxEVT_REPAINT_ALL_VIEWS, MainFrame::OnRepaintAllViews)
	EVT_TOOL(myID_CURSOR_TOOLBAR, MainFrame::OnChangeObjectMode)
	EVT_TOOL(myID_TRANSLATE_TOOLBAR, MainFrame::OnChangeObjectMode)
	EVT_TOOL(myID_ROTATE_TOOLBAR, MainFrame::OnChangeObjectMode)
	EVT_TOOL(myID_SCALE_TOOLBAR, MainFrame::OnChangeObjectMode)
	EVT_TOOL(myID_CHANGEPIVOT_TOOLBAR, MainFrame::OnChangeObjectMode)
	EVT_TOOL(myID_SELECTPOINTS_TOOLBAR, MainFrame::OnChangeObjectMode)
	EVT_MENU(MainFrame::MENU_FILE_NEW, MainFrame::OnMenuFileNew)
	EVT_MENU(MainFrame::MENU_FILE_OPEN, MainFrame::OnMenuFileOpen)
	EVT_MENU(MainFrame::MENU_FILE_IMPORT, MainFrame::OnMenuFileImport)
	EVT_MENU(MainFrame::MENU_FILE_SAVE, MainFrame::OnMenuFileSave)
	EVT_MENU(MainFrame::MENU_FILE_SAVEAS, MainFrame::OnMenuFileSaveAs)
	EVT_MENU(MainFrame::MENU_FILE_SAVE_TREE, MainFrame::OnMenuFileSaveTree)
	EVT_MENU(MainFrame::MENU_FILE_LOAD_TREE, MainFrame::OnMenuFileLoadTree)
	EVT_MENU(MainFrame::MENU_FILE_LOAD_DEBUGGING_INFO, MainFrame::OnMenuFileLoadDebuggingInfo)
    EVT_MENU(MainFrame::MENU_FILE_LOAD_SAMPLES, MainFrame::OnMenuFileLoadSamples)
//	EVT_MENU(MainFrame::MENU_FILE_SAVE_CAMERA_KEYS, MainFrame::OnSaveCameraKeys)
//	EVT_MENU(MainFrame::MENU_FILE_LOAD_CAMERA_KEYS, MainFrame::OnLoadCameraKeys)
	EVT_MENU(MainFrame::MENU_FILE_EXPORT_IMAGES, MainFrame::OnExportImageSequence)
	EVT_MENU(MainFrame::MENU_FILE_EXPORT_IMAGE, MainFrame::OnExportImage)
	EVT_MENU(MainFrame::MENU_FILE_EXPORT_TIMELAPSE, MainFrame::OnExportTimelapseImage)
	EVT_MENU(MainFrame::MENU_FILE_EXPORT_RIB, MainFrame::OnExportRib)
	EVT_MENU(MainFrame::MENU_FILE_EXPORT_MAYA_ASCII, MainFrame::OnExportMayaAscii)
	EVT_MENU(MainFrame::MENU_FILE_DUMP_GENERATION_IMAGES, MainFrame::OnMenuFileDumpGenerationImages)
#ifdef DRAW_LINES
	EVT_MENU(MainFrame::MENU_FILE_DUMP_PATHS_IMAGES, MainFrame::OnMenuFileDumpPathsImages)
#endif
	EVT_MENU(MainFrame::MENU_FILE_EXIT, MainFrame::OnExit)
	EVT_MENU(MainFrame::MENU_EDIT_UNDO, MainFrame::OnUndo)
	EVT_MENU(MainFrame::MENU_EDIT_REDO, MainFrame::OnRedo)
	EVT_MENU(MainFrame::MENU_EDIT_COPY, MainFrame::OnCopy)
	EVT_MENU(MainFrame::MENU_EDIT_PASTE, MainFrame::OnPaste)
	EVT_MENU(MainFrame::MENU_EDIT_DUPLICATE, MainFrame::OnDuplicate)
	EVT_MENU(MainFrame::MENU_EDIT_HIDE, MainFrame::OnHide)
	EVT_MENU(MainFrame::MENU_EDIT_SHOW_ALL, MainFrame::OnShowAll)
	EVT_MENU(MainFrame::MENU_EDIT_GROUP, MainFrame::OnGroup)
	EVT_MENU(MainFrame::MENU_EDIT_UNGROUP, MainFrame::OnUngroup)
	EVT_MENU(MainFrame::MENU_EDIT_CENTER_PIVOT, MainFrame::OnCenterPivot)
	EVT_MENU(MainFrame::MENU_EDIT_FIX_CENTER_OF_MASS, MainFrame::OnFixCenterOfMass)
	EVT_MENU(MainFrame::MENU_EDIT_FUSE, MainFrame::OnFuse)
	EVT_MENU(MainFrame::MENU_EDIT_UNFUSE, MainFrame::OnUnfuse)
	EVT_MENU(MainFrame::MENU_EDIT_FIT_BOXES, MainFrame::OnFitBoxes)
	EVT_MENU(MainFrame::MENU_EDIT_FIT_CONVEX_HULLS, MainFrame::OnFitConvexHulls)
	EVT_MENU(MainFrame::MENU_EDIT_CREATE_CONVEX_HULL, MainFrame::OnCreateConvexHull)
#ifdef WINDOWS_MEDIA
	EVT_MENU(MainFrame::MENU_EDIT_SETTINGS, MainFrame::OnEditSettings)
#endif // WINDOWS_MEDIA
	EVT_MENU(MainFrame::MENU_CAMERA_PROPERTIES, MainFrame::OnCameraProperties)
	EVT_MENU(MainFrame::MENU_RENDER_GEOMETRY_BOXES, MainFrame::OnSetRenderOption)
	EVT_MENU(MainFrame::MENU_RENDER_GEOMETRY_POLYS, MainFrame::OnSetRenderOption)
	EVT_MENU(MainFrame::MENU_RENDER_GEOMETRY_WIREFRAME, MainFrame::OnSetRenderOption)
	EVT_MENU(MainFrame::MENU_RENDER_GEOMETRY_POINTS, MainFrame::OnSetRenderOption)
	EVT_MENU(MainFrame::MENU_RENDER_WIREFRAME_ON_SHADED, MainFrame::OnSetRenderOption)
	EVT_MENU(MainFrame::MENU_RENDER_INERTIA_TENSOR, MainFrame::OnSetRenderOption)
	EVT_MENU(MainFrame::MENU_RENDER_CONVEX_HULLS, MainFrame::OnSetRenderOption)
	EVT_MENU(MainFrame::MENU_RENDER_CONTACT_POINTS, MainFrame::OnSetRenderOption)
	EVT_MENU(MainFrame::MENU_RENDER_SHOW_ALL_PATHS, MainFrame::OnSetRenderShowAllPaths)
	EVT_MENU(MainFrame::MENU_RENDER_SHADOWS, MainFrame::OnSetRenderOption)
	EVT_MENU(MainFrame::MENU_RENDER_GRID, MainFrame::OnSetRenderOption)
	EVT_MENU(MainFrame::MENU_RENDER_BACKGROUND_WHITE, MainFrame::OnNeedsRedraw)
	EVT_MENU(MainFrame::MENU_RENDER_BACKGROUND_BLACK, MainFrame::OnNeedsRedraw)
	EVT_MENU(MainFrame::MENU_RENDER_BACKGROUND_LIGHT_GRAY, MainFrame::OnNeedsRedraw)
	EVT_MENU(MainFrame::MENU_RENDER_BACKGROUND_DARK_GRAY, MainFrame::OnNeedsRedraw)
	EVT_MENU(MainFrame::MENU_CREATE_BOX, MainFrame::OnCreateBox)
	EVT_MENU(MainFrame::MENU_CREATE_SPHERE, MainFrame::OnCreateSphere)
	EVT_MENU(MainFrame::MENU_CREATE_SPHERES, MainFrame::OnCreateSpheres)
	EVT_MENU(MainFrame::MENU_CREATE_CAPPED_CYLINDER, MainFrame::OnCreateCappedCylinder)
	EVT_MENU(MainFrame::MENU_CREATE_CYLINDER, MainFrame::OnCreateCylinder)
	EVT_MENU(MainFrame::MENU_CREATE_CONE, MainFrame::OnCreateCone)
	EVT_MENU(MainFrame::MENU_CREATE_OBJ, MainFrame::OnCreateObj)
	EVT_MENU(MainFrame::MENU_CREATE_PLANE, MainFrame::OnCreatePlane)
	EVT_MENU(MainFrame::MENU_CREATE_ORTHO_CAMERA, MainFrame::OnCreateOrthoCam)
	EVT_MENU(MainFrame::MENU_CREATE_PERSPECTIVE_CAMERA, MainFrame::OnCreatePerspectiveCam)
	EVT_MENU(MainFrame::MENU_CREATE_BALL_JOINT, MainFrame::OnCreateBallJoint)
	EVT_MENU(MainFrame::MENU_CREATE_SLIDER_JOINT, MainFrame::OnCreateSliderJoint)
	EVT_MENU(MainFrame::MENU_CREATE_HINGE_JOINT, MainFrame::OnCreateHingeJoint)
	EVT_MENU(MainFrame::MENU_CREATE_UNIVERSAL_JOINT, MainFrame::OnCreateUniversalJoint)
	EVT_MENU(MainFrame::MENU_CREATE_HINGE2_JOINT, MainFrame::OnCreateHinge2Joint)
//	EVT_MENU(MainFrame::MENU_CREATE_ROTMOTOR1_JOINT, MainFrame::OnCreateRotMotor1Joint)
	EVT_MENU(MainFrame::MENU_CREATE_PLANE_JOINT, MainFrame::OnCreatePlaneJoint)
	EVT_MENU(MainFrame::MENU_SIMULATE_RUN, MainFrame::OnMenuSimulateRun)
	EVT_MENU(MainFrame::MENU_SIMULATE_SAMPLING, MainFrame::OnMenuSetSampling)
	EVT_MENU(MainFrame::MENU_SIMULATE_UPDATE_PATH, MainFrame::OnMenuUpdateSample)
#ifdef SOUND
	EVT_MENU(MainFrame::MENU_SOUND_COMPUTE_ALL_MODES, MainFrame::OnSoundComputeAllModes)
	EVT_MENU(MainFrame::MENU_SOUND_COMPUTE_OBJECT_MODES, MainFrame::OnSoundComputeObjectModes)
	EVT_MENU(MainFrame::MENU_SOUND_STOP_COMPUTING_MODES, MainFrame::OnSoundStopComputingModes)
#endif // SOUND
	EVT_MENU(MainFrame::MENU_HELP_ABOUT, MainFrame::OnHelpAbout)
	EVT_MENU(MainFrame::MENU_WINDOW_SIM_SORTER, MainFrame::OnWindowSimSorter)
	EVT_MENU(MainFrame::MENU_WINDOW_HIERARCHY_VIEW, MainFrame::OnWindowHierarchyView)
	EVT_MENU(MainFrame::MENU_WINDOW_PROPERTIES, MainFrame::OnWindowProperties)
	EVT_BUTTON(MainFrame::BUTTON_PLAY, MainFrame::OnPlay)
	EVT_COMMAND_SCROLL(MainFrame::myID_TIMESLIDER, MainFrame::OnScrollTime)
	EVT_KEY_DOWN(MainFrame::OnKey)
	EVT_KEY_UP(MainFrame::OnKeyUp)
	EVT_CLOSE(MainFrame::OnClose)
	EVT_CHOICE(CHOICE_MODE, MainFrame::OnSetMode)
END_EVENT_TABLE()

template <typename T>
void GenericRelease( T* object )
{
	if( object != NULL )
		object->Release();
}

#ifdef WINDOWS_MEDIA
bool checkVideoCodec( WM_MEDIA_TYPE* pType )
{
	return true;
}


bool checkAudioCodec( WM_MEDIA_TYPE* pType )
{
	WAVEFORMATEX* pWave = NULL;

	// Check that the format data is present.
	if(pType->cbFormat >= sizeof(WAVEFORMATEX))
	{
		pWave = (WAVEFORMATEX*) pType->pbFormat;
	}
	else
	{
		throw Exception( "Unexpected Error retrieving audio codecs" );
	}

	if( pWave->nChannels != 2 ) // require stereo
		return false;

	if( pWave->nSamplesPerSec != 44100 ) // require CD-quality
		return false;

	if( pWave->wBitsPerSample < 16 ) // require at least 16 bits per sample
		return false;

	// Since audio/video synchronization required, check the number
	//  of packets per second (Bps / BlockAlign). The bit rate is
	//  greater than 3200 bps, this value must be 5. 
	//  Otherwise this value is 3.
	// This is an ASF requirement.
	// see http://msdn2.microsoft.com/en-us/library/aa384570.aspx
	if((pWave->nAvgBytesPerSec / pWave->nBlockAlign) >= 
		((pWave->nAvgBytesPerSec >= 4000) ? 5.0 : 3.0))
	{
		// it's okay
	}
	else
	{
		return false;
	}

	return true;
}
#endif

struct SamplingStruct
{
	MainFrame* frame;
	std::string serverName;
	boost::shared_ptr<Event> samplingEvent;
};

DWORD WINAPI SamplingThreadProc( LPVOID lpParam )
{
	SamplingStruct* pData = reinterpret_cast<SamplingStruct*>( lpParam );

	assert( pData->frame );
#ifdef SOCKET_SAMPLING
	if( pData->serverName.empty() )
#endif
		pData->frame->sampleTree( pData->samplingEvent );
#ifdef SOCKET_SAMPLING
	else
		pData->frame->sampleTreeRemote( pData->serverName, pData->samplingEvent );
#endif


	delete pData;
	return 0;
}

MainFrame::MainFrame(wxFrame *frame, const wxString& title, const wxPoint& pos,
    const wxSize& size, long style)
	:	wxFrame(frame, -1, title, pos, size, style), 
		treeCounter_(1),
		playBitmap_( play_xpm ),
		pauseBitmap_( pause_xpm ),
		mode_( MODE_EDIT ),
		currentPath_( 0 ),
		historyComplete_( true ),
		numMachines_( 0 ),
		wireframeOnSelected_( true )
#ifdef SOUND
		,
		prevHighestAudio_( 1.0 ),
		computingModes_( false ),
		stopComputingModes_( true )
#endif
{
    AllocConsole();

    {
        int hCrt = _open_osfhandle(
                 (long) GetStdHandle(STD_OUTPUT_HANDLE),
                 _O_TEXT
              );
        FILE* hf = _fdopen( hCrt, "w" );
        *stdout = *hf;
        setvbuf( stdout, NULL, _IONBF, 0 );
    }

    {
        int hCrt = _open_osfhandle(
                 (long) GetStdHandle(STD_ERROR_HANDLE),
                 _O_TEXT
              );
        FILE* hf = _fdopen( hCrt, "w" );
        *stdout = *hf;
        setvbuf( stdout, NULL, _IONBF, 0 );
    }
    
    {
		//wxLog *logger=new wxLogStream(&std::cerr);
		//wxLog::SetActiveTarget(logger);

		logWindow_ = new wxLogWindow(this, "Error log");
		//wxLog::GetActiveTarget()->SetVerbose();
	}

	// default to set:
	historyComplete_.notify_all();

#ifdef USE_NOVODEX
	try
	{
		physicsSDK_ = instantiateSDK(0, new wxNxOutputStream);
	}
	catch( CreationException& )
	{
		wxLogMessage( "Unable to instantiate PhysX SDK; disabling support." );
	}
#endif

	aboutDialog_ = new AboutDialog( this );

#ifdef WIN32
	SYSTEM_INFO info;
	GetSystemInfo( &info );
	numberOfProcessors = /* info.dwNumberOfProcessors */ 1;
#endif

#ifdef WINDOWS_MEDIA
	// Windows Media junk
	{
		IWMProfileManager* profileManager;
		HRESULT res = WMCreateProfileManager( &profileManager );
		if( FAILED(res) )
			throw Exception( "Unable to create profile manager." );

		profileManager_.reset( profileManager, &GenericRelease<IWMProfileManager> );
	}

	boost::tie(audioCodecDescriptions_, audioStreamConfigurations_) = 
		findAvailableCodecs( WMMEDIATYPE_Audio, checkAudioCodec, false );
	boost::tie(videoCodecDescriptions_, videoStreamConfigurations_) = 
		findAvailableCodecs( WMMEDIATYPE_Video, checkVideoCodec, true );

	this->selectedVideoCodec_ = 0;
	this->selectedAudioCodec_ = 0;
#endif

	computeHistoryDone_ = false;
	computeHistoryThread_.reset( 
		new boost::thread(
			boost::bind(&MainFrame::computeHistory, this) ) );

#ifdef LOGGING
	{
		std::string logFilename;
		for( size_t i = 0; ; ++i )
		{
			std::ostringstream oss;
			//oss << "g:/users/cdtwigg/planning/logs/log" << std::setw(5) << std::setfill('0') << i << ".txt";
			oss << std::setw(5) << std::setfill('0') << i << ".txt";
			boost::filesystem::path p( oss.str() );
			if( !boost::filesystem::exists(p) )
			{
				logFilename = oss.str();
				break;
			}
		}

		this->logger_.reset( new Logger(logFilename) );
	}
#endif

	{
		this->sampling_ = true;
		this->pauseSampling_ = false;
		int numThreads = numberOfProcessors;

#ifdef SOCKET_SAMPLING
		std::vector<std::string> hostnames(1, "flip.graphics.cs.cmu.edu");
		{
			boost::filesystem::path path( "c:\\users\\cdtwigg\\planning\\servers" );
			//boost::filesystem::path path( "\\\\volley.graphics.cs.cmu.edu\\cdtwigg\\planning\\servers" );
			boost::filesystem::directory_iterator dirItr( path );
			boost::filesystem::directory_iterator endDir;
			for( ; dirItr != endDir; ++dirItr )
			{
				this->samplingEvents_.push_back( boost::shared_ptr<Event>(new Event(false)) );
				hostnames.push_back( dirItr->leaf() );
			}
		}

		numThreads += hostnames.size();
#endif

	scene_.reset( new Scene );

	menuBar_ = new wxMenuBar();

	{
		fileMenu_ = new wxMenu();
		fileMenu_->Append(MENU_FILE_NEW, wxT("New"));
		fileMenu_->Append(MENU_FILE_OPEN, wxT("Open..."));
		fileMenu_->Append(MENU_FILE_IMPORT, wxT("Import..."));
		fileMenu_->Append(MENU_FILE_SAVE, wxT("Save"));
		fileMenu_->Append(MENU_FILE_SAVEAS, wxT("Save As..."));
		fileMenu_->AppendSeparator();

		fileMenu_->Append(MENU_FILE_LOAD_TREE, wxT("Load tree...") );
		fileMenu_->Append(MENU_FILE_SAVE_TREE, wxT("Save tree...") );
		fileMenu_->Append(MENU_FILE_LOAD_SAMPLES, wxT("Load samples...") );
		fileMenu_->Append(MENU_FILE_LOAD_DEBUGGING_INFO, wxT("Load debugging info...") );
		fileMenu_->AppendSeparator();

		/*
		fileMenu_->Append(MENU_FILE_LOAD_CAMERA_KEYS, wxT("Load camera animation..."));
		fileMenu_->Append(MENU_FILE_SAVE_CAMERA_KEYS, wxT("Save camera animation..."));
		fileMenu_->AppendSeparator();
		*/

		fileMenu_->Append(MENU_FILE_EXPORT_RIB, wxT("Write .rib file..."));
		fileMenu_->Append(MENU_FILE_EXPORT_MAYA_ASCII, wxT("Write .ma file..."));
		fileMenu_->Append(MENU_FILE_EXPORT_IMAGE, wxT("Write image..."));
		fileMenu_->Append(MENU_FILE_EXPORT_IMAGES, wxT("Write image sequence..."));
		fileMenu_->Append(MENU_FILE_EXPORT_TIMELAPSE, wxT("Write timelapse image..."));
		fileMenu_->Append(MENU_FILE_DUMP_GENERATION_IMAGES, wxT("Dump dataset generation images..."));
#ifdef DRAW_LINES
		fileMenu_->Append(MENU_FILE_DUMP_PATHS_IMAGES, wxT("Dump paths images..."));
#endif
		fileMenu_->AppendSeparator();

		fileMenu_->Append(MENU_FILE_EXIT, wxT("E&xit"));
	}

	{
		editMenu_ = new wxMenu();
		editMenu_->Append(MENU_EDIT_UNDO, wxT("Undo\tCTRL+z"));
		editMenu_->Append(MENU_EDIT_REDO, wxT("Redo\tCTRL+y"));
		editMenu_->AppendSeparator();

		editMenu_->Append(MENU_EDIT_COPY, wxT("Copy\tCTRL+c"));
		editMenu_->Append(MENU_EDIT_PASTE, wxT("Paste\tCTRL+v"));
		editMenu_->Append(MENU_EDIT_DUPLICATE, wxT("Duplicate\tCTRL+d"));
		editMenu_->AppendSeparator();

		editMenu_->Append(MENU_EDIT_HIDE, wxT("Hide"));
		editMenu_->Append(MENU_EDIT_SHOW_ALL, wxT("Show all"));
		editMenu_->AppendSeparator();

		editMenu_->Append(MENU_EDIT_GROUP, wxT("Group\tCTRL+g"));
		editMenu_->Append(MENU_EDIT_UNGROUP, wxT("Ungroup\tCTRL+u"));
		editMenu_->Append(MENU_EDIT_CENTER_PIVOT, wxT("Center pivot"));
		editMenu_->Append(MENU_EDIT_FIX_CENTER_OF_MASS, wxT("Fix center of mass"));
		editMenu_->AppendSeparator();

		editMenu_->Append(MENU_EDIT_FUSE, wxT("Fuse"));
		editMenu_->Append(MENU_EDIT_UNFUSE, wxT("Unfuse"));
		editMenu_->AppendSeparator();

		editMenu_->Append(MENU_EDIT_FIT_CONVEX_HULLS, wxT("Fit convex hulls"));
		editMenu_->Append(MENU_EDIT_CREATE_CONVEX_HULL, wxT("Create convex hull from points"));
		editMenu_->AppendSeparator();

		{
			editSelectByMenu_ = new wxMenu();
			editSelectByMenu_->AppendRadioItem(MENU_EDIT_SELECT_BY_OBJECT, wxT("Object"));
			editSelectByMenu_->AppendRadioItem(MENU_EDIT_SELECT_BY_HIERARCHY, wxT("Hierarchy"));
			editSelectByMenu_->AppendRadioItem(MENU_EDIT_SELECT_BY_SIMULATED, wxT("Simulated"));
			editSelectByMenu_->Check(MENU_EDIT_SELECT_BY_OBJECT, true);
			editMenu_->Append(MENU_EDIT_SELECT_BY, wxT("Select by"), editSelectByMenu_);
		}

#ifdef WINDOWS_MEDIA
		editMenu_->Append(MENU_EDIT_SETTINGS, wxT("Settings..."));
#endif // WINDOWS_MEDIA
	}

	{
		createMenu_ = new wxMenu();
		createMenu_->Append(MENU_CREATE_BOX, wxT("Box"));
		createMenu_->Append(MENU_CREATE_PLANE, wxT("Plane"));
		createMenu_->Append(MENU_CREATE_SPHERE, wxT("Sphere"));
		createMenu_->Append(MENU_CREATE_CAPPED_CYLINDER, wxT("Capsule"));
		createMenu_->Append(MENU_CREATE_CYLINDER, wxT("Cylinder"));
		createMenu_->Append(MENU_CREATE_CONE, wxT("Cone"));
		createMenu_->Append(MENU_CREATE_OBJ, wxT(".obj geometry"));
		createMenu_->Append(MENU_CREATE_SPHERES, wxT(".sph file"));
		createMenu_->AppendSeparator();

		createMenu_->Append(MENU_CREATE_BALL_JOINT, wxT("Ball joint"));
		createMenu_->Append(MENU_CREATE_HINGE_JOINT, wxT("Hinge joint"));
		createMenu_->Append(MENU_CREATE_SLIDER_JOINT, wxT("Slider joint"));
		createMenu_->Append(MENU_CREATE_UNIVERSAL_JOINT, wxT("Universal joint"));
		createMenu_->Append(MENU_CREATE_HINGE2_JOINT, wxT("Hinge-2 (e.g. suspension)"));
		createMenu_->Append(MENU_CREATE_PLANE_JOINT, wxT("Plane joint"));
//		createMenu_->Append(MENU_CREATE_ROTMOTOR1_JOINT, wxT("Rotational motor (1-axis)"));
		createMenu_->AppendSeparator();

		createMenu_->Append(MENU_CREATE_PERSPECTIVE_CAMERA, wxT("Perspective camera"));
		createMenu_->Append(MENU_CREATE_ORTHO_CAMERA, wxT("Orthographic camera"));
	}
    
	{
		cameraMenu_ = new wxMenu();
		{
			cameraMenu_->Append( MENU_CAMERA_SET_KEY, wxT("Set key") );
			cameraMenu_->Append( MENU_CAMERA_ERASE_KEY, wxT("Delete key") );
			cameraMenu_->Append( MENU_CAMERA_NEXT_KEY, wxT("Next key") );
			cameraMenu_->Append( MENU_CAMERA_PREV_KEY, wxT("Previous key") );
			cameraMenu_->Append( MENU_CAMERA_PROPERTIES, wxT("Properties...") );
			cameraMenu_->AppendSeparator();
		}
	}

	{
		renderMenu_ = new wxMenu();

		{
			renderGeometryMenu_ = new wxMenu();
			renderGeometryMenu_->AppendRadioItem(MENU_RENDER_GEOMETRY_BOXES, wxT("Boxes"));
			renderGeometryMenu_->AppendRadioItem(MENU_RENDER_GEOMETRY_POLYS, wxT("Polys"));
			renderGeometryMenu_->AppendRadioItem(MENU_RENDER_GEOMETRY_WIREFRAME, wxT("Wireframe"));
			renderGeometryMenu_->AppendRadioItem(MENU_RENDER_GEOMETRY_POINTS, wxT("Points"));
			renderGeometryMenu_->Check(MENU_RENDER_GEOMETRY_POLYS, TRUE);
		}

		renderMenu_->Append(MENU_RENDER_GEOMETRY, wxT("Geometry"), renderGeometryMenu_);
		renderMenu_->AppendCheckItem(MENU_RENDER_WIREFRAME_ON_SHADED, wxT("Wireframe on shaded"));
		renderMenu_->Check(MENU_RENDER_WIREFRAME_ON_SHADED, TRUE);

		renderMenu_->AppendCheckItem(MENU_RENDER_INERTIA_TENSOR, wxT("Inertia tensor"));
		renderMenu_->AppendCheckItem(MENU_RENDER_CONVEX_HULLS, wxT("Show convex hulls for meshes"));
		renderMenu_->AppendCheckItem(MENU_RENDER_CONTACT_POINTS, wxT("Show contact points"));
		renderMenu_->AppendCheckItem(MENU_RENDER_SHOW_ALL_PATHS, wxT("Show all paths"));

		renderMenu_->AppendCheckItem(MENU_RENDER_GRID, wxT("Grid"));
		renderMenu_->Check(MENU_RENDER_GRID, TRUE);

		renderMenu_->AppendCheckItem(MENU_RENDER_SHADOWS, wxT("Shadows"));
		renderMenu_->Check(MENU_RENDER_SHADOWS, TRUE);

		{
			renderBackgroundMenu_ = new wxMenu();
			renderBackgroundMenu_->AppendRadioItem(MENU_RENDER_BACKGROUND_WHITE, wxT("White"));
			renderBackgroundMenu_->AppendRadioItem(MENU_RENDER_BACKGROUND_LIGHT_GRAY, wxT("Light Gray"));
			renderBackgroundMenu_->AppendRadioItem(MENU_RENDER_BACKGROUND_DARK_GRAY, wxT("Dark Gray"));
			renderBackgroundMenu_->AppendRadioItem(MENU_RENDER_BACKGROUND_BLACK, wxT("Black"));
			renderBackgroundMenu_->Check(MENU_RENDER_BACKGROUND_LIGHT_GRAY, TRUE);
		}
		renderMenu_->Append(MENU_RENDER_BACKGROUND, wxT("Background"), renderBackgroundMenu_);
	}

	{
		simulateMenu_ = new wxMenu();
		simulateMenu_->AppendCheckItem( MENU_SIMULATE_RUN, wxT("Run simulation"));
		simulateMenu_->AppendCheckItem( MENU_SIMULATE_SAMPLING, wxT("Sampling"));
		simulateMenu_->Check( MENU_SIMULATE_SAMPLING, !pauseSampling_ );
		{
			simulatorMenu_ = new wxMenu();
			simulatorMenu_->AppendRadioItem( MENU_SIMULATE_SIMULATOR_ODE_NATIVE, wxT("Open Dynamics Engine (native collision detection)") );
			simulatorMenu_->AppendRadioItem( MENU_SIMULATE_SIMULATOR_ODE_PENALTY, wxT("Simple penalty model") );
			simulatorMenu_->Check(MENU_SIMULATE_SIMULATOR_ODE_NATIVE, TRUE);
#ifdef USE_BULLET
			simulatorMenu_->AppendRadioItem( MENU_SIMULATE_SIMULATOR_ODE_BULLET, wxT("Open Dynamics Engine (Bullet collision detection)") );
			simulatorMenu_->AppendRadioItem( MENU_SIMULATE_SIMULATOR_BULLET, wxT("Bullet Physics Library") );
//			simulatorMenu_->Check(MENU_SIMULATE_SIMULATOR_ODE_BULLET, TRUE);
#endif
#ifdef USE_NEWTON
			simulatorMenu_->AppendRadioItem( MENU_SIMULATE_SIMULATOR_NEWTON, wxT("Newton Game Dynamics") );
			simulatorMenu_->Check(MENU_SIMULATE_SIMULATOR_NEWTON, TRUE);
#endif
#ifdef USE_NOVODEX
			simulatorMenu_->AppendRadioItem( MENU_SIMULATE_SIMULATOR_NOVODEX, wxT("Ageia PhysX API") );
			if( physicsSDK_ )
				simulatorMenu_->Check(MENU_SIMULATE_SIMULATOR_NOVODEX, TRUE);
			else
			{
				simulatorMenu_->Enable(MENU_SIMULATE_SIMULATOR_NOVODEX, FALSE);
				simulatorMenu_->Check(MENU_SIMULATE_SIMULATOR_NOVODEX, FALSE);
			}
#endif
		}
		simulateMenu_->Append(MENU_SIMULATE_SIMULATOR, wxT("Simulator"), simulatorMenu_);

		simulateMenu_->Append( MENU_SIMULATE_UPDATE_PATH, wxT("Update sample"));
		simulateMenu_->Check(MENU_SIMULATE_RUN, FALSE);
	}

#ifdef SOUND
	{
		soundMenu_ = new wxMenu();
		soundMenu_->AppendCheckItem( MENU_SOUND_ENABLE_SOUND, wxT("Enable sound"));
		soundMenu_->Check(MENU_SOUND_ENABLE_SOUND, TRUE);
		soundMenu_->Append( MENU_SOUND_COMPUTE_ALL_MODES, wxT("Compute modes for all objects") );
		soundMenu_->Append( MENU_SOUND_COMPUTE_OBJECT_MODES, wxT("Compute modes for selected object") );
		soundMenu_->Append( MENU_SOUND_STOP_COMPUTING_MODES, wxT("Stop computing modes") );
	}
#endif // SOUND

	{
		windowMenu_ = new wxMenu();
		windowMenu_->Append( MENU_WINDOW_SIM_SORTER, wxT("Simulation sorter") );
		windowMenu_->Append( MENU_WINDOW_HIERARCHY_VIEW, wxT("View hierarchy") );
		windowMenu_->AppendCheckItem( MENU_WINDOW_PROPERTIES, wxT("Attributes") );
	}

	{
		helpMenu_ = new wxMenu();
		helpMenu_->Append( MENU_HELP_ABOUT, wxT("About Many-Worlds Browsing Demo...") );
	}
	
	menuBar_->Append(fileMenu_, wxT("&File"));
	menuBar_->Append(editMenu_, wxT("&Edit"));
	menuBar_->Append(createMenu_, wxT("&Create"));
	menuBar_->Append(cameraMenu_, wxT("Camera"));
	menuBar_->Append(renderMenu_, wxT("Render"));
	menuBar_->Append(simulateMenu_, wxT("Simulate"));
#ifdef SOUND
	menuBar_->Append(soundMenu_, wxT("Sound"));
#endif // SOUND
	menuBar_->Append(windowMenu_, wxT("Window"));
	menuBar_->Append(helpMenu_, wxT("Help"));

	SetMenuBar(menuBar_);


	 // initialize statusbar
	static const int widths[] = {-1};
	CreateStatusBar (WXSIZEOF(widths), wxST_SIZEGRIP, myID_STATUSBAR);
	SetStatusWidths (WXSIZEOF(widths), widths);

	{
		// create toolbar
		toolBar_ = wxFrame::CreateToolBar (wxTB_TEXT|wxTB_FLAT|wxTB_NO_TOOLTIPS);

	/*
		wxBitmap playImage(play_xpm);
		toolBar_->AddCheckTool(myID_PLAY, _T("Play"), playImage, wxNullBitmap );

		wxBitmap simImage(sim_xpm);
		toolBar_->AddCheckTool(myID_SIM, _T("Simulate"), simImage, wxNullBitmap );
	*/
		wxString choices[] = { wxString("Edit"), wxString("Browse") };
		modeChoice_ = new wxChoice(toolBar_, CHOICE_MODE, wxDefaultPosition, wxDefaultSize, 2, choices );
		toolBar_->AddControl( modeChoice_ );

		toolBar_->AddSeparator();

		wxBitmap redRectImage(redRect_xpm);
		wxBitmap greenRectImage(greenRect_xpm);
		toolBar_->AddRadioTool(myID_SELECT_RECT, wxT(""), selectRect_xpm, wxNullBitmap,
			wxT("Select path"),
			wxT("Select a path") );
		toolBar_->AddRadioTool(myID_GOOD_CONSTRAINT_RECT, wxT(""), greenRect_xpm, wxNullBitmap,
			wxT("Positive constraint"),
			wxT("Add a positive constraint, filtering out paths that do not pass through the region") );
		toolBar_->AddRadioTool(myID_BAD_CONSTRAINT_RECT, wxT(""), redRect_xpm, wxNullBitmap,
			wxT("Negative constraint"),
			wxT("Add a negative constraint, filtering out paths that pass through the region") );
		toolBar_->AddRadioTool(myID_REFINE_RECT, wxT(""), refine_xpm, wxNullBitmap,
			wxT("Refine simulation"), 
			wxT("Refine the simulation of the selected objects in the selected region") );
		toolBar_->AddRadioTool(myID_SKETCH_ACTION, wxT(""), arc_xpm, wxNullBitmap,
			wxT("Sketch action"), 
			wxT("Sketch next action") );

		toolBar_->AddSeparator();

		toolBar_->AddRadioTool(myID_CURSOR_TOOLBAR, wxT(""), cursor_xpm, wxNullBitmap,
			wxT("Select"), 
			wxT("Select objects"));
		toolBar_->AddRadioTool(myID_TRANSLATE_TOOLBAR, wxT(""), translate_xpm, wxNullBitmap,
			wxT("Translate"), 
			wxT("Translate objects"));
		toolBar_->AddRadioTool(myID_ROTATE_TOOLBAR, wxT(""), rotate_xpm, wxNullBitmap,
			wxT("Rotate"), 
			wxT("Rotate objects"));
		toolBar_->AddRadioTool(myID_SCALE_TOOLBAR, wxT(""), scale_xpm, wxNullBitmap,
			wxT("Scale"), 
			wxT("Scale objects"));
		toolBar_->AddRadioTool(myID_CHANGEPIVOT_TOOLBAR, wxT(""), changePivot_xpm, wxNullBitmap,
			wxT("Change pivot"), 
			wxT("Change object pivot"));
		toolBar_->AddRadioTool(myID_SELECTPOINTS_TOOLBAR, wxT(""), points_xpm, wxNullBitmap,
			wxT("Select points"), 
			wxT("Select mesh points"));

		toolBar_->Realize();
	}

	timeView_ = new TimeView(this, -1, wxPoint(100, 100), wxSize(30, 30), wxSUNKEN_BORDER);
	propertyDialog_ = new PropertyDialog( this );
#ifdef WINDOWS_MEDIA
	settingsDialog_ = new SettingsDialog( this );
#endif

	wxBoxSizer *topsizer = new wxBoxSizer( wxVERTICAL );

	viewsSizer_ = new wxFlexGridSizer( 2, 2 );
	canvases_.push_back( 
		new GLView(this, -1, wxDefaultPosition, wxDefaultSize, wxSUNKEN_BORDER) );
	this->currentCanvas_ = canvases_.front();

	for( size_t i = 1; i < 4; ++i )
	{
		canvases_.push_back( 
			new GLView(this, currentCanvas_->GetContext(), -1, 
				wxDefaultPosition, wxDefaultSize, wxSUNKEN_BORDER) );
	}

	for( size_t i = 0; i < 2; ++i )
	{
		viewsSizer_->AddGrowableRow(i, 1);
		viewsSizer_->AddGrowableCol(i, 1);
	}

	for( size_t i = 0; i < canvases_.size(); ++i )
	{
		viewsSizer_->Add( canvases_.at(i), 1, wxEXPAND );
	}


	this->propSizer_ = new wxBoxSizer(wxHORIZONTAL);
	propSizer_->Add( viewsSizer_, 1, wxEXPAND );
	propSizer_->Add( propertyDialog_, 0, wxEXPAND );

	/*
	@todo: need to be able to switch between the 4-panel view and the normal one

	grid->Add( bestView_, 1, wxEXPAND );
	grid->Add( worstView_, 1, wxEXPAND );
	grid->Add( baseView_, 1, wxEXPAND );
	*/

	topsizer->Add( propSizer_, 1, wxEXPAND );

	wxBoxSizer *timebarSizer = new wxBoxSizer( wxHORIZONTAL );
	timebarSizer->Add( timeView_, 1, 0, 0 );
	{
		playButton_ = new wxBitmapButton(this, BUTTON_PLAY, playBitmap_);
		playing_ = false;
		playButton_->Connect( wxEVT_KEY_DOWN, wxKeyEventHandler(MainFrame::OnKey), NULL, this );
		playButton_->Connect( wxEVT_KEY_UP, wxKeyEventHandler(MainFrame::OnKeyUp), NULL, this );
		timebarSizer->Add( playButton_, 0, wxEXPAND );
	}

	topsizer->Add( timebarSizer, 0, wxALL | wxEXPAND );
	SetSizer( topsizer );
//	topsizer->SetSizeHints( this );


	cameraPropertyDialog_ = new CameraPropertyDialog( this, "Camera Properties", wxDefaultPosition );
//	graphFrame_ = new GraphFrame(this, "Graph view", wxDefaultPosition, wxSize(300, 300));
//	plotFrame_ = new PlotFrame(this, "Plot view", wxDefaultPosition, wxSize(500, 500));
	bestScoresFrame_ = new BestScoresFrame(this, canvases_.front()->GetContext(), "Sorter", wxDefaultPosition, wxSize(720, 480));
	hierarchyFrame_ = new HierarchyFrame(this, "Hierarchy View", wxDefaultPosition, wxSize(300, 480));

	wxSize mySize = this->GetSize();
	wxPoint myPosition = this->GetPosition();

	setMode( MODE_EDIT, true );
	updateMenus();

	changeScene( scene_ );

		for( size_t iSamplingThread = 0; iSamplingThread < numThreads; ++iSamplingThread )
		{
			this->samplingEvents_.push_back( boost::shared_ptr<Event>(new Event(false)) );

			SamplingStruct* pData = new SamplingStruct;
			pData->frame = this;
			pData->samplingEvent = this->samplingEvents_.back();

#ifdef SOCKET_SAMPLING
			if( iSamplingThread < hostnames.size() )
				pData->serverName = hostnames[iSamplingThread];
#endif

			DWORD id;
			HANDLE h = CreateThread( 
				NULL,              // default security attributes
				0,                 // use default stack size  
				SamplingThreadProc,        // thread function 
				pData,             // argument to thread function 
				0,                 // use default creation flags 
				&id);	           // returns the thread identifier 

			if( h == NULL )
				continue;

			samplingThreadHandles_.push_back( h );
			samplingThreadIds_.push_back( id );
			SetThreadPriority( h, THREAD_PRIORITY_BELOW_NORMAL );
		}
	}

	windowMenu_->Check( MENU_WINDOW_PROPERTIES, this->arePropertiesShown() );

}

void MainFrame::log( const std::string& message )
{
#ifdef LOGGING
	logger_->log( message );
#endif
}

MainFrame::~MainFrame()
{
    FreeConsole();
}

void MainFrame::OnHelpAbout( wxCommandEvent& event )
{
	aboutDialog_->ShowModal();
}

void MainFrame::updateActiveTimeRange( std::pair<double, double> timeRange )
{
	assert( timeRange.first <= timeRange.second );

	SimulationTreePtr tree = this->currentTree_;
	if( tree )
		tree->setTimeRange( timeRange );
}

bool MainFrame::viewMaximized() const
{
	bool allVisible = true;

	for( size_t i = 0; i < canvases_.size(); ++i )
	{
		allVisible = allVisible &&
			canvases_.at(i)->IsShown();
	}

	return !allVisible;
}

void MainFrame::setCurrentView( GLView* view, bool toggleMaximized )
{
	if( this->currentCanvas_ != view )
	{
		// mark both for redraw
		this->currentCanvas_->Refresh( FALSE );
		view->Refresh( FALSE );

		this->currentCanvas_ = view;
		this->currentCanvas_->SetFocus();
		this->updateCameraMenu();
	}

	if( toggleMaximized )
	{
		if( !this->viewMaximized() )
		{
			for( size_t i = 0; i < canvases_.size(); ++i )
			{
				this->viewsSizer_->Show( canvases_.at(i), 
					canvases_.at(i) == view );
			}
		}
		else
		{
			for( size_t i = 0; i < canvases_.size(); ++i )
			{
				this->viewsSizer_->Show( canvases_.at(i), true );
			}
		}

		this->viewsSizer_->Layout();
	}
}

void MainFrame::OnLogError( wxCommandEvent& event )
{
	wxLogError( event.GetString() );
}

GLView* MainFrame::GetCanvas()
{
	assert( currentCanvas_ != 0 );
	return this->currentCanvas_;
}

const GLView* MainFrame::GetCanvas() const
{
	assert( currentCanvas_ != 0 );
	return this->currentCanvas_;
}


void MainFrame::OnIdle(wxIdleEvent& event)
{
	if( this->idleAction_ )
	{
		if( idleAction_->complete() )
		{
			idleAction_.reset();
		}
		else
		{
			idleAction_->perform();
			event.RequestMore();
			return;
		}
	}

	SimulationTreePtr tree = this->currentTree_;

	if( tree )
	{
		if( tree->checkForNewPaths() )
			this->updateViews();

		wxStatusBar* statusBar = this->GetStatusBar();

		static int lastNumMachines = 0;
		int curNumMachines = this->numMachines_.value();

		if( curNumMachines  != lastNumMachines )
		{
			std::ostringstream oss;
			oss << "Sampling on " << numMachines_.value() << " machines.";
			statusBar->SetStatusText( wxT(oss.str().c_str()) );
			lastNumMachines = curNumMachines;
		}

		static size_t counter = 0;
		if( ((counter++) % 1000) == 0 )
		{
			std::ostringstream oss;
			oss << "Computed " << tree->numPaths() << " paths on " << this->numMachines_.value() << " machines.";
			this->log( oss.str() );
		}
	}

	if( this->playing() )
	{
		this->setTime( getTime() + 1.0/30.0, true );
	}
	else if( simulation_ )
	{
        size_t nSteps = 8;
        for( size_t i = 0; i < nSteps; ++i )
    		simulation_->step( 1.0/(30.0 * static_cast<double>(nSteps)) );
		this->updateViews();
	}
	else if( updated_ )
	{
		this->updateViews();
	}
	else
	{
		event.Skip();
		return;
	}

	event.RequestMore();
}

void MainFrame::test()
{
#ifdef _DEBUG
	testOBB();
	testQuat();
	testSphereTree( scene_->camera()->camera()->projectiveTransform() 
		* scene_->camera()->camera()->modelViewTransform() );

	this->setMode( MainFrame::MODE_SOLVE );
	SimulationTreePtr tree( new SimulationTree("test", this->scene(), 30, simulationFactory()) );
#else
	errMsg( "Test only works in debug build." );
#endif
}

void MainFrame::sampleTree(boost::shared_ptr<Event> samplingEvent)
{
	RandomStream random;
	boost::shared_ptr<SimulationTree> tree;

	const DWORD baseRetryTimeout = 200;
	DWORD retryTimeout = baseRetryTimeout;
	while( this->sampling_ )
	{
		/*
#ifdef WIN32
		PERFORMANCE_INFORMATION info;
		if( GetPerformanceInfo(&info, sizeof(PERFORMANCE_INFORMATION)) )
		{
			SIZE_T pageSize = info.PageSize;
			SIZE_T physicalAvailable = info.PhysicalAvailable;
			SIZE_T totalAvailable = pageSize*physicalAvailable;

			const SIZE_T fiftyMB = 50*1024*1024;
			if( totalAvailable < fiftyMB ) // less than 100MB available
				break;
		}
#endif
		*/

		tree = this->currentTree_;

		if( !tree || this->pauseSampling_ || this->getMode() != MODE_SOLVE )
		{
			// wait for event so we're not spinlocking
			WaitForSingleObject( samplingEvent->handle(), INFINITE );
			continue;
		}

		// draw a sample
		//if( random.uniform() < 0.2 )
			tree->sample(random, false);
		//else
		//	tree->expand(random);

		//tree->expand(random);
	}

}

#ifdef SOCKET_SAMPLING
std::string GetOsErrorMessage(DWORD error)
{
    std::string rv;
    LPVOID lpMsgBuf;
    if (FormatMessage(
        FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        error,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
        (LPTSTR) &lpMsgBuf,
        0,
        NULL ))
    {
        rv.assign(reinterpret_cast<const char*>(lpMsgBuf));
    }
    else
    {
        rv.assign("FormatMessage API failed");
    }

	// @todo strip line feed here

	LocalFree(lpMsgBuf);
    return rv;
}

class SocketException
	: public std::exception
{
public:
	SocketException( const char* functionCall, int errorCode )
		: functionCall_(functionCall), errorCode_(errorCode)
	{
	}

	const char* function() const { return functionCall_; }
	int errorCode() const { return errorCode_; }

private:
	const char* functionCall_;
	int errorCode_;
};

class SocketWrapper
{
public:
	SocketWrapper( SOCKET socket )
		: socket_(socket), overlappedEvent_(true)
	{
		assert( socket_ != INVALID_SOCKET );
		memset( &overlapped_, 0, sizeof(OVERLAPPED) );
		overlapped_.hEvent = overlappedEvent_.handle();
	}

	~SocketWrapper()
	{
		// @todo implement proper graceful shutdown
		// see http://msdn.microsoft.com/library/default.asp?url=/library/en-us/winsock/winsock/shutdown_2.asp
		CancelIo( (HANDLE) this->socket_ );
		shutdown( this->socket_, SD_BOTH );
		closesocket( socket_ );
	}

	SOCKET socket() const
	{
		return this->socket_;
	}

	OVERLAPPED& overlapped()
	{
		return this->overlapped_;
	}

	void write( float value )
	{
		const unsigned int* uintVal = reinterpret_cast<const unsigned int*>( &value );
		this->write( *uintVal );
	}

	// waits for operation to complete:
	void write( unsigned int value )
	{
		BOOST_STATIC_ASSERT( sizeof(unsigned int) == 4 );
		unsigned int networkOrder = htonl( value );
		WSABUF buf;
		DWORD numberOfBytesSent;
		buf.len = sizeof( unsigned int );
		buf.buf = reinterpret_cast<char*>(&networkOrder);
		int result = WSASend(
			this->socket_,
			&buf,
			1,
			&numberOfBytesSent,
			0,
			NULL,
			NULL );

		if( result == SOCKET_ERROR )
			throw SocketException( "WSASend", WSAGetLastError() );

		if( numberOfBytesSent != buf.len )
			throw SocketException( "WSASend", 0 );
	}

	void write( std::string value )
	{
		assert( !value.empty() );

		WSABUF buf;
		buf.len = value.size();
		buf.buf = &value[0];
		DWORD numberOfBytesSent;

		int result = WSASend(
			this->socket_,
			&buf,
			1,
			&numberOfBytesSent,
			0,
			NULL,
			NULL );

		if( result == SOCKET_ERROR )
			throw SocketException( "WSASend", WSAGetLastError() );

		if( numberOfBytesSent != buf.len )
			throw SocketException( "WSASend", 0 );
	}

	template <typename T>
	void write( const std::vector<T>& value )
	{
		BOOST_STATIC_ASSERT( sizeof(T) == 4 );
		unsigned int size = value.size();
		this->write( size );

		if( 0 == size )
			return;

		std::vector<unsigned int> networkOrder( value.size() );
		const unsigned int* uintPtr = 
			reinterpret_cast<const unsigned int*>( &value[0] );
		std::transform( uintPtr, uintPtr + value.size(), 
			networkOrder.begin(), &htonl );

		WSABUF buf;
		buf.len = sizeof(unsigned int)*networkOrder.size();
		buf.buf = reinterpret_cast<char*>(&networkOrder[0]);
		DWORD numberOfBytesSent;		
		int result = WSASend(
			this->socket_,
			&buf,
			1,
			&numberOfBytesSent,
			0,
			NULL,
			NULL );

		if( result == SOCKET_ERROR )
			throw SocketException( "WSASend", WSAGetLastError() );

		if( numberOfBytesSent != buf.len )
			throw SocketException( "WSASend", 0 );
	}

private:
	SOCKET socket_;
	OVERLAPPED overlapped_;
	Event overlappedEvent_;
};

//// LookupAddress /////////////////////////////////////////////////////
// Given an address string, determine if it's a dotted-quad IP address
// or a domain address.  If the latter, ask DNS to resolve it.  In
// either case, return resolved IP address.  If we fail, we return
// INADDR_NONE.
// from http://tangentsoft.net/wskfaq/examples/basics/basic-client.cpp
u_long LookupAddress(const char* pcHost)
{
	u_long nRemoteAddr = inet_addr(pcHost);
	if (nRemoteAddr == INADDR_NONE) {
		// pcHost isn't a dotted IP, so resolve it through DNS
		hostent* pHE = gethostbyname(pcHost);
		if (pHE == 0) {
			return INADDR_NONE;
		}
		nRemoteAddr = *((u_long*)pHE->h_addr_list[0]);
	}

	return nRemoteAddr;
}

// This reads numBytes bytes from the overlapped socket
// while watching to see if eventToWatch gets raised;
// in that case, we toss an exception (tbd)
void readOverlapped( 
	char* buffer,
	size_t numBytes,
	SocketWrapper& socket,
	const Event& eventToWatch )
{
	size_t current = 0;
	while( current < numBytes )
	{
		ResetEvent(socket.overlapped().hEvent);

		WSABUF buf;
		buf.len = numBytes - current;
		buf.buf = buffer + current;
		
		DWORD bytesReceived;
		DWORD flags = 0;

		{
			int rc = WSARecv(
				socket.socket(),
				&buf,
				1,
				&bytesReceived,
				&flags,
				&(socket.overlapped()),
				NULL );
			if ( (rc == SOCKET_ERROR) && (WSA_IO_PENDING != WSAGetLastError()))
				throw SocketException( "WSARecv", WSAGetLastError() );
		}

		{
			WSAEVENT events[2];
			events[0] = eventToWatch.handle();
			events[1] = socket.overlapped().hEvent;
			DWORD rc = WaitForMultipleObjectsEx(2, events, FALSE, INFINITE, FALSE);
			if( rc == WAIT_FAILED )
				throw SocketException( "WaitForMultipleObjectsEx", GetLastError() );

			// @todo: do something more intelligent here
			if( rc == WAIT_OBJECT_0 )
				throw SocketException( "other event raised", 0 );
		}

		{
			int rc = WSAGetOverlappedResult(socket.socket(), &(socket.overlapped()), &bytesReceived, TRUE, &flags);
			if( rc == FALSE )
				throw SocketException( "WSAGetOverlappedResult", WSAGetLastError() );
		}

		if( bytesReceived == 0 )
			throw SocketException( "Connection closed", 0 );

		current += bytesReceived;
	}
}

std::vector<char> readOverlapped( 
	size_t numBytes,
	SocketWrapper& socket,
	const Event& eventToWatch )
{
	std::vector<char> result;
	result.resize( numBytes );
	if( numBytes > 0 )
		readOverlapped( &result[0], numBytes, socket, eventToWatch );
	return result;
}


unsigned int readOverlapped_uint32_noChecksum( SocketWrapper& socket,
	const Event& eventToWatch )
{
	unsigned int value;
	readOverlapped( reinterpret_cast<char*>( &value ),
		sizeof(unsigned int), socket, eventToWatch );
	return ntohl( value );
}

unsigned int readOverlapped_uint32( SocketWrapper& socket,
	const Event& eventToWatch )
{
	unsigned int value = readOverlapped_uint32_noChecksum( socket, eventToWatch );
	unsigned int checksumRemote = readOverlapped_uint32_noChecksum( socket, eventToWatch );

	unsigned int checksumLocal = checksum( value );
	if( checksumLocal != checksumRemote )
		throw SocketException( "readOverlapped_uint32", WSAGetLastError() );
	return value;
}

const size_t maxSize = 1024*1024;
unsigned int readLength( SocketWrapper& socket,
	const Event& eventToWatch )
{
	unsigned int value = readOverlapped_uint32( socket, eventToWatch );
	if( value > maxSize )
		throw SocketException( "readOverlapped_length", WSAGetLastError() );
	return value;
}


unsigned int readOverlapped_float( SocketWrapper& socket,
	const Event& eventToWatch )
{
	unsigned int value = readOverlapped_uint32( socket, eventToWatch );
	unsigned int checksumRemote = readOverlapped_uint32_noChecksum( socket, eventToWatch );
	unsigned int checksumLocal = checksum(value);
	if( checksumLocal != checksumRemote )
		throw SocketException( "readOverlapped_float", WSAGetLastError() );

	const float* floatBytes = reinterpret_cast<const float*>( &value );
	return *floatBytes;
}

std::vector<char> readOverlapped_charArray( SocketWrapper& socket, 
	const Event& eventToWatch )
{
	unsigned int length = readLength( socket, eventToWatch );
	unsigned int checksumRemote = readOverlapped_uint32_noChecksum( socket, eventToWatch );

	std::vector<char> result = readOverlapped( length, socket, eventToWatch );
	unsigned int checksumLocal = checksum(result);
	if( checksumLocal != checksumRemote )
		throw SocketException( "readOverlapped_charArray", WSAGetLastError() );

	return result;
}

std::vector<unsigned int> readOverlapped_unsignedArray( SocketWrapper& socket, 
	const Event& eventToWatch )
{
	unsigned int length = readLength( socket, eventToWatch );
	unsigned int checksumRemote = readOverlapped_uint32_noChecksum( socket, eventToWatch );

	std::vector<unsigned int> result( length );
	if( length > 0 )
	{
		readOverlapped( reinterpret_cast<char*>( &result[0] ),
			length*4, socket, eventToWatch );
		std::transform( result.begin(), result.end(), result.begin(), &ntohl );
	}

	unsigned int checksumLocal = checksum( result );
	if( checksumLocal != checksumRemote )
		throw SocketException( "readOverlapped_unsignedArray", WSAGetLastError() );

	return result;
}

std::vector<float> readOverlapped_floatArray( SocketWrapper& socket, 
	const Event& eventToWatch )
{
	unsigned int length = readLength( socket, eventToWatch );
	unsigned int checksumRemote = readOverlapped_uint32_noChecksum( socket, eventToWatch );

	std::vector<float> result( length );
	if( length > 0 )
	{
		readOverlapped( reinterpret_cast<char*>( &result[0] ),
			length*4, socket, eventToWatch );

		unsigned int* uintPtr = reinterpret_cast<unsigned int*>( &result[0] );
		std::transform( uintPtr, uintPtr + length, uintPtr, &ntohl );
	}

	unsigned int checksumLocal = checksum( result );
	if( checksumLocal != checksumRemote )
		throw SocketException( "readOverlapped_unsignedArray", WSAGetLastError() );

	return result;
}

std::vector<vl::Vec3f> readOverlapped_vec3fArray( SocketWrapper& socket, 
	const Event& eventToWatch )
{
	unsigned int length = readLength( socket, eventToWatch );
	unsigned int checksumRemote = readOverlapped_uint32_noChecksum( socket, eventToWatch );

	std::vector<vl::Vec3f> result( length );
	if( length > 0 )
	{
		readOverlapped( reinterpret_cast<char*>( &result[0] ),
			length*4*3, socket, eventToWatch );

		unsigned int* uintPtr = reinterpret_cast<unsigned int*>( &result[0] );
		std::transform( uintPtr, uintPtr + length*3, uintPtr, &ntohl );
	}

	unsigned int checksumLocal = checksum( result );
	if( checksumLocal != checksumRemote )
		throw SocketException( "readOverlapped_unsignedArray", WSAGetLastError() );

	return result;
}

class AutoIncrement
{
public:
	AutoIncrement( IncrementDecrement& incDec )
		: incDec_( incDec )
	{
		incDec_.increment();
	}

	~AutoIncrement()
	{
		incDec_.decrement();
	}

private:
	IncrementDecrement& incDec_;
};

void MainFrame::sampleTreeRemote(const std::string& address, boost::shared_ptr<Event> samplingEvent)
{
	RandomStream random;
	boost::shared_ptr<SimulationTree> tree;
    std::pair<float, float> timeRange;
	boost::shared_ptr<SocketWrapper> socket;

	{
		WSAData wsaData;
		int err = WSAStartup(MAKEWORD(1, 1), &wsaData);
		// @todo need to communicate to the main app that sampling is broken
		assert( err == 0 );
	}

	size_t samplingRate;

	const DWORD baseRetryTimeout = 200;
	DWORD retryTimeout = baseRetryTimeout;

	boost::shared_ptr<AutoIncrement> autoIncrement;
	while( this->sampling_ )
	{
		try
		{
			if( this->pauseSampling_ || tree != this->currentTree_ || (tree && tree->timeRange() != timeRange) )
			{
				tree.reset();
				socket.reset();
				autoIncrement.reset();
				retryTimeout = baseRetryTimeout;
			}

			if( !this->pauseSampling_ )
				tree = this->currentTree_;

			if( !tree )
			{
				// wait for event so we're not spinlocking
				WaitForSingleObject( samplingEvent->handle(), INFINITE );
				continue;
			}

			samplingRate = tree->samplingRate();
			if( !socket )
			{
				// first create the socket
				{
					SOCKET sock = WSASocket(
						AF_INET,
						SOCK_STREAM,
						IPPROTO_TCP,
						NULL, 
						0,
						WSA_FLAG_OVERLAPPED );
					if( sock == INVALID_SOCKET )
						throw SocketException( "WSASocket", WSAGetLastError() );

					socket.reset( new SocketWrapper( sock ) );
				}

				const int nPort = 10037;
				u_long nRemoteAddress = LookupAddress( address.c_str() );
				if (nRemoteAddress == INADDR_NONE)
					throw SocketException( "gethostbyname", WSAGetLastError() );

				in_addr Address;
				memcpy(&Address, &nRemoteAddress, sizeof(u_long)); 
				sockaddr_in sinRemote;
				sinRemote.sin_family = AF_INET;
				sinRemote.sin_addr.s_addr = nRemoteAddress;
				sinRemote.sin_port = htons(nPort);

				// now connect to the remote server
				{
					if(connect( socket->socket(), reinterpret_cast<sockaddr*>(&sinRemote), sizeof(sockaddr_in)) ==
						SOCKET_ERROR)
					{
						throw SocketException( "connect", WSAGetLastError() );
					}
				}

				// Now, send the scene over
				socket->write( samplingRate );

				std::string scene;
				{
					MechModelObjectPtr mechModel = tree->scene()->toMechModel();
					std::ostringstream oss;
					mechModel->dump( oss, std::string() );
					scene.swap( oss.str() );
				}
				socket->write( scene.size() );
				socket->write( scene );

				std::string parentSceneDesc;
				if( tree->parentScene() )
				{
					MechModelObjectPtr mechModel = tree->parentScene()->toMechModel();
					std::ostringstream oss;
					mechModel->dump( oss, std::string() );
					parentSceneDesc.swap( oss.str() );
				}
				socket->write( parentSceneDesc.size() );
				socket->write( parentSceneDesc );

				std::vector<size_t> fixedObjects = tree->fixedObjects();
				std::vector<unsigned int> fixedObjectsTmp( fixedObjects.begin(), fixedObjects.end() );
				socket->write( fixedObjectsTmp );

				float startTime = tree->startTime();
				socket->write( startTime );

				timeRange = tree->timeRange();
				socket->write( timeRange );

				std::vector<PiecewisePath> fixedObjectsHistory = tree->fixedObjectsHistory();
				socket->write( fixedObjectsHistory.size() );
				for( size_t i = 0; i < fixedObjectsHistory.size(); ++i )
				{
					std::vector<unsigned int> positionTimes;
					std::vector<float> positionCoeffs;
					std::vector<unsigned int> rotationTimes;
					std::vector<float> rotations;
					std::vector<float> rotationCoeffs;

					fixedObjectsHistory[i].dumpToArrays( 
						positionTimes,
						positionCoeffs,
						rotationTimes,
						rotations,
						rotationCoeffs );

					socket->write( positionTimes );
					socket->write( positionCoeffs );
					socket->write( rotationTimes );
					socket->write( rotations );
					socket->write( rotationCoeffs );
				}
			}

			assert( socket );

			unsigned int numPaths = readOverlapped_uint32( 
				*socket, *samplingEvent );
			std::vector<CompressedPath> paths;
			paths.reserve( numPaths );

			if( !autoIncrement )		// wait until we actually get a path back listening before this:
				autoIncrement.reset( new AutoIncrement( this->numMachines_ ) );

			std::vector< std::vector<vl::Vec3f> > cmPaths;
			cmPaths.reserve( numPaths );

			for( size_t i = 0; i < numPaths; ++i )
			{
				paths.push_back( CompressedPath() );
				paths.back().positionTimes = 
					readOverlapped_charArray( *socket, *samplingEvent );
				for( size_t i = 0; i < 3; ++i )
					paths.back().positionLinearCoeffs[i] = readOverlapped_charArray( *socket, *samplingEvent );
				paths.back().positionQuadraticCoeffs[1] = readOverlapped_charArray( *socket, *samplingEvent );

				paths.back().rotationTimes = 
					readOverlapped_charArray( *socket, *samplingEvent );
				for( size_t i = 0; i < 3; ++i )
					paths.back().rotationCoeffs[i] = readOverlapped_charArray( *socket, *samplingEvent );

				cmPaths.push_back( readOverlapped_vec3fArray( *socket, *samplingEvent ) );
			}

			// now need to read the snapshots
			std::vector<Path::TimedSnapshot> snapshots;
			unsigned int nSnaps = readOverlapped_uint32( *socket, *samplingEvent );
			snapshots.reserve( nSnaps );
			for( size_t i = 0; i < nSnaps; ++i )
			{
				float time = readOverlapped_float( *socket, *samplingEvent );
				unsigned int numStates = readOverlapped_uint32( *socket, *samplingEvent );

				snapshots.push_back( Path::TimedSnapshot(time, std::vector<RigidStaticState>(numStates)) );
				std::vector<RigidStaticState>& state = snapshots.back().second;

				readOverlapped( reinterpret_cast<char*>( &state[0] ),
					numStates*4*7, *socket, *samplingEvent );
				unsigned int* uintPtr = reinterpret_cast<unsigned int*>( &state[0] );
				std::transform( uintPtr, uintPtr + numStates*7, uintPtr, &ntohl );
			}

			std::vector<float> metricScores = 
				readOverlapped_floatArray( *socket, *samplingEvent );

			tree->addPath( paths, snapshots, cmPaths, samplingRate, metricScores );
		}
		catch( SocketException& e )
		{
			std::string errMsg = GetOsErrorMessage( e.errorCode() );

			std::ostringstream oss;
			oss << "Sockets: error in call " << e.function() << ": " << errMsg;
			this->logError( oss.str() );

			if( socket )
			{
				// need to fire and forget here
				char end[] = { "END" };
				WSABUF buf;
				buf.len = sizeof(end);
				buf.buf = end;
				DWORD numberOfBytesSent;

				int result = WSASend(
					socket->socket(),
					&buf,
					1,
					&numberOfBytesSent,
					0,
					NULL,
					NULL );
			}

			socket.reset();
			autoIncrement.reset();

			// we use this rather than sleep so that if the program exits
			//   or some other important event happens we break out of it
			WaitForSingleObject( samplingEvent->handle(), retryTimeout );
			retryTimeout *= 2;
		}
	}

	WSACleanup();
}
#endif

void MainFrame::addNodes( const SimulationTree* tree, const std::deque<const Path*>& nodes, const boost::dynamic_bitset<>& active )
{
	ConstSimulationTreePtr currentTree = this->getTree();
	if( tree != currentTree.get() )
		return;

	ReaderWriterLock::ScopedWriterLock lineCacheLock( this->lineCacheMutex_ );
	if( this->showAllPaths() )
	{
		for( size_t iPath = 0; iPath < nodes.size(); ++iPath )
		{
			this->addToLineCache( nodes[iPath], active[iPath] );
		}
	}
	else
	{
		for( size_t iPath = 0; iPath < nodes.size(); ++iPath )
		{
			if( active[iPath] )
				this->addToLineCache( nodes[iPath], true );
		}
	}
}

void MainFrame::setNodes( const SimulationTree* tree, const boost::dynamic_bitset<>& nodes )
{
	ConstSimulationTreePtr currentTree = this->getTree();
	if( tree == currentTree.get() )
		this->setSamples( currentTree, nodes );
}

void MainFrame::OnCreateConvexHull(wxCommandEvent& event)
{
	MainFrame::SelectedPoints selectedPoints = this->selectedPoints_;
	boost::shared_ptr<MeshObject> meshObj = 
		boost::dynamic_pointer_cast<MeshObject>( selectedPoints.first );
	if( !meshObj )
		return;

	ConvexHullPtr newHull = meshObj->convexHullForPoints( selectedPoints.second );
	ConvexHullList oldHulls = meshObj->getConvexHulls();
	ConvexHullList newHulls = oldHulls;
	newHulls.push_back( newHull );

	ActionPtr action( new SetHullsAction(meshObj, oldHulls, newHulls) );
	this->performAction( action );
}

void MainFrame::OnFitConvexHulls(wxCommandEvent& event)
{
	std::deque<ActionPtr> actions;

	std::deque<SceneGraphElementPtr> selected = this->getSelected();
	for( std::deque<SceneGraphElementPtr>::const_iterator selectedItr = selected.begin();
		selectedItr != selected.end(); ++selectedItr )
	{
		boost::shared_ptr<MeshObject> meshObj = 
			boost::dynamic_pointer_cast<MeshObject>( *selectedItr );
		if( !meshObj )
			continue;

		ConvexHullList hulls = meshObj->getMesh()->computeConvexHulls();

		actions.push_back( 
			ActionPtr( 
				new SetHullsAction(meshObj, meshObj->getConvexHulls(), hulls) ) );
	}

	boost::shared_ptr<Action> compoundAction( new CompoundAction(actions) );
	this->performAction( compoundAction );
}

void MainFrame::OnFitBoxes(wxCommandEvent& event)
{
	try
	{
		std::deque<SceneGraphElementPtr> boxes;

		std::deque<SceneGraphElementPtr> selected = this->getSelected();
		for( std::deque<SceneGraphElementPtr>::const_iterator selectedItr = selected.begin();
			selectedItr != selected.end(); ++selectedItr )
		{
			boost::shared_ptr<const MeshObject> meshObject = 
				boost::dynamic_pointer_cast<const MeshObject>( *selectedItr );
			if( !meshObject )
				continue;

			TriangleMeshPtr mesh = meshObject->getMesh();
			assert( mesh );

			std::vector<Triangle> triangles = mesh->triangles();
			std::vector<vl::Vec3f> positions = mesh->positions();
			std::vector<vl::Vec3d> handledNormals;
			for( size_t iTriangle = 0; iTriangle < triangles.size(); ++iTriangle )
			{
				Triangle tri = triangles[iTriangle];
				vl::Vec3d normal = vl::norm( vl::cross( 
					toVec3d(positions[ tri[1] ]) - toVec3d(positions[ tri[0] ]),
					toVec3d(positions[ tri[2] ]) - toVec3d(positions[ tri[0] ]) ) );
				if( fabs(normal[2]) > 1e-3 )
					continue;

				bool seenBefore = false;
				for( std::vector<vl::Vec3d>::const_iterator itr = handledNormals.begin();
					itr != handledNormals.end(); ++itr )
				{
					if( vl::len(vl::cross( normal, *itr )) < 1e-3 )
						seenBefore = true;
				}

				if( seenBefore )
					continue;

				handledNormals.push_back( normal );

				double angle = atan2( normal[1], normal[0] );
				vl::Vec3d rot( 0.0, 0.0, angle );
				vl::Vec3d alongAxis = vl::norm( vl::cross( normal, vl::Vec3d(0, 0, 1) ) );

				vl::Vec3d center = toVec3d(positions[ tri[0] ]);
				double maxLen = 0.0;
				double sign = 1.0;
				for( int iVert = 1; iVert < 3; ++iVert )
				{
					double diff = vl::dot( toVec3d(positions[ tri[iVert] ]) - center, alongAxis );
					if( fabs(diff) > maxLen )
					{
						maxLen = fabs(diff);
						sign = signum(diff);
					}
				}
				center += sign*0.5*maxLen*alongAxis;

				double offAxisScale = 0.4;
				vl::Vec3d scale( offAxisScale, maxLen, 1.0 );

				vl::Vec3d pos = center /*- 1.0*offAxisScale*normal*/;

				boost::shared_ptr<PhysicsObject> box( 
					new BoxObject( meshObject->name() + "_proxyGeom_" +
							boost::lexical_cast<std::string>(boxes.size()), 
						scene_->defaultMaterial() ) );
				box->setScale( toVec3f(scale) );
				box->setRotate( toVec3f(rot) );
				box->setTranslate( vl::Vec3f( pos[0], pos[1], 0.5f ) );
				box->setPivot( vl::Vec3f(0.5f, 0.0f, 0.0f) );

				boxes.push_back( box );
			}
		}

		addObjects( boxes );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating boxes: " << e.message();
		errMsg( oss.str() );
	}
}

void MainFrame::OnFuse(wxCommandEvent& event)
{
	try
	{
		std::deque<SceneGraphElementPtr> objects =
			this->getSelected();
		if( objects.empty() )
			return;		

		boost::shared_ptr<CombinedObject> group(
			new CombinedObject( 
				scene_->uniqueObjectName("fused"), scene_->defaultMaterial() ) );

		std::deque<ActionPtr> actions;
		{
			ActionPtr addObject( 
				new CreateObjectAction(group, scene_) );
			actions.push_back( addObject );
		}

		// the proper way to do this would be to implement Least Common Ancestor
		SceneGraphElementPtr commonParent = objects.front()->parent().lock();
		for( std::deque<SceneGraphElementPtr>::const_iterator obIter = objects.begin();
			obIter != objects.end(); ++obIter )
		{
			if( (*obIter)->parent().lock() != commonParent )
				commonParent = scene_->root();
		}

		ActionPtr reparentAction( new ReparentAction(group, commonParent) );
		actions.push_back( reparentAction );

		// temporarily change the parent so the xforms are right
		group->setParent( commonParent );
		std::deque<SceneGraphElementPtr> parents(objects.size(), group);
		std::deque<ActionPtr> reparentActions = this->reparent( objects, parents, true );
		std::copy( reparentActions.begin(), reparentActions.end(),
			std::back_inserter(actions) );
		group->setParent( SceneGraphElementPtr() );

		ActionPtr changeSelection( 
			new ChangeSelectionAction(this, 
				this->getSelected(), std::deque<SceneGraphElementPtr>(1, group),
				this->selectedConstraint(), ConstraintPtr() ) );
		actions.push_back( changeSelection );

		boost::shared_ptr<Action> compoundAction( new CompoundAction(actions) );
		this->performAction( compoundAction );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error fusing: " << e.message();
		errMsg( oss.str() );
	}
}

void MainFrame::OnUnfuse(wxCommandEvent& event)
{
	try
	{
		std::deque<SceneGraphElementPtr> selected = this->getSelected();

		std::deque<SceneGraphElementPtr> objects;
		std::deque<SceneGraphElementPtr> newParents;
		for( std::deque<SceneGraphElementPtr>::iterator selectedItr = selected.begin();
			selectedItr != selected.end(); ++selectedItr )
		{
			SceneGraphElementPtr parent = *selectedItr;

			SceneGraphElementPtr grandParent = (*selectedItr)->parent().lock();
			assert( grandParent );

			boost::shared_ptr<CombinedObject> combined = 
				boost::dynamic_pointer_cast<CombinedObject>( *selectedItr );
			if( !combined )
				continue;

			for( SceneGraphElement::child_iterator childItr = parent->children_begin();
				childItr != parent->children_end(); ++childItr )
			{
				objects.push_back( *childItr );
				newParents.push_back( grandParent );
			}
		}

		this->reparent( 
			objects,
			newParents,
			true );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error ungrouping: " << e.message();
		errMsg( oss.str() );
	}
}

void MainFrame::OnGroup(wxCommandEvent& event)
{
	try
	{
		std::deque<SceneGraphElementPtr> objects =
			this->getSelected();
		if( objects.empty() )
			return;		

		boost::shared_ptr<TransformGroup> group(
			new TransformGroup( 
				scene_->uniqueObjectName("group") ) );

		std::deque<ActionPtr> actions;
		{
			ActionPtr addObject( 
				new CreateObjectAction(group, scene_) );
			actions.push_back( addObject );
		}

		// the proper way to do this would be to implement Least Common Ancestor
		SceneGraphElementPtr commonParent = objects.front()->parent().lock();
		for( std::deque<SceneGraphElementPtr>::const_iterator obIter = objects.begin();
			obIter != objects.end(); ++obIter )
		{
			if( (*obIter)->parent().lock() != commonParent )
				commonParent = scene_->root();
		}

		ActionPtr reparentAction( new ReparentAction(group, commonParent) );
		actions.push_back( reparentAction );

		// temporarily change the parent so the xforms are right
		group->setParent( commonParent );
		std::deque<SceneGraphElementPtr> parents(objects.size(), group);
		std::deque<ActionPtr> reparentActions = this->reparent( objects, parents, true );
		std::copy( reparentActions.begin(), reparentActions.end(),
			std::back_inserter(actions) );
		group->setParent( SceneGraphElementPtr() );

		ActionPtr changeSelection( 
			new ChangeSelectionAction(this, 
				this->getSelected(), std::deque<SceneGraphElementPtr>(1, group),
				this->selectedConstraint(), ConstraintPtr() ) );
		actions.push_back( changeSelection );

		boost::shared_ptr<Action> compoundAction( new CompoundAction(actions) );
		this->performAction( compoundAction );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error grouping: " << e.message();
		errMsg( oss.str() );
	}
}

void MainFrame::OnUngroup(wxCommandEvent& event)
{
	try
	{
		std::deque<SceneGraphElementPtr> selected = this->getSelected();

		std::deque<SceneGraphElementPtr> objects;
		std::deque<SceneGraphElementPtr> newParents;

		std::deque<ActionPtr> deleteActions;
		for( std::deque<SceneGraphElementPtr>::iterator selectedItr = selected.begin();
			selectedItr != selected.end(); ++selectedItr )
		{
			SceneGraphElementPtr parent = *selectedItr;

			SceneGraphElementPtr grandParent = (*selectedItr)->parent().lock();
			assert( grandParent );

			for( SceneGraphElement::child_iterator childItr = parent->children_begin();
				childItr != parent->children_end(); ++childItr )
			{
				objects.push_back( *childItr );
				newParents.push_back( grandParent );
			}

			boost::shared_ptr<TransformGroup> transformGroup = 
				boost::dynamic_pointer_cast<TransformGroup>( parent );
			if( transformGroup )
			{
				ActionPtr reparentAction( new ReparentAction(transformGroup, SceneGraphElementPtr()) );
				deleteActions.push_back( reparentAction );

				ActionPtr deleteObject(
					new DeleteObjectAction( transformGroup, scene_ ) );
				deleteActions.push_back( deleteObject );
			}
		}

		std::deque<ActionPtr> actions = this->reparent( 
			objects,
			newParents,
			true );

		std::copy( deleteActions.begin(), deleteActions.end(),
			std::back_inserter(actions) );
		boost::shared_ptr<Action> compoundAction( new CompoundAction(actions) );
		this->performAction( compoundAction );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error ungrouping: " << e.message();
		errMsg( oss.str() );
	}
}

std::deque<ActionPtr> MainFrame::reparent( const std::deque<SceneGraphElementPtr>& objects,
	const std::deque<SceneGraphElementPtr>& newParents,
	bool adjustTransforms )
{
	std::deque<ActionPtr> actions;

	std::deque<SceneGraphElementPtr> newSelection;

	assert( objects.size() == newParents.size() );
	for( size_t iObject = 0; iObject < objects.size(); ++iObject )
	{
		SceneGraphElementPtr object = objects[iObject];
		SceneGraphElementPtr parent = object->parent().lock();
		SceneGraphElementPtr newParent = newParents[iObject];
		vl::Mat4f newParentTransform = newParent->transform();
		
		TransformedElementPtr transformed = 
			boost::dynamic_pointer_cast<TransformedSceneGraphElement>( object );
		if( transformed )
		{
			// need to update the transform to factor out the parent transform
			vl::Mat4f childTransform = transformed->transform();
			// @todo user a more intelligent solver here rather than inverting matrix
			vl::Mat4f differentialTransform = 
				inv(newParentTransform) * childTransform;

			std::pair< vl::Mat3f, vl::Mat3f > matPr = 
				polarDecomposition( toMat3f(differentialTransform) );

			ActionPtr rotateAction( 
				new SetRotateAction(
					transformed,
					transformed->rotate(),
					matToEuler(matPr.first) ) );

			vl::Mat3f scaleMat = matPr.second;

			// check for shear
			float epsilon = 1e-3;
			for( vl::Int i = 0; i < 3; ++i )
			{
				for( vl::Int j = 0; j < 3; ++j )
				{
					if( i == j )
						continue;

					assert( fabs(scaleMat[i][j]) < epsilon );
				}
			}

			vl::Vec3f scale( scaleMat[0][0], scaleMat[1][1], scaleMat[2][2] );
			if( !transformed->allowNonUniformScale() )
				scale[0] = scale[1] = scale[2] = (scale[0] + scale[1] + scale[2]) / 3.0f;

			if( transformed->allowUniformScale() )
			{
				ActionPtr scaleAction(
					new SetScaleAction(
						transformed, 
						transformed->scale(),
						vl::Vec3f( scale[0], scale[1], scale[2] ) ) );
				actions.push_back( scaleAction );
			}

			ActionPtr translateAction(
				new SetTranslateAction(
					transformed,
					transformed->translate(),
					vl::xform( differentialTransform, vl::Vec3f(vl::vl_0) ) ) );

			actions.push_back( rotateAction );
			actions.push_back( translateAction );
		}

		ActionPtr reparentAction(
			new ReparentAction(
				object, newParent ) );
		actions.push_back( reparentAction );

		newSelection.push_back( object );
	}

	ActionPtr changeSelection( new ChangeSelectionAction(this, 
		this->getSelected(), newSelection,
		this->selectedConstraint(), ConstraintPtr() ) );
	actions.push_back( changeSelection );
	return actions;
}

void MainFrame::OnFixCenterOfMass(wxCommandEvent& event)
{
	try
	{
		std::deque<SceneGraphElementPtr> toHandle( 1, scene_->root() );
		std::deque<ActionPtr> actions;

		while( !toHandle.empty() )
		{
			SceneGraphElementPtr current = toHandle.front();
			toHandle.pop_front();
			boost::shared_ptr<CombinedObject> combined = 
				boost::dynamic_pointer_cast<CombinedObject>( current );
			if( combined )
			{
				dMass mass = combined->getOdeMassProps();
				vl::Vec3f center = vl::xform(combined->rigidTransform(), 
					vl::Vec3f(mass.c[0], mass.c[1], mass.c[2]) );
				
				vl::Mat4f transform = combined->transform();

				SceneGraphElementPtr parent = combined->parent().lock();
				vl::Mat4f parentTransform = parent->transform();

				vl::Vec3f newCenter = vl::xform( inv( parentTransform ), center );

				ActionPtr changeTranslateAction( 
					new SetTranslateAction( combined, combined->translate(), newCenter ) );
				actions.push_back( changeTranslateAction );

				for( SceneGraphElement::child_iterator childItr = combined->children_begin();
					childItr != combined->children_end(); ++childItr )
				{
					TransformedElementPtr transformed = 
						boost::dynamic_pointer_cast<TransformedSceneGraphElement>( *childItr );
					if( !transformed )
						continue;

					vl::Vec3f childTranslate = transformed->translate();
					vl::Vec3f childCenterLocal = vl::xform( vl::inv(transform), center );
					vl::Vec3f newChildTranslate = childTranslate - childCenterLocal;

					ActionPtr changeTranslateChildAction( 
						new SetTranslateAction( transformed, childTranslate, newChildTranslate ) );
					actions.push_back( changeTranslateChildAction );
				}
			}
			else
			{
				std::copy( current->children_begin(), current->children_end(),
					std::front_inserter(toHandle) );
			}
		}

		boost::shared_ptr<Action> compoundAction( new CompoundAction(actions) );
		this->performAction( compoundAction );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error fixing center of mass: " << e.message();
		errMsg( oss.str() );
	}
}


void MainFrame::OnCenterPivot(wxCommandEvent& event)
{
	try
	{
		std::deque<ActionPtr> actions;

		std::deque<SceneGraphElementPtr> objects =
			this->getSelected();

		for( std::deque<SceneGraphElementPtr>::iterator objectIter = objects.begin();
			objectIter != objects.end(); ++objectIter )
		{
			SceneGraphElementPtr object = *objectIter;
			TransformedElementPtr transformed = 
				boost::dynamic_pointer_cast<TransformedSceneGraphElement>( object );
			boost::shared_ptr<CombinedObject> combined = 
				boost::dynamic_pointer_cast<CombinedObject>( object );

			if( !transformed ) // pivot change only applies to transformed groups
				continue;

			vl::Vec3f center;
			vl::Mat4f transform = transformed->transform();

			if( combined )
			{
				dMass mass = combined->getOdeMassProps();
				center = vl::xform(combined->rigidTransform(), 
					vl::Vec3f(mass.c[0], mass.c[1], mass.c[2]) );
			}
			else
			{
				// find the bounding box for the object plus all its children:
				BoundingBox3f bounds;
				std::deque<SceneGraphElementPtr> queue( 1, *objectIter );
				while( !queue.empty() )
				{
					SceneGraphElementPtr current = queue.front();
					bounds.expand( current->bounds() );
					queue.pop_front();
					std::copy( current->children_begin(), current->children_end(), 
						std::front_inserter(queue) );
				}

				center = 0.5f*(bounds.minimum() + bounds.maximum());
			}

			SceneGraphElementPtr parent = transformed->parent().lock();
			vl::Mat4f parentTransform = parent->transform();

			// move pivot to center
			vl::Vec3f newPivot = vl::xform( inv( transform ), center );

			vl::Vec3f prevPivot = transformed->pivot();
			vl::Vec3f newTranslate = vl::xform(inv(parentTransform), center);

			ActionPtr changePivotAction(
				new SetPivotAction( transformed, prevPivot, newPivot ) );
			ActionPtr changeTranslateAction( 
				new SetTranslateAction( transformed, transformed->translate(), newTranslate ) );
			actions.push_back( changePivotAction );
			actions.push_back( changeTranslateAction );
		}

		boost::shared_ptr<Action> compoundAction( new CompoundAction(actions) );
		this->performAction( compoundAction );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error centering pivot: " << e.message();
		errMsg( oss.str() );
	}
	
}

void MainFrame::OnShowAll(wxCommandEvent& event)
{
	std::deque<SceneGraphElementPtr> hiddenObjects;

	for( Scene::object_iter iter = scene_->begin_objects();
		iter != scene_->end_objects(); ++iter )
	{
		if( !(*iter)->visible() )
			hiddenObjects.push_back( (*iter) );
	}

	std::vector<bool> oldState( hiddenObjects.size(), false );
	std::vector<bool> newState( hiddenObjects.size(), true );
	ActionPtr action( new SetVisibleAction(hiddenObjects, oldState, newState) );
	this->performAction( action );
}


void MainFrame::OnHide(wxCommandEvent& event)
{
	std::deque<SceneGraphElementPtr> selected = this->getSelected();
	std::vector<bool> oldState( selected.size() );
	for( size_t iObject = 0; iObject < selected.size(); ++iObject )
		oldState[iObject] = selected[iObject]->visible();

	std::vector<bool> newState( selected.size(), false );
	ActionPtr action( new SetVisibleAction(selected, oldState, newState) );
	this->performAction( action );
}

void MainFrame::OnSetMode( wxCommandEvent& event )
{
	int mode = modeChoice_->GetSelection();

	if( mode == 0 )
	{
		setMode( MODE_EDIT );
	}
	else if( mode == 1 )
	{
		if( this->getMode() == MODE_SOLVE )
			return;

		// just to make extra-sure:
		setMode( MODE_SOLVE );
	}
	else
		assert( false );
}

MainFrame::SimulatorType MainFrame::simulatorType() const
{
	if( simulatorMenu_->IsChecked(MENU_SIMULATE_SIMULATOR_ODE_NATIVE) )
		return SIMULATOR_ODE_NATIVE;
	else if( simulatorMenu_->IsChecked(MENU_SIMULATE_SIMULATOR_ODE_PENALTY) )
		return SIMULATOR_ODE_PENALTY;
#ifdef USE_NOVODEX
	else if( simulatorMenu_->IsChecked(MENU_SIMULATE_SIMULATOR_NOVODEX) )
		return SIMULATOR_NOVODEX;
#endif
#ifdef USE_NEWTON
	else if( simulatorMenu_->IsChecked(MENU_SIMULATE_SIMULATOR_NEWTON) )
		return SIMULATOR_NEWTON;
#endif
#ifdef USE_BULLET
	else if( simulatorMenu_->IsChecked(MENU_SIMULATE_SIMULATOR_BULLET) )
		return SIMULATOR_BULLET;
	else if( simulatorMenu_->IsChecked(MENU_SIMULATE_SIMULATOR_ODE_BULLET) )
		return SIMULATOR_ODE_BULLET;
#endif
	else
	{
		assert( false );
		return SIMULATOR_ODE_NATIVE;
	}
}

MainFrame::MouseMode MainFrame::mouseMode() const
{
	if( this->toolBar_->GetToolState( myID_SELECT_RECT ) )
		return MODE_MOUSE_SELECT;
	else if( this->toolBar_->GetToolState( myID_GOOD_CONSTRAINT_RECT ) )
		return MODE_MOUSE_ADDITIVE_CONSTRAINT;
	else if( this->toolBar_->GetToolState( myID_BAD_CONSTRAINT_RECT ) )
		return MODE_MOUSE_SUBTRACTIVE_CONSTRAINT;
	else if( this->toolBar_->GetToolState( myID_REFINE_RECT ) )
		return MODE_MOUSE_REFINE;
	else if( this->toolBar_->GetToolState( myID_SKETCH_ACTION ) )
		return MODE_MOUSE_SKETCH_ACTION;
	
	assert( false );
	return MODE_MOUSE_SELECT;
}

void MainFrame::updateCameraMenu()
{
	// add all the cameras:
	{
		assert( scene()->begin_cameras() != scene()->end_cameras() );

		int i = 0;
		for( Scene::camera_iterator iter = scene()->begin_cameras();
			iter != scene()->end_cameras(); ++iter )
		{
			int currentId = CAMERA_MENU_START + i++;
			std::string name = (*iter)->name();

			wxMenuItem *item = cameraMenu_->FindItem(currentId);
			if(item)
			{
				item->SetText( wxT( name.c_str() ) );
			}
			else
			{
				item = cameraMenu_->AppendRadioItem( currentId, wxT( name.c_str() ) );
				//wxObjectEventFunction fn = &MainFrame::OnSetCamera;
				Connect( currentId, wxEVT_COMMAND_MENU_SELECTED, 
					wxCommandEventHandler(MainFrame::OnSetCamera) );
			}

			item->Check( *iter == currentCanvas_->cameraWrapper() );
		}

		// clear out any entries above and beyond the current total
		for( ; cameraMenu_->FindItem(CAMERA_MENU_START + i) != 0; 
			++i )
		{
			cameraMenu_->Destroy(CAMERA_MENU_START + i);
		}
	}

}

#ifdef WINDOWS_MEDIA
void MainFrame::OnEditSettings(wxCommandEvent& event)
{
	settingsDialog_->setAudioCodec( this->selectedAudioCodec_ );
	settingsDialog_->setVideoCodec( this->selectedVideoCodec_ );

	if( settingsDialog_->ShowModal() == wxID_OK )
	{
		this->selectedAudioCodec_ = settingsDialog_->getAudioCodec();
		this->selectedVideoCodec_ = settingsDialog_->getVideoCodec();
	}
}
#endif

void MainFrame::OnCameraProperties(wxCommandEvent& event)
{
	CameraWrapperPtr cam = currentCanvas_->cameraWrapper();

	// transfer settings over to dialog:
	cameraPropertyDialog_->setName( cam->name() );
	cameraPropertyDialog_->setFarZ( cam->camera()->zMax() );
	cameraPropertyDialog_->setNearZ( cam->camera()->zMin() );

	boost::shared_ptr<ProjectiveFixedYCamera> projCam =
		boost::dynamic_pointer_cast<ProjectiveFixedYCamera>( cam->camera() );
	cameraPropertyDialog_->enableFOV( projCam );
	if( projCam )
		cameraPropertyDialog_->setFOV( projCam->FOV() );

	while( cameraPropertyDialog_->ShowModal() == wxID_OK )
	{
		std::deque<ActionPtr> actions;

		std::string newName = cameraPropertyDialog_->name();

		if( newName != cam->name() )
		{
			if( scene()->camera( newName ) )
			{
				errMsg( "Name '" + newName + "' already in use." );
				continue;
			}			
			
			// need to change name
			ActionPtr action( 
				new CameraNameChangeAction(cam, this->scene(), cam->name(), newName, this) );
			actions.push_back( action );
		}

		// need an action for changing the camera properties:
		{
			boost::shared_ptr<Action> propAction(
				new CameraPropertyChangeAction(cam->camera(), cameraPropertyDialog_->nearZ(), cameraPropertyDialog_->farZ()) );
			actions.push_back( propAction );
		}

		if( projCam )
		{
			boost::shared_ptr<Action> fovAction( 
				new ProjectiveCameraChangeFOVAction(projCam, cameraPropertyDialog_->fov()) );
			actions.push_back( fovAction );
		}

		boost::shared_ptr<CompoundAction> action( new CompoundAction(actions) );
		this->performAction( action );
		return;
	}
}

void MainFrame::OnSetCamera( wxCommandEvent& event )
{
	// have to do this because for whatever reason the event won't tell
	//   me which menu item was pressed
	for( int i = 0; ; ++i )
	{
		int id = CAMERA_MENU_START+i;
		wxMenuItem* item = cameraMenu_->FindItem(id);
		if( !item )
			break;

		if( item->IsChecked() )
		{
			wxString nm = item->GetLabel();
			std::string name( nm.c_str() );
			assert( !name.empty() );
			currentCanvas_->setCamera( scene()->camera( name ) );
			return;
		}
	}

	// should never reach here
	assert( false );
}

void MainFrame::OnWindowSimSorter(wxCommandEvent& event)
{
	this->bestScoresFrame_->Show();
	this->bestScoresFrame_->Raise();
}

void MainFrame::OnWindowHierarchyView(wxCommandEvent& event)
{
	this->hierarchyFrame_->Show();
	this->hierarchyFrame_->Raise();
}

void MainFrame::setMode( MainFrame::Mode mode, bool force )
{
	if( mode == mode_ && !force )
		return;

	// might eventually not want to do this
	this->undoQueue_.clear();
	this->redoQueue_.clear();

	mode_ = mode;

	switch( mode )
	{
	case MODE_EDIT:
		{
			modeChoice_->SetSelection(0);

//			graphFrame_->Hide();
//			plotFrame_->Hide();
			bestScoresFrame_->Hide();

			this->updateMenus();
			this->setSelected( this->getSelected() );

			break;
		}
	case MODE_SOLVE:
		{
			modeChoice_->SetSelection(1);

			if( this->trees_.empty() /*|| !dirtyObjects_.empty()*/ )
			{
				std::pair<double, double> timeRange = 
					this->timeView_->activeTimeRange();

				ScenePtr newScene( new Scene(*this->scene_) );
				SimulationTreePtr tree( new SimulationTree(
					"tree" + boost::lexical_cast<std::string>(this->treeCounter_++), 
					newScene, 120, simulationFactory() ) );
				tree->setTimeRange( timeRange );
				addTree( tree );
				this->setPath( tree, 0 );
			}

			std::for_each( samplingEvents_.begin(), samplingEvents_.end(), boost::bind( &Event::notify_all, _1 ) );
			showProperties( false );

			this->updateMenus();
			break;
		}
	default:
		assert( false );
	}
}

bool MainFrame::wireframeOnSelected() const
{
	return wireframeOnSelected_;
}

void MainFrame::OnPlay(wxCommandEvent& event)
{
	if( playing_ )
	{
		playButton_->SetBitmapLabel( playBitmap_ );
		playing_ = false;
	}
	else
	{
		playButton_->SetBitmapLabel( pauseBitmap_ );
		playing_ = true;
	}

	playButton_->Refresh( FALSE );
	std::for_each( canvases_.begin(), canvases_.end(), 
		boost::bind( &GLView::Refresh, _1, false, (const wxRect*) NULL ) );
}

bool MainFrame::playing() const
{
	return playing_;
}

void MainFrame::removeTree( SimulationTreePtr tree )
{
	{
		boost::mutex::scoped_lock lock( this->treesMutex_ );

		for( std::deque<CallbackRegistrarPtr>::iterator iter = this->callbackRegistrars_.begin();
			iter != callbackRegistrars_.end(); )
		{
			if( (*iter)->tree() == tree )
				iter = callbackRegistrars_.erase(iter);
			else
				++iter;
		}

        tree->setTimeRange( std::pair<float, float>(0.0, 0.0) );
		trees_.erase( std::remove( trees_.begin(), trees_.end(), tree ) );
	}

	this->bestScoresFrame_->removePanel( tree );
}

void MainFrame::addTree( SimulationTreePtr tree )
{
	{
		boost::mutex::scoped_lock lock( this->treesMutex_ );
		this->trees_.push_back( tree );
		this->callbackRegistrars_.push_back(
			CallbackRegistrarPtr( new SimulationNodeCallbackRegistrar(tree, this) ) );
		std::for_each( samplingEvents_.begin(), samplingEvents_.end(), boost::bind( &Event::notify_all, _1 ) );
	}

	this->bestScoresFrame_->newPanel( tree, tree->name() );
}

void MainFrame::OnMenuUpdateSample( wxCommandEvent& event )
{
	if( this->getPath() == 0 )
	{
		errMsg( "No path selected." );
		return;
	}

	this->startSubtree( 
		this->getTree(),
		this->piecewisePaths(),
		this->piecewisePathsFrameRate(),
		0.0f,
		std::deque<ConstPhysicsObjectPtr>(),
		false);
}

void MainFrame::startSubtree( 
	SimulationTreePtr parent,
	const std::vector<PiecewisePath>& paths, 
	size_t pathsFrameRate,
	float startTime,
	const std::deque<ConstPhysicsObjectPtr>& active,
	bool backwards )
{
	{
		std::ostringstream oss;
		oss << "Refining tree " << parent->name() << " at time " 
			<< startTime << "; " << active.size() << "objects.";
		this->log( oss.str() );
	}

    if( paths.empty() )
        return;

	std::pair<double, double> timeRange = 
		this->timeView_->activeTimeRange();
	ScenePtr newScene( new Scene(*this->scene_) );
	SimulationTreePtr tree( new
		SimulationTree( 
			"subtree",
			newScene,
			parent, 
			active, 
			paths, 
			pathsFrameRate,
			startTime,
            backwards) );
	tree->setTimeRange( timeRange );

	std::deque< ActionPtr > actions;
	actions.push_back( ActionPtr( 
		new AddTreeAction(this, tree) ) );
	actions.push_back( ActionPtr( 
		new SetPathSelectionAction( 
			SetPathSelectionAction::PathWithTree(tree, 0),
			SetPathSelectionAction::PathWithTree(this->getTree(), this->getPath()),
			this ) ) );

	boost::shared_ptr<Action> action( new CompoundAction(actions) );
	this->performAction( action );
}

// We need to be able to do things like refine, etc. using just the timeline
// This is a function the timeline can call to report that such an event has
// happened
void MainFrame::applyTimeSelection( float startTime, float endTime )
{
    if( this->getMode() == MainFrame::MODE_SOLVE
        && this->mouseMode() == MainFrame::MODE_MOUSE_REFINE )
    {
        bool backwards = (startTime > endTime);

        SimulationTreePtr currentTree = this->getTree();
		std::deque<SceneGraphElementPtr> selected = this->getSelected();
		std::vector<ConstPhysicsObjectPtr> dynamicObjects = currentTree->dynamicObjects();
		std::deque<ConstPhysicsObjectPtr> refinedObjects;

		for( std::vector<ConstPhysicsObjectPtr>::const_iterator itr = dynamicObjects.begin();
			itr != dynamicObjects.end(); ++itr )
		{
			ConstSceneGraphElementPtr current = *itr;
			while( current )
			{
				if( std::find_if( selected.begin(), selected.end(), 
					NameEquals(current->name()) ) != selected.end() )
				{
					refinedObjects.push_back( *itr );
					break;
				}

				current = current->parent().lock();
			}
		}

        this->startSubtree(
		    currentTree,
		    this->piecewisePaths(),
		    this->piecewisePathsFrameRate(),
		    startTime,
		    std::deque<ConstPhysicsObjectPtr>( dynamicObjects.begin(), dynamicObjects.end() ),
		    backwards);
    }
}

#ifdef WINDOWS_MEDIA
std::vector<std::string> MainFrame::audioCodecs() const
{
	return this->audioCodecDescriptions_;
}

std::vector<std::string> MainFrame::videoCodecs() const
{
	return this->videoCodecDescriptions_;
}

std::pair< std::vector<std::string>,
	std::vector< boost::shared_ptr<IWMStreamConfig> > >
		MainFrame::findAvailableCodecs( REFGUID codecType, bool (*checkType) (WM_MEDIA_TYPE*), bool requireVBR )
{
	boost::shared_ptr<IWMProfileManager> profileManager = this->profileManager_;
	assert( profileManager );

	// Get the codec information interface.
	boost::shared_ptr<IWMCodecInfo3> codecInfo;
	{
		IWMCodecInfo3* pCodecInfo = NULL;
		HRESULT res = profileManager->QueryInterface(IID_IWMCodecInfo3, (void**)&pCodecInfo);
		if( FAILED(res) )
			throw Exception( "Unable to load query interface." );

		codecInfo.reset( pCodecInfo, &GenericRelease<IWMCodecInfo3> );
	}

	// Get the number of codecs for which there is information.
	DWORD cEntries = 0;
	{
		HRESULT res = codecInfo->GetCodecInfoCount(codecType, &cEntries);
		if( FAILED(res) )
			throw Exception( "Error retrieving codecs" );
	}

	std::vector<std::string> codecDescriptions;
	std::vector< boost::shared_ptr<IWMStreamConfig> > streamConfigs;

	// Find the codecIndex of the codec corresponding to the requested subytpe.
	for(DWORD codecIndex = 0; codecIndex < cEntries; codecIndex++)
	{
		/*
		if( requireVBR )
		{
			DWORD numPasses = 1;
			HRESULT res = codecInfo->SetCodecEnumerationSetting( codecType,
				codecIndex,
				g_wszNumPasses,
				WMT_TYPE_DWORD,
				reinterpret_cast<const BYTE *>( &numPasses ),
				sizeof(DWORD) );
			if( FAILED(res) )
				throw Exception( "Error retrieving codecs" );

			BOOL vbr = TRUE;
			res = codecInfo->SetCodecEnumerationSetting( codecType,
				codecIndex,
				g_wszVBREnabled,
				WMT_TYPE_BOOL,
				reinterpret_cast<const BYTE *>( &vbr ),
				sizeof(BOOL) );
			if( FAILED(res) )
				throw Exception( "Error retrieving codecs" );
		}
		*/

		wchar_string codecName;
		DWORD codecNameLength;
		{
			HRESULT res = codecInfo->GetCodecName(codecType, codecIndex, NULL, &codecNameLength);
			if( FAILED(res) )
				throw Exception( "Error retrieving codecs" );

			boost::scoped_array<WCHAR> codecNameArr( new WCHAR[codecNameLength] );
			res = codecInfo->GetCodecName(codecType, codecIndex, codecNameArr.get(), &codecNameLength);
			if( FAILED(res) )
				throw Exception( "Error retrieving codecs" );

			codecName = wchar_string( codecNameArr.get() );
		}

		// Get the number of formats supported for the codec.
		DWORD numFormats = 0;
		{
			HRESULT res = codecInfo->GetCodecFormatCount(codecType, codecIndex, &numFormats);
			if( FAILED(res) )
				throw Exception( "Error retrieving codecs" );
		}

		for(DWORD formatIndex = 0; formatIndex < numFormats; ++formatIndex)
		{
			// Get the format description for each codec. 
			// first, call to get the string length
			boost::shared_ptr<IWMStreamConfig> streamConfig;
			DWORD formatDescLength;
			wchar_string formatDesc;
			{
				HRESULT res = codecInfo->GetCodecFormatDesc(
					codecType, codecIndex, formatIndex, NULL, NULL, &formatDescLength);
				if( FAILED(res) )
					throw Exception( "Error retrieving codecs" );
			}

			{
				std::vector<WCHAR> formatDescArr( formatDescLength );
				IWMStreamConfig* pStreamConfig;
				HRESULT res = codecInfo->GetCodecFormatDesc(
					codecType, codecIndex, formatIndex, &pStreamConfig, &formatDescArr[0], &formatDescLength);
				if( FAILED(res) )
					throw Exception( "Error retrieving codecs" );

				streamConfig.reset( pStreamConfig, &GenericRelease<IWMStreamConfig> );

				formatDesc = wchar_string( &formatDescArr[0] );
			}

			// Get the media properties interface.
			boost::shared_ptr<IWMMediaProps> props;
			{
				IWMMediaProps* pProps;
				HRESULT res = streamConfig->QueryInterface(IID_IWMMediaProps, (void**)&pProps);
				if( FAILED(res) )
					throw Exception( "Error retrieving codecs" );

				props.reset( pProps, &GenericRelease<IWMMediaProps> );
			}

			// Get the size required for the media type structure.
			DWORD cbType;
			{
				HRESULT res = props->GetMediaType(NULL, &cbType);
				if( FAILED(res) )
					throw Exception( "Error retrieving codecs" );
			}

			std::vector<BYTE> typeBytes( cbType );
			WM_MEDIA_TYPE* pType = reinterpret_cast<WM_MEDIA_TYPE*>( &typeBytes[0] );
			{
				HRESULT res = props->GetMediaType(pType, &cbType);
				if( FAILED(res) )
					throw Exception( "Error retrieving codecs" );
			}

			if( !checkType( pType ) )
				continue;

			std::ostringstream oss;
			oss << toString(codecName);
			if( !formatDesc.empty() )
				oss << " (" << toString(formatDesc) << ")";
			std::string fullName = oss.str();
			codecDescriptions.push_back( fullName );

			streamConfigs.push_back( streamConfig );
		}
	}

	return std::make_pair( codecDescriptions, streamConfigs );
}
#endif // WINDOWS_MEDIA

void MainFrame::clearTrees()
{
	this->log( "Cleared all trees." );

	boost::mutex::scoped_lock lock( this->treesMutex_ );
	this->bestScoresFrame_->clearPanels();
	this->callbackRegistrars_.clear();
	this->trees_.clear();
	this->setPath( SimulationTreePtr(), 0 );
}


void MainFrame::setPath( SimulationTreePtr tree, const Path* path )
{
	{
		boost::mutex::scoped_lock lock( this->currentPathMutex_ );
		if( tree == this->currentTree_ && (path == this->currentPath_) )
			return;

		this->currentPath_ = path;
		if( tree != this->currentTree_ )
		{
			{
				std::string currentName; if( currentTree_ ) currentName = currentTree_->name();
				std::string newName; if( tree ) currentName = tree->name();
				this->log( "Changed current tree from " + currentName + " to " + newName );
			}

			this->currentTree_ = tree;

			std::for_each( samplingEvents_.begin(), samplingEvents_.end(), boost::bind( &Event::notify_all, _1 ) );
			this->updateDirtyObjects();

			if( tree )
			{
				std::pair<double, double> timeRange = 
					this->timeView_->activeTimeRange();
	            tree->setTimeRange( timeRange );

				this->setSamples( tree, tree->activePaths() );
			}
			else
			{
				this->clearSamples();
			}
		}
	}

	{
		boost::mutex::scoped_lock lock( this->computeHistoryPendingMutex_ );
		this->computeHistoryPending_.push_back( 
			std::make_pair(tree, path) );
		this->computeHistoryPendingCondition_.notify_all();

		while( !this->computeHistoryPending_.empty() )
			this->computeHistoryPendingCondition_.wait( lock );
	}

	this->bestScoresFrame_->setPath( tree, path );
//	this->plotFrame_->changeSelected();
	this->updateSelectedPath();
	//this->setTime( 0.0, true );
}

void MainFrame::computeHistory()
{
	RandomStream random;

	SimulationTreePtr currentTree = this->currentTree_;
	while( !computeHistoryDone_ )
	{
		const Path* currentPath = 0;

		{
			boost::mutex::scoped_lock lock( this->computeHistoryPendingMutex_ );
			while( computeHistoryPending_.empty() && !computeHistoryDone_ )
				computeHistoryPendingCondition_.wait( lock );

			while( !this->computeHistoryPending_.empty() )
			{
				ReaderWriterLock::ScopedWriterLock lock( this->historyLock_ );

				currentTree = this->computeHistoryPending_.front().first;
				currentPath = this->computeHistoryPending_.front().second;
				this->computeHistoryPending_.pop_front();

				this->piecewisePaths_.clear();
				this->stateHistory_.clear();

				if( currentPath == 0 )
					continue;

				historyComplete_.reset();

				this->piecewisePathsFrameRate_ = currentPath->frameRate();

				std::vector<Simulation::State> states = currentPath->samples();
				for( std::vector<Simulation::State>::const_iterator iter = states.begin();
					iter != states.end(); ++iter )
				{
					stateHistory_.push_back( StateWithNode(*iter, 0) );
				}
                assert( std::is_sorted(stateHistory_.begin(), stateHistory_.end()) );
			}
		}

		computeHistoryPendingCondition_.notify_all();

		if( currentPath )
		{
			std::vector<PiecewisePath> piecewisePaths
				= currentTree->decode( currentPath );

			ReaderWriterLock::ScopedWriterLock lock( this->historyLock_ );
			this->piecewisePaths_.swap( piecewisePaths );
			historyComplete_.notify_all();
		}
	}
}

void MainFrame::changeSelected()
{
	updateMenus();
	if( arePropertiesShown() )
		this->propertyDialog_->updateSelection();

	this->hierarchyFrame_->setSelected( selected_ );
	this->bestScoresFrame_->updateSelected();

//	if( plotFrame_->IsShown() )
//		plotFrame_->changeSelected();

	this->updateViews();
	timeView_->Refresh();
}

void MainFrame::sceneChanged()
{
	this->updateDirtyObjects();

	dirty_ = true;
	if( arePropertiesShown() )
		this->propertyDialog_->updateValues();

	timeView_->Refresh();

	this->updateViews();

	// updates the MouseAction
	this->setSelected( this->getSelected() );
	updateMenus();
}

void MainFrame::updateDirtyObjects()
{
	dirtyObjects_.clear();

	SimulationTreePtr tree = this->currentTree_;
	if( !tree )
		return;

	ConstScenePtr curScene = this->scene_;
	ConstScenePtr treeScene = tree->scene();
	std::deque<ConstSceneGraphElementPtr> changed = 
		symmetricDiffScene( *curScene, *treeScene );
	for( std::deque<ConstSceneGraphElementPtr>::const_iterator itr = changed.begin();
		itr != changed.end(); ++itr )
		dirtyObjects_.insert( itr->get() );
}

ConstSimulationTreePtr MainFrame::getTree() const
{
	boost::mutex::scoped_lock lock( this->currentPathMutex_ );
	SimulationTreePtr result( this->currentTree_ );
	return result;
}

SimulationTreePtr MainFrame::getTree()
{
	boost::mutex::scoped_lock lock( this->currentPathMutex_ );
	SimulationTreePtr result( this->currentTree_ );
	return result;
}

const Path* MainFrame::getPath() const
{
	boost::mutex::scoped_lock lock( this->currentPathMutex_ );
	const Path* result = this->currentPath_;
	return result;
}

boost::shared_ptr<Scene> MainFrame::scene()
{
	return scene_;
}

boost::shared_ptr<const Scene> MainFrame::scene() const
{
	return scene_;
}

Simulation::State MainFrame::state() const
{
	if( simulation_ )
		return simulation_->getSimulationState();

	ConstSimulationTreePtr tree = this->getTree();
	ReaderWriterLock::ScopedReaderLock historyLock( this->historyLock_ );

	if( currentPath_ && tree )
	{
		assert( !this->stateHistory_.empty() || !this->piecewisePaths_.empty() );

		float curtime = this->getTime();
		if( this->piecewisePaths_.empty() )
		{
			std::vector<RigidDynamicState> tmp; // for some reason gcc wants this
			const Simulation::State state( tmp, curtime );
			typedef std::deque<MainFrame::StateWithNode>::const_iterator Iter;
			Iter iter = std::lower_bound( stateHistory_.begin(), stateHistory_.end(), 
				MainFrame::StateWithNode(state, (const SimulationNode*) 0), 
				Compare1st<MainFrame::StateWithNode>() );
			if( iter != stateHistory_.begin() )
				--iter;

			return iter->first;
		}
		else
		{
			float iFrame = curtime 
				* boost::numeric_cast<float>( this->piecewisePathsFrameRate_ );
			std::vector< RigidDynamicState > states;
			states.reserve( piecewisePaths_.size() );

			for( size_t iPath = 0; iPath < piecewisePaths_.size(); ++iPath )
			{
				states.push_back( 
					RigidDynamicState(
						piecewisePaths_[iPath].position( iFrame ),
						piecewisePaths_[iPath].rotation( iFrame ),
						piecewisePaths_[iPath].linearVelocity( iFrame ),
						piecewisePaths_[iPath].angularVelocity( iFrame ) ) );
			}

			return Simulation::State( states, curtime );
		}
	}

	return Simulation::State();
}

std::vector<Simulation::ContactPoint> MainFrame::contactPoints() const
{
	if( this->simulation_ )
		return simulation_->contactPoints();
	else
		return std::vector<Simulation::ContactPoint>();
}

void MainFrame::OnCopy(wxCommandEvent& event)
{
	copy( this->getSelected() );
}

void MainFrame::OnPaste(wxCommandEvent& event)
{
	paste();
}

void MainFrame::OnDuplicate(wxCommandEvent& event)
{
	duplicate( this->getSelected() );
}

void MainFrame::OnUndo(wxCommandEvent& event)
{
	undo();
}

void MainFrame::OnRedo(wxCommandEvent& event)
{
	redo();
}

void MainFrame::updateMenus()
{
	editMenu_->Enable( MENU_EDIT_UNDO, !undoQueue_.empty() );
	editMenu_->Enable( MENU_EDIT_REDO, !redoQueue_.empty()  );

	editMenu_->Enable( MENU_EDIT_PASTE, !clipboard_.empty() );
	editMenu_->Enable( MENU_EDIT_DUPLICATE, !this->getSelected().empty() );
	editMenu_->Enable( MENU_EDIT_COPY, !this->getSelected().empty() );

	menuBar_->EnableTop(0, true);
	menuBar_->EnableTop(1, true);
	menuBar_->EnableTop(2, this->getMode() == MainFrame::MODE_EDIT );
	menuBar_->EnableTop(3, true);
	menuBar_->EnableTop(4, true);
	//menuBar_->EnableTop(5, this->getMode() == MainFrame::MODE_EDIT );
	menuBar_->EnableTop(5, true );
	menuBar_->EnableTop(6, true);

	fileMenu_->Enable( MENU_FILE_LOAD_TREE, true );
	fileMenu_->Enable( MENU_FILE_SAVE_TREE, this->getMode() == MainFrame::MODE_SOLVE );

	windowMenu_->Enable( MENU_WINDOW_SIM_SORTER, this->getMode() == MainFrame::MODE_SOLVE );

	simulateMenu_->Enable( MENU_SIMULATE_SAMPLING, this->getMode() == MainFrame::MODE_SOLVE );
	simulateMenu_->Enable( MENU_SIMULATE_RUN, this->getMode() == MainFrame::MODE_EDIT );
	simulateMenu_->Enable( MENU_SIMULATE_UPDATE_PATH, 
		this->getMode() == MainFrame::MODE_SOLVE && !dirtyObjects_.empty() );

#ifdef SOUND
	soundMenu_->Enable( MENU_SOUND_ENABLE_SOUND, 
		this->simulatorType() == SIMULATOR_ODE_NATIVE ||
		this->simulatorType() == SIMULATOR_ODE_BULLET );
	soundMenu_->Enable( MENU_SOUND_COMPUTE_ALL_MODES, !computingModes_ );
	soundMenu_->Enable( MENU_SOUND_COMPUTE_OBJECT_MODES, !computingModes_ );
	soundMenu_->Enable( MENU_SOUND_STOP_COMPUTING_MODES, computingModes_ );
#endif // SOUND

	toolBar_->EnableTool(myID_GOOD_CONSTRAINT_RECT, this->getMode() == MainFrame::MODE_SOLVE);
	toolBar_->EnableTool(myID_BAD_CONSTRAINT_RECT, this->getMode() == MainFrame::MODE_SOLVE);
	toolBar_->EnableTool(myID_REFINE_RECT, this->getMode() == MainFrame::MODE_SOLVE
		&& dirtyObjects_.empty() );
	toolBar_->EnableTool(myID_SKETCH_ACTION, this->getMode() == MainFrame::MODE_SOLVE );

	toolBar_->EnableTool(myID_CURSOR_TOOLBAR, this->getMode() == MainFrame::MODE_EDIT);
	toolBar_->EnableTool(myID_TRANSLATE_TOOLBAR, this->getMode() == MainFrame::MODE_EDIT);
	toolBar_->EnableTool(myID_ROTATE_TOOLBAR, this->getMode() == MainFrame::MODE_EDIT);
	toolBar_->EnableTool(myID_SCALE_TOOLBAR, this->getMode() == MainFrame::MODE_EDIT);
	toolBar_->EnableTool(myID_CHANGEPIVOT_TOOLBAR, this->getMode() == MainFrame::MODE_EDIT);
	toolBar_->EnableTool(myID_SELECTPOINTS_TOOLBAR, this->getMode() == MainFrame::MODE_EDIT);
}

void MainFrame::copy( const std::deque<SceneGraphElementPtr>& objects )
{
	clipboard_ = objects;
	updateMenus();
}

std::deque<SceneGraphElementPtr> MainFrame::cloneObjects( const std::deque<SceneGraphElementPtr>& objects )
{
	// @todo: make this not run in quadratic time
	std::deque< SceneGraphElementPtr > result;

	STLEXT hash_set<std::string> usedNames;
	SceneGraphElementMap oldToNewMap;
	for( std::deque< SceneGraphElementPtr >::const_iterator iter = objects.begin();
		iter != objects.end(); ++iter )
	{
		result.push_back( this->cloneObject( *iter, usedNames, oldToNewMap ) );
	}

	// now need to fix all the joints
	std::deque< SceneGraphElementPtr > toHandle = result;
	while( !toHandle.empty() )
	{
		SceneGraphElementPtr current = toHandle.front();
		toHandle.pop_front();

		boost::shared_ptr<Joint> joint = 
			boost::dynamic_pointer_cast<Joint>( current );
		if( !joint )
		{
			std::copy( current->children_begin(), current->children_end(),
				std::front_inserter(toHandle) );
			continue;
		}

		std::vector<ConstPhysicsObjectPtr> objects = joint->objects();
		std::vector<PhysicsObjectPtr> newObjects; newObjects.reserve( objects.size() );
		for( std::vector<ConstPhysicsObjectPtr>::const_iterator objectItr = objects.begin();
			objectItr != objects.end(); ++objectItr )
		{
			SceneGraphElementMap::const_iterator mapItr = 
				oldToNewMap.find( objectItr->get() );
			if( mapItr == oldToNewMap.end() )
			{
				// this joint connects something not dup'd, so we will delete it
				break;
			}

			PhysicsObjectPtr physicsObject = 
				boost::dynamic_pointer_cast<PhysicsObject>(mapItr->second);
			assert( physicsObject );
			newObjects.push_back( physicsObject );
		}

		if( newObjects.size() == objects.size() )
		{
			assert( newObjects.size() == 2 );
			joint->setObjects( PhysicsObjectPair(newObjects.at(0), newObjects.at(1)) );

			std::copy( current->children_begin(), current->children_end(),
				std::front_inserter(toHandle) );
		}
		else
		{
			// need to delete this joint; we will also have to toss its
			//   children since we don't know what to do with them
			//   (replace it with a basic group?)
			std::deque< SceneGraphElementPtr >::iterator itr = 
				std::find( result.begin(), result.end(), current );
			if( itr != result.end() )
				result.erase( itr );

			current->parent().lock()->removeChild( current );
		}
	}

	return result;
}

SceneGraphElementPtr MainFrame::cloneObject( ConstSceneGraphElementPtr object, 
	STLEXT hash_set<std::string>& used, SceneGraphElementMap& oldToNewMap )
{
	std::string name = object->name();
	std::string::size_type pos = name.size();
	while( pos > 1 && isdigit(name[pos-1]) )
		--pos;

	name = name.substr(0, pos);
	for( int i = 1; ; ++i )
	{
		std::ostringstream oss;
		oss << name << i;
		std::string newName = oss.str();

		if( scene_->object(oss.str()) )
			continue;

		if( used.find( newName ) != used.end() )
			continue;

		name = newName;
		break;
	}

	SceneGraphElementPtr result = object->clone();
	result->setName( name );
	used.insert( name );

	std::pair<SceneGraphElementMap::iterator, bool> inserted = 
		oldToNewMap.insert( SceneGraphElementMap::value_type(object.get(), result) );
	// should not have any duplicate entries!
	assert( inserted.second );

	for( SceneGraphElement::const_child_iterator childItr = object->children_begin();
		childItr != object->children_end(); ++childItr )
	{
		SceneGraphElementPtr newChild = this->cloneObject( *childItr, used, oldToNewMap );
		newChild->setParent( result );
		result->addChild( newChild );
	}

	return result;
}

void MainFrame::paste()
{
	addObjects( cloneObjects(clipboard_) );
}

void MainFrame::duplicate( const std::deque<SceneGraphElementPtr>& objects )
{
	addObjects( cloneObjects(objects) );
}

void MainFrame::OnCreateBox(wxCommandEvent& event)
{
	try
	{
		PhysicsObjectPtr boxObject(
			new BoxObject( 
				scene_->uniqueObjectName("box"),
				scene_->defaultMaterial() ) );
		addObject( boxObject );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating box: " << e.message();
		errMsg( oss.str() );
	}
}

void MainFrame::OnCreateCylinder(wxCommandEvent& event)
{
	try
	{
		PhysicsObjectPtr cylinderObject(
			new CylinderObject( 
				scene_->uniqueObjectName("cylinder"),
				scene_->defaultMaterial() ) );
		addObject( cylinderObject );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating cylinder: " << e.message();
		errMsg( oss.str() );
	}
}

void MainFrame::OnCreateCone(wxCommandEvent& event)
{
	try
	{
		PhysicsObjectPtr coneObject(
			new ConeObject( 
				scene_->uniqueObjectName("cone"),
				scene_->defaultMaterial() ) );
		addObject( coneObject );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating cone: " << e.message();
		errMsg( oss.str() );
	}
}


void MainFrame::OnCreateCappedCylinder(wxCommandEvent& event)
{
	try
	{
		PhysicsObjectPtr cylinderObject(
			new CappedCylinderObject( 
				scene_->uniqueObjectName("cappedCylinder"),
				scene_->defaultMaterial() ) );
		addObject( cylinderObject );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating capsule: " << e.message();
		errMsg( oss.str() );
	}
}

void MainFrame::OnCreateSpheres(wxCommandEvent& event)
{
	std::string filename;
	try
	{
		boost::shared_ptr<wxFileDialog> dlg(
			new wxFileDialog(
				this, 
				"Load .obj geometry",
				"", 
				"", 
				".sph files (*.sph)|*.sph|All files(*.*)|*.*",
				wxFD_OPEN, 
				wxDefaultPosition), 
			boost::mem_fn(&wxFileDialog::Destroy) );

		if ( dlg->ShowModal() == wxID_OK )
		{
			filename = dlg->GetPath().GetData();

			boost::filesystem::path path( filename );
			std::string name = path.leaf();
			std::string::size_type dotPos = name.rfind(".");
			if( dotPos != std::string::npos )
				name = name.substr(0, dotPos);

			SphereWrapperTree tree( filename );
			std::vector<SphereWrapperTree::Sphere> spheres = tree.level(4);

			size_t iSphere = 1;
			std::deque< SceneGraphElementPtr > newSpheres;
			for( std::vector<SphereWrapperTree::Sphere>::const_iterator itr = spheres.begin();
				itr != spheres.end(); ++itr )
			{
				std::string sphereName = name;
				sphereName += "_sphere";
				sphereName += boost::lexical_cast<std::string>( iSphere++ );

				PhysicsObjectPtr sphereObject(
					new SphereObject( sphereName,
						scene_->defaultMaterial() ) );
				sphereObject->setTranslate( itr->center );
				float scale = itr->radius;
				sphereObject->setScale( vl::Vec3f(scale, scale, scale) );
				newSpheres.push_back( sphereObject );
			}

			addObjects( newSpheres );
		}
	}
	catch( ParserException& e )
	{
		std::ostringstream oss;
		oss << "Error loading '" << filename << "': " << e;
		errMsg( oss.str() );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error loading '" << filename << "': " << e.message();
		errMsg( oss.str() );
	}
}

void MainFrame::OnCreateSphere(wxCommandEvent& event)
{
	try
	{
		PhysicsObjectPtr sphereObject(
			new SphereObject( 
				scene_->uniqueObjectName("sphere"),
				scene_->defaultMaterial() ) );
		addObject( sphereObject );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating sphere: " << e.message();
		errMsg( oss.str() );
	}
}



PhysicsObjectPair MainFrame::getJointObjects() const
{
	std::deque<SceneGraphElementPtr> selected = this->getSelected();
	if( selected.size() != 2 )
		throw Exception( "Must select exactly 2 objects to create joint." );

	PhysicsObjectPtr first = boost::dynamic_pointer_cast<PhysicsObject>(selected[0]);
	PhysicsObjectPtr second = boost::dynamic_pointer_cast<PhysicsObject>(selected[1]);
	if( !first || !second )
		throw Exception( "Joints must be applied to leaf objects" );

	return PhysicsObjectPair( first, second );
}

void MainFrame::OnCreatePlaneJoint(wxCommandEvent& event)
{
    try
    {
	    std::deque<SceneGraphElementPtr> selected = this->getSelected();

        std::deque<PhysicsObjectPtr> physicsObjects;
        for( std::deque<SceneGraphElementPtr>::const_iterator itr = selected.begin();
            itr != selected.end(); ++itr )
        {
            PhysicsObjectPtr physOb = boost::dynamic_pointer_cast<PhysicsObject>( *itr );
            if( physOb )
                physicsObjects.push_back( physOb );
        }

        if( physicsObjects.empty() )
            throw Exception( "No physics objects selected." );

        std::deque<SceneGraphElementPtr> joints;
        for( std::deque<PhysicsObjectPtr>::const_iterator physObItr = physicsObjects.begin();
            physObItr != physicsObjects.end(); ++physObItr )
        {
		    boost::shared_ptr<Joint> joint( 
			    new PlaneJoint( scene_->uniqueObjectName("ballJoint") ) );
            joint->setObjects( *physObItr );
            joints.push_back( joint );
        }
        this->addObjects( joints );
    }
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating plane joint: " << e.message();
		errMsg( oss.str() );
	}
}

void MainFrame::OnCreateBallJoint(wxCommandEvent& event)
{
	try
	{
		boost::shared_ptr<Joint> joint( 
			new BallAndSocketJoint( scene_->uniqueObjectName("ballJoint") ) );
		joint->setObjects( this->getJointObjects() );
		addObject( joint );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating ball joint: " << e.message();
		errMsg( oss.str() );
	}
}

void MainFrame::OnCreateSliderJoint(wxCommandEvent& event)
{
	try
	{
		boost::shared_ptr<Joint> joint( 
			new SliderJoint( scene_->uniqueObjectName("sliderJoint") ) );
		joint->setObjects( this->getJointObjects() );
		addObject( joint );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating ball joint: " << e.message();
		errMsg( oss.str() );
	}
}


void MainFrame::OnCreateUniversalJoint(wxCommandEvent& event)
{
	try
	{
		boost::shared_ptr<Joint> joint( 
			new UniversalJoint( scene_->uniqueObjectName("ballJoint") ) );
		joint->setObjects( this->getJointObjects() );
		addObject( joint );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating universal joint: " << e.message();
		errMsg( oss.str() );
	}
}

void MainFrame::OnCreateHinge2Joint(wxCommandEvent& event)
{
	try
	{
		boost::shared_ptr<Joint> joint( 
			new Hinge2Joint( scene_->uniqueObjectName("hinge2Joint") ) );
		joint->setObjects( this->getJointObjects() );
		addObject( joint );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating hinge-2 joint: " << e.message();
		errMsg( oss.str() );
	}
}

void MainFrame::OnCreateHingeJoint(wxCommandEvent& event)
{
	try
	{
		boost::shared_ptr<Joint> joint( 
			new HingeJoint( scene_->uniqueObjectName("hingeJoint") ) );
		joint->setObjects( this->getJointObjects() );
		addObject( joint );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating hinge joint: " << e.message();
		errMsg( oss.str() );
	}
}

/*
void MainFrame::OnCreateRotMotor1Joint(wxCommandEvent& event)
{
	try
	{
		boost::shared_ptr<Joint> joint( 
			new OneAxisRotationalMotor( scene_->uniqueObjectName("hinge2Joint") ) );
		joint->setObjects( this->getJointObjects() );
		addObject( joint );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating 1-axis rotational motor joint: " << e.message();
		errMsg( oss.str() );
	}
}
*/

void MainFrame::OnCreateOrthoCam(wxCommandEvent& event)
{
	try
	{
		boost::shared_ptr<Camera> cam( new OrthoCamera );
		CameraWrapperPtr wrapper( 
			new CameraWrapper(scene_->uniqueCameraName("orthoCam"), cam) );
        addCamera( wrapper );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating camera: " << e.message();
		errMsg( oss.str() );
	}
}

void MainFrame::OnCreatePerspectiveCam(wxCommandEvent& event)
{
	try
	{
		boost::shared_ptr<Camera> cam( new ProjectiveFixedYCamera );
		CameraWrapperPtr wrapper( 
			new CameraWrapper(scene_->uniqueCameraName("perspectiveCam"), cam) );
        addCamera( wrapper );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating camera: " << e.message();
		errMsg( oss.str() );
	}
}

void MainFrame::OnCreatePlane(wxCommandEvent& event)
{
	try
	{
		PhysicsObjectPtr planeObject(
			new PlaneObject( scene_->uniqueObjectName("plane"), scene_->defaultMaterial() ) );
		addObject( planeObject );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating plane: " << e.message();
		errMsg( oss.str() );
	}
}

void MainFrame::OnCreateObj(wxCommandEvent& event)
{
	try
	{
		boost::shared_ptr<wxFileDialog> dlg(
			new wxFileDialog(
				this, 
				"Load .obj geometry",
				"", 
				"", 
				".obj files (*.obj)|*.obj|All files(*.*)|*.*",
				wxFD_OPEN, 
				wxDefaultPosition), 
			boost::mem_fn(&wxFileDialog::Destroy) );

		if ( dlg->ShowModal() == wxID_OK )
		{
			std::string filename(dlg->GetPath().GetData());

			boost::filesystem::path path( filename );
			std::string name = path.leaf();
			std::string::size_type dotPos = name.rfind(".");
			if( dotPos != std::string::npos )
				name = name.substr(0, dotPos);

			TriangleMeshPtr mesh = scene_->triMesh( filename );
			PhysicsObjectPtr meshObj(
				new MeshObject( 
					scene_->uniqueObjectName(name), 
					mesh, 
					scene_->defaultMaterial() ) );
			addObject( meshObj );
		}
	}
	catch( ParserException& e )
	{
		std::ostringstream oss;
		oss << "Error creating mesh: " << e;
		errMsg( oss.str() );
	}
	catch( Exception& e )
	{
		std::ostringstream oss;
		oss << "Error creating mesh: " << e.message();
		errMsg( oss.str() );
	}
}

void MainFrame::showProperties( bool show )
{
	if( show )
		this->propertyDialog_->updateSelection();

	this->propSizer_->Show( this->propertyDialog_, show );
	this->propSizer_->Layout();

	windowMenu_->Check(MENU_WINDOW_PROPERTIES, show);
	this->currentCanvas_->SetFocus();
}

bool MainFrame::arePropertiesShown() const
{
	return propertyDialog_->IsShown();
}

void MainFrame::OnWindowProperties(wxCommandEvent& event)
{
	this->showProperties( windowMenu_->IsChecked(MENU_WINDOW_PROPERTIES) );
}

void MainFrame::performAction( ActionPtr action )
{
	if( !action )
		return;

	redoQueue_.clear();

	undoQueue_.push_back( action );
	undoQueue_.back()->doIt();

	// don't really need to set scene changed in this case:
	if( action->changesScene() )
		sceneChanged();

	if( action->changesHierarchy() )
		this->hierarchyFrame_->setScene( scene_, selected_ );

	this->updateViews();
}

float MainFrame::getTime() const
{
	return time_;
}

void MainFrame::setTime( float time, bool moving )
{
	this->time_ = time;

	if( moving )
	{
		for( Scene::camera_iterator iter = scene()->begin_cameras();
			iter != scene()->end_cameras(); ++iter )
		{
			(*iter)->setTime( time );
		}
	}

	Refresh(FALSE);
	Update();
}

void MainFrame::undo()
{
	if( !undoQueue_.empty() )
	{
		ActionPtr action = undoQueue_.back();
		action->undoIt();
		redoQueue_.push_back( action );
		undoQueue_.pop_back();

		if( action->changesScene() )
			sceneChanged();

		if( action->changesHierarchy() )
			this->hierarchyFrame_->setScene( scene_, selected_ );

		this->updateViews();
	}
}

void MainFrame::redo()
{
	if( !redoQueue_.empty() )
	{
		ActionPtr action = redoQueue_.back();
		action->doIt();
		undoQueue_.push_back( action );
		redoQueue_.pop_back();

		if( action->changesScene() )
			sceneChanged();

		if( action->changesHierarchy() )
			this->hierarchyFrame_->setScene( scene_, selected_ );

		this->updateViews();
	}
}

std::vector<PiecewisePath> MainFrame::piecewisePaths() const
{
	WaitForSingleObject( historyComplete_.handle(), INFINITE );
	ReaderWriterLock::ScopedReaderLock piecewisePath( this->historyLock_ );
	return this->piecewisePaths_;
}

size_t MainFrame::piecewisePathsFrameRate() const
{
	ReaderWriterLock::ScopedReaderLock piecewisePath( this->historyLock_ );
	return this->piecewisePathsFrameRate_;
}

std::deque<Simulation::State> MainFrame::stateHistory() const
{
	WaitForSingleObject( historyComplete_.handle(), INFINITE );
	ReaderWriterLock::ScopedReaderLock piecewisePath( this->historyLock_ );

	ConstSimulationTreePtr tree = this->getTree();
	if( !tree )
		return std::deque<Simulation::State>();

	assert( this->stateHistory_.empty() || this->piecewisePaths_.empty() );

	if( !this->stateHistory_.empty() )
	{
		ReaderWriterLock::ScopedReaderLock historyLock( this->historyLock_ );
		std::deque<Simulation::State> result;
		std::transform( this->stateHistory_.begin(), this->stateHistory_.end(),
			std::back_inserter(result), std::select1st<MainFrame::StateWithNode>() );

		return result;
	}
	else
	{
		assert( false );
		return std::deque<Simulation::State>();
	}
}

void MainFrame::addObject( SceneGraphElementPtr object )
{
	addObjects( std::deque<SceneGraphElementPtr>(1, object) );
}

void MainFrame::addCamera( CameraWrapperPtr camera )
{
	ActionPtr addCamera(
		new CreateCameraAction( 
			camera,
            scene_,
			this) );
	performAction( addCamera );
	std::for_each( canvases_.begin(), canvases_.end(), 
		boost::bind( &GLView::clearCache, _1 ) );
}

std::deque<SceneGraphElementPtr> allChildren( const std::deque<SceneGraphElementPtr>& objects )
{
	std::deque<SceneGraphElementPtr> result;
	{
		std::deque<SceneGraphElementPtr> queue( objects );
		while( !queue.empty() )
		{
			SceneGraphElementPtr current = queue.front();
			queue.pop_front();
			result.push_back( current );
			std::copy( current->children_begin(), current->children_end(), 
				std::front_inserter(queue) );
		}
	}

	return result;
}

void MainFrame::addObjects( const std::deque<SceneGraphElementPtr>& objects )
{
	SceneGraphElementPtr root = scene_->root();

	std::deque<ActionPtr> actions;
	{
		// need to tease apart the whole hierarchy
		std::deque<SceneGraphElementPtr> allObjects = allChildren( objects );
		
		ActionPtr addObject( 
			new CreateObjectAction(allObjects, scene_) );
		actions.push_back( addObject );

		for( std::deque<SceneGraphElementPtr>::const_iterator objectIter = objects.begin();
			objectIter != objects.end(); ++objectIter )
		{
			ActionPtr reparentAction( new ReparentAction(*objectIter, root) );
			actions.push_back( reparentAction );
		}
	}

	ActionPtr changeSelection( new ChangeSelectionAction(this, this->getSelected(), objects,
				this->selectedConstraint(), ConstraintPtr() ) );
	actions.push_back( changeSelection );

	boost::shared_ptr<Action> compoundAction( new CompoundAction(actions) );
	this->performAction( compoundAction );

	std::for_each( canvases_.begin(), canvases_.end(), 
		boost::bind( &GLView::clearCache, _1 ) );
}

void MainFrame::deleteObjects( const std::deque<SceneGraphElementPtr>& objects )
{
	std::deque<ActionPtr> actions;

	ActionPtr changeSelection( new ChangeSelectionAction(this, this->getSelected(), 
		std::deque<SceneGraphElementPtr>(), this->selectedConstraint(), ConstraintPtr() ) );
	actions.push_back( changeSelection );

	// @todo remove all joints involved with these objects

	{
		for( std::deque<SceneGraphElementPtr>::const_iterator objectIter = objects.begin();
			objectIter != objects.end(); ++objectIter )
		{
			ActionPtr reparentAction( new ReparentAction(*objectIter, SceneGraphElementPtr()) );
			actions.push_back( reparentAction );
		}

		std::deque<SceneGraphElementPtr> allObjects = allChildren( objects );

		ActionPtr deleteObjects(
			new DeleteObjectAction( allObjects, scene_ ) );
		actions.push_back( deleteObjects );
	}

	boost::shared_ptr<Action> compoundAction( new CompoundAction(actions) );
	this->performAction( compoundAction );

	std::for_each( canvases_.begin(), canvases_.end(), 
		boost::bind( &GLView::clearCache, _1 ) );
}

#ifdef SOUND
void MainFrame::OnSoundComputeAllModes( wxCommandEvent& event )
{
	assert( !computingModes_ );
	std::deque<PhysicsObjectPtr> physicsObjects;

	std::deque<SceneGraphElementPtr> toHandle( 1, scene_->root() );
	while( !toHandle.empty() )
	{
		SceneGraphElementPtr current = toHandle.front();
		toHandle.pop_front();

		PhysicsObjectPtr po = boost::dynamic_pointer_cast<PhysicsObject>( current );
		if( po )
			physicsObjects.push_back( po );

		if( !boost::dynamic_pointer_cast<CombinedObject>(current) )
		{
			std::copy( current->children_begin(), current->children_end(), 
				std::front_inserter(toHandle) );
		}
	}

	computingModes_ = true;
	stopComputingModes_.reset();
	computeModesThread_.reset( new boost::thread(
		boost::bind(&MainFrame::findModes, this, physicsObjects, true) ) );

	this->updateMenus();
}

void MainFrame::OnSoundComputeObjectModes(wxCommandEvent& event)
{
	assert( !computingModes_ );
	std::deque<PhysicsObjectPtr> physicsObjects;

	std::deque<SceneGraphElementPtr> selected = this->getSelected();
	for( std::deque<SceneGraphElementPtr>::const_iterator selectedItr = selected.begin();
		selectedItr != selected.end(); ++selectedItr )
	{
		boost::shared_ptr<PhysicsObject> po = 
			boost::dynamic_pointer_cast<PhysicsObject>( *selectedItr );
		if( po )
		{
			SceneGraphElementPtr current = po;
			while( current )
			{
				boost::shared_ptr<CombinedObject> combined = 
					boost::dynamic_pointer_cast<CombinedObject>( current );
				if( combined )
					po = combined;

				current = current->parent().lock();
			}

			physicsObjects.push_back( po );
		}
	}

	std::sort( physicsObjects.begin(), physicsObjects.end() );
	physicsObjects.erase( 
		std::unique(physicsObjects.begin(), physicsObjects.end()), 
		physicsObjects.end() );

	computingModes_ = true;
	stopComputingModes_.reset();
	computeModesThread_.reset( new boost::thread(
		boost::bind(&MainFrame::findModes, this, physicsObjects, true) ) );

	this->updateMenus();
}

void MainFrame::OnFinishedComputingModes( wxCommandEvent& event )
{
	computingModes_ = false;
	this->updateMenus();
}

void MainFrame::findModes( 
	const std::deque<PhysicsObjectPtr>& physicsObjects, bool compute )
{
	// note there's a threading bug here in that it is possible to 
	// change the physicsObjects in question even while we're computing modes which 
	// could lead to incorrect modes in come cases.  The correct way to fix this is
	// to dup the objects in question before we do this.  Slot this as a "todo."
	for( std::deque<PhysicsObjectPtr>::const_iterator itr = physicsObjects.begin();
		itr != physicsObjects.end(); ++itr )
	{
		AudioHash hash = hashForAudio( *itr );

		/*
		{
			// debugging the hashing.
			std::ofstream ofs ( ((*itr)->name() + ".audioHash").c_str() );
			ofs << hash.first;
			ofs << "\nShear:";
			for( size_t i = 0; i < hash.second.size(); ++i )
				ofs << " " << hash.second[i];
			ofs << "\n";
		}
		*/

		{
			ReaderWriterLock::ScopedReaderLock lock( descriptionToSquashingCubesModesLock_ );
			DescriptionToModesMap::const_iterator mapItr = 
				descriptionToSquashingCubesModes_.find( hash );
			if( mapItr != descriptionToSquashingCubesModes_.end() )
			{
				(*itr)->setSquashingCubesModes( mapItr->second );
				continue;
			}
		}

		if( !compute )
			continue;

		if( WaitForSingleObject(stopComputingModes_.handle(), 0) == WAIT_OBJECT_0 )
			break;

		try
		{
			boost::shared_ptr<SquashingCubesModes> modes = 
				planning::computeModes( *itr, this->stopComputingModes_ );

			if( modes )
			{
				ReaderWriterLock::ScopedWriterLock lock( descriptionToSquashingCubesModesLock_ );
				descriptionToSquashingCubesModes_.insert( 
					DescriptionToModesMap::value_type(hash, modes) );
			}
		}
		catch( ParserException& e )
		{
			std::ostringstream oss;
			oss << e;
			this->logError( oss.str() );
		}
		catch( Exception& e )
		{
			this->logError( e.message() );
		}
	}

	// post event saying we are finished
	{
		wxCommandEvent event(wxEVT_FINISHED_COMPUTING_MODES, this->GetId());
		event.SetEventObject( this );
		::wxPostEvent( this, event );
	}
}
#endif // SOUND

// it is only safe to post errors to the log window in the main message handling loop,
//   so if we want to log them from somewhere else we will need to create a message
//   and pass it along
void MainFrame::logError( const std::string& error )
{
	wxCommandEvent event(wxEVT_LOG_ERROR, this->GetId());
	event.SetEventObject( this );
	event.SetString( wxT(error.c_str()) );
	::wxPostEvent( this, event );
}

#ifdef SOUND
void MainFrame::OnSoundStopComputingModes(wxCommandEvent& event)
{
	if( !computingModes_ )
	{
		updateMenus();
		return;
	}

	assert( computeModesThread_ );
	stopComputingModes_.notify_all();
	computeModesThread_->join();
	computingModes_ = false;

	this->updateMenus();
}
#endif // SOUND

void MainFrame::deleteConstraint( ConstraintPtr constraint )
{
	if( !constraint || !this->getTree() )
		return;

	this->log( "Deleted constraint" );
	std::deque<ActionPtr> actions;

	ActionPtr changeSelection( new ChangeSelectionAction(this, this->getSelected(), 
		this->getSelected(), this->selectedConstraint(), ConstraintPtr() ) );
	actions.push_back( changeSelection );

	ActionPtr deleteConstraint( new RemoveConstraintAction(constraint, this->getTree()) );
	actions.push_back( deleteConstraint );

	boost::shared_ptr<Action> compoundAction( new CompoundAction(actions) );
	this->performAction( compoundAction );
}

void MainFrame::loadAny( const std::string& filename )
{
	std::string::size_type dotPos = filename.rfind(".");
	std::string extension;
	if( dotPos != std::string::npos )
	{
		extension = filename.substr( dotPos );
	}

	if( extension == ".tree" )
	{
		loadTree( filename );
	}
	else
	{
		load( filename );
#ifdef SOUND
		loadModesFile( filename );
#endif
	}
}

void MainFrame::load( const std::string& filename )
{
	try
	{
		changeDirectory( filename );
		boost::shared_ptr<Scene> scene = 
			parseSimulation( filename );

		changeScene( scene );
		sceneFilename_ = filename;
		this->log( "Loaded scene: " + sceneFilename_ );
	}
	catch( ParserException& e )
	{
		std::ostringstream oss;
		oss << "Error loading scene: \n" << e;
		this->errMsg( oss.str() );
	}
	catch( IOException& e )
	{
		std::ostringstream oss;
		oss << "Error loading scene: \n" << e.message();
		this->errMsg( oss.str() );
	}
}

#ifdef SOUND
void MainFrame::loadModesFile( const std::string origFilename )
{
	boost::filesystem::path p( origFilename );
	boost::filesystem::path soundDataPath( p.branch_path() / (boost::filesystem::basename(p) + ".dat") );
	std::ifstream ifs( soundDataPath.native_file_string().c_str(), std::ios::binary );
	if( !ifs )
		return;

	ifs.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );

	std::ifstream::pos_type len = filelen( ifs );
	try
	{
		while( ifs.tellg() < len )
		{
			AudioHash desc;
			readArray( desc.first, ifs );
			ifs.read( reinterpret_cast<char*>(&desc.second[0]), 
				sizeof(float) * 6 );

			BoundingBox3f bounds;
			ifs.read( reinterpret_cast<char*>( &bounds ), sizeof( BoundingBox3f ) );

			boost::shared_ptr<SquashingCubesModes> modes( new SquashingCubesModes(
				BoundedVoxelGrid( VoxelGrid(ifs), bounds ) ) );

			{
				std::vector<SmallIndexTriple> vertices;
				readArray( vertices, ifs );
				modes->vertices.resize( vertices.size() );
				std::copy( vertices.begin(), vertices.end(), modes->vertices.begin() );
			}

			readArray( modes->frequencies, ifs );
			readMatrix( modes->modes, ifs );

			// todo do some sanity checking here
			ReaderWriterLock::ScopedWriterLock lock( this->descriptionToSquashingCubesModesLock_ );
			this->descriptionToSquashingCubesModes_.insert(
				DescriptionToModesMap::value_type( desc, modes ) );
		}
	}
	catch( std::ifstream::failure& e )
	{
		std::ostringstream oss;
		oss << "Failure reading from file '" << soundDataPath.native_file_string() << "': " << e.what();
		logError( oss.str() );
		return;
	}

	std::deque< boost::shared_ptr<PhysicsObject> > physicsObjects;
	std::deque<SceneGraphElementPtr> toHandle( 1, scene_->root() );
	while( !toHandle.empty() )
	{
		SceneGraphElementPtr current = toHandle.front();
		toHandle.pop_front();

		PhysicsObjectPtr po = boost::dynamic_pointer_cast<PhysicsObject>( current );
		if( po )
			physicsObjects.push_back( po );

		if( !boost::dynamic_pointer_cast<CombinedObject>(current) )
		{
			std::copy( current->children_begin(), current->children_end(), 
				std::front_inserter(toHandle) );
		}
	}

	this->findModes( physicsObjects, false );
}
#endif // SOUND

bool MainFrame::wireframeOnShaded() const
{
	return renderMenu_->IsChecked(MENU_RENDER_WIREFRAME_ON_SHADED);
}

bool MainFrame::visualizeInertiaTensor() const
{
	return renderMenu_->IsChecked(MENU_RENDER_INERTIA_TENSOR);
}

bool MainFrame::visualizeConvexHulls() const
{
	return renderMenu_->IsChecked(MENU_RENDER_CONVEX_HULLS);
}

bool MainFrame::visualizeContactPoints() const
{
	return renderMenu_->IsChecked(MENU_RENDER_CONTACT_POINTS);
}

bool MainFrame::showAllPaths() const
{
	return renderMenu_->IsChecked(MENU_RENDER_SHOW_ALL_PATHS);
}

bool MainFrame::hasShadows() const
{
	return renderMenu_->IsChecked(MENU_RENDER_SHADOWS);
}

bool MainFrame::renderGrid() const
{
	return renderMenu_->IsChecked(MENU_RENDER_GRID);
}

void MainFrame::OnMenuSetSampling(wxCommandEvent& event)
{
	if( simulateMenu_->IsChecked( MainFrame::MENU_SIMULATE_SAMPLING ) )
		this->pauseSampling_ = false;
	else
		this->pauseSampling_ = true;

	std::for_each( samplingEvents_.begin(), samplingEvents_.end(), boost::bind( &Event::notify_all, _1 ) );
}


void MainFrame::OnMenuSimulateRun(wxCommandEvent& event)
{
	if( simulateMenu_->IsChecked( MainFrame::MENU_SIMULATE_RUN ) )
		this->startSimulation();
	else
		this->stopSimulation();
}

void MainFrame::stopSimulating()
{
	simulateMenu_->Check( MainFrame::MENU_SIMULATE_RUN, false );
	this->stopSimulation();
}

MainFrame::RenderStyle MainFrame::renderStyle() const
{
	if( renderGeometryMenu_->IsChecked(MENU_RENDER_GEOMETRY_POINTS) )
		return MainFrame::STYLE_POINTS;
	else if( renderGeometryMenu_->IsChecked(MENU_RENDER_GEOMETRY_WIREFRAME) )
		return MainFrame::STYLE_WIREFRAME;
	else if( renderGeometryMenu_->IsChecked(MENU_RENDER_GEOMETRY_BOXES) )
		return MainFrame::STYLE_BOXES;
	else
		return MainFrame::STYLE_POLYS;
}

vl::Vec3f MainFrame::backgroundColor() const
{
	switch( backgroundStyle() )
	{
	case MainFrame::BACKGROUND_WHITE:
		return vl::Vec3f( 1.0, 1.0, 1.0 );
	case MainFrame::BACKGROUND_LIGHT_GRAY:
		return vl::Vec3f(60./255., 60./255., 60./255. );
	case MainFrame::BACKGROUND_DARK_GRAY:
		return vl::Vec3f(30./255., 30./255., 30./255. );
	case MainFrame::BACKGROUND_BLACK:
		return vl::Vec3f( 0.0, 0.0, 0.0 );
	}

	assert( false );
	return vl::Vec3f( 0.0, 0.0, 0.0 ); 
}

MainFrame::BackgroundStyle MainFrame::backgroundStyle() const
{
	if( renderBackgroundMenu_->IsChecked(MENU_RENDER_BACKGROUND_WHITE) )
		return MainFrame::BACKGROUND_WHITE;
	if( renderBackgroundMenu_->IsChecked(MENU_RENDER_BACKGROUND_LIGHT_GRAY) )
		return MainFrame::BACKGROUND_LIGHT_GRAY;
	if( renderBackgroundMenu_->IsChecked(MENU_RENDER_BACKGROUND_DARK_GRAY) )
		return MainFrame::BACKGROUND_DARK_GRAY;
	else
		return MainFrame::BACKGROUND_BLACK;
}

void MainFrame::setBackgroundStyle( MainFrame::BackgroundStyle style )
{
	switch( style )
	{
	case MainFrame::BACKGROUND_WHITE:
		renderBackgroundMenu_->Check( MENU_RENDER_BACKGROUND_WHITE, TRUE );
		break;
	case MainFrame::BACKGROUND_LIGHT_GRAY:
		renderBackgroundMenu_->Check( MENU_RENDER_BACKGROUND_LIGHT_GRAY, TRUE );
		break;
	case MainFrame::BACKGROUND_DARK_GRAY:
		renderBackgroundMenu_->Check( MENU_RENDER_BACKGROUND_DARK_GRAY, TRUE );
		break;
	case MainFrame::BACKGROUND_BLACK:
		renderBackgroundMenu_->Check( MENU_RENDER_BACKGROUND_BLACK, TRUE );
		break;
	};
}

void MainFrame::OnExit(wxCommandEvent& event)
{
    Close();
}

void MainFrame::OnClose(wxCloseEvent& event)
{
	if(checkDirty() || !event.CanVeto())
	{
		this->preDestroy();
		Destroy();
	}
}

void MainFrame::preDestroy()
{
	computeHistoryDone_ = true;
	this->computeHistoryPendingCondition_.notify_all();
	this->computeHistoryThread_->join();

	// need to stop all sampling threads first, otherwise we'll see
	//   access violations:
	std::for_each( trees_.begin(), trees_.end(), 
        boost::bind(&SimulationTree::setTimeRange, _1, std::make_pair(0.0f, 0.0f)) );

	this->sampling_ = false;

	// tell sampling threads to stop
	std::for_each( samplingEvents_.begin(), samplingEvents_.end(), boost::bind( &Event::notify_all, _1 ) );
	// wait for sampling threads
	WaitForMultipleObjects(this->samplingThreadHandles_.size(), &samplingThreadHandles_[0], TRUE, INFINITE);

	// close all thread handles
	std::for_each( samplingThreadHandles_.begin(), samplingThreadHandles_.end(), CloseHandle );
}

void MainFrame::OnSetRenderShowAllPaths(wxCommandEvent& event)
{
	SimulationTreePtr tree = this->getTree();
	if( tree )
	{
		this->setSamples( tree, tree->activePaths() );
		bestScoresFrame_->updateShowAllPaths();
	}

	this->updateViews();
}

void MainFrame::OnSetRenderOption(wxCommandEvent& event)
{
	renderMenu_->Enable(MENU_RENDER_WIREFRAME_ON_SHADED, 
		renderGeometryMenu_->IsChecked(MENU_RENDER_GEOMETRY_POLYS) );
	this->updateViews();
}

void MainFrame::OnNeedsRedraw(wxCommandEvent& event)
{
	std::for_each( canvases_.begin(), canvases_.end(), 
		boost::bind( &GLView::clearCache, _1 ) );
}

void MainFrame::OnRepaintAllViews(wxCommandEvent& event)
{
	// send Update message to all those that are visible
	if( this->viewMaximized() )
		this->currentCanvas_ ->Update();
	else
		std::for_each( canvases_.begin(), canvases_.end(), 
			boost::bind( &GLView::Update, _1 ) );
}

void MainFrame::changeScene( boost::shared_ptr<Scene> newScene )
{
	stopSimulating();
	this->setMode( MODE_EDIT );
	this->clearTrees();

	undoQueue_.clear();
	redoQueue_.clear();
	dirty_ = false;

	updateViews();
	timeView_->Refresh( TRUE );
	sceneFilename_.clear();
	this->setTime( 0.0f, true );

	this->setSelected( std::deque<SceneGraphElementPtr>() );

	scene_ = newScene;

	this->hierarchyFrame_->setScene( scene_, selected_ );

	// @todo do something intelligent to get different views here
	std::deque<CameraWrapperPtr> cameras( 
		scene_->begin_cameras(), scene_->end_cameras() );

	for( size_t iCanvas = 0; iCanvas < this->canvases_.size(); ++iCanvas )
	{
		canvases_.at(iCanvas)->setCamera( cameras.at( iCanvas % cameras.size() ) );
	}

	{
		BoundingBox3f bounds = scene_->bounds();
		float maxAxis = 0.0;
		for( size_t i = 0; i < 3; ++i )
			maxAxis = std::max( maxAxis, bounds.maximum()[i] - bounds.minimum()[i] );
//		if( maxAxis > 1e-4 )
//			scene_->camera()->camera()->setScale( 2.0 / maxAxis );

		if( maxAxis > 1e-4 )
			scene_->camera()->camera()->setTranslationalSensitivity( 0.1*maxAxis );
	}

	this->updateCameraMenu();

	if( arePropertiesShown() )
		this->propertyDialog_->updateValues();

	timeView_->Refresh();
	this->updateViews();

	updateMenus();
}

bool MainFrame::save()
{
	if( this->sceneFilename_.empty() )
		return saveAs();
	else
		return save( sceneFilename_ );
}

bool MainFrame::save( const std::string& filename )
{
	try
	{
		{
			MechModelObjectPtr mechModel = scene_->toMechModel();
			std::ofstream ofs( filename.c_str() );
			if( !ofs )
			{
				errMsg( "Unable to open file '" + filename + "' for writing." );
				return false;
			}

			mechModel->dump( ofs, std::string() );
			dirty_ = false;
			this->sceneFilename_ = filename;
		}

#ifdef SOUND
		{
			boost::filesystem::path p( filename );
			boost::filesystem::path soundDataPath( p.branch_path() / (boost::filesystem::basename(p) + ".dat") );
			std::ofstream ofs( soundDataPath.native_file_string().c_str(), std::ios::binary );
			if( ofs )
			{
				ReaderWriterLock::ScopedReaderLock lock( this->descriptionToSquashingCubesModesLock_ );
				for( DescriptionToModesMap::const_iterator itr = descriptionToSquashingCubesModes_.begin();
					itr != descriptionToSquashingCubesModes_.end(); ++itr )
				{
					dumpArray( itr->first.first, ofs );
					ofs.write( reinterpret_cast<const char*>(&itr->first.second[0]), 
						sizeof(float) * 6 );
					std::vector<SmallIndexTriple> vertices( 
						itr->second->vertices.begin(), itr->second->vertices.end() );
					BoundingBox3f bounds = itr->second->voxels.bounds();
					ofs.write( reinterpret_cast<const char*>( &bounds ), sizeof( BoundingBox3f ) );
					itr->second->voxels.voxelGrid().dump( ofs );

					dumpArray( vertices, ofs );
					dumpArray( itr->second->frequencies, ofs );
					dumpMatrix( itr->second->modes, ofs );
				}
			}
		}
#endif // SOUND

		this->log( "Saved scene: " + sceneFilename_ );

		return true;
	}
	catch( IOException& e )
	{
		errMsg( "Unable to save to file '" + filename + "': " + e.message() );
		return false;
	}
	catch( boost::filesystem::filesystem_error& )
	{
		errMsg( "Unable to save to file '" + filename + "': path error" );
		return false;
	}
}

bool MainFrame::saveAs()
{
	do
	{
		boost::shared_ptr<wxFileDialog> dlg(
			new wxFileDialog(
				this, 
				"Save scene",
				"", 
				"", 
				".plan files (*.plan)|*.plan|All files(*.*)|*.*",
				wxFD_SAVE | wxFD_OVERWRITE_PROMPT | wxFD_CHANGE_DIR, 
				wxDefaultPosition), 
			boost::mem_fn(&wxFileDialog::Destroy) );

		if( dlg->ShowModal() != wxID_OK )
			return false;

		std::string filename(dlg->GetPath().GetData());
		std::string::size_type dotPos = filename.rfind(".");
		if( dotPos == std::string::npos )
			filename += ".plan";

		if( save(filename) )
			return true;
	} while( 1 );
}

bool MainFrame::checkDirty()
{
	if( !dirty_ )
		return true;

	boost::shared_ptr<wxMessageDialog> dlg( new wxMessageDialog(this,
		"Scene has changed since last save.  Save changes made?", 
		"Save changes?",
		wxYES_NO | wxCANCEL ),
		boost::mem_fn(&wxMessageDialog::Destroy) );
	switch( dlg.get()->ShowModal() )
	{
	case wxID_CANCEL:
		return false;
	case wxID_YES:
		return this->save();
	case wxID_NO:
		return true;
	default:
		assert( false );
		return false;
	}
}

void MainFrame::OnMenuFileNew(wxCommandEvent& event)
{
	if( checkDirty() )
	{
		changeScene( boost::shared_ptr<Scene>( new Scene ) );
		this->log( "New scene." );
	}
}

void MainFrame::OnMenuFileImport(wxCommandEvent& event)
{
	boost::shared_ptr<wxFileDialog> dlg( 
		new wxFileDialog(
			this, 
			"Import a plan file",
			"", 
			"", 
			"Plan file(*.plan)|*.plan|All files(*.*)|*.*",
			wxFD_OPEN),
		boost::mem_fn(&wxFileDialog::Destroy) );

	if ( dlg->ShowModal() == wxID_OK )
	{
		try
		{
			std::string filename( dlg->GetPath().GetData() );
			changeDirectory( filename );
			boost::shared_ptr<Scene> scene = 
				parseSimulation( filename );

			// don't copy the cameras for now

			std::deque<ActionPtr> actions;

			{
				STLEXT hash_set<std::string> usedNames;
				// it's okay to reuse materials in the new scene because
				//   we are just going to be toasting it anyway.  Need to
				//   fix the names though:
				for( size_t iMaterial = 0; iMaterial < scene->materialCount(); ++iMaterial )
				{
					MaterialPtr mat = scene->material(iMaterial);
					if( scene_->material( mat->name() ) || usedNames.find( mat->name() ) != usedNames.end() )
					{
						std::string prefix = mat->name();
						// need to rename
						for( size_t i = 0; 
							scene_->material(mat->name()) || 
							(usedNames.find( mat->name() ) != usedNames.end()); ++i )
						{
							std::string newName = prefix;
							newName.push_back( '_' );
							newName.append( boost::lexical_cast<std::string>( i ) );
							mat->setName(newName);
						}
					}

					usedNames.insert( mat->name() );

					ActionPtr createMat( new CreateMaterialAction(mat, this->scene_) );
					actions.push_back( createMat );
				}
			}

			{
				STLEXT hash_set< std::string > usedNames;
				// same for the scene objects, just need to run through them all and 
				//   check that they do not have duplicate names:
				for( Scene::object_iter obItr = scene->begin_objects(); 
					obItr != scene->end_objects(); ++obItr )
				{
					if( *obItr == scene->root() )
						continue;

					if( usedNames.find( (*obItr)->name() ) != usedNames.end() || 
						this->scene_->object( (*obItr)->name() ) )
					{
						std::string prefix = (*obItr)->name();
						// need to rename
						for( size_t i = 0; 
							scene_->object((*obItr)->name()) || 
							(usedNames.find( (*obItr)->name() ) != usedNames.end()); ++i )
						{
							std::string newName = prefix;
							newName.push_back( '_' );
							newName.append( boost::lexical_cast<std::string>( i ) );
							(*obItr)->setName(newName);
						}
					}

					usedNames.insert( (*obItr)->name() );

					ActionPtr createOb( new CreateObjectAction(*obItr, this->scene_) );
					actions.push_back( createOb );
				}

				SceneGraphElementPtr root = scene->root();
				for( SceneGraphElement::child_iterator childItr = root->children_begin();
					childItr != root->children_end(); ++childItr )
				{
					ActionPtr action( new ReparentAction( *childItr, scene_->root() ) );
					actions.push_back( action );
				}
			}

			ActionPtr action( new CompoundAction(actions) );
			this->performAction( action );
		}
		catch( ParserException& e )
		{
			std::ostringstream oss;
			oss << "Error loading scene: \n" << e;
			this->errMsg( oss.str() );
		}
		catch( IOException& e )
		{
			std::ostringstream oss;
			oss << "Error loading scene: \n" << e.message();
			this->errMsg( oss.str() );
		}

		/*
		try
		{
			m_canvas->loadSkins( std::string(dlg->GetPath().GetData()) );
			updateMenus();
		}
		catch( IOException& e )
		{
			this->errMsg( e.message() );
		}
		*/
	}
}

void MainFrame::OnMenuFileLoadDebuggingInfo(wxCommandEvent& event)
{
	boost::shared_ptr<wxFileDialog> dlg( 
		new wxFileDialog(
			this, 
			"Open debugging info",
			"", 
			"", 
			"Debugging info(*.bin)|*.bin|All files(*.*)|*.*",
			wxFD_OPEN),
		boost::mem_fn(&wxFileDialog::Destroy) );
	if ( dlg->ShowModal() == wxID_OK )
	{
        std::ifstream ifs( dlg->GetPath().GetData(), std::ios::binary );
        float time;
        ifs.read( reinterpret_cast<char*>( &time ), sizeof( float ) );
        std::vector<RigidDynamicState> dynamicStates;
        readArray( dynamicStates, ifs );
        Simulation::State state( dynamicStates, time );

	    std::ifstream::pos_type len = filelen( ifs );
        std::vector< Simulation::ContactPoint > contactPoints;
		while( ifs.tellg() < len )
        {
            vl::Vec3f pos;
            ifs.read( reinterpret_cast<char*>( pos.Ref() ), sizeof( vl::Vec3f ) );

            vl::Vec3f normal;
            ifs.read( reinterpret_cast<char*>( normal.Ref() ), sizeof( vl::Vec3f ) );

            contactPoints.push_back( Simulation::ContactPoint(pos, normal, -1, -1) );
        }

        this->simulation_.reset( new DummySimulation(this->scene_, state, contactPoints) );
    }
}

void MainFrame::OnMenuFileOpen(wxCommandEvent& event)
{
	if( !checkDirty() )
		return;

	boost::shared_ptr<wxFileDialog> dlg( 
		new wxFileDialog(
			this, 
			"Open a plan file",
			"", 
			"", 
			"Plan file(*.plan)|*.plan|All files(*.*)|*.*",
			wxFD_OPEN),
		boost::mem_fn(&wxFileDialog::Destroy) );
	if ( dlg->ShowModal() == wxID_OK )
	{
		load( dlg->GetPath().GetData() );

#ifdef SOUND
		{
			ReaderWriterLock::ScopedWriterLock lock( this->descriptionToSquashingCubesModesLock_ );
			this->descriptionToSquashingCubesModes_.clear();
		}
		loadModesFile( dlg->GetPath().GetData() );
#endif // SOUND

		/*
		try
		{
			m_canvas->loadSkins( std::string(dlg->GetPath().GetData()) );
			updateMenus();
		}
		catch( IOException& e )
		{
			this->errMsg( e.message() );
		}
		*/
	}
}

void MainFrame::OnMenuFileSave(wxCommandEvent& event)
{
	save();
}

void MainFrame::OnMenuFileSaveAs(wxCommandEvent& event)
{
	saveAs();
}

void writePiecewisePath( const PiecewisePath& value, std::ofstream& ofs )
{
	std::vector<int> positionTimes;
	std::vector<float> positionCoeffs;
	std::vector<int> rotationTimes;
	std::vector<float> rotations;
	std::vector<float> rotationCoeffs;

	value.dumpToArrays( 
		positionTimes,
		positionCoeffs,
		rotationTimes,
		rotations,
		rotationCoeffs );

	dumpArray( positionTimes, ofs );
	dumpArray( positionCoeffs, ofs );
	dumpArray( rotationTimes, ofs );
	dumpArray( rotations, ofs );
	dumpArray( rotationCoeffs, ofs );
}

PiecewisePath readPiecewisePath( std::ifstream& ifs )
{
	std::vector<int> positionTimes;
	std::vector<float> positionCoeffs;
	std::vector<int> rotationTimes;
	std::vector<float> rotations;
	std::vector<float> rotationCoeffs;

	readArray( positionTimes, ifs );
	readArray( positionCoeffs, ifs );
	readArray( rotationTimes, ifs );
	readArray( rotations, ifs );
	readArray( rotationCoeffs, ifs );

	return PiecewisePath(
		positionTimes,
		positionCoeffs,
		rotationTimes,
		rotations,
		rotationCoeffs );
}


ScenePtr readScene( std::ifstream& ifs )
{
	std::string str;
	readArray( str, ifs );
	if( str.empty() )
		return ScenePtr();

	// @todo fix this plz!!
#ifdef WIN32
	const char* tmpFile = "c:\\temp\\tmpScene.txt";
#else
	const char* tmpFile = "/tmp/tmpFile.txt";
#endif

	{
		std::ofstream ofs( tmpFile );
		assert( ofs );
		ofs << str;
	}
	return parseSimulation( tmpFile );
}

void writeScene( boost::shared_ptr<const Scene> scene, std::ofstream& ofs )
{
	std::string str;
	if( scene )
	{
		MechModelObjectPtr mechModel = scene->toMechModel();
		std::ostringstream oss;
		mechModel->dump( oss, std::string() );
		str = oss.str();
	}
	dumpArray( str, ofs );
}

template <typename T>
T readPrimitive( std::ifstream& ofs )
{
	T result;
	ofs.read( reinterpret_cast<char*>( &result ), sizeof(T) );
	return result;
}


template <typename T>
void writePrimitive( T value, std::ofstream& ofs )
{
	ofs.write( reinterpret_cast<const char*>( &value ), sizeof(T) );
}

void MainFrame::OnMenuFileLoadTree(wxCommandEvent& event)
{
	if( !checkDirty() )
		return;

	boost::shared_ptr<wxFileDialog> dlg( 
		new wxFileDialog(
			this, 
			"Load tree",
			"", 
			"", 
			"Tree files(*.tree)|*.tree|All files(*.*)|*.*",
			wxFD_OPEN, 
			wxDefaultPosition),
		boost::mem_fn(&wxFileDialog::Destroy) );

	if ( dlg.get()->ShowModal() == wxID_OK )
	{
		std::string filename = dlg.get()->GetPath().GetData();
		this->loadTree( filename );
	}
}

void MainFrame::OnMenuFileLoadSamples(wxCommandEvent& event)
{
	SimulationTreePtr tree = this->getTree();
	if( !tree )
		return;

	boost::shared_ptr<wxFileDialog> dlg( 
		new wxFileDialog(
			this, 
			"Load samples",
			"", 
			"", 
			"Samples files(*.samples)|*.samples|All files(*.*)|*.*",
			wxFD_OPEN | wxFD_MULTIPLE, 
			wxDefaultPosition),
		boost::mem_fn(&wxFileDialog::Destroy) );

	if ( dlg.get()->ShowModal() == wxID_OK )
	{
		try
		{
			wxArrayString paths;
			dlg.get()->GetPaths(paths);
			for( size_t i = 0; i < paths.Count(); ++i )
			{
				std::string filename( paths[i].GetData() );
				this->loadSamples( filename );
			}
		}
		catch( Exception& e )
		{
			errMsg( std::string("Error loading samples: ") + e.message() );
		}
	}
}

void MainFrame::loadSamples( const std::string& filename )
{
	SimulationTreePtr tree = this->getTree();
	if( !tree )
		throw IOException( "No tree currently loaded" );
	
	tree->read( filename );
}

void MainFrame::loadTree( const std::string& filename )
{
	changeDirectory( filename );
	std::ifstream ifs( filename.c_str(), std::ios::binary );
	ifs.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );

	if( !ifs )
	{
		errMsg( "Unable to open file '" + filename + "' for reading" );
		return;
	}

	try
	{
		size_t samplingRate = readPrimitive<size_t>( ifs );
		ScenePtr scene = readScene( ifs );
		ScenePtr parentScene = readScene( ifs );

		std::vector<size_t> fixedObjects;
		readArray( fixedObjects, ifs );

		float startTime = readPrimitive<float>( ifs );
		std::vector<PiecewisePath> fixedObjectsHistory;
		{
			size_t numPaths = readPrimitive<size_t>( ifs );
			for( size_t i = 0; i < numPaths; ++i )
				fixedObjectsHistory.push_back( readPiecewisePath(ifs) );
		}

		SimulationTreePtr tree;
		if( parentScene )
		{
			SimulationTreePtr parentTree( new SimulationTree("tree", parentScene, samplingRate, simulationFactory()) );
			std::sort( fixedObjects.begin(), fixedObjects.end() );
			std::vector<ConstPhysicsObjectPtr> dynamicObjects = 
				parentTree->dynamicObjects();
			std::deque<ConstPhysicsObjectPtr> activeObjectPtrs;
			for( size_t i = 0; i < dynamicObjects.size(); ++i )
			{
				if( !std::binary_search(fixedObjects.begin(), fixedObjects.end(), i) )
					activeObjectPtrs.push_back( dynamicObjects.at(i) );
			}

			tree.reset( new SimulationTree(
				"childTree",
				scene,
				parentTree,
				activeObjectPtrs,
				fixedObjectsHistory,
				samplingRate,
				startTime ) );
		}
		else
		{
			assert( fixedObjects.empty() );
			tree.reset( new SimulationTree(
				"tree",
				scene,
				samplingRate,
				simulationFactory() ) );
		}

		std::vector<ConstPhysicsObjectPtr> dynamicObjects = 
			tree->dynamicObjects();
		{
			size_t numPaths = readPrimitive<size_t>(ifs);
			for( size_t i = 0; i < numPaths; ++i )
			{
				Path p( ifs, dynamicObjects.size(), tree->samplingRate(), filename );
				tree->addPath( p );
			}
		}

        this->pauseSampling_ = true;
		changeScene( ScenePtr( new Scene( *tree->scene() ) ) );
		this->addTree( tree );
		dirtyObjects_.clear();
		this->setMode( MODE_SOLVE );
		this->setPath( tree, 0 );

		this->log( "Loaded tree: " + filename );
	}
	catch( std::ifstream::failure& e )
	{
		std::ostringstream oss;
		oss << "Failure reading from file '" << filename << "': " << e.what();
		errMsg( oss.str() );
		return;
	}
}

void MainFrame::OnMenuFileSaveTree(wxCommandEvent& event)
{
	boost::shared_ptr<const SimulationTree> tree = this->getTree();
	if( !tree )
		return;

	boost::shared_ptr<wxFileDialog> dlg( 
		new wxFileDialog(
			this, 
			"Write tree",
			"", 
			"", 
			"Tree files(*.tree)|*.tree|All files(*.*)|*.*",
			wxFD_SAVE, 
			wxDefaultPosition),
		boost::mem_fn(&wxFileDialog::Destroy) );

	if ( dlg.get()->ShowModal() == wxID_OK )
	{
		std::string filename = dlg.get()->GetPath().GetData();
        std::ofstream ofs( filename.c_str(), std::ios::binary );
		if( !ofs )
		{
			errMsg( "Unable to open file '" + filename + "' for writing" );
			return;
		}

		writePrimitive<size_t>( tree->samplingRate(), ofs );
		writeScene( tree->scene(), ofs );
		writeScene( tree->parentScene(), ofs );

		std::vector<size_t> fixedObjects = tree->fixedObjects();
		dumpArray( fixedObjects, ofs );

		float startTime = tree->startTime();
		writePrimitive<float>( startTime, ofs );

		std::vector<PiecewisePath> fixedObjectsHistory = 
			tree->fixedObjectsHistory();
		writePrimitive<size_t>( fixedObjectsHistory.size(), ofs );
		for( size_t i = 0; i < fixedObjectsHistory.size(); ++i )
		{
			writePiecewisePath( fixedObjectsHistory.at(i), ofs );
		}

		size_t numPaths = tree->numPaths();
		writePrimitive<size_t>( numPaths, ofs );
		for( size_t i = 0; i < numPaths; ++i )
			tree->path(i)->write( ofs );

		this->log( "Saved tree: " + filename );
	}
}

void MainFrame::OnExportMayaAscii(wxCommandEvent& event)
{
	boost::shared_ptr<wxFileDialog> dlg( 
		new wxFileDialog(
			this, 
			"Export Maya Ascii file",
			"", 
			"", 
			"Maya Ascii files(*.ma)|*.ma|All files(*.*)|*.*",
			wxFD_SAVE, 
			wxDefaultPosition),
		boost::mem_fn(&wxFileDialog::Destroy) );

	if ( dlg.get()->ShowModal() == wxID_OK )
	{
		SimulationTreePtr tree;
		{
			boost::mutex::scoped_lock lock( currentPathMutex_ );
			tree = this->currentTree_;
		}

		std::string filename = dlg.get()->GetPath().GetData();

		// @todo fix this
		std::vector<ConstPhysicsObjectPtr> dynamicObjects;
		if( tree )
			dynamicObjects = tree->dynamicObjects();

		scene_->writeMayaAsciiFile( filename, this->piecewisePaths(), this->piecewisePathsFrameRate(), dynamicObjects );

		this->log( "Dump as Maya Ascii: " + filename );
	}
}

void MainFrame::OnExportRib(wxCommandEvent& event)
{
	boost::shared_ptr<wxFileDialog> dlg( 
		new wxFileDialog(
			this, 
			"Export Renderman .rib file",
			"", 
			"", 
			"RIB files(*.metarib)|*.metarib|All files(*.*)|*.*",
			wxFD_SAVE, 
			wxDefaultPosition),
		boost::mem_fn(&wxFileDialog::Destroy) );

	if ( dlg.get()->ShowModal() == wxID_OK )
	{
		SimulationTreePtr tree;
		{
			boost::mutex::scoped_lock lock( currentPathMutex_ );
			tree = this->currentTree_;
		}

		std::string filename = dlg.get()->GetPath().GetData();

		// @todo fix this
		std::vector<ConstPhysicsObjectPtr> dynamicObjects;
		if( tree )
			dynamicObjects = tree->dynamicObjects();
		std::vector<ConstSceneGraphElementPtr> dynObjects( dynamicObjects.begin(), dynamicObjects.end() );

		scene_->writeRibFile( filename, this->stateHistory(), dynObjects );

		this->log( "Dump as RIB: " + filename );
	}
}

void MainFrame::OnExportImage(wxCommandEvent& event)
{
	boost::shared_ptr<wxFileDialog> dlg( 
		new wxFileDialog(
			this, 
			"Export image",
			"", 
			"", 
			"PNG Files(*.png)|*.png|JPEG Files(*.jpg)|*.jpg|All files(*.*)|*.*",
			wxFD_SAVE, 
			wxDefaultPosition),
		boost::mem_fn(&wxFileDialog::Destroy) );

	if ( dlg.get()->ShowModal() == wxID_OK )
	{
		boost::shared_ptr<wxTextEntryDialog> textDlg;

		int width = 720;
		textDlg.reset( new wxTextEntryDialog(this, 
			"Width of image", "Width of image to export"),
			boost::mem_fn(&wxTextEntryDialog::Destroy) );
		do
		{
			std::ostringstream oss;
			oss << width;
			textDlg.get()->SetValue( wxString(oss.str().c_str()) );

			if( textDlg.get()->ShowModal() != wxID_OK )
				return;
			width = atoi( textDlg.get()->GetValue().GetData() );
		} while( width <= 0 || width >= 20000 );

		int height = 480;
		textDlg.reset( new wxTextEntryDialog(this, 
			"Height of image", "Height of image to export"),
			boost::mem_fn(&wxTextEntryDialog::Destroy) );
		{
			std::ostringstream oss;
			oss << height;
			textDlg.get()->SetValue( wxString(oss.str().c_str()) );

			if( textDlg.get()->ShowModal() != wxID_OK )
				return;
			textDlg.get()->Destroy();
			height = atoi( textDlg.get()->GetValue().GetData() );
		} while( height <= 0 || height >= 20000 );

		std::string filename = dlg.get()->GetPath().GetData();
		currentCanvas_->dumpFrame( width, height, filename );
		currentCanvas_->Refresh(FALSE);
	}
}

void MainFrame::OnExportTimelapseImage(wxCommandEvent& event)
{
	SimulationTreePtr tree = this->currentTree_;
	const Path* path = currentPath_;
	if( path == 0 )
	{
		errMsg( "No path selected." );
		return;
	}

	std::deque<SceneGraphElementPtr> backgroundElements;
	std::deque<SceneGraphElementPtr> foregroundElements;
	std::deque<SceneGraphElementPtr> dontRenderElements;

	std::deque<SceneGraphElementPtr> backupSelection;
	backupSelection.swap( this->selected_ );

	{
		ODENativeSimulation sim( this->scene_ );
		std::vector<ConstPhysicsObjectPtr> dynamicObjects = sim.dynamicObjects();
		std::sort( dynamicObjects.begin(), dynamicObjects.end() );

		std::deque<SceneGraphElementPtr> toHandle( 1, scene_->root() );
		while( !toHandle.empty() )
		{
			SceneGraphElementPtr current = toHandle.front();
			toHandle.pop_front();

			if( !current->visible() )
				continue;

			if( boost::dynamic_pointer_cast<const Joint>( current ) )
				dontRenderElements.push_back( current );
			else if( boost::dynamic_pointer_cast<const PhysicsObject>(current) )
			{
				if( std::binary_search(dynamicObjects.begin(), dynamicObjects.end(), current) )
					foregroundElements.push_back( current );
				else
					backgroundElements.push_back( current );
			}
			else
			{
				std::copy( current->children_begin(), current->children_end(), 
					std::front_inserter(toHandle) );
			}
		}
	}

	boost::shared_ptr<wxFileDialog> dlg( 
		new wxFileDialog(
			this, 
			"Export timelapse image",
			"", 
			"", 
			"EXR file (*.exr)|*.exr|All files(*.*)|*.*",
			wxFD_SAVE, 
			wxDefaultPosition),
		boost::mem_fn(&wxFileDialog::Destroy) );

	if ( dlg.get()->ShowModal() == wxID_OK )
	{
		std::string outFilename = dlg.get()->GetPath().GetData();

		boost::shared_ptr<wxTextEntryDialog> textDlg;

		int width = 720;
		textDlg.reset( new wxTextEntryDialog(this, 
			"Width of image", "Width of image to export"),
			boost::mem_fn(&wxTextEntryDialog::Destroy) );
		do
		{
			std::ostringstream oss;
			oss << width;
			textDlg.get()->SetValue( wxString(oss.str().c_str()) );

			if( textDlg.get()->ShowModal() != wxID_OK )
				return;
			width = atoi( textDlg.get()->GetValue().GetData() );
		} while( width <= 0 || width >= 20000 );

		int height = 480;
		textDlg.reset( new wxTextEntryDialog(this, 
			"Height of image", "Height of image to export"),
			boost::mem_fn(&wxTextEntryDialog::Destroy) );
		{
			std::ostringstream oss;
			oss << height;
			textDlg.get()->SetValue( wxString(oss.str().c_str()) );

			if( textDlg.get()->ShowModal() != wxID_OK )
				return;
			textDlg.get()->Destroy();
			height = atoi( textDlg.get()->GetValue().GetData() );
		} while( height <= 0 || height >= 20000 );

		// Background needs to be black
		MainFrame::BackgroundStyle oldStyle = this->backgroundStyle();
		this->setBackgroundStyle( BACKGROUND_BLACK );

		float oldTime = time_;
		bool oldWireframeOnShaded = this->wireframeOnShaded();
		bool oldRenderGrid = renderMenu_->IsChecked(MENU_RENDER_GRID);
		renderMenu_->Check(MENU_RENDER_WIREFRAME_ON_SHADED, false);
		renderMenu_->Check(MENU_RENDER_GRID, false);

		try
		{
			const size_t rescale = 1;

			// Write out background plate
			std::string backgroundPlateFilename = outFilename + ".background.png";
			std::for_each( backgroundElements.begin(), backgroundElements.end(),
				boost::bind( &SceneGraphElement::setVisible, _1, true ) );
			std::for_each( foregroundElements.begin(), foregroundElements.end(),
				boost::bind( &SceneGraphElement::setVisible, _1, false ) );

			this->currentCanvas_->dumpFrame( rescale*width, rescale*height, backgroundPlateFilename );


			// now write the timelapse
			std::for_each( dontRenderElements.begin(), dontRenderElements.end(),
				boost::bind( &SceneGraphElement::setVisible, _1, false ) );
			std::for_each( backgroundElements.begin(), backgroundElements.end(),
				boost::bind( &SceneGraphElement::setVisible, _1, false ) );
			std::for_each( foregroundElements.begin(), foregroundElements.end(),
				boost::bind( &SceneGraphElement::setVisible, _1, true ) );

			using namespace Imf;

			std::vector<BYTE> buffer( (rescale*width)*(rescale*height)*3 );
			std::vector<size_t> summedBuffer( width*height*3 );

			std::pair<double, double> timeRange = 
				this->timeView_->activeTimeRange();

			this->time_ = timeRange.first;
			this->currentCanvas_->dumpFrame( rescale*width, rescale*height, outFilename + ".start.png" );

			this->time_ = timeRange.second;
			this->currentCanvas_->dumpFrame( rescale*width, rescale*height, outFilename + ".end.png" );

			size_t curFrame = 0;
			double fps = 360.0;
			double lastDump = 0.0;
			for (double curTime = timeRange.first; curTime < timeRange.second; curTime += 1.0/fps)
			{
				++curFrame;

				this->time_ = curTime;

				// so the user has feedback
				this->timeView_->Refresh( true );
				this->timeView_->Update();

				this->currentCanvas_->dumpFrame( rescale*width, rescale*height, &buffer[0] );
				for( unsigned int i = 0; i < 3*width*height; ++i )
					summedBuffer[i] += buffer[i];
				/*
				wxImage image( rescale*width, rescale*height, &buffer[0], true );
				image.Rescale( width, height, wxIMAGE_QUALITY_HIGH );

				for( int x = 0; x < width; ++x )
				{
					for( int y = 0; y < height; ++y )
					{
						summedBuffer[3*(y*width + x)+0] += image.GetRed(x, y);
						summedBuffer[3*(y*width + x)+1] += image.GetGreen(x, y);
						summedBuffer[3*(y*width + x)+2] += image.GetBlue(x, y);
					}
				}
				*/

				if( (curTime + 1.0/fps) >= timeRange.second || (curTime - lastDump) > 1.0 )
				{
					lastDump = curTime;

					std::vector<Rgba> pixels( width*height );
					for( size_t i = 0; i < width*height; ++i )
					{
						pixels[i].r = half( static_cast<double>(summedBuffer[3*i + 2]) / (10.0*fps) );
						pixels[i].g = half( static_cast<double>(summedBuffer[3*i + 1]) / (10.0*fps) );
						pixels[i].b = half( static_cast<double>(summedBuffer[3*i + 0]) / (10.0*fps) );
						pixels[i].a = half(1.0f);
					}

					try
					{
						Imf::RgbaOutputFile outputFile( outFilename.c_str(),
							width,
							height,
							WRITE_RGBA );
						outputFile.setFrameBuffer( &pixels[0], 1, width );
						outputFile.writePixels( height );
					}
					catch( std::exception& e )
					{
						errMsg( "Error writing to " + outFilename + ": "
							+ e.what() );
						goto restoreScreen;
					}
				}
			}

			boost::shared_ptr<wxMessageDialog> 
				dlg( new wxMessageDialog(this, "Image export", "Image export complete.", wxOK ),
				boost::mem_fn(&wxMessageDialog::Destroy) );
			dlg.get()->ShowModal();
		}
		catch( IOException& e )
		{
			errMsg( e.message() );
		}

restoreScreen:
		time_ = oldTime;

		backupSelection.swap( this->selected_ );

		this->setBackgroundStyle( oldStyle );
		renderMenu_->Check(MENU_RENDER_WIREFRAME_ON_SHADED, oldWireframeOnShaded);
		renderMenu_->Check(MENU_RENDER_GRID, oldRenderGrid);

		std::for_each( dontRenderElements.begin(), dontRenderElements.end(),
			boost::bind( &SceneGraphElement::setVisible, _1, true ) );
		std::for_each( backgroundElements.begin(), backgroundElements.end(),
			boost::bind( &SceneGraphElement::setVisible, _1, true ) );
		std::for_each( foregroundElements.begin(), foregroundElements.end(),
			boost::bind( &SceneGraphElement::setVisible, _1, true ) );

		currentCanvas_->Refresh(FALSE);
	}
}


// this is straight from the Microsoft example code
// http://msdn2.microsoft.com/en-us/library/aa391082.aspx
#ifndef SAFE_RELEASE
#define SAFE_RELEASE(x) \
   if(x != NULL)        \
   {                    \
      x->Release();     \
      x = NULL;         \
   }
#endif

#ifndef SAFE_ARRAY_DELETE
#define SAFE_ARRAY_DELETE(x) \
   if(x != NULL)             \
   {                         \
      delete[] x;            \
      x = NULL;              \
   }
#endif

#ifndef GOTO_EXIT_IF_FAILED
#define GOTO_EXIT_IF_FAILED(hr) if(FAILED(hr)) goto Exit;
#endif

boost::shared_ptr<IWMInputMediaProps> FindInputFormat(IWMWriter* pWriter, 
                       DWORD dwInput,
                       GUID guidSubType)
{
    HRESULT hr       = S_OK;

    // Find the number of formats supported by this input.
    DWORD   cFormats = 0;
    hr = pWriter->GetInputFormatCount(dwInput, &cFormats);
	if( FAILED(hr) )
		throw Exception( "Unable to get input format count" );

    // Loop through all of the supported formats.
    for (DWORD formatIndex = 0; formatIndex < cFormats; formatIndex++)
    {
        // Get the input media properties for the input format.
		boost::shared_ptr<IWMInputMediaProps> mediaProps;
		{
			IWMInputMediaProps* pProps = NULL;
			hr = pWriter->GetInputFormat(dwInput, formatIndex, &pProps);
			if( FAILED(hr) )
				throw Exception( "Unable to get input format" );
			mediaProps.reset( pProps, &GenericRelease<IWMInputMediaProps> );
		}

        // Get the size of the media type structure.
		DWORD   cbSize   = 0;
        hr = mediaProps->GetMediaType(NULL, &cbSize);
		if( FAILED(hr) )
			throw Exception( "Unable to get format input" );

		std::vector<BYTE> bytes( cbSize );
		WM_MEDIA_TYPE* pType = reinterpret_cast<WM_MEDIA_TYPE*>( &bytes[0] );
        // Get the media type structure.
        hr = mediaProps->GetMediaType(pType, &cbSize);
		if( FAILED(hr) )
			throw Exception( "Unable to get format input" );

        if(pType->subtype == guidSubType)
			return mediaProps;
    }

	return boost::shared_ptr<IWMInputMediaProps>();
}

void runMediaPlayer( const std::string& inputFile )
{
	STARTUPINFO si;
	PROCESS_INFORMATION pi;
	ZeroMemory( &si, sizeof(si) );
	si.cb = sizeof(si);
	ZeroMemory( &pi, sizeof(pi) );

	// todo make this somehow user-configurable
	std::string wmplayerLoc( "c:\\Program Files\\Windows Media Player\\wmplayer.exe" );
	std::string commandLine( wmplayerLoc );
	commandLine += " " + inputFile;

	std::vector<char> commandLineVec( commandLine.begin(), commandLine.end() );
	commandLineVec.push_back( 0 );

	if( !CreateProcess(
		wmplayerLoc.c_str(),
		&commandLineVec[0],
		NULL,                // Process handle not inheritable. 
		NULL,                // Thread handle not inheritable. 
		FALSE,               // Set handle inheritance to FALSE. 
		IDLE_PRIORITY_CLASS, // Set to very low priority
		NULL,                // Use parent's environment block. 
		NULL,                // Use parent's starting directory. 
		&si,                 // Pointer to STARTUPINFO structure.
		&pi )                // Pointer to PROCESS_INFORMATION structure.
		)
	{
		std::ostringstream oss;
		oss << "Unable to execute process '" << &commandLineVec[0] << "'." << std::endl;
		oss << "CreateProcess returned error: " << GetLastError() << std::endl;
		throw CreationException( oss.str() );
	}

	// Close process and thread handles. 
	CloseHandle( pi.hProcess );
	CloseHandle( pi.hThread );
}

#ifdef SOUND
class RecordSoundIdleAction
	: public IdleAction
{
public:
	RecordSoundIdleAction(
		MainFrame* mainFrame,
		size_t audioFrameRate )
		:	mainFrame_(mainFrame),
			audioFrameRate_(audioFrameRate),
			origLastAudioValue_(0.0), 
			dcRemovedLastAudioValue_(0.0),
			currentAudioFrame_(0),
			maxSeenAudioPressure_(0.0),
			interrupted_(false),
			mouseCapture_(mainFrame)
	{
		mainFrame->SetFocus();

		this->timeRange_ = mainFrame->timeView_->activeTimeRange();

		this->oldSimulation_ = mainFrame->simulation_;
		if( mainFrame->mode_ == MainFrame::MODE_EDIT )
		{
			boost::shared_ptr<SimulationFactory> factory = mainFrame->simulationFactory();
			mainFrame->simulation_ = factory->construct( mainFrame->scene() );
		}

		this->oldTime_ = mainFrame->time_;
	}

	void perform()
	{
		if( timeRange_.first + boost::numeric_cast<double>(currentAudioFrame_) 
			/ boost::numeric_cast<double>(audioFrameRate_) > timeRange_.second )
		{
			this->interrupted_ = true;
			return;
		}

		boost::shared_ptr<Simulation> simulation = mainFrame_->simulation_;
		boost::shared_ptr<ODESimulation> odeSimulation = 
			boost::dynamic_pointer_cast<ODESimulation>( simulation );

		// for the leaky integrator:
		double leakyIntegratorConstant = 0.99;

		const size_t subTimestepsPerFrame = 32;
		const size_t videoFrameRate = 30;
		const size_t timestepsPerSecond = videoFrameRate*subTimestepsPerFrame;
		const double stepSize = 1.0 / boost::numeric_cast<double>(timestepsPerSecond);

		if( simulation )
		{
			for( size_t j = 0; j < subTimestepsPerFrame; ++j )
			{
				Simulation::State prevState = simulation->getSimulationState();
				simulation->step( stepSize );

				const double endTime = 
					boost::numeric_cast<double>(currentAudioFrame_) / 
					boost::numeric_cast<double>(audioFrameRate_) + stepSize;
				size_t goalAudioFrame = 
					boost::numeric_cast<size_t>( endTime * boost::numeric_cast<double>(audioFrameRate_) );
				assert( currentAudioFrame_ < goalAudioFrame );
				size_t numAudioFrames = goalAudioFrame - currentAudioFrame_;
				std::vector<double> audioFrames;
				if( odeSimulation )
					audioFrames = odeSimulation->stepAudio( prevState, stepSize, numAudioFrames );
				else
					audioFrames = std::vector<double>( numAudioFrames, 0.0 );

				size_t prevSize = allAudioFrames_.size();
				allAudioFrames_.resize( prevSize + audioFrames.size() );

				// apply leaky integrator
				for( size_t i = 0; i < audioFrames.size(); ++i )
				{
					double curValue = audioFrames[i];
					double dcRemovedValue = leakyIntegratorConstant*dcRemovedLastAudioValue_ + 
						curValue - origLastAudioValue_;
					allAudioFrames_[prevSize + i] = dcRemovedValue;
					dcRemovedLastAudioValue_ = dcRemovedValue;
					origLastAudioValue_ = curValue;
					maxSeenAudioPressure_ = std::max( 
						fabs(maxSeenAudioPressure_), dcRemovedValue );
				}

				currentAudioFrame_ = goalAudioFrame;
			}
		}

		mainFrame_->time_ = timeRange_.first + boost::numeric_cast<double>(currentAudioFrame_) / 
			boost::numeric_cast<double>(audioFrameRate_);
		mainFrame_->timeView_->Refresh(FALSE);
		mainFrame_->timeView_->Update();
	}

	bool complete() const
	{
		return interrupted_;
	}

	void interrupt()
	{
		this->interrupted_ = true;
	}

	~RecordSoundIdleAction()
	{
		mainFrame_->simulation_ = this->oldSimulation_;
		mainFrame_->time_ = this->oldTime_;
		mainFrame_->prevHighestAudio_ = 
			std::max(this->maxSeenAudioPressure_, 0.5*mainFrame_->prevHighestAudio_);

		std::vector<short> quantizedFrames( 2*allAudioFrames_.size() );
		const double scale = 32767.0 / maxSeenAudioPressure_;
		for( size_t i = 0; i < allAudioFrames_.size(); ++i )
		{
			short val = static_cast<short>( scale * allAudioFrames_[i] );
			quantizedFrames[2*i+0] = val;
			quantizedFrames[2*i+1] = val;
		}

		createWavFile( "test.wav", quantizedFrames );
	}

private:
	SimulationPtr oldSimulation_;
	float oldTime_;

	const size_t audioFrameRate_;
	MainFrame* mainFrame_;

	std::deque<double> allAudioFrames_;

	double origLastAudioValue_;
	double dcRemovedLastAudioValue_;
	double maxSeenAudioPressure_;

	size_t currentAudioFrame_;
	bool interrupted_;

	// need to capture the mouse so that we don't lose focus
	MouseCapture mouseCapture_;

	std::pair<double, double> timeRange_;
};

void MainFrame::recordSound()
{
	if( this->idleAction_ )
		return;

	const size_t audioFrameRate = 44100;
	this->idleAction_.reset( 
		new RecordSoundIdleAction( 
			this, 
			44100 ) );
}
#endif // SOUND

IdleAction::~IdleAction()
{
}

#ifdef WINDOWS_MEDIA
class RecordMovieIdleAction
	: public IdleAction
{
public:
	RecordMovieIdleAction( 
		MainFrame* mainFrame,
		boost::shared_ptr<IWMStreamConfig> videoStreamConfig,
		boost::shared_ptr<IWMStreamConfig> audioStreamConfig,
		size_t width,
		size_t height,
		size_t videoBitrate,
		size_t videoFrameRate,
		size_t audioFrameRate )
		:	origLastAudioValue_(0.0), 
			dcRemovedLastAudioValue_(0.0),
			currentAudioFrame_(0),
			currentVideoFrame_(0),
			maxSeenAudioPressure_(0.0),
			interrupted_(false),
			mouseCapture_(mainFrame)
	{
		/*
   // Get the default control word.
   int cw = _controlfp( 0,0 );

   // Set the exception masks OFF, turn exceptions on.
   cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_ZERODIVIDE|EM_DENORMAL);

   // Set the control word.
   _controlfp( cw, MCW_EM );
   */

		mainFrame->SetFocus();

		this->width_ = width;
		this->height_ = height;
		this->audioFrameRate_ = audioFrameRate;
		this->videoFrameRate_ = videoFrameRate;
		this->mainFrame_ = mainFrame;

		this->timeRange_ = mainFrame->timeView_->activeTimeRange();
		if( !mainFrame->piecewisePaths_.empty() )
		{
			const double frameRate = mainFrame_->piecewisePathsFrameRate_;
			// round to frame rate of paths to improve accuracy:
			timeRange_.first = floor(timeRange_.first * frameRate) / frameRate;
			timeRange_.second = ceil(timeRange_.second * frameRate) / frameRate;
		}

		HRESULT hr;

		// first step is to create a profile
		boost::shared_ptr<IWMProfile> profile;
		{
			IWMProfile* ppProfile;
			hr = mainFrame->profileManager_->CreateEmptyProfile(WMT_VER_9_0, &ppProfile);
			if( FAILED(hr) )
				throw Exception( "Unable to create profile" );

			profile.reset( ppProfile, &GenericRelease<IWMProfile> );
		}

		// add the video stream
		wchar_string videoConnectionName = toWideString( std::string("video") );
		{
			// video
			boost::shared_ptr<IWMStreamConfig> streamConfig = videoStreamConfig;

			streamConfig->SetBitrate( videoBitrate );
			streamConfig->SetBufferWindow( -1 );
			streamConfig->SetStreamNumber( 1 );

			std::vector<WCHAR> connectionNameVec( videoConnectionName.begin(), videoConnectionName.end() );
			connectionNameVec.push_back( 0 );
			streamConfig->SetConnectionName( &connectionNameVec[0] );

			// Get the media properties interface.
			boost::shared_ptr<IWMMediaProps> props;
			{
				IWMMediaProps* pProps;
				hr = streamConfig->QueryInterface(IID_IWMMediaProps, (void**)&pProps);
				if( FAILED(hr) )
					throw Exception( "Error retrieving codecs" );

				props.reset( pProps, &GenericRelease<IWMMediaProps> );
			}

			// Get the size required for the media type structure.
			DWORD cbType;
			hr = props->GetMediaType(NULL, &cbType);
			if( FAILED(hr) )
				throw Exception( "Error retrieving codecs" );

			// get the media type structure
			std::vector<BYTE> typeBytes( cbType );
			WM_MEDIA_TYPE* pType = reinterpret_cast<WM_MEDIA_TYPE*>( &typeBytes[0] );
			hr = props->GetMediaType(pType, &cbType);
			if( FAILED(hr) )
				throw Exception( "Error retrieving codecs" );

			WMVIDEOINFOHEADER* pVidHeader = NULL;

			// Check that the format data is present.
			if(pType->cbFormat >= sizeof(WMVIDEOINFOHEADER))
				pVidHeader = (WMVIDEOINFOHEADER*) pType->pbFormat;
			else
				throw Exception( "Unexpected Error retrieving video codecs" );

			// modify the profile with the correct size, bitrate information
			pVidHeader->rcSource.left = 0;
			pVidHeader->rcSource.top = 0;
			pVidHeader->rcSource.right = this->width_;
			pVidHeader->rcSource.bottom = this->height_;

			pVidHeader->rcTarget = pVidHeader->rcSource;

			pVidHeader->dwBitRate = videoBitrate;

			pVidHeader->AvgTimePerFrame = 
				boost::numeric_cast<LONGLONG>( 1.0e7 / 
					boost::numeric_cast<double>( this->videoFrameRate_ ) );

			pVidHeader->bmiHeader.biWidth = this->width_;
			pVidHeader->bmiHeader.biHeight = this->height_;

			// now, write the media type back to the properties
			{
				HRESULT hr = props->SetMediaType( pType );
				if( FAILED(hr) )
					throw Exception( "Unable to create profile" );
			}

			// We're finished with this stream, add it to the profile
			{
				HRESULT hr = profile->AddStream( streamConfig.get() );
				if( FAILED(hr) )
					throw Exception( "Unable to create profile" );
			}
		}

		// add the audio stream
		wchar_string audioConnectionName = toWideString( std::string("audio") );
		{
			boost::shared_ptr<IWMStreamConfig> streamConfig = audioStreamConfig;

			std::vector<WCHAR> connectionNameVec( audioConnectionName.begin(), audioConnectionName.end() );
			connectionNameVec.push_back( 0 );
			streamConfig->SetConnectionName( &connectionNameVec[0] );
			streamConfig->SetStreamNumber( 2 );

			// We're finished with this stream, add it to the profile
			{
				HRESULT hr = profile->AddStream( streamConfig.get() );
				if( FAILED( hr ) )
				{
					if( E_INVALIDARG == hr )
						throw Exception( "AddStream failed: The pConfig parameter is NULL" );
					else if( E_OUTOFMEMORY == hr )
						throw Exception( "AddStream failed: There is not enough available memory" );
					else if( NS_E_INVALID_STREAM == hr )
						throw Exception( "AddStream failed: The stream is not valid, possibly because it does not have a valid stream number" );
					else if( E_FAIL == hr )
						throw Exception( "AddStream failed: The method failed for an unspecified reason" );
					else
						throw Exception( "AddStream failed: Unexpected error" );
				}
			}
		}

		// Okay, now, start writing the .asf file.
		// first step is to create the writer object
		{
			IWMWriter* pWriter;
			hr = WMCreateWriter( NULL, &pWriter );
			if( FAILED(hr) )
				throw Exception( "Unable to create profile" );
			asfWriter_.reset( pWriter, &GenericRelease<IWMWriter> );
		}

		hr = asfWriter_->SetProfile( profile.get() );
		if( FAILED(hr) )
			throw Exception( "Unable to create profile" );

		wchar_string outputFilename = toWideString( std::string("test.asf") );
		hr = asfWriter_->SetOutputFilename( outputFilename.c_str() );
		if( FAILED(hr) )
			throw Exception( "Unable to create profile" );

		// we need to map the inputs back to audio and video
		// unfortunately, the MSDN docs specify that there is no
		//  guarantee that they will retain the profile order,
		//  so a linear search is necessary
		{
			DWORD cInputs;
			hr = asfWriter_->GetInputCount(&cInputs);

			videoConnectionId_ = cInputs;
			audioConnectionId_ = cInputs;

			for( DWORD iInput = 0; iInput < cInputs; ++iInput )
			{
				boost::shared_ptr<IWMInputMediaProps> props;
				{
					IWMInputMediaProps* pProps = NULL;
					hr = asfWriter_->GetInputProps(iInput, &pProps);
					if( FAILED(hr) )
						throw Exception( "Unable to get input properties" );

					props.reset( pProps, &GenericRelease<IWMInputMediaProps> );
				}

				WORD cchName = 0;
				hr = props->GetConnectionName(NULL, &cchName);
				if( FAILED(hr) )
					throw Exception( "Unable to get connection name" );

				std::vector<WCHAR> connNameVec( cchName );
				hr = props->GetConnectionName(&connNameVec[0], &cchName);
				if( FAILED(hr) )
					throw Exception( "Unable to get connection name" );
				if(cchName == 0)
					continue;
				wchar_string connName( &connNameVec[0] );

				if( connName == videoConnectionName )
					videoConnectionId_ = iInput;
				else if( connName == audioConnectionName )
					audioConnectionId_ = iInput;
			}

			if( videoConnectionId_ >= cInputs )
				throw Exception( "Video connection not found" );
			if( audioConnectionId_ >= cInputs )
				throw Exception( "Audio connection not found" );
		}

		// need to set the video to accept RGB 24
		{
			boost::shared_ptr<IWMInputMediaProps> props =
				FindInputFormat(asfWriter_.get(), videoConnectionId_, WMMEDIASUBTYPE_RGB24);
			if( !props )
				throw Exception( "Unable to find RGB24 input format" );

			DWORD cbSize = 0;
			hr = props->GetMediaType(NULL, &cbSize);
			if( FAILED(hr) )
				throw Exception( "Unable to get media type" );

			std::vector<BYTE> bytes( cbSize );
			WM_MEDIA_TYPE* pType = reinterpret_cast<WM_MEDIA_TYPE*>( &bytes[0] );
			hr = props->GetMediaType(pType, &cbSize);
			if( FAILED(hr) )
				throw Exception( "Unable to get media type" );

			// Adjust the format to match source images.
			// First set pointers to the video structures.
			WMVIDEOINFOHEADER* pVidHdr = reinterpret_cast<WMVIDEOINFOHEADER*>(pType->pbFormat);
			BITMAPINFOHEADER*  pBMHdr  = &(pVidHdr->bmiHeader);
			pBMHdr->biWidth  = this->width_;
			pBMHdr->biHeight = this->height_;

			// Stride = (width * bytes/pixel), rounded to the next DWORD boundary
			this->stride_ = (width * (pBMHdr->biBitCount / 8) + 3) & ~3;

			// Image size = stride * height. 
			pBMHdr->biSizeImage = this->height_ * this->stride_;

			// Apply the adjusted type to the video input. 
			hr = props->SetMediaType(pType);
			if( FAILED(hr) )
				throw Exception( "Unable to get media type" );

			hr = asfWriter_->SetInputProps(videoConnectionId_, props.get());
			if( FAILED(hr) )
				throw Exception( "Unable to get media type" );
		}

		// set audio to uncompressed PCM
		{
			boost::shared_ptr<IWMInputMediaProps> props =
				FindInputFormat(asfWriter_.get(), audioConnectionId_, WMMEDIASUBTYPE_PCM);
			if( !props )
				throw Exception( "Unable to find PCM input format" );

			DWORD cbSize = 0;
			hr = props->GetMediaType(NULL, &cbSize);
			if( FAILED(hr) )
				throw Exception( "Unable to get media type" );

			std::vector<BYTE> bytes( cbSize );
			WM_MEDIA_TYPE* pType = reinterpret_cast<WM_MEDIA_TYPE*>( &bytes[0] );
			hr = props->GetMediaType(pType, &cbSize);
			if( FAILED(hr) )
				throw Exception( "Unable to get media type" );

			WAVEFORMATEX* pAudioHdr = reinterpret_cast<WAVEFORMATEX*>(pType->pbFormat);
			assert( pAudioHdr->wFormatTag == 0x0001 ); // PCM audio
			pAudioHdr->nChannels = 2;
			pAudioHdr->nSamplesPerSec = this->audioFrameRate_;
			pAudioHdr->wBitsPerSample = 16;
			pAudioHdr->nBlockAlign = pAudioHdr->nChannels*(pAudioHdr->wBitsPerSample / 8);
			pAudioHdr->cbSize = 0;
			pAudioHdr->nAvgBytesPerSec = pAudioHdr->nSamplesPerSec*pAudioHdr->nBlockAlign;

			// not sure about this one:
			pType->lSampleSize = pAudioHdr->nBlockAlign;

			// Apply the adjusted type to the audio input. 
			hr = props->SetMediaType(pType);
			if( FAILED(hr) )
				throw Exception( "Unable to get media type" );

			hr = asfWriter_->SetInputProps(audioConnectionId_, props.get());
			if( FAILED(hr) )
				throw Exception( "Error in SetInputProps" );
		}


		this->oldSimulation_ = mainFrame->simulation_;
		if( !mainFrame_->piecewisePaths_.empty() )
			this->simulation_.reset( new ODEDefaultSimulationType(mainFrame->scene()) );
		else
		{
			boost::shared_ptr<SimulationFactory> factory = mainFrame->simulationFactory();
			this->simulation_ = factory->construct( mainFrame->scene() );
		}

		if( mainFrame->mode_ == MainFrame::MODE_EDIT )
			mainFrame->simulation_ = this->simulation_;

		this->oldTime_ = mainFrame->time_;
		this->maxAudioPressure_ = mainFrame->prevHighestAudio_;

		// start writing
		hr = asfWriter_->BeginWriting();
		if( FAILED( hr ) )
		{
			if( NS_E_VIDEO_CODEC_NOT_INSTALLED == hr )
				throw Exception( "BeginWriting failed: Video Codec not installed" );
			if( NS_E_AUDIO_CODEC_NOT_INSTALLED == hr )
				throw Exception( "BeginWriting failed: Audio Codec not installed" );
			else if( NS_E_INVALID_OUTPUT_FORMAT == hr )
				throw Exception( "BeginWriting failed: Invalid Output Format" );
			else if( NS_E_VIDEO_CODEC_ERROR == hr )
				throw Exception( "BeginWriting failed: An unexpected error occurred with the video codec" );
			else if( NS_E_AUDIO_CODEC_ERROR == hr )
				throw Exception( "BeginWriting failed: An unexpected error occurred with the audio codec" );
			else
				throw Exception( "BeginWriting failed: Error" );
		}
	}

	void perform()
	{
		if( timeRange_.first + boost::numeric_cast<double>(currentVideoFrame_) / 
			boost::numeric_cast<double>(videoFrameRate_) > timeRange_.second )
		{
			this->interrupted_ = true;
			return;
		}

		boost::shared_ptr<Simulation> simulation = this->simulation_;
		boost::shared_ptr<ODESimulation> odeSimulation = 
			boost::dynamic_pointer_cast<ODESimulation>( simulation );

		// for the leaky integrator:
		double leakyIntegratorConstant = 0.99;

		glPixelStorei( GL_PACK_ALIGNMENT, 8 );
		glPixelStorei( GL_PACK_SWAP_BYTES, GL_FALSE );

		HRESULT hr;

		size_t finalGoalAudioFrame = boost::numeric_cast<double>(audioFrameRate_)
			* boost::numeric_cast<double>(currentVideoFrame_+1)
				/ boost::numeric_cast<double>(videoFrameRate_);
		size_t startAudioFrame = currentAudioFrame_;
		std::vector<double> allAudioFrames( finalGoalAudioFrame - startAudioFrame, 0.0 );

		if( mainFrame_->mode_ == MainFrame::MODE_SOLVE && !mainFrame_->piecewisePaths_.empty() )
		{
			assert( odeSimulation );

			assert( mainFrame_->piecewisePathsFrameRate_ % videoFrameRate_ == 0 );
			const size_t subTimestepsPerFrame = mainFrame_->piecewisePathsFrameRate_ / videoFrameRate_;

			for( size_t j = 0; j < subTimestepsPerFrame; ++j )
			{
				double startTime = timeRange_.first + 
					boost::numeric_cast<double>(currentVideoFrame_*subTimestepsPerFrame + j) 
					/ boost::numeric_cast<double>(videoFrameRate_*subTimestepsPerFrame);
				double endTime = timeRange_.first + 
					boost::numeric_cast<double>(currentVideoFrame_*subTimestepsPerFrame + j + 1) 
					/ boost::numeric_cast<double>(videoFrameRate_*subTimestepsPerFrame);

				mainFrame_->time_ = startTime;
				Simulation::State prevState = mainFrame_->state();

				mainFrame_->time_ = endTime;
				Simulation::State nextState = mainFrame_->state();

				size_t goalAudioFrame = 
					boost::numeric_cast<size_t>( endTime * boost::numeric_cast<double>(audioFrameRate_) );
				assert( currentAudioFrame_ < goalAudioFrame );
				assert( goalAudioFrame <= finalGoalAudioFrame );
				size_t numAudioFrames = goalAudioFrame - currentAudioFrame_;
#ifdef SOUND
				std::vector<double> audioFrames = odeSimulation->stepAudio( 
					prevState, nextState, mainFrame_->piecewisePathsFrameRate_,
					0.005, numAudioFrames );
#else
				std::vector<double> audioFrames( numAudioFrames, 0.0 );
#endif // SOUND

				std::copy( audioFrames.begin(), audioFrames.end(),
					allAudioFrames.begin() + (currentAudioFrame_ - startAudioFrame) );
				currentAudioFrame_ += audioFrames.size();
			}
		}
		else if( simulation )
		{
			const size_t subTimestepsPerFrame = 32;
			const double stepSize = 1.0 / boost::numeric_cast<double>(subTimestepsPerFrame*videoFrameRate_);

			for( size_t j = 0; j < subTimestepsPerFrame; ++j )
			{
				Simulation::State prevState = simulation->getSimulationState();
				simulation->step( stepSize );

				const double endTime = 
					boost::numeric_cast<double>(currentVideoFrame_*subTimestepsPerFrame + j + 1) * stepSize;
				size_t goalAudioFrame = 
					boost::numeric_cast<size_t>( endTime * boost::numeric_cast<double>(audioFrameRate_) );
				assert( currentAudioFrame_ < goalAudioFrame );
				size_t numAudioFrames = goalAudioFrame - currentAudioFrame_;
				std::vector<double> audioFrames;
#ifdef SOUND
				if( odeSimulation )
					audioFrames = odeSimulation->stepAudio( prevState, stepSize, numAudioFrames );
				else
#endif
					audioFrames = std::vector<double>( numAudioFrames, 0.0 );

				std::copy( audioFrames.begin(), audioFrames.end(),
					allAudioFrames.begin() + (currentAudioFrame_ - startAudioFrame) );
				currentAudioFrame_ += audioFrames.size();
			}
		}

		{
			// now, need to encode in a chunk for windows media
			size_t sampleSize = allAudioFrames.size()*4;
			boost::shared_ptr<INSSBuffer> buffer;
			{
				INSSBuffer* bufPtr = 0;
				hr = asfWriter_->AllocateSample( sampleSize, &bufPtr );
				if( FAILED(hr) )
					throw Exception( "Error allocating buffer" );

				buffer.reset( bufPtr, &GenericRelease<INSSBuffer> );
			}

			BYTE* ppdwBuffer;
			hr = buffer->GetBuffer( &ppdwBuffer );
			if( FAILED(hr) )
				throw Exception( "Error retrieving buffer" );

			assert( sizeof(short) == 2 );
			short* values = reinterpret_cast<short*>( ppdwBuffer );

			// apply leaky integrator
			for( size_t i = 0; i < allAudioFrames.size(); ++i )
			{
				const double origValue = allAudioFrames[i];
				double dcRemoved = leakyIntegratorConstant*dcRemovedLastAudioValue_ + 
					(origValue - origLastAudioValue_);
				if( fabs(dcRemoved) < 1e-10 )
					dcRemoved = 0.0;

				allAudioFrames[i] = dcRemoved;
				dcRemovedLastAudioValue_ = dcRemoved;
				origLastAudioValue_ = origValue;

				maxSeenAudioPressure_ = std::max( maxSeenAudioPressure_, 
					fabs(dcRemoved) );
				maxAudioPressure_ = std::max( maxAudioPressure_, maxSeenAudioPressure_ );

				double renormalized = std::max( std::min( 32767.0,
					32767.0 * dcRemoved / maxAudioPressure_ ), -32767.0 );
				// for now, same in left and right channel
				short value = static_cast<short>( renormalized );
				values[ 2*i+0 ] = value;
				values[ 2*i+1 ] = value;
			}

			buffer->SetLength( sampleSize );
			DWORD frameTime = 
				boost::numeric_cast<double>(currentAudioFrame_)
					/ (boost::numeric_cast<double>(audioFrameRate_)*1.0e-7);

			hr = asfWriter_->WriteSample( audioConnectionId_,
				frameTime,
				0,
				buffer.get() );
			if( FAILED(hr) )
				throw Exception( "Error writing sample" );
		}

		boost::shared_ptr<INSSBuffer> buffer;
		size_t sampleSize = this->stride_*this->height_;
		{
			INSSBuffer* bufPtr = 0;
			hr = asfWriter_->AllocateSample( sampleSize, &bufPtr );
			if( FAILED(hr) )
				throw Exception( "Error allocating buffer" );

			buffer.reset( bufPtr, &GenericRelease<INSSBuffer> );
		}

		BYTE* ppdwBuffer;
		hr = buffer->GetBuffer( &ppdwBuffer );
		if( FAILED(hr) )
			throw Exception( "Error retrieving buffer" );

		mainFrame_->time_ = timeRange_.first + boost::numeric_cast<double>(currentVideoFrame_) 
			/ boost::numeric_cast<double>(videoFrameRate_);
		mainFrame_->timeView_->Refresh(FALSE);
		mainFrame_->timeView_->Update();

		mainFrame_->currentCanvas_->dumpFrame( this->width_, this->height_, ppdwBuffer );

		buffer->SetLength( sampleSize );
		DWORD frameTime = 
			boost::numeric_cast<double>(currentVideoFrame_)
				/ (boost::numeric_cast<double>(videoFrameRate_)*1.0e-7);
		hr = asfWriter_->WriteSample( videoConnectionId_,
			frameTime,
			0,
			buffer.get() );
		if( FAILED(hr) )
			throw Exception( "Error writing sample" );

		++currentVideoFrame_;
	}

	bool complete() const
	{
		return interrupted_;
	}

	void interrupt()
	{
		this->interrupted_ = true;
	}

	~RecordMovieIdleAction()
	{
		mainFrame_->simulation_ = this->oldSimulation_;
		mainFrame_->time_ = this->oldTime_;
		mainFrame_->prevHighestAudio_ = 
			std::max(this->maxSeenAudioPressure_, 0.5*mainFrame_->prevHighestAudio_);

		glPixelStorei( GL_PACK_ALIGNMENT, 1 );

		// done writing asf file
		HRESULT hr = asfWriter_->EndWriting();
		if( FAILED(hr) )
			throw Exception( "EndWriting failed" ); 

		try
		{
			runMediaPlayer( "test.asf" );
		}
		catch( Exception& )
		{
			// ignore; if we can't run media player it's not
			//   the end of the world
		}

		mainFrame_->currentCanvas_->Refresh(FALSE);
	}

private:
	MainFrame* mainFrame_;
	boost::shared_ptr<IWMWriter> asfWriter_;
	SimulationPtr oldSimulation_;
	float oldTime_;

	double origLastAudioValue_;
	double dcRemovedLastAudioValue_;
	size_t currentAudioFrame_;
	size_t currentVideoFrame_;
	double maxSeenAudioPressure_;
	double maxAudioPressure_;
	bool interrupted_;

	size_t width_;
	size_t height_;
	size_t stride_;

	size_t videoFrameRate_;
	size_t audioFrameRate_;

	DWORD videoConnectionId_;
	DWORD audioConnectionId_;

	// need to capture the mouse so that we don't lose focus
	MouseCapture mouseCapture_;

	std::pair<double, double> timeRange_;

	boost::shared_ptr<Simulation> simulation_;
};

void MainFrame::recordMovie()
{
	if( this->idleAction_ )
		return;

	// todo make these user-specifiable
	size_t width = 360;
	size_t height = 240;
	size_t bitRate = 4000000;

	this->idleAction_.reset( 
		new RecordMovieIdleAction( 
			this, 
			videoStreamConfigurations_.at(selectedVideoCodec_),
			audioStreamConfigurations_.at(selectedAudioCodec_),
			width,
			height,
			bitRate,
			30,
			44100 ) );
}
#endif // WINDOWS_MEDIA

void MainFrame::OnExportImageSequence(wxCommandEvent& event)
{
	boost::shared_ptr<wxFileDialog> dlg( 
		new wxFileDialog(
			this, 
			"Export image sequence",
			"", 
			"", 
			"PNG Files(*.png)|*.png|JPEG Files(*.jpg)|*.jpg|All files(*.*)|*.*",
			wxFD_SAVE, 
			wxDefaultPosition),
		boost::mem_fn(&wxFileDialog::Destroy) );

	if ( dlg.get()->ShowModal() == wxID_OK )
	{
		boost::shared_ptr<wxTextEntryDialog> textDlg;

		int width = 720;
		textDlg.reset( new wxTextEntryDialog(this, 
			"Width of image", "Width of image to export"),
			boost::mem_fn(&wxTextEntryDialog::Destroy) );
		do
		{
			std::ostringstream oss;
			oss << width;
			textDlg.get()->SetValue( wxString(oss.str().c_str()) );

			if( textDlg.get()->ShowModal() != wxID_OK )
				return;
			width = atoi( textDlg.get()->GetValue().GetData() );
		} while( width <= 0 || width >= 20000 );

		int height = 480;
		textDlg.reset( new wxTextEntryDialog(this, 
			"Height of image", "Height of image to export"),
			boost::mem_fn(&wxTextEntryDialog::Destroy) );
		{
			std::ostringstream oss;
			oss << height;
			textDlg.get()->SetValue( wxString(oss.str().c_str()) );

			if( textDlg.get()->ShowModal() != wxID_OK )
				return;
			textDlg.get()->Destroy();
			height = atoi( textDlg.get()->GetValue().GetData() );
		} while( height <= 0 || height >= 20000 );

		int frames = 0;
		textDlg.reset( new wxTextEntryDialog(this, 
			"Number of frames", "Number of frames to export"),
			boost::mem_fn(&wxTextEntryDialog::Destroy) );
		do
		{
			std::ostringstream oss;
			oss << frames;
			textDlg.get()->SetValue( wxString(oss.str().c_str()) );

			if( textDlg.get()->ShowModal() != wxID_OK )
				return;
			frames = atoi( textDlg.get()->GetValue().GetData() );
		} while( frames <= 0 );

		float oldTime = time_;

		try
		{
			std::string filename = dlg.get()->GetPath().GetData();
			std::string::size_type lastPeriod = filename.find_last_of( '.' );
			if( lastPeriod == std::string::npos )
				throw IOException( std::string("Invalid filename: ") + filename );

			for (int iFrame = 0; iFrame < frames; ++iFrame)
			{
				std::ostringstream number;
				number << ".";
				number.width(4);
				number.fill('0');
				number << iFrame;

				this->time_ = boost::numeric_cast<double>(iFrame) / 30.0;
				std::string currentFrameFilename = filename;
				currentFrameFilename.insert( lastPeriod, "." + number.str() );
				this->currentCanvas_->dumpFrame( width, height, currentFrameFilename );
			}
			currentCanvas_->Refresh(FALSE);

			boost::shared_ptr<wxMessageDialog> 
				dlg( new wxMessageDialog(this, "Image export", "Image export complete.", wxOK ),
				boost::mem_fn(&wxMessageDialog::Destroy) );
			dlg.get()->ShowModal();
		}
		catch( IOException& e )
		{
			errMsg( e.message() );
		}
		time_ = oldTime;
	}
}

void MainFrame::OnMenuFileDumpGenerationImages(wxCommandEvent& event)
{
	boost::shared_ptr<SimulationTree> tree = 
		this->getTree();
	if( !tree )
		return;

	boost::shared_ptr<wxFileDialog> dlg( 
		new wxFileDialog(
			this, 
			"Export image sequence",
			"", 
			"", 
			"PNG Files(*.png)|*.png|JPEG Files(*.jpg)|*.jpg|All files(*.*)|*.*",
			wxFD_SAVE, 
			wxDefaultPosition),
		boost::mem_fn(&wxFileDialog::Destroy) );

	if ( dlg.get()->ShowModal() == wxID_OK )
	{
		boost::shared_ptr<wxTextEntryDialog> textDlg;

		int width = 720;
		textDlg.reset( new wxTextEntryDialog(this, 
			"Width of image", "Width of image to export"),
			boost::mem_fn(&wxTextEntryDialog::Destroy) );
		do
		{
			std::ostringstream oss;
			oss << width;
			textDlg.get()->SetValue( wxString(oss.str().c_str()) );

			if( textDlg.get()->ShowModal() != wxID_OK )
				return;
			width = atoi( textDlg.get()->GetValue().GetData() );
		} while( width <= 0 || width >= 20000 );

		int height = 480;
		textDlg.reset( new wxTextEntryDialog(this, 
			"Height of image", "Height of image to export"),
			boost::mem_fn(&wxTextEntryDialog::Destroy) );
		{
			std::ostringstream oss;
			oss << height;
			textDlg.get()->SetValue( wxString(oss.str().c_str()) );

			if( textDlg.get()->ShowModal() != wxID_OK )
				return;
			textDlg.get()->Destroy();
			height = atoi( textDlg.get()->GetValue().GetData() );
		} while( height <= 0 || height >= 20000 );

		std::string filename = dlg.get()->GetPath().GetData();
		std::string::size_type lastPeriod = filename.find_last_of( '.' );
		if( lastPeriod == std::string::npos )
			throw IOException( std::string("Invalid filename: ") + filename );

		const Path* backupPath = this->getPath();
		float oldTime = time_;

		size_t iGlobalFrame = 0;
		boost::dynamic_bitset<> activePaths = tree->activePaths();
		if( activePaths.empty() )
			return;

		// todo back these up
		lineCache_.clear();
		lineCacheColors_.clear();
		lineCacheStarts_.clear();
		lineCacheCounts_.clear();
		selectedLineCache_.clear();

		{
			const Path* path = tree->path( 0 );
			lineCache_.resize( path->treeCount() );
			lineCacheColors_.resize( lineCache_.size() );
			lineCacheStarts_.resize( lineCache_.size() );
			lineCacheCounts_.resize( lineCache_.size() );
			selectedLineCache_.resize( lineCache_.size() );
		}


		std::deque<SceneGraphElementPtr> prevSelected = this->getSelected();
		boost::mutex::scoped_lock lock( this->selectedLock_ );

		bool prevWireframeOnSelected = wireframeOnSelected_;
		wireframeOnSelected_ = false;

		try
		{
			const size_t frameRate = tree->path( 0 )->frameRate();
			size_t skip = frameRate / 30;

			for( size_t iPath = 0; iPath < activePaths.size(); ++iPath )
			{
				if( iPath >= 2 )
					skip = frameRate * 2;

				if( !activePaths[iPath] )
					continue;

				// block the computeHistoryThread from doing anything

				std::vector<PiecewisePath> piecewisePaths;
				const Path* path = tree->path( iPath );
				{
					{
						boost::mutex::scoped_lock lock( this->currentPathMutex_ );
						this->currentPath_ = path;
					}

					{
						ReaderWriterLock::ScopedWriterLock lock( this->historyLock_ );
						piecewisePaths = tree->decode( path );
						this->piecewisePaths_ = piecewisePaths;
						this->piecewisePathsFrameRate_ = path->frameRate();
					}

					stateHistory_.clear();
					std::vector<Simulation::State> states = path->samples();
					for( std::vector<Simulation::State>::const_iterator iter = states.begin();
						iter != states.end(); ++iter )
					{
						stateHistory_.push_back( StateWithNode(*iter, 0) );
					}
				}

				selectedLineCache_.clear();
				selectedLineCache_.resize( lineCache_.size() );

                PiecewisePath::FrameType endFrame = boost::numeric::bounds<PiecewisePath::FrameType>::lowest();
				PiecewisePath::FrameType startFrame = boost::numeric::bounds<PiecewisePath::FrameType>::highest();
				for( std::vector<PiecewisePath>::const_iterator iter = piecewisePaths.begin();
					iter != piecewisePaths.end(); ++iter )
				{
					endFrame = std::max( endFrame, iter->endFrame() );
					startFrame = std::min( startFrame, iter->startFrame() );
				}

				const float scale = this->lineCacheStream_.uniform( 0.4f, 1.0f );
				vl::Vec3f lineColor(0.0f/255.0f, 108.0f/255.0f, 255.0f/255.0f);


				for( size_t iObject = 0; iObject < lineCache_.size(); ++iObject )
				{
					lineCacheStarts_.at(iObject).push_back( lineCache_.at(iObject).size() );
					lineCacheCounts_.at(iObject).push_back( 0 );
				}

				size_t prevFrame = 0;
				for( size_t iFrame = startFrame; iFrame < endFrame; iFrame += skip )
				{
					double time = boost::numeric_cast<double>(iFrame)
						/ boost::numeric_cast<double>(frameRate);
					this->time_ = time;

					for( size_t curFrame = prevFrame; curFrame <= iFrame; ++curFrame )
					{
						for( size_t iObject = 0; iObject < path->treeCount(); ++iObject )
						{
							lineCache_.at(iObject).push_back( 
								piecewisePaths.at(iObject).position(curFrame) );
							lineCacheColors_.at(iObject).push_back( scale*lineColor );
							lineCacheCounts_.at(iObject).back() += 1;
						}
					}
					prevFrame = iFrame + 1;

					std::ostringstream number;
					number << ".";
					number.width(4);
					number.fill('0');
					number << iGlobalFrame++;

					std::string currentFrameFilename = filename;
					currentFrameFilename.insert( lastPeriod, number.str() );
					this->currentCanvas_->dumpFrame( width, height, currentFrameFilename );
				}
			}
			this->setPath( tree, backupPath );
			time_ = oldTime;

			currentCanvas_->Refresh(FALSE);

			boost::shared_ptr<wxMessageDialog> 
				dlg( new wxMessageDialog(this, "Image export", "Image export complete.", wxOK ),
				boost::mem_fn(&wxMessageDialog::Destroy) );
			dlg.get()->ShowModal();
		}
		catch( IOException& e )
		{
			errMsg( e.message() );
		}

		wireframeOnSelected_ = prevWireframeOnSelected;
	}
}

#ifdef DRAW_LINES
class PathWriter
	: boost::noncopyable
{
public:
	PathWriter( boost::shared_ptr<cairo_surface_t> surface, const vl::Mat4f& transform )
		: transform_(transform), surface_( surface )
	{
		width_ = cairo_image_surface_get_width( surface_.get() );
		height_ = cairo_image_surface_get_height( surface_.get() );

		cr_ = cairo_create( surface_.get() );

		/*
		cairo_set_font_size (cr_, 10.0);

		drawCoordinateText( 10.0, 10.0 );
		drawCoordinateText( 10.0, height_ - 10.0 );
		drawCoordinateText( width_ - 50.0, height_ - 10.0 );
		drawCoordinateText( width_ - 50.0, 10.0 );
		drawCoordinateText( 0.5*width_, 0.5*height_ );
		*/

		assert( cairo_status(cr_) == CAIRO_STATUS_SUCCESS );
	}

	~PathWriter()
	{
		cairo_destroy(cr_);
	}

	void writePath( const std::vector<vl::Vec3f>& path, const vl::Vec4f& color, float thickness )
	{
		assert( thickness > 0.0f );

		if( path.empty() )
			return;

		cairo_new_path( cr_ );

		vl::Vec2d firstPoint = toLocalCoords( path.front() );
		cairo_move_to (cr_, firstPoint[0], firstPoint[1]);
		for( std::vector<vl::Vec3f>::const_iterator iter = path.begin() + 1;
			iter != path.end(); ++iter )
		{
			vl::Vec2d curPoint = toLocalCoords( *iter );
			cairo_line_to (cr_, curPoint[0], curPoint[1]);
		}

		cairo_set_line_width( cr_, thickness );
		cairo_set_source_rgba(cr_, color[0], color[1], color[2], color[3]);

		cairo_stroke (cr_);
		assert( cairo_status(cr_) == CAIRO_STATUS_SUCCESS );
	}

private:
	void drawCoordinateText( double x, double y )
	{
		std::ostringstream oss;
		oss << "(" << x << ", " << y << ")";
		std::string str = oss.str();

		cairo_move_to (cr_, x, y);
		cairo_show_text (cr_, str.c_str());
	}

	vl::Vec2d toLocalCoords( const vl::Vec3f& point )
	{
		vl::Vec3f transformed = vl::xform( transform_, point );
		vl::Vec2d result( (0.5*transformed[0] + 0.5)*static_cast<double>(width_),
			(0.5 - 0.5*transformed[1])*static_cast<double>(height_) );
		return result;
	}
	
	vl::Mat4f transform_;

	int width_;
	int height_;

	boost::shared_ptr<cairo_surface_t> surface_;
	cairo_t* cr_;
};

void MainFrame::OnMenuFileDumpPathsImages(wxCommandEvent& event)
{
	std::string pathsFilename;
	{
		boost::shared_ptr<wxFileDialog> dlg( 
			new wxFileDialog(
				this, 
				"Load paths file",
				"", 
				"", 
				"bin Files(*.bin)|*.bin|All files(*.*)|*.*",
				wxFD_OPEN, 
				wxDefaultPosition),
			boost::mem_fn(&wxFileDialog::Destroy) );
		if ( dlg.get()->ShowModal() != wxID_OK )
			return;

		pathsFilename = dlg.get()->GetPath().GetData();
	}

	std::string imagesFilename;
	{
		boost::shared_ptr<wxFileDialog> dlg( 
			new wxFileDialog(
				this, 
				"Export image sequence",
				"", 
				"", 
				"PNG Files(*.png)|*.png|JPEG Files(*.jpg)|*.jpg|All files(*.*)|*.*",
				wxFD_SAVE, 
				wxDefaultPosition),
			boost::mem_fn(&wxFileDialog::Destroy) );
		if ( dlg.get()->ShowModal() != wxID_OK )
			return;

		imagesFilename = dlg.get()->GetPath().GetData();
	}

	int width = 720;
	int height = 480;

	{
		boost::shared_ptr<wxTextEntryDialog> textDlg;

		textDlg.reset( new wxTextEntryDialog(this, 
			"Width of image", "Width of image to export"),
			boost::mem_fn(&wxTextEntryDialog::Destroy) );
		do
		{
			std::ostringstream oss;
			oss << width;
			textDlg.get()->SetValue( wxString(oss.str().c_str()) );

			if( textDlg.get()->ShowModal() != wxID_OK )
				return;
			width = atoi( textDlg.get()->GetValue().GetData() );
		} while( width <= 0 || width >= 20000 );

		textDlg.reset( new wxTextEntryDialog(this, 
			"Height of image", "Height of image to export"),
			boost::mem_fn(&wxTextEntryDialog::Destroy) );
		{
			std::ostringstream oss;
			oss << height;
			textDlg.get()->SetValue( wxString(oss.str().c_str()) );

			if( textDlg.get()->ShowModal() != wxID_OK )
				return;
			textDlg.get()->Destroy();
			height = atoi( textDlg.get()->GetValue().GetData() );
		} while( height <= 0 || height >= 20000 );
	}

	boost::shared_ptr<SimulationFactory> factory = simulationFactory();
	boost::shared_ptr<Simulation> oldSimulation = simulation_;
	simulation_ = factory->construct( this->scene() );

	try
	{
		std::string::size_type lastPeriod = imagesFilename.find_last_of( '.' );
		if( lastPeriod == std::string::npos )
		{
			errMsg( "Invalid filename: " + imagesFilename );
			return;
		}

		std::string prefix = imagesFilename.substr(0, lastPeriod);
		if( prefix.empty() )
		{
			errMsg( "Invalid filename: " + imagesFilename );
			return;
		}

		// should also set the correct ratio:
		this->currentCanvas_->dumpFrame( width, height, prefix + ".background.png" );
		const Camera& camera = this->currentCanvas_->camera();
		vl::Mat4f transformMatrix = camera.projectiveTransform() * camera.modelViewTransform();

		std::ifstream ifs( pathsFilename.c_str(), std::ios::binary );
		if( !ifs )
			throw IOException( std::string("Unable to open file '") + pathsFilename + "' for reading." );

		std::ifstream::pos_type fileSize = filelen( ifs );

		typedef std::vector<vl::Vec3f> Path;
		typedef std::deque< Path > PathList;

		cairo_format_t format = CAIRO_FORMAT_ARGB32;

		typedef std::pair<PathList, size_t> PathWithTime;
		std::list< PathWithTime > activePaths;

		typedef std::pair< Simulation::State, size_t > StateWithTime;
		std::list< StateWithTime > activeStates;

		size_t iFrame = 0;

		size_t lengthBetweenPaths = 2;
		size_t lengthOfFade = 10;
		double fadePerFrame = 0.1;

		float thicknesses[] =  { 0.5, 1.0, 2.0, 3.0, 2.0, 1.5, 1.0, 1.0, 1.0, 1.0 };
		float colorWeights[] = { 0.8, 0.9, 1.0, 1.0, 1.0, 0.9, 0.7, 0.5, 0.3, 0.1 };
		assert( sizeof(thicknesses) == lengthOfFade * sizeof(float) );
		assert( sizeof(colorWeights) == lengthOfFade * sizeof(float) );

		vl::Vec3f initialColor( 1.0f, 1.0f, 0.5f );
		vl::Vec3f fadedColor( 0.9f, 0.9f, 0.9f );

		while( !activePaths.empty() || ifs.tellg() < fileSize )
		{
			if( (activePaths.empty() || activePaths.back().second > lengthBetweenPaths)
				&& ifs.tellg() < fileSize )
			{
				int numPaths;
				ifs.read( reinterpret_cast<char*>( &numPaths ), sizeof(int) );

				// introduce a new path
				activePaths.push_back( PathWithTime(PathList(), 0) );
				assert( numPaths > 0 );
				PathList& currentPathList = activePaths.back().first;

				for( int i = 0; i < numPaths; ++i )
				{
					currentPathList.push_back( std::vector<vl::Vec3f>() );
					readArray( currentPathList.back(), ifs );
				}

				int numSnaps;
				ifs.read( reinterpret_cast<char*>( &numSnaps ), sizeof(int) );

				for( int iSnap = 0; iSnap < numSnaps; ++iSnap )
				{
					float time;
					ifs.read( reinterpret_cast<char*>( &time ), sizeof(float) );

					std::vector<RigidStaticState> states;
					readArray(states, ifs);

					std::vector<RigidDynamicState> dynStates( states.begin(), states.end() );
					StateWithTime state( Simulation::State(dynStates, time), 0 );

					if( iSnap == 0 || (iSnap+1) == numSnaps /*|| iSnap == numSnaps/2*/ )
						activeStates.push_back( state );
				}

			}

			// start with background plate
			const int stride = ((width*4) + (0xF - 1)) & ~0xF;
			assert( stride >= (width*4) );
			std::vector<unsigned char> imageData( stride*height, 0 );

			boost::shared_ptr<cairo_surface_t> currentSurface(
				cairo_image_surface_create_for_data(&imageData[0], format, width, height, stride),
				&cairo_surface_destroy );


			for( unsigned int j = 0; j < height; ++j )
				for( unsigned int i = 0; i < width; ++i )
					imageData[ j*stride + i*4 + 3 ] = 255;

			{
				for( std::list< StateWithTime >::iterator stateItr = activeStates.begin();
					stateItr != activeStates.end(); )
				{
					size_t curTime = stateItr->second;
					std::vector<Simulation::State> states;
					std::list< StateWithTime >::iterator endOfRange = stateItr;
					while( endOfRange != activeStates.end() && endOfRange->second == curTime )
					{
						endOfRange->second += 1;
						states.push_back( endOfRange->first );
						++endOfRange;
					}

					if( curTime >= lengthOfFade )
						stateItr = activeStates.erase( stateItr, endOfRange );
					else
					{
						wxImage image = currentCanvas_->dumpFrame( width, height, states );

						/*
						std::ostringstream oss;
						oss << prefix << ".objects." << std::setw(4) << std::setfill('0') << iFrame << ".png";
						std::string filename = oss.str();
						if( !image.SaveFile( wxT(filename.c_str()) ) )
							throw IOException( "Unable to write to image '" + filename + "'." );
							*/

						unsigned char* curData = image.GetData();
						float colorWeight = 0.8*colorWeights[curTime];
						for( unsigned int j = 0; j < height; ++j )
						{
							for( unsigned int i = 0; i < width; ++i )
							{
								for( unsigned int component = 0; component < 3; ++component )
								{
									unsigned char& origDataRef = imageData[ j*stride + i*4 + (2 - component) ];
									unsigned char& curDataRef = curData[ (j*width + i)*3 + component ];
									origDataRef = std::max( origDataRef, 
										static_cast<unsigned char>( colorWeight*static_cast<float>(curDataRef) ) );
								}
							}
						}

						stateItr = endOfRange;
					}
				}
			}

			{
				PathWriter pathWriter( currentSurface, transformMatrix );

				for( std::list< PathWithTime >::iterator currentPathItr = activePaths.begin();
					currentPathItr != activePaths.end(); )
				{
					if( currentPathItr->second >= lengthOfFade )
					{
						size_t numSteps = currentPathItr->second - lengthOfFade;
						double alpha = 1.0 - static_cast<double>( numSteps )*fadePerFrame;
						if( alpha <= 0.0 )
						{
							currentPathItr = activePaths.erase( currentPathItr );
						}
						else
						{
							float thickness = thicknesses[lengthOfFade - 1];
							const PathList& paths = currentPathItr->first;
							std::for_each( paths.begin(), paths.end(), 
								boost::bind( &PathWriter::writePath, boost::ref(pathWriter), _1, vl::Vec4f(fadedColor, alpha), thickness ) );

							++currentPathItr->second;
							++currentPathItr;
						}
					}
					else
					{
						float colorWeight = colorWeights[currentPathItr->second];
						float thickness = thicknesses[currentPathItr->second];
						vl::Vec3f color = colorWeight*initialColor + (1.0f - colorWeight)*fadedColor;
						const PathList& paths = currentPathItr->first;
						std::for_each( paths.begin(), paths.end(), 
							boost::bind( &PathWriter::writePath, boost::ref(pathWriter), _1, vl::Vec4f(color, 1.0), thickness ) );

						++currentPathItr->second;
						++currentPathItr;
					}
				}
			}

			{
				std::ostringstream oss;
				oss << prefix << "." << std::setw(4) << std::setfill('0') << iFrame << ".png";
				std::string outFilename = oss.str();
				cairo_surface_write_to_png( currentSurface.get(), outFilename.c_str() );
			}

			++iFrame;
		}
	}
	catch( IOException& e )
	{
		errMsg( e.message() );
	}

	simulation_ = oldSimulation;
}
#endif

void clearMenu( wxMenu* menu )
{
	while( menu->GetMenuItemCount() > 0 )
	{
		wxMenuItem* item = menu->FindItemByPosition(0);
		menu->Remove( item );
	}
}

void MainFrame::OnScrollTime(wxScrollEvent& event)
{
//	assert( !simulating() );
	updateViews();
}

void MainFrame::OnKey(wxKeyEvent& event)
{
	if( idleAction_ )
	{
		idleAction_->interrupt();
		return;
	}

	currentCanvas_->OnKey( event );
}

void MainFrame::OnKeyUp(wxKeyEvent& event)
{
	currentCanvas_->OnKeyUp( event );
}

void MainFrame::errMsg(const std::string& message)
{
	this->log( "Caught error: " + message );
	planning::errMsg( message, this );
}

void MainFrame::updateViews()
{
	this->updated_ = false;
	std::for_each( canvases_.begin(), canvases_.end(), 
		boost::bind( &GLView::Refresh, _1, false, (const wxRect*) NULL ) );

	wxCommandEvent event(wxEVT_REPAINT_ALL_VIEWS, this->GetId());
	event.SetEventObject( this );
	::wxPostEvent( this, event );
}

bool MainFrame::isModifiedFromTree( ConstSceneGraphElementPtr object ) const
{
	DirtyObjectSet::const_iterator setItr = 
		dirtyObjects_.find( object.get() );
	return (setItr != dirtyObjects_.end());
}

void MainFrame::OnChangeObjectMode(wxCommandEvent& event)
{
	this->setSelected( selected_ );
	this->updateViews();
}

void MainFrame::setObjectMode( ObjectMode mode )
{
	ObjectMode prevMode = this->getObjectMode();

	switch( mode )
	{
	case MODE_OBJECT_SELECT:
		this->toolBar_->ToggleTool( myID_CURSOR_TOOLBAR, true );
		break;
	case MODE_OBJECT_TRANSLATE:
		this->toolBar_->ToggleTool( myID_TRANSLATE_TOOLBAR, true );
		break;
	case MODE_OBJECT_CHANGEPIVOT:
		this->toolBar_->ToggleTool( myID_CHANGEPIVOT_TOOLBAR, true );
		break;
	case MODE_OBJECT_ROTATE:
		this->toolBar_->ToggleTool( myID_ROTATE_TOOLBAR, true );
		break;
	case MODE_OBJECT_SCALE:
		this->toolBar_->ToggleTool( myID_SCALE_TOOLBAR, true );
		break;
	case MODE_OBJECT_SELECTPOINTS:
		this->toolBar_->ToggleTool( myID_SELECTPOINTS_TOOLBAR, true );
		break;
	}

	if( prevMode != mode )
		this->setSelected( selected_ );
}

MainFrame::ObjectMode MainFrame::getObjectMode() const
{
	if( this->toolBar_->GetToolState( myID_CURSOR_TOOLBAR ) )
		return MODE_OBJECT_SELECT;
	else if( this->toolBar_->GetToolState( myID_TRANSLATE_TOOLBAR ) )
		return MODE_OBJECT_TRANSLATE;
	else if( this->toolBar_->GetToolState( myID_CHANGEPIVOT_TOOLBAR ) )
		return MODE_OBJECT_CHANGEPIVOT;
	else if( this->toolBar_->GetToolState( myID_SELECTPOINTS_TOOLBAR ) )
		return MODE_OBJECT_SELECTPOINTS;
	else if( this->toolBar_->GetToolState( myID_ROTATE_TOOLBAR ) )
		return MODE_OBJECT_ROTATE;
	else if( this->toolBar_->GetToolState( myID_SCALE_TOOLBAR ) )
		return MODE_OBJECT_SCALE;
	else
	{
		assert( false );
		return MODE_OBJECT_SELECT;
	}
}

std::deque<SceneGraphElementPtr> MainFrame::getSelected() const
{
	boost::mutex::scoped_lock lock( this->selectedLock_ );
	std::deque<SceneGraphElementPtr> result = this->selected_;
	return result;
}

ConstraintPtr MainFrame::selectedConstraint() const
{
	return this->selectedConstraint_;
}

void MainFrame::setSelectedConstraint( ConstraintPtr constraint )
{
	bool changed = selectedConstraint_ != constraint;
	this->selectedConstraint_ = constraint;
	if( changed )
		this->changeSelected();
}

class HashOrdering
{
	typedef STLEXT hash_map<const SceneGraphElement*, size_t, hash_ptr<const SceneGraphElement*> > HashType;
	const HashType& hash_;
public:
	HashOrdering( const HashType& hash )
		: hash_(hash) {}

	bool operator()( SceneGraphElementPtr left, SceneGraphElementPtr right ) const
	{
		HashType::const_iterator firstIter = 
			hash_.find( left.get() );
		assert( firstIter != hash_.end() );
		
		HashType::const_iterator secondIter = 
			hash_.find( right.get() );
		assert( secondIter != hash_.end() );

		std::less<size_t> comparator;
		return comparator( firstIter->second, secondIter->second );
	}
};

std::deque<SceneGraphElementPtr> cleanSelected( const std::deque<SceneGraphElementPtr>& selected )
{
	STLEXT hash_map<const SceneGraphElement*, size_t, hash_ptr<const SceneGraphElement*> > orderMap;
	for( size_t i = 0; i < selected.size(); ++i )
	{
		orderMap.insert( std::make_pair(selected.at(i).get(), i) );
	}

	std::deque<SceneGraphElementPtr> result;
	// remove anything whose parent is already in the set.
	for( std::deque<SceneGraphElementPtr>::const_iterator selectedIter = selected.begin();
		selectedIter != selected.end(); ++selectedIter )
	{
		SceneGraphElementPtr current = (*selectedIter)->parent().lock();
		while( current && orderMap.find(current.get()) == orderMap.end() )
			current = current->parent().lock();

		if( !current )
			result.push_back( *selectedIter );
	}

	// remove dups
	std::sort( result.begin(), result.end() );
	result.erase( std::unique(result.begin(), result.end()), result.end() );

	std::sort( result.begin(), result.end(), 
		HashOrdering(orderMap) );
	return result;
}


void MainFrame::upHierarchy()
{
	std::deque<SceneGraphElementPtr> newSelected;
	for( std::deque<SceneGraphElementPtr>::const_iterator objectIter = selected_.begin();
		objectIter != selected_.end(); ++objectIter )
	{
		SceneGraphElementPtr parent = (*objectIter)->parent().lock();

		// parent should never be 0 as everything is a child of the root:
		assert( parent );

		if( parent == scene_->root() )
			newSelected.push_back( *objectIter );
		else
			newSelected.push_back( parent );
	}

	ActionPtr changeSelectionAction( 
		new ChangeSelectionAction(
			this,
			selected_,
			cleanSelected(newSelected),
			selectedConstraint(), ConstraintPtr() ) );
	this->performAction( changeSelectionAction );
}

void MainFrame::downHierarchy()
{
	std::deque<SceneGraphElementPtr> newSelected;
	for( std::deque<SceneGraphElementPtr>::const_iterator objectIter = selected_.begin();
		objectIter != selected_.end(); ++objectIter )
	{
		if( (*objectIter)->children_begin() == (*objectIter)->children_end() )
			newSelected.push_back( *objectIter );
		else
			newSelected.push_back( *((*objectIter)->children_begin()) );
	}

	ActionPtr changeSelectionAction( 
		new ChangeSelectionAction(
			this,
			selected_,
			cleanSelected(newSelected),
			selectedConstraint(), ConstraintPtr() ) );
	this->performAction( changeSelectionAction );
}


std::deque<SceneGraphElementPtr> MainFrame::convertSelected( 
	const std::deque<SceneGraphElementPtr>& selected )
{
	std::deque<SceneGraphElementPtr> result;
	for( std::deque<SceneGraphElementPtr>::const_iterator itr = selected.begin();
		itr != selected.end(); ++itr )
	{
		SceneGraphElementPtr current = *itr;
		if( this->editSelectByMenu_->IsChecked(MENU_EDIT_SELECT_BY_HIERARCHY) )
		{
			while( current->parent().lock() != scene_->root() )
				current = current->parent().lock();

			result.push_back( current );
		}
		else if( this->editSelectByMenu_->IsChecked(MENU_EDIT_SELECT_BY_SIMULATED) )
		{
			SceneGraphElementPtr tmp = current;
			while( tmp->parent().lock() )
			{
				tmp = tmp->parent().lock();
				if( boost::dynamic_pointer_cast<CombinedObject>( tmp ) )
					current = tmp;
			}

			if( boost::dynamic_pointer_cast<PhysicsObject>( current ) )
				result.push_back( current );
		}
		else
		{
			result.push_back( current );
		}
	}

	// we don't worry about duplicates since the caller is going to call "cleanSelected" at some point
	return result;
}

void MainFrame::setSelected( std::deque<SceneGraphElementPtr> selected )
{
	{
		boost::mutex::scoped_lock lock( this->selectedLock_ );
		selected_.swap( selected );
	}

	std::vector<TransformedElementPtr> transformedSelected;
	transformedSelected.reserve( selected_.size() );

	for( std::deque<SceneGraphElementPtr>::const_iterator iter = selected_.begin();
		iter != selected_.end(); ++iter )
	{
		TransformedElementPtr transformed = 
			boost::dynamic_pointer_cast<TransformedSceneGraphElement>( *iter );
		if( transformed )
			transformedSelected.push_back( transformed );
	}

	mouseAction_.reset();
	if( !transformedSelected.empty() && this->mode_ == MODE_EDIT )
	{
		switch( this->getObjectMode() )
		{
		case MODE_OBJECT_SELECT:
		case MODE_OBJECT_SELECTPOINTS:
			break;
		case MODE_OBJECT_TRANSLATE:
			mouseAction_.reset( new TranslateMousingAction(transformedSelected, false) );
			break;
		case MODE_OBJECT_CHANGEPIVOT:
			mouseAction_.reset( new TranslateMousingAction(transformedSelected, true) );
			break;
		case MODE_OBJECT_ROTATE:
			mouseAction_.reset( new RotateMousingAction(transformedSelected) );
			break;
		case MODE_OBJECT_SCALE:
			mouseAction_.reset( new ScaleMousingAction(transformedSelected) );
			break;
		}
	}

	if( selected_ == selected )
		return;

	this->changeSelected();
}

void MainFrame::setSelectedPoints( const MainFrame::SelectedPoints& selected )
{
	this->selectedPoints_ = selected;
}

std::pair<SceneGraphElementPtr, std::deque<unsigned int> > MainFrame::getSelectedPoints() const
{
	return this->selectedPoints_;
}

void MainFrame::updateSelectedPath()
{
	const Path* path = this->getPath();
	std::deque< LineCache > lineCache;
	if( path )
	{
		lineCache.resize( path->treeCount() );
		for( size_t iObject = 0; iObject < path->treeCount(); ++iObject )
		{
			std::vector<vl::Vec3f> points;
			path->obbTree(iObject).points( points );
			for( size_t j = 1; j < points.size(); ++j )
			{
				lineCache[iObject].push_back( points[j-1] );
				lineCache[iObject].push_back( points[j] );
			}
		}
	}

	{
		ReaderWriterLock::ScopedWriterLock lock( this->selectedLineCacheMutex_ );
		this->selectedLineCache_.swap( lineCache );
	}

	this->Refresh(FALSE);
}

void MainFrame::clearSamples()
{
	ReaderWriterLock::ScopedWriterLock lineCacheLock( this->lineCacheMutex_ );
	this->lineCache_.clear();
	this->lineCacheColors_.clear();
	this->lineCacheStarts_.clear();
	this->lineCacheCounts_.clear();

	// start out with the path roots from the tree
	SimulationTreePtr tree = this->getTree();
	if( tree )
	{
		size_t numDynOb = tree->dynamicObjects().size();
		const vl::Vec3f lineColor(0.3f, 0.3f, 0.3f);

		lineCache_.resize( numDynOb );
		lineCacheColors_.resize( lineCache_.size() );
		lineCacheStarts_.resize( lineCache_.size() );
		lineCacheCounts_.resize( lineCache_.size() );

		for( size_t iObject = 0; iObject < numDynOb; ++iObject )
		{
			const size_t prevSize = this->lineCache_[iObject].size();

			std::vector<vl::Vec3f> root = tree->pathRoot(iObject);
			if( root.empty() )
				continue;

			std::copy( root.begin(), root.end(), 
				std::back_inserter( this->lineCache_[iObject] ) );

			const size_t newSize = lineCache_[iObject].size();
			lineCacheColors_.at(iObject).resize( lineCache_[iObject].size(), lineColor );

			if( (lineCache_[iObject].size() - prevSize) == 0 )
				continue;

			lineCacheStarts_.at(iObject).push_back( prevSize );
			lineCacheCounts_.at(iObject).push_back( lineCache_[iObject].size() - prevSize );
		}
	}
}

void MainFrame::setSamples( ConstSimulationTreePtr tree, const boost::dynamic_bitset<>& samples )
{
	this->clearSamples();
	ReaderWriterLock::ScopedWriterLock lineCacheLock( this->lineCacheMutex_ );

	bool showAllPaths = this->showAllPaths();

	for( size_t i = 0; i < samples.size(); ++i )
	{
		bool active = samples[i];
		if( !active && !showAllPaths )
			continue;

		const Path* path = tree->path( i );

		this->addToLineCache( path, active );
	}
}

void MainFrame::addToLineCache( const Path* path, bool active )
{
	const float scale = this->lineCacheStream_.uniform( 0.4f, 1.0f );
	vl::Vec3f lineColor(0.0f/255.0f, 108.0f/255.0f, 255.0f/255.0f);
	if( !active )
		lineColor = vl::Vec3f( 0.3, 0.3, 0.3 );

	lineCache_.resize( path->treeCount() );
	lineCacheColors_.resize( lineCache_.size() );
	lineCacheStarts_.resize( lineCache_.size() );
	lineCacheCounts_.resize( lineCache_.size() );

	assert( path->treeCount() == lineCache_.size() );

	for( size_t iObject = 0; iObject < path->treeCount(); ++iObject )
	{
		const size_t prevSize = this->lineCache_[iObject].size();
		path->obbTree( iObject ).points( lineCache_[iObject] );

		/*
#ifdef _DEBUG
		{
			std::vector<vl::Vec3f> points1 = path->tree( iObject ).points();
			std::vector<vl::Vec3f> points2;
			path->obbTree( iObject ).points( points2 );
			assert( points1 == points2 );
		}
#endif
		*/

		const size_t newSize = lineCache_[iObject].size();
		lineCacheColors_.at(iObject).resize( lineCache_[iObject].size(), scale * lineColor );

		if( (lineCache_[iObject].size() - prevSize) == 0 )
			continue;

		lineCacheStarts_.at(iObject).push_back( prevSize );
		lineCacheCounts_.at(iObject).push_back( lineCache_[iObject].size() - prevSize );
	}
}

boost::shared_ptr<SimulationFactory> MainFrame::simulationFactory() const
{
	switch( this->simulatorType() )
	{
	case SIMULATOR_NOVODEX:
#ifdef USE_NOVODEX
		return boost::shared_ptr<SimulationFactory>(
			new NovodexSimulationFactory(this->physicsSDK_) );
#endif

	case SIMULATOR_NEWTON:
#ifdef USE_NEWTON
		return boost::shared_ptr<SimulationFactory>(
			new SimulationFactoryInst<NewtonSimulation> );
#endif

	case SIMULATOR_BULLET:
#ifdef USE_BULLET
		return boost::shared_ptr<SimulationFactory>(
			new SimulationFactoryInst<BulletSimulation> );
#endif

	case SIMULATOR_ODE_BULLET:
#ifdef USE_BULLET
		return boost::shared_ptr<SimulationFactory>(
			new SimulationFactoryInst<ODEBulletSimulation> );
#endif

    case SIMULATOR_ODE_PENALTY:
		return boost::shared_ptr<SimulationFactory>(
			new SimulationFactoryInst<PenaltyForceSimulation> );

    case SIMULATOR_ODE_NATIVE:
	default:
		return boost::shared_ptr<SimulationFactory>(
			new SimulationFactoryInst<ODENativeSimulation> );
	}
}

void MainFrame::startSimulation()
{
	if( !simulation_ )
	{
		boost::shared_ptr<SimulationFactory> factory = simulationFactory();
		simulation_ = factory->construct( this->scene() );
	}
}


const std::deque<MainFrame::LineCache>& MainFrame::lineCache() const
{
	return this->lineCache_;
}

const std::deque<MainFrame::LineCache>& MainFrame::lineCacheColors() const
{
	return this->lineCacheColors_;
}

const std::deque<MainFrame::LineCacheIndex>& MainFrame::lineCacheStarts() const
{
	return this->lineCacheStarts_;
}

const std::deque<MainFrame::LineCacheCount>& MainFrame::lineCacheCounts() const
{
	return this->lineCacheCounts_;
}

const std::deque<MainFrame::LineCache>& MainFrame::selectedLineCache() const
{
	return this->selectedLineCache_;
}

ReaderWriterLock& MainFrame::lineCacheMutex() const
{
	return this->lineCacheMutex_;
}

ReaderWriterLock& MainFrame::selectedLineCacheMutex() const
{
	return this->selectedLineCacheMutex_;
}

boost::shared_ptr<MousingAction> MainFrame::mouseAction()
{
	if( simulation_ )
		return boost::shared_ptr<MousingAction>();
	else
		return this->mouseAction_;
}

boost::shared_ptr<const MousingAction> MainFrame::mouseAction() const
{
	if( simulation_ )
		return boost::shared_ptr<const MousingAction>();
	else
		return this->mouseAction_;
}

void MainFrame::stopSimulation()
{
	if( simulation_ )
	{
		simulation_.reset();
		this->Refresh(FALSE);
	}
}

bool MainFrame::isSimulating() const
{
	return this->simulation_;
}


MousingAction::MousingAction( const std::vector<TransformedElementPtr>& objects )
	: objects_(objects), screenRadius_(100.0), moving_(false)
{
	assert( !objects_.empty() );

	vl::Vec3f center = vl::xform( objects_.back()->transform(), vl::Vec3f(objects_.back()->pivot()) );
	axisCenter_ = toVec3d( center );
}

bool MousingAction::pressMouse( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter )
{
	this->moving_ = false;
	if( select( screenPos, converter ) )
	{
		moving_ = true;
		return true;
	}
	else
		return false;
}

void MousingAction::moveMouse( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter )
{
	if( !moving_ )
		return;

	move( screenPos, converter );
}

void MousingAction::releaseMouse( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter )
{
	if( !moving_ )
		return;

	move( screenPos, converter );
	moving_ = false;
}

bool MousingAction::moving() const
{
	return moving_;
}

MousingAction::~MousingAction()
{
}

double MousingAction::getAxisScale(const ScreenSpaceConverter& converter) const
{
	vl::Vec3d screenSpace = converter.toScreenSpace( vl::vl_0 );
	screenSpace[1] += screenRadius_;
	vl::Vec3d offset = converter.toWorldSpace( screenSpace );
	return len( offset );
}

ActionPtr MousingAction::action()
{
	ActionPtr result = getAction();
	return result;
}

unsigned int MousingAction::skipAxis(const ScreenSpaceConverter& converter) const
{
	// this is the minimum distance in screen space between the root and a particular
	// axis endpoint
	const double minScreenDist = 5.0;

	unsigned int result = 10;
	{
		vl::Vec2d screenRoot = vl::strip( converter.toScreenSpace( axisCenter_ ) );

		vl::Vec3f rotation = objects_.back()->rotate();
		const Quaternion<float> quat = eulerToQuat(rotation);
		const vl::Mat3d mat = quat.toRotMatd();

		vl::Mat4d transform( vl::vl_1 );
		for( vl::Int i = 0; i < 3; ++i )
			for( vl::Int j = 0; j < 3; ++j )
				transform[i][j] = mat[i][j];

		for( unsigned int iAxis = 0; iAxis < 3; ++iAxis )
		{
			vl::Vec3d axis( vl::vl_0 );
			axis[iAxis] = 100.0f;
			vl::Vec2d screenPosAxis = vl::strip( converter.toScreenSpace( axisCenter_ + xform(transform, axis)) );
			if( vl::len( screenPosAxis - screenRoot ) < minScreenDist )
				result = iAxis;
		}
	}

	return result;
}

vl::Vec3d nearestPoint( vl::Vec2d point, vl::Vec3d center, const ScreenSpaceConverter& converter )
{
	vl::Vec3d nearVec = converter.toWorldSpace( vl::Vec3d( point, -1.0 ) );
	vl::Vec3d farVec = converter.toWorldSpace( vl::Vec3d( point, 1.0 ) );

	Wm4::Line3d viewSegment;
	viewSegment.Origin = toMagicVec( nearVec );
	viewSegment.Direction = toMagicVec( vl::norm(farVec - nearVec) );

	Wm4::DistVector3Line3d dist( toMagicVec( center ), viewSegment );
	dist.GetSquared();
	Wm4::Vector3d result = viewSegment.Origin + dist.GetLineParameter()*viewSegment.Direction;
	return vl::Vec3d( result[0], result[1], result[2] );
}


TranslateMousingAction::TranslateMousingAction( const std::vector<TransformedElementPtr>& objects, bool changePivot )
	: MousingAction( objects ), changePivot_(changePivot), translate_(vl::vl_0)
{
	this->backupTranslates_.reserve( objects.size() );
	std::transform( objects.begin(), objects.end(),
		std::back_inserter(backupTranslates_),
		boost::bind( &TransformedSceneGraphElement::translate, _1 ) );

	this->backupPivots_.reserve( objects.size() );
	std::transform( objects.begin(), objects.end(),
		std::back_inserter(backupPivots_),
		boost::bind( &TransformedSceneGraphElement::pivot, _1 ) );
}

TranslateMousingAction::~TranslateMousingAction()
{
}

vl::Vec3d TranslateMousingAction::projectPosition( const vl::Vec2d& pos, const ScreenSpaceConverter& converter ) const
{
	if( this->clickedAxis_ < 3 )
	{
		vl::Vec3d nearVec = converter.toWorldSpace( vl::Vec3d( pos, -1.0 ) );
		vl::Vec3d farVec = converter.toWorldSpace( vl::Vec3d( pos, 1.0 ) );

		Wm4::Ray3d viewSegment;
		viewSegment.Origin = toMagicVec( nearVec );
		viewSegment.Direction = toMagicVec( vl::norm(farVec - nearVec) );

		Wm4::Line3d axisLine;
		axisLine.Origin = toMagicVec( this->axisCenter_ );
		vl::Vec3d axisDir( vl::vl_0 );
		axisDir[ this->clickedAxis_ ] = 1.0;
		axisLine.Direction = toMagicVec( axisDir );

		Wm4::DistLine3Ray3d dist( axisLine, viewSegment );
		dist.GetSquared();
		Wm4::Vector3d result = dist.GetClosestPoint0();
		return vl::Vec3d( result[0], result[1], result[2] );
	}
	else
	{
		return nearestPoint( pos, this->axisCenter_, converter );
	}
}

bool TranslateMousingAction::select( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter )
{
	double axisScale = this->getAxisScale( converter );

	vl::Vec2d screenRoot = vl::strip( converter.toScreenSpace( axisCenter_ ) );
	unsigned int skipAxis = this->skipAxis( converter );

	if( len(screenPos - screenRoot) < 15 )
	{
		clickedPos_ = screenPos;
		clickedAxis_ = 10;
		return true;
	}

	double minDist = std::numeric_limits<double>::max();
	unsigned int minAxis = 0;
	for( unsigned int iAxis = 0; iAxis < 3; ++iAxis )
	{
		vl::Vec3d end = axisCenter_;
		end[iAxis] += axisScale;
		// project both ends of the axis into screen space:
		vl::Vec2d screenEnd = strip( converter.toScreenSpace( end ) );

		if( iAxis == skipAxis )
			continue;

		double dist = segmentDistance( screenPos, screenRoot, screenEnd );

		if( dist < minDist )
		{
			minDist = dist;
			minAxis = iAxis;
		} 
	}

	if( minDist < 10 )
	{
		clickedPos_ = screenPos;
		clickedAxis_ = minAxis;
		return true;
	}

	return false;
}

void TranslateMousingAction::move( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter )
{
	vl::Vec3d startPos = projectPosition( clickedPos_, converter );
	vl::Vec3d endPos = projectPosition( screenPos, converter );
	vl::Vec3d translate = endPos - startPos;

	translate_ = vl::vl_0;
	if( clickedAxis_ < 3 )
		translate_[clickedAxis_] = translate[clickedAxis_];
	else
		translate_ = translate;

	for( unsigned int i = 0; i < objects_.size(); ++i )
	{
		TransformedElementPtr object = objects_.at(i);
		SceneGraphElementPtr parent = object->parent().lock();

		if( changePivot_ )
		{
			object->setPivot( this->backupPivots_.at(i), true );
			object->setTranslate( this->backupTranslates_.at(i), true );

			vl::Mat4f transform = object->transform();
			vl::Mat4f parentTransform = parent->transform();

			vl::Vec3f origPivotWorld = vl::xform( transform, object->pivot() );
			vl::Vec3f newPivotWorld = origPivotWorld + toVec3f(translate_);

			vl::Vec3f newPivot = vl::xform( inv(transform), newPivotWorld ) - object->pivot();
			newPivot += this->backupPivots_.at(i);
			vl::Vec3f newTranslate = vl::xform( inv(parentTransform), newPivotWorld );

			object->setPivot( newPivot );
			object->setTranslate( newTranslate );

			static int iTest = 0;
			iTest++;
			if( iTest > 50 )
			{
				std::cout << "Test." << std::endl;
			}
		}
		else
		{
			object->setTranslate( 
				this->backupTranslates_.at(i) + 
					vl::xform( inv(toMat3f(parent->transform())), toVec3f(translate_) ), true );
		}
	}
}

ActionPtr TranslateMousingAction::getAction() const
{
	if( len(translate_) < 1e-5 )
		return ActionPtr();

	std::vector<vl::Vec3f> newTranslations( objects_.size() );
	std::vector<vl::Vec3f> newPivots( this->backupPivots_ );

	for( unsigned int i = 0; i < objects_.size(); ++i )
	{
		TransformedElementPtr object = objects_.at(i);
		SceneGraphElementPtr parent = object->parent().lock();

		if( changePivot_ )
		{
			object->setPivot( this->backupPivots_.at(i), true );
			object->setTranslate( this->backupTranslates_.at(i), true );

			vl::Mat4f transform = object->transform();
			vl::Mat4f parentTransform = parent->transform();

			vl::Vec3f origPivotWorld = vl::xform( transform, object->pivot() );
			vl::Vec3f newPivotWorld = origPivotWorld + toVec3f(translate_);

			vl::Vec3f newPivot = vl::xform( inv(transform), newPivotWorld ) - object->pivot();
			newPivot += this->backupPivots_.at(i);
			vl::Vec3f newTranslate = vl::xform( inv(parentTransform), newPivotWorld );

			newTranslations.at(i) = newTranslate;
			newPivots.at(i) = newPivot;
		}
		else
		{
			newTranslations.at(i) = this->backupTranslates_.at(i) + 
				vl::xform( inv(toMat3f(parent->transform())), toVec3f(translate_) );
		}

	}

	std::deque<ActionPtr> actions;
	ActionPtr translate( 
		new SetTranslateAction(objects_, this->backupTranslates_, newTranslations) );
	ActionPtr pivot( 
		new SetPivotAction(objects_, this->backupPivots_, newPivots) );
	actions.push_back( translate );
	actions.push_back( pivot );

	ActionPtr result( new CompoundAction( actions ) );
	return result;
}

void TranslateMousingAction::render(const ScreenSpaceConverter& converter) const
{
	unsigned int skipAxis = this->skipAxis(converter);

	glClear( GL_DEPTH_BUFFER_BIT);
	glPolygonMode( GL_FRONT, GL_FILL );

	vl::Vec3d newCenter = translate_ + axisCenter_;
	{
		GLMatrixStackHandler pushMatrix;
		glTranslatef( newCenter[0], newCenter[1], newCenter[2] );
		double axisScale = this->getAxisScale(converter);

		for( unsigned int iAxis = 0; iAxis < 3; ++iAxis )
		{
			if(iAxis == skipAxis)
				continue;

			if( this->moving() && this->clickedAxis_ != iAxis )
				continue;

			vl::Vec4f color( 0.0, 0.0, 0.0, 1.0 );
			color[iAxis] = 1.0;

			vl::Vec3 orientation(vl::vl_0);
			orientation[iAxis] = 1.0;

			if( this->changePivot_ )
			{
				glDisable( GL_LIGHTING );
				glColor4fv( color.Ref() );
				GLActionHandler lines( GL_LINES );
				glVertex3dv( vl::Vec3d(0.0, 0.0, 0.0).Ref() );
				glVertex3dv( (axisScale*orientation).Ref() );
			}
			else
			{
				drawShadedArrow( vl::Vec3d(0.0, 0.0, 0.0), orientation, axisScale, color, false );
			}
		}
	}

	{
		glDisable( GL_LIGHTING );
		vl::Vec4f yellow( 1.0, 1.0, 0.0, 1.0 );
		glColor4fv( yellow.Ref() );

		vl::Vec2d center = strip( converter.toScreenSpace(newCenter) );
		vl::Vec3d right = nearestPoint( center + vl::Vec2d( 0.1*screenRadius_, 0.0 ), newCenter, converter ) - newCenter;
		vl::Vec3d up = nearestPoint( center + vl::Vec2d( 0.0, 0.1*screenRadius_ ), newCenter, converter ) - newCenter;

		GLDisableHandler depthTest( GL_DEPTH_TEST );
		GLActionHandler lineStrip( GL_LINE_STRIP );
		glVertex3dv( ( newCenter - right - up).Ref() );
		glVertex3dv( ( newCenter - right + up).Ref() );
		glVertex3dv( ( newCenter + right + up).Ref() );
		glVertex3dv( ( newCenter + right - up).Ref() );
		glVertex3dv( ( newCenter - right - up).Ref() );
	}
}

RotateMousingAction::RotateMousingAction( const std::vector<TransformedElementPtr>& objects )
	: MousingAction( objects ), theta_(0.0)
{
	this->backupRotates_.reserve( objects.size() );
	std::transform( objects.begin(), objects.end(),
		std::back_inserter(backupRotates_),
		boost::bind( &TransformedSceneGraphElement::rotate, _1 ) );
	newRotates_ = backupRotates_;
}

RotateMousingAction::~RotateMousingAction()
{
}

double RotateMousingAction::projectPosition( const vl::Vec2d& screenPos, 
											unsigned int iAxis, 
											const vl::Mat3d& rotMat, 
											vl::Vec3d& pointOnCircle,
											vl::Vec3d& pointOnLine,
											const ScreenSpaceConverter& converter )
{
	double axisScale = this->getAxisScale( converter );

	vl::Vec3d nearVec = converter.toWorldSpace( vl::Vec3d(screenPos, -1.0) );
	vl::Vec3d farVec = converter.toWorldSpace( vl::Vec3d(screenPos, 1.0) );

	// now, find the distance between that intersection point and the given circle
    vl::Vec3d N(vl::vl_0);
	N[iAxis] = 1.0;

	vl::Vec3d x( vl::vl_0 );
	x[(iAxis+1)%3] = 1.0;

	vl::Vec3d y( vl::vl_0 );
	y[(iAxis+2)%3] = 1.0;

	Wm4::Circle3d circle;
	circle.N = toMagicVec( rotMat*N ); // rotation matrix, so inv(trans(N)) == N
	circle.U = toMagicVec( rotMat*x );
	circle.V = toMagicVec( rotMat*y );
	circle.Center = toMagicVec(axisCenter_);
	circle.Radius = axisScale;

	// first, find the intersection with the sphere of rotations
	Ray3d ray( nearVec, farVec - nearVec );
	double t;
	if( intersectSphere( ray, t, axisCenter_, axisScale ) )
	{
		pointOnLine = ray.at(t);
		Wm4::Vector3d intersectPos = toMagicVec( pointOnLine );

		Wm4::DistVector3Circle3d wmlDist( intersectPos, circle );
		double dist = wmlDist.Get();
		Wm4::Vector3d circleClosest = wmlDist.GetClosestPoint1();

		pointOnCircle = vl::Vec3d( circleClosest[0], circleClosest[1], circleClosest[2] );

		return dist;
	}
	else
	{
		// if it doesn't intersect sphere, the nearest point is going to automatically
		//   be the edge of the circle
		Wm4::Line3d viewSegment;
		viewSegment.Origin = toMagicVec( nearVec );
		viewSegment.Direction = toMagicVec( vl::norm(farVec - nearVec) );

		Wm4::DistLine3Circle3d wmlDist( viewSegment, circle );
		double dist = wmlDist.Get();
		Wm4::Vector3d lineClosest = wmlDist.GetClosestPoint0();
		Wm4::Vector3d circleClosest = wmlDist.GetClosestPoint1();

		pointOnLine = vl::Vec3d( lineClosest[0], lineClosest[1], lineClosest[2] );
		pointOnCircle = vl::Vec3d( circleClosest[0], circleClosest[1], circleClosest[2] );

		return dist;
	}
}

bool RotateMousingAction::select( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter )
{
	vl::Mat3d rotMat = eulerToQuat( this->backupRotates_.back() ).toRotMatd();

	double minDist = std::numeric_limits<double>::max();
	unsigned int minAxis = 0;
	for( unsigned int iAxis = 0; iAxis < 3; ++iAxis )
	{
		vl::Vec3d pointOnLine, pointOnCircle;
		projectPosition( screenPos, iAxis, rotMat, pointOnCircle, pointOnLine, converter );
		vl::Vec2d screenCirclePt = vl::strip( converter.toScreenSpace(pointOnCircle) );
		vl::Vec2d screenLinePt = vl::strip( converter.toScreenSpace(pointOnLine) );

		double dist = vl::len( screenLinePt - screenCirclePt );
		if( dist < minDist )
		{
			minDist = dist;
			minAxis = iAxis;
		} 
	}

	if( minDist < 10 )
	{
		clickedPos_ = screenPos;
		clickedAxis_ = minAxis;
		return true;
	}

	return false;
}

void RotateMousingAction::move( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter )
{
	vl::Mat3d rotMat = eulerToQuat( this->backupRotates_.back() ).toRotMatd();

	vl::Vec3d pointOnLine;
	vl::Vec3d newPointOnCircle, origPointOnCircle;

	projectPosition( clickedPos_, clickedAxis_, rotMat, origPointOnCircle, pointOnLine, converter );
	projectPosition( screenPos, clickedAxis_, rotMat, newPointOnCircle, pointOnLine, converter );

	vl::Vec3d axis = vl::vl_0;
	axis[ clickedAxis_ ] = 1.0;

	newPointOnCircle = vl::norm( newPointOnCircle - axisCenter_ );
	origPointOnCircle = vl::norm( origPointOnCircle - axisCenter_ );

	theta_ = acos( dot(newPointOnCircle, origPointOnCircle) );
	vl::Vec3d crossProd = cross( origPointOnCircle, newPointOnCircle );
	if( dot(crossProd, rotMat*axis) < 0.0 )
		theta_ = -theta_;

	Quaternion<double> rotation( theta_ * axis );
	assert( !isBad(rotation) );

	for( unsigned int i = 0; i < objects_.size(); ++i )
	{
		Quaternion<float> rot = eulerToQuat( this->backupRotates_.at(i) );
		Quaternion<float> changedRot = rot * rotation;
		newRotates_.at(i) = quatToEuler(changedRot);
		objects_.at(i)->setRotate( newRotates_.at(i), true );
	}
}

void RotateMousingAction::render(const ScreenSpaceConverter& converter) const
{
	{
		glClear( GL_DEPTH_BUFFER_BIT);

		GLMatrixStackHandler pushMatrix;
		glTranslatef( axisCenter_[0], axisCenter_[1], axisCenter_[2] );
		double axisScale = this->getAxisScale(converter);

		{
			Quaternion<float> quat = eulerToQuat( this->backupRotates_.back() );
			vl::Mat3f mat = quat.toRotMatf();

			float glmat[16];
			for( vl::Int i = 0; i < 3; ++i )
				for( vl::Int j = 0; j < 3; ++j )
					glmat[j*4 + i] = mat[i][j];

			vl::Vec3f t = this->objects_.back()->translate();
			glmat[12] = 0.0;
			glmat[13] = 0.0;
			glmat[14] = 0.0;
			glmat[3] = glmat[7] = glmat[11] = 0.0f;
			glmat[15] = 1.0f;
			glMultMatrixf( glmat );
		}

		static std::pair<bool, GLuint> displayList
			= std::make_pair(false, 0);
		if( !displayList.first )
		{
			GLuint dispList = glGenLists(1);
			glNewList(dispList, GL_COMPILE);
			const int divisions = 12;
			GLUquadricObj* gluq;
			gluq = gluNewQuadric();
			gluQuadricDrawStyle( gluq, GLU_FILL );
			gluQuadricTexture( gluq, GL_TRUE );
			gluSphere(gluq, 0.99, divisions, divisions);
			gluDeleteQuadric( gluq );
			glEndList();

			displayList = std::make_pair( true, dispList );
		}


		glDisable( GL_LIGHTING );
		glLineWidth( 1.0 );
		glScalef( axisScale, axisScale, axisScale );

		{
			glColorMask( GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE );
			glDepthMask( GL_TRUE );
			GLEnableHandler polygonOffset(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset (1., 1.);
			glCallList( displayList.second );
#ifdef _DEBUG
			checkGLError( "rotateSphere" );
#endif
			glColorMask( GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE );
		}

		for( unsigned int iAxis = 0; iAxis < 3; ++iAxis )
		{
			vl::Vec3f x( vl::vl_0 );
			x[(iAxis+1)%3] = 1.0;

			vl::Vec3f y( vl::vl_0 );
			y[(iAxis+2)%3] = 1.0;

			vl::Vec4f color( 0.0, 0.0, 0.0, 1.0 );
			color[iAxis] = 1.0;
			glColor4fv( color.Ref() );

			GLActionHandler lineStrip( GL_LINE_STRIP );

			const unsigned int nPoints = 30;
			for( unsigned int i = 0; i <= nPoints; ++i )
			{
				double theta = 2*M_PI*boost::numeric_cast<double>(i) 
					/ boost::numeric_cast<double>(nPoints);
				glVertex3fv( (sin(theta)*x + cos(theta)*y).Ref() );
			}
		}
	}

	{
		glDisable(GL_DEPTH_TEST);
		vl::Vec3d screenCenter = converter.toScreenSpace( axisCenter_ );
		vl::Vec3d x = converter.toWorldSpace( 
			vl::Vec3d(screenCenter[0] + screenRadius_, screenCenter[1], screenCenter[2]) ) - axisCenter_;
		vl::Vec3d y = converter.toWorldSpace( 
			vl::Vec3d(screenCenter[0], screenCenter[1] + screenRadius_, screenCenter[2]) ) - axisCenter_;

		glDisable( GL_LIGHTING );
		vl::Vec4f color( 0.7, 0.7, 0.7, 1.0 );
		glColor4fv( color.Ref() );
		GLEnableHandler lineStipple( GL_LINE_STIPPLE );
		glLineStipple( 1, 0xCCCC );

		GLActionHandler lineStrip( GL_LINE_LOOP );

		const unsigned int nPoints = 30;
		for( unsigned int i = 0; i < nPoints; ++i )
		{
			double theta = 2*M_PI*boost::numeric_cast<double>(i) / boost::numeric_cast<double>(nPoints);
			glVertex3dv( (axisCenter_ + sin(theta)*x + cos(theta)*y).Ref() );
		}
		glEnable(GL_DEPTH_TEST);
	}

}

ActionPtr RotateMousingAction::getAction() const
{
	if( fabs(theta_) < 1e-5 )
		return ActionPtr();

	ActionPtr result( 
		new RotateAction(objects_, backupRotates_, newRotates_) );
	return result;
}

ScaleMousingAction::ScaleMousingAction( const std::vector<TransformedElementPtr>& objects )
	: MousingAction( objects ), offsetScale_(1.0, 1.0, 1.0)
{
	this->backupScales_.reserve( objects.size() );
	std::transform( objects.begin(), objects.end(),
		std::back_inserter(backupScales_),
		boost::bind( &TransformedSceneGraphElement::scale, _1 ) );

	this->allowNonUniformScale_ = true;
	this->allowUniformScale_ = true;

	// need to see if anyone in the tree doesn't allow nonuniform scale:
	for( std::vector<TransformedElementPtr>::const_iterator objectIter = objects.begin();
		objectIter != objects.end(); ++objectIter )
	{
		TransformedElementPtr top = *objectIter;

		std::deque<SceneGraphElementPtr> queue( 1, top );

		while( !queue.empty() )
		{
			SceneGraphElementPtr current = queue.front();
			queue.pop_front();

			this->allowUniformScale_ = this->allowUniformScale_ && 
				current->allowUniformScale();
			this->allowNonUniformScale_ = this->allowNonUniformScale_ && 
				current->allowNonUniformScale();

			// a shear is applied if we apply a nonuniform scale to 
			//   the _parent_ and the _child_ has a rotation applied
			// Since we can't know whether we are going to rotate the child
			//   later, we block all scales at the parent level unless the
			//   child allows shear.
			for( SceneGraphElement::child_iterator iter = current->children_begin(); 
				iter != current->children_end(); ++iter )
			{
				this->allowNonUniformScale_ = this->allowNonUniformScale_ &&
					(*iter)->allowShear();
			}

			// dfs
			std::copy( current->children_begin(), current->children_end(),
				std::front_inserter( queue ) );
		}
	}
}

ScaleMousingAction::~ScaleMousingAction()
{
}

bool ScaleMousingAction::select( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter )
{
	double axisScale = this->getAxisScale( converter );

	vl::Vec2d screenRoot = vl::strip( converter.toScreenSpace( axisCenter_ ) );
	unsigned int skipAxis = this->skipAxis( converter );

	double minDist = len(screenPos - screenRoot);
	unsigned int minAxis = 10;

	const vl::Vec3f rotation = objects_.back()->rotate();
	const Quaternion<float> quat = eulerToQuat(rotation);
	const vl::Mat3d mat = quat.toRotMatd();

	vl::Mat4d transform( vl::vl_1 );
	for( vl::Int i = 0; i < 3; ++i )
		for( vl::Int j = 0; j < 3; ++j )
			transform[i][j] = mat[i][j];

	for( unsigned int iAxis = 0; iAxis < 3; ++iAxis )
	{
		vl::Vec3d axis( vl::vl_0 );
		axis[iAxis] = axisScale;
		vl::Vec2d screenPosAxis = vl::strip( converter.toScreenSpace( axisCenter_ + xform(transform, axis)) );
		if( iAxis == skipAxis )
			continue;

		if( len(screenPos - screenPosAxis) < minDist )
		{
			minDist = len(screenPos - screenPosAxis);
			minAxis = iAxis;
		}
	}

	if( minDist < 15 )
	{
		clickedPos_ = screenPos;
		clickedAxis_ = minAxis;
		return true;
	}

	return false;
}

vl::Vec3d ScaleMousingAction::projectPosition( const vl::Vec2d& pos, const ScreenSpaceConverter& converter ) const
{
	if( this->clickedAxis_ < 3 )
	{
		const vl::Vec3f rotation = objects_.back()->rotate();
		const Quaternion<float> quat = eulerToQuat(rotation);
		const vl::Mat3d mat = quat.toRotMatd();
		vl::Mat4d transform( vl::vl_1 );
		for( vl::Int i = 0; i < 3; ++i )
			for( vl::Int j = 0; j < 3; ++j )
				transform[i][j] = mat[i][j];

		vl::Vec3d nearVec = converter.toWorldSpace( vl::Vec3d( pos, -1.0 ) );
		vl::Vec3d farVec = converter.toWorldSpace( vl::Vec3d( pos, 1.0 ) );

		Wm4::Segment3d viewSegment;
		vl::Vec3d dir = farVec - nearVec;
		viewSegment.Origin = toMagicVec( nearVec + 0.5*dir );
		double dirLen = vl::len( dir );
		dir /= dirLen;
		viewSegment.Direction = toMagicVec( dir );
		viewSegment.Extent = 0.5*dirLen;

		Wm4::Line3d axisLine;
		axisLine.Origin = toMagicVec( this->axisCenter_ );
		vl::Vec3d axisDir( vl::vl_0 );
		axisDir[ this->clickedAxis_ ] = 1.0;
		axisLine.Direction = toMagicVec( xform(transform, axisDir) );

		Wm4::DistLine3Segment3d dist( axisLine, viewSegment );
		dist.GetSquared();
		return this->axisCenter_ + dist.GetLineParameter()*axisDir;
	}
	else
	{
		return nearestPoint( pos, this->axisCenter_, converter );
	}
}

void ScaleMousingAction::move( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter )
{
	double axisScale = this->getAxisScale( converter );

	unsigned int skipAxis = this->skipAxis( converter );

	vl::Vec3d startPos = projectPosition( clickedPos_, converter );
	vl::Vec3d endPos = projectPosition( screenPos, converter );
	vl::Vec3f translation = toVec3f(endPos - startPos);

	if( clickedAxis_ < 3 && this->allowNonUniformScale_ )
		this->offsetScale_[clickedAxis_] = 1.0 + translation[clickedAxis_] / axisScale;
	else if( clickedAxis_ >= 3 && this->allowUniformScale_ )
	{
		double translateAxis = 0;
		if( skipAxis == 0 )
			translateAxis = 2;

		this->offsetScale_ = vl::Vec3f(vl::vl_1) + translation[translateAxis] / axisScale * vl::Vec3f(vl::vl_1);
	}

	for( unsigned int i = 0; i < 3; ++i )
		offsetScale_[i] = std::max<float>( std::numeric_limits<float>::epsilon(), offsetScale_[i] );

	for( unsigned int i = 0; i < objects_.size(); ++i )
        objects_[i]->setScale( backupScales_[i] * offsetScale_, true );

}

ActionPtr ScaleMousingAction::getAction() const
{
	ActionPtr result(
		new ScaleAction( this->objects_, this->backupScales_, this->offsetScale_ ) );
	return result;
}

void ScaleMousingAction::render(const ScreenSpaceConverter& converter) const
{
	unsigned int skipAxis = this->skipAxis(converter);

	vl::Vec2d screenRoot = vl::strip( converter.toScreenSpace( axisCenter_ ) );

	glPolygonMode( GL_FRONT, GL_FILL );
	glClear( GL_DEPTH_BUFFER_BIT);

	GLMatrixStackHandler pushMatrix;
	glTranslatef( axisCenter_[0], axisCenter_[1], axisCenter_[2] );
	double axisScale = this->getAxisScale(converter);

	vl::Vec4f grayColor( 0.6f, 0.6f, 0.6f, 1.0f );

	{
		const vl::Vec3f rotation = objects_.back()->rotate();
		const Quaternion<float> quat = eulerToQuat(rotation);
		const vl::Mat3d mat = quat.toRotMatd();

		float glmat[16];
		for( vl::Int i = 0; i < 3; ++i )
			for( vl::Int j = 0; j < 3; ++j )
				glmat[j*4 + i] = mat[i][j];
		vl::Vec3f t = objects_.back()->translate();
		glmat[12] = 0.0;
		glmat[13] = 0.0;
		glmat[14] = 0.0;
		glmat[3] = glmat[7] = glmat[11] = 0.0f;
		glmat[15] = 1.0f;
		glMultMatrixf( glmat );
	}

	glScalef( axisScale, axisScale, axisScale );
	{
		glDisable( GL_LIGHTING );
		vl::Vec4f color( 0.7, 0.7, 0.7, 1.0 );
		glColor4fv( color.Ref() );
		GLEnableHandler lineStipple( GL_LINE_STIPPLE );
		glLineStipple( 1, 0xCCCC );

		GLActionHandler lines( GL_LINES );
		for( unsigned int iAxis = 0; iAxis < 3; ++iAxis )
		{
			if( iAxis == skipAxis )
				continue;

			vl::Vec3f x( vl::vl_0 );
			x[iAxis] = 1.0;

			glVertex3fv( vl::Vec3f(0.0, 0.0, 0.0).Ref() );
			glVertex3fv( x.Ref() );
		}

	}

	{
		glEnable( GL_LIGHTING );

		for( unsigned int iAxis = 0; iAxis < 3; ++iAxis )
		{
			if( iAxis == skipAxis )
				continue;

			vl::Vec4f color( 0.0, 0.0, 0.0, 1.0 );

			if( this->allowNonUniformScale_ )
				color[iAxis] = 1.0;
			else
				color = grayColor;

			glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, color.Ref() );
			GLMatrixStackHandler pushMatrix;
			vl::Vec3f x( vl::vl_0 );
			x[iAxis] = 1.0;
			glTranslatef(x[0], x[1], x[2]);
			glutSolidCube(0.1);
		}

		vl::Vec4f color( 1.0, 1.0, 0.0, 1.0 );
		if( !this->allowUniformScale_ )
			color = grayColor;

		glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, color.Ref() );
		glutSolidCube(0.1);
	}
}

BEGIN_EVENT_TABLE(CameraPropertyDialog, wxDialog)
	EVT_TEXT(FOCAL_LENGTH_BOX, CameraPropertyDialog::OnFocalLengthChange)
	EVT_BUTTON(wxID_OK, CameraPropertyDialog::OnOK)
END_EVENT_TABLE()

CameraPropertyDialog::CameraPropertyDialog(wxFrame *parent, const wxString& title, const wxPoint& pos, long style)
	: wxDialog( parent, -1, title, pos, wxDefaultSize, style | wxSTAY_ON_TOP )
{
	wxBoxSizer *topsizer = new wxBoxSizer( wxVERTICAL );

	wxFlexGridSizer* propsizer = new wxFlexGridSizer( 2, 2, 5, 5 );

	propsizer->Add( new wxStaticText(this, -1, "Name"), 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL);
	nameBox_ = new wxTextCtrl(this, NAME_BOX);
	propsizer->Add( nameBox_, 0, wxEXPAND );

	propsizer->Add( new wxStaticText(this, -1, "Focal length (mm)"), 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL);
	focalLengthBox_ = new wxTextCtrl(this, FOCAL_LENGTH_BOX);
	propsizer->Add( focalLengthBox_, 0, wxEXPAND  );

	propsizer->Add( new wxStaticText(this, -1, "FOV"), 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL);
	fovBox_ = new wxTextCtrl(this, FOV_BOX);
	fovBox_->Enable( false );
	propsizer->Add( fovBox_, 0, wxEXPAND  );

	propsizer->Add( new wxStaticText(this, -1, "Near z"), 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL);
	nearZBox_ = new wxTextCtrl(this, NEAR_Z_BOX);
	propsizer->Add( nearZBox_, 0, wxEXPAND );

	propsizer->Add( new wxStaticText(this, -1, "Far z"), 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL);
	farZBox_ = new wxTextCtrl(this, FAR_Z_BOX);
	propsizer->Add( farZBox_, 0, wxEXPAND );

	wxBoxSizer *button_sizer = new wxBoxSizer( wxHORIZONTAL );
	button_sizer->Add(
		new wxButton( this, wxID_OK, "OK" ),
		0,           // make horizontally unstretchable
		wxALL,       // make border all around (implicit top alignment)
		10 );        // set border width to 10
	button_sizer->Add(
		new wxButton( this, wxID_CANCEL, "Cancel" ),
		0,           // make horizontally unstretchable
		wxALL,       // make border all around (implicit top alignment)
		10 );        // set border width to 10

	int borderWidth = 5;
	topsizer->Add(
		propsizer, 
		0,
		wxEXPAND | wxALL,
		borderWidth );
	topsizer->Add(
		button_sizer,
		0,
		wxALIGN_CENTER );
    
	SetSizer( topsizer );
	topsizer->SetSizeHints( this );
}

std::string CameraPropertyDialog::name() const
{
	return std::string(nameBox_->GetValue().c_str());
}

void CameraPropertyDialog::OnFocalLengthChange(wxCommandEvent& event)
{
	using boost::lexical_cast;
	using boost::bad_lexical_cast;

	std::string text = event.GetString().c_str();
	if( text.empty() )
		return;

	try
	{

		double focalLength = boost::lexical_cast<double>(text);
		if( focalLength > 0.0 )
		{
			double fovDeg = this->fovDeg(focalLength);
			std::string text = boost::lexical_cast<std::string>( fovDeg );
			fovBox_->SetValue( wxT(text.c_str()) );
			return;
		}
	}
	catch(boost::bad_lexical_cast &)
	{
	}

	fovBox_->SetValue( wxT("") );
}

double CameraPropertyDialog::fovDeg( double focalLength )
{
	// vertical fov:
	const double frameSize = 24;

	double fovRad = 2.0 * atan2( frameSize, focalLength*2.0 );
	return fovRad * (180.0 / M_PI);
}

void CameraPropertyDialog::OnOK(wxCommandEvent& event)
{
	std::string nameStr = nameBox_->GetValue().c_str();
	if( nameStr.empty() )
	{
		errMsg( "Name cannot be empty", this );
		return;
	}

	if( focalLengthBox_->IsEnabled() )
	{
		std::string focalLengthStr = focalLengthBox_->GetValue().c_str();
		try
		{
			double focalLength = boost::lexical_cast<double>(focalLengthStr);
			if( focalLength <= 0.0 )
			{
				errMsg( "Focal length must be positive.", this );
				return;
			}

			double fovDeg = this->fovDeg(focalLength);
			fov_ = fovDeg;
		}
		catch(boost::bad_lexical_cast &)
		{
			errMsg( "Invalid focal length: '" + focalLengthStr + "'", this );
			return;
		}
	}

	std::string nearZString = nearZBox_->GetValue().c_str();
	float nearZ;
	try
	{
		nearZ = boost::lexical_cast<float>(nearZString);
		if( nearZ <= 0.0 )
		{
			errMsg( "Near z must be positive.", this );
			return;
		}
	}
	catch( boost::bad_lexical_cast &)
	{
		errMsg( "Invalid near Z: '" + nearZString + "'", this );
		return;
	}

	std::string farZString = farZBox_->GetValue().c_str();
	float farZ;
	try
	{
		farZ = boost::lexical_cast<float>( farZString );
		if( farZ <= nearZ )
		{
			errMsg( "Far z must be strictly greater than near Z.", this );
			return;
		}
	}
	catch( boost::bad_lexical_cast &)
	{
		errMsg( "Invalid far Z: '" + farZString + "'", this );
		return;
	}

	this->nearZ_ = nearZ;
	this->farZ_ = farZ;

	EndModal(wxID_OK);
}

void CameraPropertyDialog::setNearZ( float value )
{
	std::string val = boost::lexical_cast<std::string>(value);
	nearZBox_->SetValue( wxT(val.c_str()) );
}

void CameraPropertyDialog::setFarZ( float value )
{
	std::string val = boost::lexical_cast<std::string>(value);
	farZBox_->SetValue( wxT(val.c_str()) );
}

void CameraPropertyDialog::enableFOV( bool enable )
{
	focalLengthBox_->Enable( enable );
}

void CameraPropertyDialog::setFOV( double value )
{
	double fovRad = value * (M_PI/180.0);
	const double frameSize = 24;
	double focalLength = frameSize / (2.0*tan(0.5*fovRad));

	std::string val = boost::lexical_cast<std::string>(focalLength);
	this->focalLengthBox_->SetValue( wxT(val.c_str()) );
}

void CameraPropertyDialog::setName( const std::string& value )
{
	nameBox_->SetValue( wxT(value.c_str()) );
}

// material operators to abstract away particular properties
class MaterialDensityOperator
	: public MaterialOperator
{
public:
	float get( const Material& material ) const { return material.density(); }
	void set( Material& material, float value ) const { material.setDensity(value); }
};

class MaterialStaticFrictionOperator
	: public MaterialOperator
{
public:
	float get( const Material& material ) const { return material.staticFriction(); }
	void set( Material& material, float value ) const { material.setStaticFriction(value); }
};

class MaterialDynamicFrictionOperator
	: public MaterialOperator
{
public:
	float get( const Material& material ) const { return material.dynamicFriction(); }
	void set( Material& material, float value ) const { material.setDynamicFriction(value); }
};

class MaterialRestitutionOperator
	: public MaterialOperator
{
public:
	float get( const Material& material ) const { return material.restitution(); }
	void set( Material& material, float value ) const { material.setRestitution(value); }
};

class MaterialYoungsModulusOperator
	: public MaterialOperator
{
public:
	float get( const Material& material ) const { return material.youngsModulus(); }
	void set( Material& material, float value ) const { material.setYoungsModulus(value); }
};

class MaterialPoissonRatioOperator
	: public MaterialOperator
{
public:
	float get( const Material& material ) const { return material.poissonRatio(); }
	void set( Material& material, float value ) const { material.setPoissonRatio(value); }
};

class MaterialRaleighAlphaOperator
	: public MaterialOperator
{
public:
	float get( const Material& material ) const { return material.raleighAlpha(); }
	void set( Material& material, float value ) const { material.setRaleighAlpha(value); }
};

class MaterialRaleighBetaOperator
	: public MaterialOperator
{
public:
	float get( const Material& material ) const { return material.raleighBeta(); }
	void set( Material& material, float value ) const { material.setRaleighBeta(value); }
};


BEGIN_EVENT_TABLE(PropertyDialog, wxNotebook)
	EVT_TEXT_ENTER(PropertyDialog::OBJECT_NAME_BOX, PropertyDialog::OnChangeObjectName)
	EVT_TEXT_ENTER(PropertyDialog::MATERIAL_NAME_BOX, PropertyDialog::OnChangeMaterialName)
	EVT_BUTTON(PropertyDialog::CREATE_MATERIAL_BUTTON, PropertyDialog::OnCreateMaterial)
	EVT_BUTTON(PropertyDialog::DUPLICATE_MATERIAL_BUTTON, PropertyDialog::OnDuplicateMaterial)
	EVT_BUTTON(PropertyDialog::DELETE_MATERIAL_BUTTON, PropertyDialog::OnDeleteMaterial)
	EVT_CHOICE(PropertyDialog::MATERIAL_CHOICE, PropertyDialog::OnChooseMaterial)
	EVT_MENU(PropertyDialog::RANDOMIZE_PROPERTY, PropertyDialog::OnRandomizeProperty)
	EVT_MENU(PropertyDialog::CLEAR_RANDOMIZED_PROPERTY, PropertyDialog::OnUnrandomizeProperty)
	EVT_CONTEXT_MENU(PropertyDialog::OnContextMenu)
END_EVENT_TABLE()

PropertyDialog::PropertyDialog(wxFrame *parent)
	: wxNotebook( parent, -1, wxDefaultPosition, wxDefaultSize, wxNB_TOP ),
	defaultTextBoxBackground_( *wxWHITE ),
	randomizedTextBoxBackground_( 72, 159, 70 )
{
	//colourDialog_ = new wxColourDialog(this);

	randomizedPropDialog_ = new RandomizedPropertyDialog(this, "Randomize value", wxDefaultPosition);


	{
		propertyPopup_.reset( new wxMenu );
		propertyPopup_->Append( RANDOMIZE_PROPERTY, "Randomize property..." );
		propertyPopup_->Append( CLEAR_RANDOMIZED_PROPERTY, "Remove randomization" );
	}

	lastId_ = END + 1;

	textHeight_ = 20;
	textWidth_ = 60;
	textSpacing_ = 5;
	textLeft_ = 100;

	unsigned int labelWidth = textLeft_ - 2*textSpacing_;
	unsigned int totalWidth = 2*textSpacing_ + 3*textWidth_;

	// @todo replace all this with Sizers?
	// first the global properties
	{
		generalProperties_ = new wxPanel(this);
		//propsSizer->AddGrowableCol( 0, 1 );

		objectNameBox_ = new wxTextCtrl( generalProperties_, OBJECT_NAME_BOX, "", wxDefaultPosition, wxDefaultSize, wxTE_PROCESS_ENTER );
		objectNameLabel_ = 
			new wxStaticText( generalProperties_, -1, "Name", wxDefaultPosition, wxSize(textWidth_, textHeight_), wxALIGN_RIGHT );

		allTextBoxes_.push_back( objectNameBox_ );

	/*
		addTextBoxes< StateValueWrapper<PositionStateOperator> >( "Translate", top, panel, "position" );
		addTextBoxes< StateValueWrapper<OrientationStateOperator> >( "Rotate", top, panel, "orientation" );
		addTextBoxes< StateValueWrapper<LinearVelocityStateOperator> >( "Linear velocity", top, panel, "linearVelocity" );
		addTextBoxes< StateValueWrapper<AngularVelocityStateOperator> >( "Angular velocity", top, panel, "angularVelocity" );
		addTextBoxes< ScaleValueWrapper >( "Scale", top, panel, std::string() );

		staticObjectCheck_ = new wxCheckBox( panel, STATIC_CHECK, "Static", wxPoint( textLeft_, top ), wxSize( totalWidth, textHeight_ ), wxCHK_3STATE );
		top += textHeight_;
	*/

		materialChoice_ = new wxChoice( generalProperties_, MATERIAL_CHOICE, wxDefaultPosition, wxDefaultSize );

		this->AddPage( generalProperties_, wxT("General"), true );
	}

	// now the material properties
	{
		wxPanel* panel = new wxPanel(this);

		int top = 30;

		materialNameLabel_ = 
			new wxStaticText( generalProperties_, -1, "Material", wxDefaultPosition, wxSize(textWidth_, textHeight_), wxALIGN_RIGHT );
		materialNameBox_ = new wxTextCtrl( panel, MATERIAL_NAME_BOX, "", wxPoint( textLeft_, top ), wxSize( totalWidth, textHeight_ ), wxTE_PROCESS_ENTER );
		top += textHeight_ + textSpacing_;
		allTextBoxes_.push_back( materialNameBox_ );

		addBox( new PositiveValidator, new MaterialDensityOperator, top, "Density", panel );
		addBox( new PositiveValidator, new MaterialStaticFrictionOperator, top, "Static friction", panel );
		addBox( new PositiveValidator, new MaterialDynamicFrictionOperator, top, "Dynamic friction", panel );
		addBox( new RangeValidator(0.0, 1.0), new MaterialRestitutionOperator, top, "Restitution", panel );
		addBox( new PositiveValidator, new MaterialYoungsModulusOperator, top, "Young's modulus", panel );
		addBox( new PositiveValidator, new MaterialPoissonRatioOperator, top, "Poisson ratio", panel );
		addBox( new PositiveValidator, new MaterialRaleighAlphaOperator, top, "Raleigh alpha", panel );
		addBox( new PositiveValidator, new MaterialRaleighBetaOperator, top, "Raleigh beta", panel );

		unsigned int labelWidth = textLeft_ - 2*textSpacing_;
		wxStaticText* label = 
			new wxStaticText( panel, -1, wxT("Color"), wxPoint(textSpacing_, top), wxSize(labelWidth, textHeight_), wxALIGN_RIGHT );
		colorPanel_ = new wxPanel( panel, COLOR_PANEL, wxPoint(textLeft_, top), wxSize( textWidth_, textHeight_ ) );
		top += textHeight_ + textSpacing_;

		createMaterialButton_    = new wxButton( panel, CREATE_MATERIAL_BUTTON,          "New", wxPoint( textLeft_ + 0*(textSpacing_ + textWidth_), top ), wxSize( textWidth_, textHeight_ ) );
		duplicateMaterialButton_ = new wxButton( panel, DUPLICATE_MATERIAL_BUTTON, "Duplicate", wxPoint( textLeft_ + 1*(textSpacing_ + textWidth_), top ), wxSize( textWidth_, textHeight_ ) );
		deleteMaterialButton_    = new wxButton( panel, DELETE_MATERIAL_BUTTON,       "Delete", wxPoint( textLeft_ + 2*(textSpacing_ + textWidth_), top ), wxSize( textWidth_, textHeight_ ) );

		this->AddPage( panel, wxT("Material"), true );
	}

	updateSelection();
}

/*
template <typename ValueWrapperType>
void PropertyDialog::addTextBoxes( 
	const std::string& name, 
	int& top,
	wxWindow* parent,
	const std::string& plugName )
{
	unsigned int labelWidth = textLeft_ - 2*textSpacing_;
	wxStaticText* label = 
		new wxStaticText( parent, -1, wxT(name.c_str()), wxPoint(textSpacing_, top), wxSize(labelWidth, textHeight_), wxALIGN_RIGHT );

	PlugPtr plug;
	if( !plugName.empty() )
		plug = PlugRegistry::instance().plug( plugName );

	for( unsigned int i = 0; i < 3; ++i )
	{
		TextBox box;
		box.plug = plug;
		unsigned int left = textLeft_ + textSpacing_*i + textWidth_*i;
		box.wxBox = new wxTextCtrl( parent, lastId_++, "", wxPoint(left, top), wxSize(textWidth_, textHeight_), wxTE_PROCESS_ENTER );
		Connect( box.wxBox->GetId(), wxEVT_COMMAND_TEXT_ENTER, 
			(wxObjectEventFunction) &PropertyDialog::OnSetObjectProperty );
		box.wrapper.reset( new ValueWrapperType(i) );
		textValueWrappers_.insert( std::make_pair( box.wxBox->GetId(), box ) );
		allTextBoxes_.push_back( box.wxBox );
	}

	top += textHeight_ + textSpacing_;
}
*/

void PropertyDialog::addBox( Validator* validator, MaterialOperator* matOp, int& top, std::string name, wxWindow* parent )
{
	unsigned int totalWidth = 2*textSpacing_ + 3*textWidth_;
	unsigned int labelWidth = textLeft_ - 2*textSpacing_;
	wxStaticText* label = 
		new wxStaticText( parent, -1, wxT(name.c_str()), wxPoint(textSpacing_, top), wxSize(labelWidth, textHeight_), wxALIGN_RIGHT );


	MaterialBox box;
	box.materialOperator.reset( matOp );
	box.validator.reset( validator );

	box.wxBox = new wxTextCtrl( parent, lastId_++, "", wxPoint(textLeft_, top), wxSize(totalWidth, textHeight_), wxTE_PROCESS_ENTER );
	Connect( box.wxBox->GetId(), wxEVT_COMMAND_TEXT_ENTER, 
		(wxObjectEventFunction) &PropertyDialog::OnSetMaterialProperty );
	allTextBoxes_.push_back( box.wxBox );
	materialBoxWrappers_.insert( std::make_pair( box.wxBox->GetId(), box ) );

	top += textHeight_ + textSpacing_;
}

std::string toFixedNum( float value )
{
	std::ostringstream oss;
	oss.precision(5);
	oss << std::fixed << value;
	return oss.str();
}

vl::Vec3f wxColorToVec3f( const wxColour& color )
{
	return vl::Vec3f( 
		boost::numeric_cast<float>(color.Red()) / 255.0f,
		boost::numeric_cast<float>(color.Green()) / 255.0f,
		boost::numeric_cast<float>(color.Blue()) / 255.0f );
}

wxColour vec3fToWxColor( const vl::Vec3f& color )
{
	wxColour result( 
		boost::numeric_cast<unsigned char>( color[0] * 255.0f ),
		boost::numeric_cast<unsigned char>( color[1] * 255.0f ),
		boost::numeric_cast<unsigned char>( color[2] * 255.0f ) );
	return result;
}


std::deque<PhysicsObjectPtr> PropertyDialog::physicsObjects( const std::deque<SceneGraphElementPtr>& selected )
{
	std::deque<PhysicsObjectPtr> result;

	// should maybe descend the tree here?
	for( std::deque<SceneGraphElementPtr>::const_iterator selectedItr = selected.begin();
		selectedItr != selected.end(); ++selectedItr )
	{
		PhysicsObjectPtr object = boost::dynamic_pointer_cast<PhysicsObject>( *selectedItr );
		if( object )
			result.push_back( object );
	}

	return result;
}

void PropertyDialog::updateSelection()
{
	unsigned int totalWidth = 2*textSpacing_ + 3*textWidth_;
	unsigned int labelWidth = textLeft_ - 2*textSpacing_;

	this->Freeze();

	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());
	std::deque<SceneGraphElementPtr> selected = parent->getSelected();

	if( selected.empty() )
	{
		objectNameBox_->Hide();
		objectNameLabel_->Hide();
		materialChoice_->Hide();
		materialNameLabel_->Hide();

		setEnable( false );
		return;
	}

	{
		this->currentObject_ = selected.back();

		int top = 30;

		objectNameLabel_->SetPosition( wxPoint(textSpacing_, top) );
		objectNameLabel_->Show();
		objectNameBox_->SetPosition( wxPoint(textLeft_, top) );
		objectNameBox_->Show();
		top += textHeight_ + textSpacing_;

		materialNameLabel_->SetPosition( wxPoint(textSpacing_, top) );
		materialNameLabel_->Show();
		materialChoice_->SetPosition( wxPoint(textLeft_, top) );
		materialChoice_->Show();
		top += textHeight_ + textSpacing_;

		{
			std::vector<std::string> floatAttributes = currentObject_->floatAttributes();
			// hide all boxes that aren't currently being used
			for( std::vector<TextBox>::iterator textBoxItr = attributeBoxes_.begin();
				textBoxItr != attributeBoxes_.end(); ++textBoxItr )
			{
				if( std::find( floatAttributes.begin(), floatAttributes.end(), 
					textBoxItr->attributeName ) != floatAttributes.end() )
				{
					textBoxItr->wxBox->Show( true );
					textBoxItr->label->Show( true );
				}
				else
				{
					textBoxItr->wxBox->Show( false );
					textBoxItr->label->Show( false );
				}
			}

			std::string currentName;
			unsigned int left = textLeft_;
			for( std::vector<std::string>::const_iterator attItr = floatAttributes.begin();
				attItr != floatAttributes.end(); ++attItr )
			{
				std::string::size_type dotPos = attItr->rfind(".");
				if( dotPos != std::string::npos )
				{
					std::string tmpName = attItr->substr(0, dotPos);
					if( currentName != tmpName )
					{
						currentName = tmpName;
						left = textLeft_;
						top += textHeight_ + textSpacing_;
					}
				}
				else
				{
					left = textLeft_;
					currentName = *attItr;
					top += textHeight_ + textSpacing_;
				}

				// reuse if possible
				std::vector<TextBox>::iterator boxItr = 
					std::find( this->attributeBoxes_.begin(), this->attributeBoxes_.end(),
						TextBox(*attItr, 0) );
				if( boxItr == this->attributeBoxes_.end() )
				{
					// create a text box
					wxStaticText* label =
							new wxStaticText( generalProperties_, -1, 
								wxT(currentName.c_str()), 
								wxDefaultPosition, 
								wxSize(labelWidth, textHeight_), 
								wxALIGN_RIGHT );

					wxTextCtrl* box = 
						new wxTextCtrl( generalProperties_, lastId_++, "", 
							wxDefaultPosition, wxSize(textWidth_, textHeight_), wxTE_PROCESS_ENTER );
					boxItr = this->attributeBoxes_.insert( this->attributeBoxes_.end(),
						TextBox(*attItr, box) );
					boxItr->label = label;
				}

				// set positions
				if( left == textLeft_ )
				{
					boxItr->label->SetPosition( wxPoint(textSpacing_, top) );
					boxItr->label->Show();
				}
				else
				{
					boxItr->label->Hide();
				}

				boxItr->wxBox->SetPosition( wxPoint(left, top) );
				boxItr->wxBox->SetValue( toFixedNum(currentObject_->getAttribute(*attItr)) );

				if( this->currentObject_->isRandomized( boxItr->attributeName ) )
					boxItr->wxBox->SetBackgroundColour( randomizedTextBoxBackground_ );
				else
					boxItr->wxBox->SetBackgroundColour( defaultTextBoxBackground_ );

				boxItr->dirty = false;
				boxItr->wxBox->Show();
				Connect( boxItr->wxBox->GetId(), wxEVT_KILL_FOCUS, 
					(wxObjectEventFunction) &PropertyDialog::OnLeaveObjectProperty );
				Connect( boxItr->wxBox->GetId(), wxEVT_COMMAND_TEXT_ENTER, 
					(wxObjectEventFunction) &PropertyDialog::OnSetObjectProperty );
				Connect( boxItr->wxBox->GetId(), wxEVT_COMMAND_TEXT_UPDATED, 
					(wxObjectEventFunction) &PropertyDialog::OnChangeObjectProperty );
				left += textSpacing_ + textWidth_;
			}
		}

		{
			std::vector<std::string> boolAttributes = currentObject_->boolAttributes();
			// hide all boxes that aren't currently being used
			for( std::vector<CheckBox>::iterator checkBoxItr = boolAttributeBoxes_.begin();
				checkBoxItr != boolAttributeBoxes_.end(); ++checkBoxItr )
			{
				if( std::find( boolAttributes.begin(), boolAttributes.end(), 
					checkBoxItr->attributeName ) != boolAttributes.end() )
					checkBoxItr->wxBox->Show( true );
				else
					checkBoxItr->wxBox->Show( false );
			}

			top += textHeight_ + textSpacing_;
			for( std::vector<std::string>::const_iterator attItr = boolAttributes.begin();
				attItr != boolAttributes.end(); ++attItr )
			{
				// reuse if possible
				std::vector<CheckBox>::iterator boxItr = 
					std::find( this->boolAttributeBoxes_.begin(), this->boolAttributeBoxes_.end(),
						CheckBox(*attItr, 0) );
				if( boxItr == this->boolAttributeBoxes_.end() )
				{
					wxCheckBox* box = new wxCheckBox( 
						generalProperties_, lastId_++, wxT(attItr->c_str()), 
						wxDefaultPosition, wxSize( totalWidth, textHeight_ ), wxCHK_3STATE );
					boxItr = this->boolAttributeBoxes_.insert( this->boolAttributeBoxes_.end(),
						CheckBox(*attItr, box) );
				}

				boxItr->wxBox->SetPosition( wxPoint(textLeft_, top) );

				bool value = currentObject_->getBoolAttribute(*attItr);
				boxItr->wxBox->SetValue( value );
				for( std::deque<SceneGraphElementPtr>::const_iterator selectedItr = selected.begin();
					selectedItr != selected.end(); ++selectedItr )
				{
					std::vector<std::string> atts = (*selectedItr)->boolAttributes();
					if( std::find( atts.begin(), atts.end(), *attItr ) != atts.end() 
						&& (*selectedItr)->getBoolAttribute(*attItr) != value )
					{
						boxItr->wxBox->Set3StateValue( wxCHK_UNDETERMINED );
						break;
					}
				}
				boxItr->wxBox->Show();
				Connect( boxItr->wxBox->GetId(), wxEVT_COMMAND_CHECKBOX_CLICKED, 
					(wxObjectEventFunction) &PropertyDialog::OnSetCheckBoxProperty );
				top += textHeight_ + textSpacing_;
			}
		}
		
		/*
		for( std::vector<AttributePtr>::const_iterator iter = attributes.begin();
			iter != attributes.end(); ++iter )
		{
			wxStaticText* text = 
				new wxStaticText( generalProperties_, -1, 
					(*iter)->name().c_str(), 
					wxDefaultPosition, 
					wxSize(textWidth_, textHeight_), 
					wxALIGN_RIGHT );
			allGeneralProperties_.push_back( text );
			propsSizer->Add( 
				text,
				1,
				wxALIGN_RIGHT | wxALIGN_BOTTOM | wxEXPAND );

			boost::shared_ptr< TypedAttribute<vl::Vec3f> > vecAtt = 
				boost::dynamic_pointer_cast< TypedAttribute<vl::Vec3f> >( *iter );
			if( vecAtt )
			{
				wxGridSizer *vecSizer = new wxGridSizer( 1, 3, textSpacing_, 0 );

				vl::Vec3f value = vecAtt->get();
				for( size_t iComponent = 0; iComponent < 3; ++iComponent )
				{
					wxTextCtrl* box = new wxTextCtrl( 
						generalProperties_, 
						currentId++, 
						boost::lexical_cast<std::string>( value[iComponent] ).c_str(), 
						wxDefaultPosition, 
						wxDefaultSize, 
						wxTE_PROCESS_ENTER );
					vecSizer->Add( box, 1, wxALIGN_BOTTOM | wxEXPAND );

					allTextBoxes_.push_back( box );
					allGeneralProperties_.push_back( box );
				}

				propsSizer->Add( vecSizer, 1, wxEXPAND );
			}
		}
		*/
	}


	setEnable( true );

	if( selected.size() == 1 )
	{
		objectNameBox_->Enable( true );
		objectNameBox_->SetValue( wxT(selected.front()->name().c_str()) );
	}
	else
	{
		objectNameBox_->SetValue("");
		objectNameBox_->Enable( false );
	}

	this->updateValues();

	this->Thaw();
}

void PropertyDialog::updateValues()
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());

	boost::shared_ptr<Scene> scene = parent->scene();
	std::deque<SceneGraphElementPtr> selected = parent->getSelected();

	if( selected.size() == 1 )
	{
		objectNameBox_->Enable( true );
		objectNameBox_->SetValue( wxT(selected.front()->name().c_str()) );
	}
	else
	{
		objectNameBox_->SetValue("");
		objectNameBox_->Enable( false );
	}

	if( this->currentObject_ )
	{
		for( std::vector<TextBox>::iterator boxIter = attributeBoxes_.begin();
			boxIter != attributeBoxes_.end(); ++boxIter )
		{
			update( *boxIter );
		}

		for( std::vector<CheckBox>::iterator boxIter = boolAttributeBoxes_.begin();
			boxIter != boolAttributeBoxes_.end(); ++boxIter )
		{
			if( !boxIter->wxBox->IsShown() )
				continue;

			bool value = currentObject_->getBoolAttribute( boxIter->attributeName );
			boxIter->wxBox->SetValue( value );
			for( std::deque<SceneGraphElementPtr>::const_iterator selectedItr = selected.begin();
				selectedItr != selected.end(); ++selectedItr )
			{
				std::vector<std::string> atts = (*selectedItr)->boolAttributes();
				if( std::find( atts.begin(), atts.end(), boxIter->attributeName ) != atts.end() 
					&& (*selectedItr)->getBoolAttribute(boxIter->attributeName) != value )
				{
					boxIter->wxBox->Set3StateValue( wxCHK_UNDETERMINED );
					break;
				}
			}
		}
	}

	if( selected.empty() )
		return;

	// decide whather all objects have picked the same material:
	{
		std::deque<PhysicsObjectPtr> objects = this->physicsObjects( selected );
		if( objects.empty() )
		{
			currentMaterial_.reset();
		}
		else
		{
			currentMaterial_ = objects.front()->material();
			for( std::deque<PhysicsObjectPtr>::iterator objectItr = objects.begin();
				objectItr != objects.end(); ++objectItr )
			{
				if( (*objectItr)->material() != currentMaterial_ )
					currentMaterial_.reset();
			}
		}
	}

	/*
	if( selected.empty() )
		staticObjectCheck_->Set3StateValue( wxCHK_UNDETERMINED );
	else
	{
		staticObjectCheck_->SetValue( selected.front()->isStatic() );
		for( std::deque<PhysicsObjectPtr>::iterator objectItr = selected.begin();
			objectItr != selected.end(); ++objectItr )
		{
			if( (*objectItr)->isStatic() != staticObjectCheck_->GetValue() )
				staticObjectCheck_->Set3StateValue( wxCHK_UNDETERMINED );
		}
	}
	*/

	materialChoice_->Clear();
	int selection = wxNOT_FOUND;
	for( size_t iMaterial = 0; iMaterial < scene->materialCount(); ++iMaterial )
	{
		MaterialPtr mat = scene->material( iMaterial );
		int id = materialChoice_->Append( wxT(mat->name().c_str()) );
		if( mat == currentMaterial_ )
			selection = id;
	}

	materialChoice_->SetSelection(selection);

	if( currentMaterial_ )
	{
		// update the material portion of the dialog.
		for( MaterialBoxMap::iterator matItr = materialBoxWrappers_.begin();
			matItr != materialBoxWrappers_.end(); ++matItr )
		{
			matItr->second.wxBox->SetValue( 
				wxT(toFixedNum(matItr->second.materialOperator->get(*currentMaterial_))).c_str() );
		}

		materialNameBox_->SetValue( wxT(currentMaterial_->name()).c_str() );
		this->colorPanel_->SetBackgroundColour( vec3fToWxColor(currentMaterial_->color()) );
	}
}

void PropertyDialog::update( const PropertyDialog::TextBox& box )
{
	if( !this->currentObject_ )
	{
		box.wxBox->SetValue( "" );
		box.wxBox->SetBackgroundColour( defaultTextBoxBackground_ );
		return;
	}

	if( this->currentObject_->isRandomized( box.attributeName ) )
		box.wxBox->SetBackgroundColour( randomizedTextBoxBackground_ );
	else
		box.wxBox->SetBackgroundColour( defaultTextBoxBackground_ );
	box.wxBox->Refresh(FALSE);

	if( box.wxBox->IsShown() )
	{
		float value = currentObject_->getAttribute( box.attributeName );
		box.wxBox->SetValue( wxT(toFixedNum(value).c_str()) );
	}
}

void PropertyDialog::setEnable( bool enabled )
{
	if( !enabled )
	{
		std::for_each( allTextBoxes_.begin(), allTextBoxes_.end(), 
			boost::bind( &wxTextCtrl::SetValue, _1, "" ) );
	}

	std::for_each( allTextBoxes_.begin(), allTextBoxes_.end(), 
		boost::bind( &wxWindow::Enable, _1, enabled ) );

	materialChoice_->Enable( enabled );

	createMaterialButton_->Enable( enabled );
	deleteMaterialButton_->Enable( enabled );
	duplicateMaterialButton_->Enable( enabled );
}

void PropertyDialog::OnContextMenu( wxContextMenuEvent& event )
{
	if( !currentObject_ )
	{
		event.Skip();
		return;
	}
	
	wxObject* object = event.GetEventObject();
	wxTextCtrl* textCtrl = dynamic_cast<wxTextCtrl*>( object );
	if( textCtrl == 0 )
	{
		event.Skip();
		return;
	}

	std::vector<TextBox>::const_iterator itr = std::find_if(
		attributeBoxes_.begin(), attributeBoxes_.end(), 
		TextBoxEquals(textCtrl) );
	if( itr == attributeBoxes_.end() )
	{
		event.Skip();
		return;
	}

	if( this->currentObject_->isRandomized( itr->attributeName ) )
	{
		propertyPopup_->Enable( RANDOMIZE_PROPERTY, true );
		propertyPopup_->Enable( CLEAR_RANDOMIZED_PROPERTY, true );
	}
	else
	{
		propertyPopup_->Enable( RANDOMIZE_PROPERTY, true );
		propertyPopup_->Enable( CLEAR_RANDOMIZED_PROPERTY, false );
	}

	this->currentProperty_ = itr->attributeName;
	bool result = this->PopupMenu( propertyPopup_.get() );
}

void PropertyDialog::OnRandomizeProperty( wxCommandEvent& event )
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());
	assert( this->currentObject_ );
	assert( !this->currentProperty_.empty() );

	std::ostringstream oss;
	oss << "Randomize property: " << this->currentProperty_;
	randomizedPropDialog_->SetTitle( wxT( oss.str().c_str() ) );

    if( this->currentObject_->isRandomized( this->currentProperty_ ) )
	{
		RandomizedProperty randProp = 
			this->currentObject_->getRandomized( this->currentProperty_ );
		randomizedPropDialog_->setMinValue( randProp.low );
		randomizedPropDialog_->setMaxValue( randProp.high );
	}
	else
	{
		float value = currentObject_->getAttribute( this->currentProperty_ );
		randomizedPropDialog_->setMinValue( value );
		randomizedPropDialog_->setMaxValue( value );
	}

	if( randomizedPropDialog_->ShowModal() == wxID_OK )
	{
		RandomizedProperty newProp( this->currentProperty_ );
		newProp.low = std::min( randomizedPropDialog_->getMinValue(), randomizedPropDialog_->getMaxValue() );
		newProp.high = std::max( randomizedPropDialog_->getMinValue(), randomizedPropDialog_->getMaxValue() );

		ActionPtr action( new RandomizePropertyAction( this->currentObject_, newProp ) );
		parent->performAction( action );
	}
}

void PropertyDialog::OnUnrandomizeProperty( wxCommandEvent& event )
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());
	assert( this->currentObject_ );
	assert( !this->currentProperty_.empty() );
	ActionPtr action( new UnrandomizePropertyAction( this->currentObject_, this->currentProperty_ ) );
	parent->performAction( action );
}

void PropertyDialog::OnChangeMaterialName( wxCommandEvent& event )
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());

	boost::shared_ptr<Scene> scene = parent->scene();
	wxTextCtrl* control = dynamic_cast<wxTextCtrl*>( event.GetEventObject() );
	std::string newName( control->GetValue().c_str() );
	if( newName.empty() || newName == currentMaterial_->name() )
	{
		// just return to old name.
	}
	else if( scene->material( newName ) )
	{
		std::ostringstream oss;
		oss << "'" << newName << "' already in use.";
		errMsg( oss.str(), this );
		return;
	}
	else
	{
		// new name is okay
		boost::shared_ptr<Action> nameChangeAction(
			new MaterialNameChangeAction( currentMaterial_, scene, currentMaterial_->name(), newName ) );
		parent->performAction( nameChangeAction );
	}

	control->SetValue( wxT(currentMaterial_->name().c_str()) );
	control->SetSelection(-1, -1);
}

void PropertyDialog::OnChangeObjectName( wxCommandEvent& event )
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());

	boost::shared_ptr<Scene> scene = parent->scene();
	assert( parent->getSelected().size() == 1 );
	SceneGraphElementPtr selected = parent->getSelected().front();

	wxTextCtrl* control = dynamic_cast<wxTextCtrl*>( event.GetEventObject() );
	std::string newName( control->GetValue().c_str() );
	if( newName.empty() || selected->name() == newName )
	{
		// just return to old name.
	}
	else if( scene->object( newName ) )
	{
		std::ostringstream oss;
		oss << "'" << newName << "' already in use.";
		errMsg( oss.str(), this );
		return;
	}
	else
	{
		// new name is okay
		boost::shared_ptr<Action> nameChangeAction(
			new NameChangeAction( selected, scene, selected->name(), newName ) );
		parent->performAction( nameChangeAction );
	}

	control->SetValue( wxT(selected->name().c_str()) );
	control->SetSelection(-1, -1);
}

PropertyDialog::TextBox& PropertyDialog::boxForId( wxObject* object )
{
	wxTextCtrl* control = dynamic_cast<wxTextCtrl*>( object );
	assert( control != 0 );

	for( std::vector<TextBox>::iterator itr = this->attributeBoxes_.begin();
		itr != this->attributeBoxes_.end(); ++itr )
	{
		if( itr->wxBox == control )
			return *itr;
	}

	assert( false );
	return this->attributeBoxes_.front();
}

PropertyDialog::CheckBox& PropertyDialog::checkBoxForId( wxObject* object )
{
	wxCheckBox* control = dynamic_cast<wxCheckBox*>( object );
	assert( control != 0 );

	for( std::vector<CheckBox>::iterator itr = this->boolAttributeBoxes_.begin();
		itr != this->boolAttributeBoxes_.end(); ++itr )
	{
		if( itr->wxBox == control )
			return *itr;
	}

	assert( false );
	return this->boolAttributeBoxes_.front();
}


void PropertyDialog::OnChangeObjectProperty( wxCommandEvent& event )
{
	TextBox& box = this->boxForId( event.GetEventObject() );
	box.dirty = true;
}

void PropertyDialog::OnLeaveObjectProperty( wxFocusEvent& event )
{
	TextBox& box = this->boxForId( event.GetEventObject() );
	if( !box.dirty )
	{
		event.Skip();
		return; // do nothing
	}

	if( !currentObject_ )
	{
		event.Skip();
		return;
	}

	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());

	std::istringstream iss( box.wxBox->GetValue().c_str() );
	float newVal;
	iss >> newVal;

	ActionPtr action( new SetFloatAttributeAction(
		box.attributeName,
		currentObject_->getAttribute( box.attributeName ),
		newVal, 
		currentObject_ ) );

	parent->performAction( action );
	update( box );
}

void PropertyDialog::OnSetCheckBoxProperty( wxCommandEvent& event )
{
	if( !currentObject_ )
	{
		event.Skip();
		return;
	}

	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());
	CheckBox& box = this->checkBoxForId( event.GetEventObject() );

	if( box.wxBox->Get3StateValue() == wxCHK_UNDETERMINED )
	{
		event.Skip();
		return;
	}

	std::deque<SceneGraphElementPtr> selected = parent->getSelected();
	bool newVal = box.wxBox->GetValue();
	std::deque<ActionPtr> actions;
	for( std::deque<SceneGraphElementPtr>::const_iterator selectedItr = selected.begin();
		selectedItr != selected.end(); ++selectedItr )
	{
		std::vector<std::string> atts = (*selectedItr)->boolAttributes();
		if( std::find( atts.begin(), atts.end(), box.attributeName ) != atts.end() 
			&& (*selectedItr)->getBoolAttribute(box.attributeName) != newVal )
		{
			ActionPtr action( new SetBoolAttributeAction(
				box.attributeName,
				!newVal,
				newVal, 
				*selectedItr ) );
			actions.push_back( action );
		}
	}

	ActionPtr action( new CompoundAction(actions) );
	parent->performAction( action );
}

void PropertyDialog::OnSetObjectProperty( wxCommandEvent& event )
{
	if( !currentObject_ )
	{
		event.Skip();
		return;
	}

	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());
	TextBox& box = this->boxForId( event.GetEventObject() );

	std::istringstream iss( box.wxBox->GetValue().c_str() );
	float newVal;
	iss >> newVal;

	ActionPtr action( new SetFloatAttributeAction(
		box.attributeName,
		currentObject_->getAttribute( box.attributeName ),
		newVal, 
		currentObject_ ) );

	parent->performAction( action );

	// update the control and select all the text in it.
	update( box );
	box.wxBox->SetSelection(-1, -1);
}

void PropertyDialog::OnSetMaterialProperty( wxCommandEvent& event )
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());

	boost::shared_ptr<Scene> scene = parent->scene();
	std::deque<SceneGraphElementPtr> selected = parent->getSelected();
	if( !currentMaterial_ )
	{
		event.Skip();
		return;
	}

	wxTextCtrl* control = dynamic_cast<wxTextCtrl*>( event.GetEventObject() );
	int id = control->GetId();

	std::istringstream iss( control->GetValue().c_str() );
	float newVal;
	iss >> newVal;

	MaterialBoxMap::iterator itr = materialBoxWrappers_.find( id );
	assert( itr != materialBoxWrappers_.end() );

	std::string err = itr->second.validator->validate(newVal);
	if( ! err.empty() )
	{
		errMsg( err, this );
	}
	else
	{
		Material oldMat = *currentMaterial_;
		Material newMat = oldMat;
		itr->second.materialOperator->set( newMat, newVal );
		boost::shared_ptr<Action> newAction( 
			new MaterialPropertyAction( currentMaterial_, oldMat, newMat ) );
		parent->performAction( newAction );
	}

	control->SetValue( 
		wxT(toFixedNum(itr->second.materialOperator->get(*currentMaterial_)).c_str()) );
	control->SetSelection(-1, -1);
}

void PropertyDialog::OnChooseMaterial( wxCommandEvent& event )
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());

	boost::shared_ptr<Scene> scene = parent->scene();
	std::deque<PhysicsObjectPtr> objects = this->physicsObjects( parent->getSelected() );

	if( objects.empty() )
	{
		event.Skip();
		return;
	}

	int choice = event.GetInt();
	MaterialPtr newMat = scene->material( choice );

	// assign the new material to all the selected objects:
	std::deque< MaterialPtr > oldMaterials;
	for( std::deque<PhysicsObjectPtr>::const_iterator objectItr = objects.begin();
		objectItr != objects.end(); ++objectItr )
	{
		oldMaterials.push_back( (*objectItr)->material() );
	}

	ActionPtr action( new ChangeMaterialAction( objects, oldMaterials, newMat ) );
	parent->performAction( action );
}

void PropertyDialog::OnCreateMaterial( wxCommandEvent& event )
{
	createMaterial( Material("Default") );
}

void PropertyDialog::createMaterial( Material baseMat )
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());

	boost::shared_ptr<Scene> scene = parent->scene();
	std::deque<PhysicsObjectPtr> selected = this->physicsObjects( parent->getSelected() );

	// first, need to find name
	std::string prefix = "material";
	std::string name;

	for( unsigned int i = 0; name.empty() || scene->material(name); ++i )
	{
		std::ostringstream oss;
		oss << prefix << i;
		name = oss.str();
	}

	// use default properties
	boost::shared_ptr<Material> newMat( new Material(baseMat) );
	newMat->setName( name );

	std::deque< ActionPtr > actions;
	actions.push_back( ActionPtr(new CreateMaterialAction( newMat, scene ) ) );

	// assign the new material to all the selected objects:
	std::deque< MaterialPtr > oldMaterials;
	for( std::deque<PhysicsObjectPtr>::const_iterator selectedItr = selected.begin();
		selectedItr != selected.end(); ++selectedItr )
	{
		oldMaterials.push_back( (*selectedItr)->material() );
	}
	
	actions.push_back( ActionPtr( new ChangeMaterialAction( selected, oldMaterials, newMat ) ) ); 

	ActionPtr action( new CompoundAction(actions) );
	parent->performAction( action );
}


void PropertyDialog::OnDeleteMaterial( wxCommandEvent& event )
{
	MainFrame* parent = dynamic_cast<MainFrame*>(this->GetParent());

	if( !currentMaterial_ || currentMaterial_ == parent->scene()->defaultMaterial() )
	{
		event.Skip();
		return;
	}

	boost::shared_ptr<Scene> scene = parent->scene();

	std::deque< ActionPtr > actions;
	actions.push_back( ActionPtr(new DeleteMaterialAction( currentMaterial_, scene ) ) );

	// which objects lost their material?
	std::deque<PhysicsObjectPtr> cutOff;
	for( Scene::object_iter iter = scene->begin_objects(); iter != scene->end_objects(); ++iter )
	{
		PhysicsObjectPtr obj = boost::dynamic_pointer_cast<PhysicsObject>( *iter );
		if( !obj )
			continue;

		if( obj->material() == currentMaterial_ )
			cutOff.push_back( obj );
	}

	std::deque<MaterialPtr> oldMaterials( cutOff.size(), currentMaterial_ );
	actions.push_back( 
		ActionPtr( 
			new ChangeMaterialAction( 
				cutOff, 
				oldMaterials, 
				parent->scene()->defaultMaterial() ) ) ); 

	ActionPtr action( new CompoundAction(actions) );
	parent->performAction( action );


}

void PropertyDialog::OnDuplicateMaterial( wxCommandEvent& event )
{
	if( currentMaterial_ )
		createMaterial( *currentMaterial_ );
}

std::string PositiveValidator::validate( float value ) const
{
	if( value <= 0.0 )
		return std::string( "Value must be positive." );

	return std::string();
}

std::string RangeValidator::validate( float value ) const
{
	if( value < this->min_ || value > this->max_ )
	{
		std::ostringstream oss;
		oss << "Value must lie within the range [" << this->min_ << ", " << this->max_ << "].";
		return oss.str();
	}

	return std::string();
}

BEGIN_EVENT_TABLE(RandomizedPropertyDialog, wxDialog)
	EVT_BUTTON(wxID_OK, RandomizedPropertyDialog::OnOK)
END_EVENT_TABLE()

RandomizedPropertyDialog::RandomizedPropertyDialog(wxWindow* parent, const wxString& title, const wxPoint& pos, long style)
	: wxDialog( parent, -1, title, pos, wxDefaultSize, style | wxSTAY_ON_TOP )
{
	wxBoxSizer *topsizer = new wxBoxSizer( wxVERTICAL );

	wxFlexGridSizer* propsizer = new wxFlexGridSizer( 2, 2, 5, 5 );

	propsizer->Add( new wxStaticText(this, -1, "Minimum value"), 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL);
	minValueBox_ = new wxTextCtrl(this, MIN_VALUE_BOX);
	propsizer->Add( minValueBox_, 0, wxEXPAND );

	propsizer->Add( new wxStaticText(this, -1, "Maximum value"), 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL);
	maxValueBox_ = new wxTextCtrl(this, MAX_VALUE_BOX);
	propsizer->Add( maxValueBox_, 0, wxEXPAND );

	wxBoxSizer *button_sizer = new wxBoxSizer( wxHORIZONTAL );
	button_sizer->Add(
		new wxButton( this, wxID_OK, "OK" ),
		0,           // make horizontally unstretchable
		wxALL,       // make border all around (implicit top alignment)
		10 );        // set border width to 10
	button_sizer->Add(
		new wxButton( this, wxID_CANCEL, "Cancel" ),
		0,           // make horizontally unstretchable
		wxALL,       // make border all around (implicit top alignment)
		10 );        // set border width to 10

	int borderWidth = 5;
	topsizer->Add(
		propsizer, 
		0,
		wxEXPAND | wxALL,
		borderWidth );
	topsizer->Add(
		button_sizer,
		0,
		wxALIGN_CENTER );
    
	SetSizer( topsizer );
	topsizer->SetSizeHints( this );
}

void RandomizedPropertyDialog::OnOK(wxCommandEvent& event)
{
	std::string minValueStr = minValueBox_->GetValue().c_str();
	std::string maxValueStr = maxValueBox_->GetValue().c_str();
	this->minValue_ = boost::lexical_cast<float>(minValueStr);
	this->maxValue_ = boost::lexical_cast<float>(maxValueStr);

	EndModal(wxID_OK);
}

void RandomizedPropertyDialog::setMinValue( float value )
{
	std::string val = toFixedNum(value);
	minValueBox_->SetValue( wxT(val.c_str()) );
}


void RandomizedPropertyDialog::setMaxValue( float value )
{
	std::string val = toFixedNum(value);
	maxValueBox_->SetValue( wxT(val.c_str()) );
}

#ifdef WINDOWS_MEDIA
SettingsDialog::SettingsDialog(wxFrame* parent)
	: wxPropertySheetDialog( parent, -1, "Settings" )
{
	MainFrame* mainFrame = dynamic_cast<MainFrame*>( parent );

	CreateButtons(wxOK|wxCANCEL|wxHELP);
	videoOptionsPanel_ = new wxPanel(GetBookCtrl());

	// first, compute all audio codec options
	const size_t border = 5;
	wxFlexGridSizer* videoOptionsSizer = new wxFlexGridSizer(2, 2, border, border);

	std::vector<std::string> videoCodecs = mainFrame->videoCodecs();
	std::vector< wxString > videoCodecStrings;
	videoCodecStrings.reserve( videoCodecs.size() );
	for( size_t iCodec = 0; iCodec < videoCodecs.size(); ++iCodec )
		videoCodecStrings.push_back( wxString( wxT(videoCodecs[iCodec].c_str()) ) );

	videoOptionsSizer->Add( 
		new wxStaticText(videoOptionsPanel_, -1, "Video codec", wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT),
		0, wxEXPAND | wxBORDER | wxALIGN_CENTER_VERTICAL | wxALIGN_RIGHT, border );
	wxString choices[1];
	videoCodecChoice_ = new wxChoice( videoOptionsPanel_, -1, wxDefaultPosition, wxDefaultSize, 
		videoCodecStrings.size(), &videoCodecStrings[0] );
	videoOptionsSizer->Add( videoCodecChoice_, 1, wxEXPAND | wxBORDER, border );

	std::vector<std::string> audioCodecs = mainFrame->audioCodecs();
	std::vector< wxString > audioCodecStrings;
	audioCodecStrings.reserve( audioCodecs.size() );
	for( size_t iCodec = 0; iCodec < audioCodecs.size(); ++iCodec )
		audioCodecStrings.push_back( wxString( wxT(audioCodecs[iCodec].c_str()) ) );

	videoOptionsSizer->Add( 
		new wxStaticText(videoOptionsPanel_, -1, "Audio codec", wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT),
		0, wxEXPAND | wxBORDER | wxALIGN_CENTER_VERTICAL | wxALIGN_RIGHT, border );
	audioCodecChoice_ = new wxChoice( videoOptionsPanel_, -1, wxDefaultPosition, wxDefaultSize, 
		audioCodecStrings.size(), &audioCodecStrings[0] );
	videoOptionsSizer->Add( audioCodecChoice_, 1, wxEXPAND | wxBORDER, border );

	videoOptionsPanel_->SetSizer( videoOptionsSizer );

	GetBookCtrl()->AddPage(videoOptionsPanel_, wxT("Codec settings"));

	LayoutDialog();
}

void SettingsDialog::setVideoCodec( size_t iCodec )
{
	this->videoCodecChoice_->SetSelection( iCodec );
}

void SettingsDialog::setAudioCodec( size_t iCodec )
{
	this->audioCodecChoice_->SetSelection( iCodec );
}

size_t SettingsDialog::getVideoCodec() const
{
	return this->videoCodecChoice_->GetSelection();
}

size_t SettingsDialog::getAudioCodec() const
{
	return this->audioCodecChoice_->GetSelection();
}
#endif // WINDOWS_MEDIA

BEGIN_EVENT_TABLE(BitmapControl, wxControl)
    EVT_PAINT(BitmapControl::OnPaint)
    EVT_SIZE(BitmapControl::OnSize)
END_EVENT_TABLE()

BitmapControl::BitmapControl( wxWindow* parent, wxWindowID id, 
	const wxPoint& pos, 
	const wxSize& size,
	long style )
	: wxControl( parent, id, pos, size, style )
{
	currentBitmap_ = 0;
}

void BitmapControl::OnPaint( wxPaintEvent& event )
{
	// must always be here
	wxPaintDC dc(this);

	if( currentBitmap_ )
		dc.DrawBitmap( *currentBitmap_, 0, 0, true );
}

void BitmapControl::OnSize( wxSizeEvent& event )
{
	if ( GetAutoLayout() )
		Layout();
}


void BitmapControl::addBitmap( wxBitmap* bitmap )
{
	wxSize size = this->GetSize();
	wxSize newSize( std::max(size.GetWidth(), bitmap->GetWidth()),
		std::max( size.GetHeight(), bitmap->GetHeight() ) );
	this->SetMinSize( newSize );

	bitmaps_.push_back( bitmap );
	if( currentBitmap_ == 0 )
		setBitmap( 0 );
}

void BitmapControl::setBitmap( size_t iBitmap )
{
	currentBitmap_ = bitmaps_.at( iBitmap );
	Refresh( FALSE );
}



BEGIN_EVENT_TABLE(AboutDialog, wxDialog)
END_EVENT_TABLE()


wxBitmap* bitmapForArray( const unsigned char* image, size_t size, wxBitmapType type )
{
	wxMemoryInputStream stream( image, size );
	wxImage wxi( stream, type );
	return new wxBitmap( wxi, 24 );
}

AboutDialog::AboutDialog(wxWindow* parent)
	: wxDialog(parent, -1, wxT("About Many-Worlds Browsing"), wxDefaultPosition, wxSize(400, 500))
{
	wxBoxSizer* mainSizer = new wxBoxSizer( wxVERTICAL );

	BitmapControl* bitmapControl = new BitmapControl( this, -1, wxDefaultPosition, wxDefaultSize, wxNO_BORDER );
	bitmapControl->addBitmap( bitmapForArray(splash_jpg, sizeof(splash_jpg), wxBITMAP_TYPE_JPEG) );
	mainSizer->Add(
		bitmapControl,
		0,
		0,
		0 );

	text_ = new wxTextCtrl(
		this,
		TEXT,
		"",
		wxDefaultPosition, 
		wxSize(-1, 300),
		wxTE_MULTILINE | wxTE_READONLY | wxTE_RICH | wxTE_AUTO_URL | wxTE_LEFT );

	mainSizer->Add( 
		text_,
		1,
		wxEXPAND | wxALL,
		5 );

	wxColour myColor = this->GetBackgroundColour();
	text_->SetBackgroundColour(myColor);

	*text_
		<< "Many-Worlds Browsing demo v1.0\n"
		<< "Copyright 2007, Christopher D. Twigg <cdtwigg@cs.cmu.edu> and Doug L. James <djames@cs.cmu.edu>\n"
		<< "For details, please see http://graphics.cs.cmu.edu/projects/mwb/"
		<< "\n\n----\n"
		<< "This program contains code from the following libraries:\n"
		
		<< "\n"
		<< "Open Dynamics Engine\n"
		<< "Copyright 2000-2007, Russell Smith et al.\n"
		<< "http://www.ode.org/\n"

		<< "\n"
		<< "Bullet Physics Library\n"
		<< "Copyright 2005-2006, Erwin Coumans et al.\n"
		<< "http://www.continuousphysics.com/Bullet/\n"

		<< "\n"
		<< "Wild Magic Engine\n"
		<< "Copyright 2007, David Eberly\n"

		<< "\n"
		<< "Convex decompoosition library\n"
		<< "Copyright 2006, John W. Ratcliff\n"

		<< "\n"
		<< "OpenEXR\n"
		<< "Copyright Industrial Light & Magic\n"

		<< "\n"
		<< "GLEW: The OpenGL Extension Wrangler Library\n"
		<< "Copyright 2002-2005, Milan Ikits <milan ikits[]ieee org>\n"
		<< "Copyright 2002-2005, Marcelo E. Magallon <mmagallo[]debian org>\n"
		<< "Copyright 2002, Lev Povalahev\n"
		<< "All rights reserved.\n"
		<< "http://glew.sourceforge.net/\n"

		<< "\n"
		<< "wxWidgets\n"
		<< "Copyright 1998-2005, Julian Smart, Robert Roebling et al.\n"
		<< "http://www.wxwidgets.org/\n"

		<< "\n"
		<< "Boost C++ Libraries\n"
		<< "http://www.boost.org/\n"

		<< "\n"
		<< "VL: Vector Library\n"
		<< "Copyright 1995-2000, Andrew Willmott\n"
		<< "http://www.cs.cmu.edu/~ajw/doc/vl.html\n";

	*text_ << "\n\nTHIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n";

	text_->SetInsertionPoint(0);

    mainSizer->Add( 
		CreateStdDialogButtonSizer(wxOK),
		0,
		wxALL | wxALIGN_CENTER,
		5 );
	SetSizer( mainSizer );      // use the sizer for layout
	mainSizer->SetSizeHints( this );   // set size hints to honour minimum size
}


} // namespace planning


