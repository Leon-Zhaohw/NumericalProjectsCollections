#ifdef WIN32
#pragma once
#endif

#include "novodexWrappers.h"
#include "physicsFwd.h"
#include "action.h"
#include "simulationTree.h"
#include "threads.h"
#include "sound.h"

#include "twigg/keyable.h"
#include "twigg/camera.h"

#include <wx/defs.h>
#include <wx/app.h>
#include <wx/menu.h>
#include <wx/dcclient.h>
#include <wx/log.h>


#include <wx/glcanvas.h>
#include <wx/propdlg.h>
#include <wx/choice.h>
#include <wx/timer.h>
#include <wx/colordlg.h>
#include <wx/sizer.h>
#include <wx/bmpbuttn.h>
#include <wx/panel.h>
#include <wx/textctrl.h>
#include <wx/notebook.h>

#include <stack>
#include <bitset>


BEGIN_DECLARE_EVENT_TYPES()
#ifdef SOUND
	DECLARE_EVENT_TYPE(wxEVT_FINISHED_COMPUTING_MODES, -1)
#endif
	DECLARE_EVENT_TYPE(wxEVT_LOG_ERROR, -1)
	DECLARE_EVENT_TYPE(wxEVT_REPAINT_ALL_VIEWS, -1)
END_DECLARE_EVENT_TYPES()

namespace planning {

extern size_t numberOfProcessors;

void errMsg(const std::string& message, wxWindow* parent );


class GLView;
class SubView;

class TimeView;
class PropertyDialog;
class SettingsDialog;
class CameraPropertyDialog;
//class GraphFrame;
//class PlotFrame;
class BestScoresFrame;
class HierarchyFrame;
class AboutDialog;

class SimulationNode;
class SimulationTree;
class MousingAction;

class MouseCapture
{
public:
	MouseCapture( wxWindow* window );
	~MouseCapture();

private:
	wxWindow* window_;
};


class IdleAction
{
public:
	virtual ~IdleAction();
	virtual void perform() = 0;
	virtual bool complete() const = 0;
	virtual void interrupt() = 0;
};

#ifdef LOGGING
class Logger
{
public:
	Logger( const std::string& filename );
	~Logger();
	void log( const std::string& message );

private:
	std::string filename_;
	std::ofstream ofs_;
	DWORD startTime_;
	boost::mutex lock_;
};
#endif

class UsefulBits
{
public:
	bool currentSelected() const
	{
		return bits_.test(3);
	}

	bool selected() const
	{
		return bits_.test(0);
	}

	bool involvedInJoint() const
	{
		return bits_.test(1);
	}

	bool involvedInConstraint() const
	{
		return bits_.test(4);
	}

	bool changedFromTree() const
	{
		return bits_.test(2);
	}

	void setCurrentSelected( bool value = true )
	{
		bits_.set(3, value);
	}

	void setSelected( bool value = true )
	{
		bits_.set(0, value);
	}

	void setInvolvedInJoint( bool value = true )
	{
		bits_.set(1, value);
	}

	void setInvolvedInConstraint( bool value = true )
	{
		bits_.set(4, value);
	}

	void setChangedFromTree( bool value = true )
	{
		bits_.set(2, value);
	}

private:
	std::bitset<5> bits_;
};


class MainFrame
	:	public wxFrame, 
		public SimulationNodeCallback
{
public:
	enum RenderStyle
	{
		STYLE_BOXES,
		STYLE_POLYS,
		STYLE_WIREFRAME,
		STYLE_POINTS,
	};

	enum BackgroundStyle
	{
		BACKGROUND_WHITE,
		BACKGROUND_LIGHT_GRAY,
		BACKGROUND_DARK_GRAY,
		BACKGROUND_BLACK,
	};

	enum Mode
	{
		MODE_EDIT,
		MODE_SOLVE
	};

    MainFrame(wxFrame *frame, const wxString& title, const wxPoint& pos, const wxSize& size,
        long style = wxDEFAULT_FRAME_STYLE);
	~MainFrame();
	void preDestroy();

	void addNodes( const SimulationTree* tree, const std::deque<const Path*>& nodes, const boost::dynamic_bitset<>& active );
	void setNodes( const SimulationTree* tree, const boost::dynamic_bitset<>& nodes );

	void test();

	void OnIdle(wxIdleEvent& event);

	void OnMenuFileOpen(wxCommandEvent& event);
	void OnMenuFileImport(wxCommandEvent& event);
	void OnMenuFileNew(wxCommandEvent& event);
	void OnMenuFileSave(wxCommandEvent& event);
	void OnMenuFileSaveAs(wxCommandEvent& event);
	void OnMenuFileLoadTree(wxCommandEvent& event);
	void OnMenuFileLoadSamples(wxCommandEvent& event);
	void OnMenuFileSaveTree(wxCommandEvent& event);
    void OnMenuFileLoadDebuggingInfo(wxCommandEvent& event);

	void OnMenuFileDumpGenerationImages(wxCommandEvent& event);
#ifdef DRAW_LINES
	void OnMenuFileDumpPathsImages(wxCommandEvent& event);
#endif

	void OnMenuSetSampling(wxCommandEvent& event);
	void OnMenuUpdateSample(wxCommandEvent& event);

#ifdef SOUND
	void OnFinishedComputingModes(wxCommandEvent& event);
#endif
	void OnLogError(wxCommandEvent& event);

	void OnExit(wxCommandEvent& event);
    void OnExportImageSequence(wxCommandEvent& event);
    void OnExportImage(wxCommandEvent& event);
    void OnExportTimelapseImage(wxCommandEvent& event);

	void OnExportRib(wxCommandEvent& event);
	void OnExportMayaAscii(wxCommandEvent& event);

	void OnCopy(wxCommandEvent& event);
	void OnPaste(wxCommandEvent& event);
	void OnDuplicate(wxCommandEvent& event);
	void OnUndo(wxCommandEvent& event);
	void OnRedo(wxCommandEvent& event);

	void OnHelpAbout(wxCommandEvent& event);

#ifdef WINDOWS_MEDIA
	void OnEditSettings(wxCommandEvent& event);
#endif

	void OnShowAll(wxCommandEvent& event);
	void OnHide(wxCommandEvent& event);

	void OnGroup(wxCommandEvent& event);
	void OnUngroup(wxCommandEvent& event);

	void OnCenterPivot(wxCommandEvent& event);
	void OnFixCenterOfMass(wxCommandEvent& event);

	void OnFuse(wxCommandEvent& event);
	void OnUnfuse(wxCommandEvent& event);
	void OnFitBoxes(wxCommandEvent& event);
	void OnFitConvexHulls(wxCommandEvent& event);
	void OnCreateConvexHull(wxCommandEvent& event);

	void OnSetCamera(wxCommandEvent& event);
	void OnCameraProperties(wxCommandEvent& event);

	void OnPlay(wxCommandEvent& event);

	void updateMenus();

	void OnClose(wxCloseEvent& event);

    GLView *GetCanvas();
    const GLView *GetCanvas() const;

    void OnSize(wxSizeEvent& event);

	void OnSetRenderOption(wxCommandEvent& event);
	void OnSetRenderShowAllPaths(wxCommandEvent& event);
	void OnNeedsRedraw(wxCommandEvent& event);
	void OnRepaintAllViews(wxCommandEvent& event);

	void OnCreateBox(wxCommandEvent& event);
	void OnCreateSphere(wxCommandEvent& event);
	void OnCreateSpheres(wxCommandEvent& event);
	void OnCreateCappedCylinder(wxCommandEvent& event);
	void OnCreateCylinder(wxCommandEvent& event);
	void OnCreateCone(wxCommandEvent& event);
	void OnCreateObj(wxCommandEvent& event);
	void OnCreatePlane(wxCommandEvent& event);
	void OnCreateOrthoCam(wxCommandEvent& event);
	void OnCreatePerspectiveCam(wxCommandEvent& event);
	void OnCreateBallJoint(wxCommandEvent& event);
	void OnCreateSliderJoint(wxCommandEvent& event);
	void OnCreateHingeJoint(wxCommandEvent& event);
	void OnCreateUniversalJoint(wxCommandEvent& event);
	void OnCreateHinge2Joint(wxCommandEvent& event);
//	void OnCreateRotMotor1Joint(wxCommandEvent& event);
    void OnCreatePlaneJoint(wxCommandEvent& event);

	void OnLoadCameraKeys(wxCommandEvent& event);
	void OnSaveCameraKeys(wxCommandEvent& event);
	void OnMenuSimulateRun(wxCommandEvent& event);
	void OnSolve(wxCommandEvent& event);

	void OnSoundComputeAllModes(wxCommandEvent& event);
	void OnSoundComputeObjectModes(wxCommandEvent& event);
	void OnSoundStopComputingModes(wxCommandEvent& event);

	void OnScrollTime(wxScrollEvent& event);
	void OnSetMode(wxCommandEvent& event);

	void setMode( MainFrame::Mode mode, bool force = false );
	MainFrame::Mode getMode() const { return this->mode_; }

	void OnWindowSimSorter(wxCommandEvent& event);
	void OnWindowHierarchyView(wxCommandEvent& event);
	void OnWindowProperties(wxCommandEvent& event);

	void OnKey(wxKeyEvent& event);
	void OnKeyUp(wxKeyEvent& event);

	void errMsg(const std::string& message);

	bool wireframeOnShaded() const;
	bool showAllPaths() const;
	bool visualizeInertiaTensor() const;
	bool visualizeConvexHulls() const;
	bool visualizeContactPoints() const;
	RenderStyle renderStyle() const;
	BackgroundStyle backgroundStyle() const;
	void setBackgroundStyle( BackgroundStyle style );
	vl::Vec3f backgroundColor() const;
	bool renderGrid() const;

	bool hasShadows() const;
	void load( const std::string& filename );
	void loadAny( const std::string& filename );
	void loadTree( const std::string& filename );
	void loadSamples( const std::string& filename );
#ifdef SOUND
	void loadModesFile( const std::string origFilename );
#endif
	void stopSimulating();
	bool isSimulating() const;

//	NxPhysicsSDKPtr physicsSDK() const;

	void undo();
	void redo();

	void log( const std::string& message );

	void setTime( float time, bool move );
	float getTime() const;

	void addCamera( CameraWrapperPtr camera );
	void addObject( SceneGraphElementPtr object );
	void addObjects( const std::deque<SceneGraphElementPtr>& object );
	void deleteObjects( const std::deque<SceneGraphElementPtr>& objects );
	void deleteConstraint( ConstraintPtr constraint );

	void performAction( ActionPtr action );

	void copy( const std::deque<SceneGraphElementPtr>& objects );
	void paste();
	void duplicate( const std::deque<SceneGraphElementPtr>& objects );

	void showProperties( bool show = true );
	bool arePropertiesShown() const;

#ifdef WINDOWS_MEDIA
	void recordMovie();
#endif // WINDOWS_MEDIA
#ifdef SOUND
	void recordSound();
#endif // SOUND

	enum ObjectMode
	{
		MODE_OBJECT_SELECT,
		MODE_OBJECT_SELECTPOINTS,
		MODE_OBJECT_TRANSLATE,
		MODE_OBJECT_CHANGEPIVOT,
		MODE_OBJECT_ROTATE,
		MODE_OBJECT_SCALE
	};
	void setObjectMode( ObjectMode mode );
	ObjectMode getObjectMode() const;
	void OnChangeObjectMode(wxCommandEvent& event);

	boost::shared_ptr<Scene> scene();
	boost::shared_ptr<const Scene> scene() const;

	void changeSelected();
	void sceneChanged();

	void logError( const std::string& error );

	void setPath( SimulationTreePtr tree, const Path* path );

	void setSamples( ConstSimulationTreePtr tree, const boost::dynamic_bitset<>& samples );

	// generate the line cache from the given list
	void generateLineCache(
		ConstSimulationTreePtr tree, 
		const std::deque<const Path*>& list,
		bool replace );
	void updateSelectedPath();

	const Path* getPath() const;
	SimulationTreePtr getTree();
	ConstSimulationTreePtr getTree() const;

	void addTree( SimulationTreePtr tree );
	void removeTree( SimulationTreePtr tree );
	void clearTrees();

	bool playing() const;
	void updateViews();

	void setCurrentView( GLView* view, bool toggleMaximized );
	bool viewMaximized() const;

	Simulation::State state() const;
	std::vector<Simulation::ContactPoint> contactPoints() const;

	template <typename Handler>
		void walkTree( Handler& handler, const Simulation::State& state ) const;

	bool isModifiedFromTree( ConstSceneGraphElementPtr object ) const;
	void addSamples( const std::deque<const SimulationNode*>& samples );

	void updateCameraMenu();

	enum MouseMode
	{
		MODE_MOUSE_SELECT,
		MODE_MOUSE_ADDITIVE_CONSTRAINT,
		MODE_MOUSE_SUBTRACTIVE_CONSTRAINT,
		MODE_MOUSE_REFINE,
        MODE_MOUSE_SKETCH_VELOCITY,
        MODE_MOUSE_SKETCH_ACTION
	};

	MouseMode mouseMode() const;

	enum SimulatorType
	{
		SIMULATOR_ODE_NATIVE,
		SIMULATOR_ODE_BULLET,
        SIMULATOR_ODE_PENALTY,
        SIMULATOR_NOVODEX,
		SIMULATOR_NEWTON,
		SIMULATOR_BULLET,
	};

	SimulatorType simulatorType() const;

	std::deque<Simulation::State> stateHistory() const;
	std::vector<PiecewisePath> piecewisePaths() const;
	size_t piecewisePathsFrameRate() const;

	void stopSimulation();
	void startSimulation();
	boost::shared_ptr<SimulationFactory> simulationFactory() const;

	void startSubtree( 
		SimulationTreePtr parent,
		const std::vector<PiecewisePath>& paths, 
		size_t pathsFrameRate,
		float startTime,
		const std::deque<ConstPhysicsObjectPtr>& active,
		bool backwards );

    // We need to be able to do things like refine, etc. using just the timeline
    // This is a function the timeline can call to report that such an event has
    // happened
    void applyTimeSelection( float startTime, float endTime );

	void sampleTree( boost::shared_ptr<Event> samplingEvent );
#ifdef SOCKET_SAMPLING
	void sampleTreeRemote(const std::string& hostname, boost::shared_ptr<Event> samplingEvent);
#endif

	void setSelected( std::deque<SceneGraphElementPtr> selected );
	std::deque<SceneGraphElementPtr> getSelected() const;

	typedef std::pair<SceneGraphElementPtr, std::deque<unsigned int> > SelectedPoints;
	void setSelectedPoints( const SelectedPoints& selectedPoints );
	SelectedPoints getSelectedPoints() const;

	// this applies the current selection operator; e.g., "select by hierarchy"
	std::deque<SceneGraphElementPtr> convertSelected( const std::deque<SceneGraphElementPtr>& selected );

	void upHierarchy();
	void downHierarchy();

	ConstraintPtr selectedConstraint() const;
	void setSelectedConstraint( ConstraintPtr constraint );

	void updateActiveTimeRange( std::pair<double, double> timeRange );

	bool wireframeOnSelected() const;

public:
	std::deque<ActionPtr> reparent( const std::deque<SceneGraphElementPtr>& elts,
		const std::deque<SceneGraphElementPtr>& newParents,
		bool adjustTransforms = true );

	boost::shared_ptr<MousingAction> mouseAction();
	boost::shared_ptr<const MousingAction> mouseAction() const;

	std::vector<size_t> objectListToDynamicObjectIds( 
		const std::deque<SceneGraphElementPtr>& objects,
		const std::vector<ConstPhysicsObjectPtr>& dynamicObjects ) const;

public:
	typedef std::vector<vl::Vec3f> LineCache;
	typedef std::vector<GLint> LineCacheIndex;
	typedef std::vector<GLsizei> LineCacheCount;

	const std::deque<LineCache>& lineCache() const;
	const std::deque<LineCache>& lineCacheColors() const;
	const std::deque<LineCacheIndex>& lineCacheStarts() const;
	const std::deque<LineCacheCount>& lineCacheCounts() const;
	ReaderWriterLock& lineCacheMutex() const;
	const std::deque<LineCache>& selectedLineCache() const;
	ReaderWriterLock& selectedLineCacheMutex() const;

private:
	void clearSamples();
	void addToLineCache( const Path* path, bool active );

	// current line cache, used by the various views
	mutable std::deque<LineCache> lineCache_;
	mutable std::deque<LineCache> lineCacheColors_;

	mutable std::deque<LineCacheIndex> lineCacheStarts_;
	mutable std::deque<LineCacheCount> lineCacheCounts_;

	mutable RandomStream lineCacheStream_;
	mutable ReaderWriterLock lineCacheMutex_;

	mutable std::deque<LineCache> selectedLineCache_;
	mutable ReaderWriterLock selectedLineCacheMutex_;


private:
	PhysicsObjectPair getJointObjects() const;

	typedef STLEXT hash_set<const SceneGraphElement*, 
		hash_ptr<const SceneGraphElement*> > DirtyObjectSet;
	DirtyObjectSet dirtyObjects_;
	void updateDirtyObjects();

	// We need to make sure the "dirty" flag gets propagated through
	//   joints, etc.
	void updateDirtyConnections( const Scene& scene );

	void addActiveSamples(const std::deque<const Path*>& samples );
	void clearActiveSamples();

	// this function is going to watch the current jobs queue and clear it
	//   out as needed.
	void watchJobsQueue();
	void computeHistory();

	void stopAllSampling();

	void changeScene( boost::shared_ptr<Scene> newScene );

	// These all return true if successful
	bool checkDirty();
	bool save();
	bool save(const std::string& filename);
	bool saveAs();

	std::deque<SceneGraphElementPtr> cloneObjects( const std::deque<SceneGraphElementPtr>& objects );
	SceneGraphElementPtr cloneObject( ConstSceneGraphElementPtr object, 
		STLEXT hash_set<std::string>& used, SceneGraphElementMap& oldToNewMap );

#ifdef LOGGING
	boost::shared_ptr<Logger> logger_;
#endif

	std::deque<GLView*> canvases_;
	GLView* currentCanvas_;

	TimeView*  timeView_;

	wxMenuBar *menuBar_;
	wxMenu *fileMenu_;
	wxMenu *editMenu_;
	wxMenu *editSelectByMenu_;
	wxMenu *createMenu_;
	wxMenu *cameraMenu_;
	wxMenu *renderMenu_;
	wxMenu *renderGeometryMenu_;
	wxMenu *renderBackgroundMenu_;
	wxMenu *simulateMenu_;
	wxMenu *simulatorMenu_;
	wxMenu *windowMenu_;
#ifdef SOUND
	wxMenu *soundMenu_;
#endif // SOUND
	wxMenu *helpMenu_;

	wxToolBar *toolBar_;
	wxChoice *modeChoice_;
	wxLogWindow* logWindow_;

	AboutDialog* aboutDialog_;

	PropertyDialog* propertyDialog_;
	SettingsDialog* settingsDialog_;
	wxBoxSizer* propSizer_;
	wxFlexGridSizer* viewsSizer_;

	CameraPropertyDialog* cameraPropertyDialog_;
	BestScoresFrame* bestScoresFrame_;
	HierarchyFrame* hierarchyFrame_;

	wxBitmap playBitmap_;
	wxBitmap pauseBitmap_;
	wxBitmapButton* playButton_;
	bool playing_;

	float time_;

	bool wireframeOnSelected_;

	enum
	{
		myID_TIMER,
		myID_STATUSBAR,
		myID_TIMESLIDER,
		myID_PLAY,
		myID_SELECT_RECT,
		myID_GOOD_CONSTRAINT_RECT,
		myID_BAD_CONSTRAINT_RECT,
		myID_REFINE_RECT,
		myID_SKETCH_ACTION,
		myID_CURSOR_TOOLBAR,
		myID_TRANSLATE_TOOLBAR,
		myID_ROTATE_TOOLBAR,
		myID_SCALE_TOOLBAR,
		myID_CHANGEPIVOT_TOOLBAR,
		myID_SELECTPOINTS_TOOLBAR,
//		myID_SIM,
		MENU_FILE,
		MENU_FILE_OPEN,
		MENU_FILE_IMPORT,
		MENU_FILE_SAVE,
		MENU_FILE_SAVEAS,
		MENU_FILE_NEW,
		MENU_FILE_LOAD_TREE,
		MENU_FILE_SAVE_TREE,
		MENU_FILE_LOAD_SAMPLES,
//		MENU_FILE_SAVE_CAMERA_KEYS,
//		MENU_FILE_LOAD_CAMERA_KEYS,
        MENU_FILE_LOAD_DEBUGGING_INFO,
		MENU_FILE_EXPORT_IMAGE,
		MENU_FILE_EXPORT_IMAGES,
		MENU_FILE_EXPORT_RIB,
		MENU_FILE_EXPORT_TIMELAPSE,
		MENU_FILE_EXPORT_MAYA_ASCII,
		MENU_FILE_DUMP_GENERATION_IMAGES,
#ifdef DRAW_LINES
		MENU_FILE_DUMP_PATHS_IMAGES,
#endif
		MENU_FILE_EXIT,
		MENU_EDIT,
		MENU_EDIT_COPY,
		MENU_EDIT_PASTE,
		MENU_EDIT_DUPLICATE,
		MENU_EDIT_UNDO,
		MENU_EDIT_REDO,
		MENU_EDIT_HIDE,
		MENU_EDIT_SHOW_ALL,
		MENU_EDIT_GROUP,
		MENU_EDIT_UNGROUP,
		MENU_EDIT_CENTER_PIVOT,
		MENU_EDIT_FIX_CENTER_OF_MASS,
		MENU_EDIT_FUSE,
		MENU_EDIT_UNFUSE,
		MENU_EDIT_FIT_BOXES,
		MENU_EDIT_FIT_CONVEX_HULLS,
		MENU_EDIT_CREATE_CONVEX_HULL,
#ifdef WINDOWS_MEDIA
		MENU_EDIT_SETTINGS,
#endif // WINDOWS_MEDIA
		MENU_EDIT_SELECT_BY,
		MENU_EDIT_SELECT_BY_OBJECT,
		MENU_EDIT_SELECT_BY_HIERARCHY,
		MENU_EDIT_SELECT_BY_SIMULATED,
		MENU_RENDER,
		MENU_RENDER_GEOMETRY,
		MENU_RENDER_GEOMETRY_BOXES,
		MENU_RENDER_GEOMETRY_POLYS,
		MENU_RENDER_GEOMETRY_WIREFRAME,
		MENU_RENDER_GEOMETRY_POINTS,
		MENU_RENDER_WIREFRAME_ON_SHADED,
		MENU_RENDER_INERTIA_TENSOR,
		MENU_RENDER_CONVEX_HULLS,
		MENU_RENDER_CONTACT_POINTS,
		MENU_RENDER_SHOW_ALL_PATHS,
		MENU_RENDER_SHADOWS,
		MENU_RENDER_GRID,
		MENU_RENDER_BACKGROUND,
		MENU_RENDER_BACKGROUND_WHITE,
		MENU_RENDER_BACKGROUND_BLACK,
		MENU_RENDER_BACKGROUND_LIGHT_GRAY,
		MENU_RENDER_BACKGROUND_DARK_GRAY,
		MANU_CAMERA,
		MENU_CAMERA_SET_KEY,
		MENU_CAMERA_ERASE_KEY,
		MENU_CAMERA_NEXT_KEY,
		MENU_CAMERA_PREV_KEY,
		MENU_CAMERA_PROPERTIES,
		MENU_CREATE,
		MENU_CREATE_BOX,
		MENU_CREATE_SPHERE,
		MENU_CREATE_SPHERES,
		MENU_CREATE_CAPPED_CYLINDER,
		MENU_CREATE_CYLINDER,
		MENU_CREATE_CONE,
		MENU_CREATE_OBJ,
		MENU_CREATE_PLANE,
		MENU_CREATE_ORTHO_CAMERA,
		MENU_CREATE_PERSPECTIVE_CAMERA,
		MENU_CREATE_BALL_JOINT,
		MENU_CREATE_HINGE_JOINT,
		MENU_CREATE_SLIDER_JOINT,
		MENU_CREATE_UNIVERSAL_JOINT,
		MENU_CREATE_HINGE2_JOINT,
		MENU_CREATE_ROTMOTOR1_JOINT,
		MENU_CREATE_PLANE_JOINT,
		MENU_SIMULATE,
		MENU_SIMULATE_RUN,
		MENU_SIMULATE_SAMPLING,
		MENU_SIMULATE_SOLVE,
		MENU_SIMULATE_UPDATE_PATH,
		MENU_SIMULATE_SIMULATOR,
		MENU_SIMULATE_SIMULATOR_ODE_NATIVE,
		MENU_SIMULATE_SIMULATOR_ODE_BULLET,
        MENU_SIMULATE_SIMULATOR_ODE_PENALTY,
        MENU_SIMULATE_SIMULATOR_NOVODEX,
		MENU_SIMULATE_SIMULATOR_NEWTON,
		MENU_SIMULATE_SIMULATOR_BULLET,
#ifdef SOUND
		MENU_SOUND_ENABLE_SOUND,
		MENU_SOUND_COMPUTE_ALL_MODES,
		MENU_SOUND_COMPUTE_OBJECT_MODES,
		MENU_SOUND_STOP_COMPUTING_MODES,
#endif // SOUND
		MENU_HELP_ABOUT,
		MENU_WINDOW_SIM_SORTER,
		MENU_WINDOW_HIERARCHY_VIEW,
		MENU_WINDOW_PROPERTIES,
		BUTTON_PLAY,
		CHOICE_MODE,
		CAMERA_MENU_START = 1000,
	};

	Mode mode_;

	std::deque< ActionPtr > undoQueue_;
	std::deque< ActionPtr > redoQueue_;

	std::deque< SceneGraphElementPtr > clipboard_;

	bool dirty_;
	boost::shared_ptr<Scene> scene_;
	std::string sceneFilename_;

#ifdef USE_NOVODEX
	NxPhysicsSDKPtr physicsSDK_;
#endif

	SimulationTreePtr currentTree_;
	const Path* currentPath_;
	mutable boost::mutex currentPathMutex_;

	std::deque<SimulationTreePtr> trees_;
	typedef boost::shared_ptr<SimulationNodeCallbackRegistrar> CallbackRegistrarPtr;
	std::deque<CallbackRegistrarPtr> callbackRegistrars_;
	boost::mutex treesMutex_;
	size_t treeCounter_;

	/**** stuff to do with computing history *****/
	typedef std::pair<Simulation::State, const SimulationNode*> StateWithNode;
	// this array will keep track of the current state history
	std::deque<StateWithNode> stateHistory_;

	std::vector<PiecewisePath> piecewisePaths_;
	size_t piecewisePathsFrameRate_;

	// when we want to compute history for a new node, we'll put it here:
	std::deque< std::pair<SimulationTreePtr, const Path*> > computeHistoryPending_;
	// mutex to lock access:
	mutable boost::mutex computeHistoryPendingMutex_;
	// condition to watch when someone puts something in computeHistoryPending_:
	mutable boost::condition computeHistoryPendingCondition_;
	mutable ReaderWriterLock historyLock_;
	mutable Event historyComplete_;

	// the thread that will do all the history computing:
	boost::scoped_ptr<boost::thread> computeHistoryThread_;
	// a bool to tell it to stop:
	bool computeHistoryDone_;

	// keep track of how many threads are currently sampling
	std::vector<HANDLE> samplingThreadHandles_;
	std::vector<DWORD> samplingThreadIds_;

	bool sampling_;
	bool pauseSampling_;
	IncrementDecrement numMachines_;

	// We will need to notify the sampling threads when either
	//    (1) the simulationTree changes, or
	//    (2) we are terminating and need to stop all sampling
	// We will use this event
	std::vector< boost::shared_ptr<Event> > samplingEvents_;

	// this particular simulation is useful for testing the current scene
	boost::shared_ptr<Simulation> simulation_;

	// current selection
	std::deque<SceneGraphElementPtr> selected_;
	mutable boost::mutex selectedLock_;

	ConstraintPtr selectedConstraint_;

	std::pair< SceneGraphElementPtr, std::deque<unsigned int> > selectedPoints_;

	boost::shared_ptr<MousingAction> mouseAction_;

	mutable bool updated_;

    DECLARE_EVENT_TABLE()

	boost::shared_ptr<IdleAction> idleAction_;

#ifdef WINDOWS_MEDIA
	friend class RecordMovieIdleAction;
	friend class RecordSoundIdleAction;

	// Windows Media stuff
private:

	boost::shared_ptr<IWMProfileManager> profileManager_;

	std::vector<std::string> audioCodecDescriptions_;
	std::vector< boost::shared_ptr<IWMStreamConfig> > audioStreamConfigurations_;

	std::vector<std::string> videoCodecDescriptions_;
	std::vector< boost::shared_ptr<IWMStreamConfig> > videoStreamConfigurations_;

	std::pair< std::vector<std::string>,
		std::vector< boost::shared_ptr<IWMStreamConfig> > >
			findAvailableCodecs( REFGUID codecType, bool (*CheckType) (WM_MEDIA_TYPE*), bool requireVBR );

	size_t selectedVideoCodec_;
	size_t selectedAudioCodec_;

	double prevHighestAudio_;

#ifdef SOUND
	typedef std::map<AudioHash, 
		boost::shared_ptr<SquashingCubesModes>, CompareAudioHash > DescriptionToModesMap;
	DescriptionToModesMap descriptionToSquashingCubesModes_;
	ReaderWriterLock descriptionToSquashingCubesModesLock_;
	bool computingModes_;
	Event stopComputingModes_;
	boost::shared_ptr<boost::thread> computeModesThread_;

	void findModes( const std::deque<PhysicsObjectPtr>& physicsObjects, bool compute );
#endif // SOUND
public:
	std::vector<std::string> audioCodecs() const;
	std::vector<std::string> videoCodecs() const;
#endif // WINDOWS_MEDIA
};

// A mousing action is something that has a UI in the 
//  window; specifically, I mean "translate", "rotate", etc.
class MousingAction
{
public:
	MousingAction( const std::vector<TransformedElementPtr>& objects );
	virtual ~MousingAction();

	bool pressMouse( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter );
	void moveMouse( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter );
	void releaseMouse( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter );

	virtual void render(const ScreenSpaceConverter& converter) const = 0;

	bool moving() const;

	ActionPtr action();

protected:
	virtual ActionPtr getAction() const = 0;


protected:
	virtual bool select( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter ) = 0;
	virtual void move( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter ) = 0;

	unsigned int skipAxis(const ScreenSpaceConverter& converter) const;
	double getAxisScale(const ScreenSpaceConverter& converter) const;

	std::vector<TransformedElementPtr> objects_;

	const double screenRadius_;
	vl::Vec3d axisCenter_;

private:
	bool moving_;
};

class TranslateMousingAction
	: public MousingAction
{
public:
	TranslateMousingAction( const std::vector<TransformedElementPtr>& objects, bool changePivot );
	~TranslateMousingAction();

	void render(const ScreenSpaceConverter& converter) const;

protected:
	bool select( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter );
	void move( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter );

	ActionPtr getAction() const;

private:
	vl::Vec3d projectPosition( const vl::Vec2d& pos, const ScreenSpaceConverter& converter ) const;

	std::vector<vl::Vec3f> backupTranslates_;
	std::vector<vl::Vec3f> backupPivots_;
	vl::Vec2d clickedPos_;
	unsigned int clickedAxis_;
	vl::Vec3d translate_;

	// for moving the pivot point.  Having a separate class would
	// be the more elegant option here, but there is a lot of shared
	// code between the two and I don't feel like refactoring up into
	// a common superclass right now.
	const bool changePivot_;
};

class RotateMousingAction
	: public MousingAction
{
public:
	RotateMousingAction( const std::vector<TransformedElementPtr>& objects );
	~RotateMousingAction();

	void render(const ScreenSpaceConverter& converter) const;

protected:
	bool select( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter );
	void move( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter );

	ActionPtr getAction() const;

private:
	// Find the point on the circle with this axis that is closest to position
	double projectPosition( const vl::Vec2d& position, unsigned int iAxis, 
		const vl::Mat3d& rotMat, vl::Vec3d& pointOnCircle, vl::Vec3d& pointOnLine, const ScreenSpaceConverter& converter );

	std::vector<vl::Vec3f> backupRotates_;
	std::vector<vl::Vec3f> newRotates_;

	unsigned int clickedAxis_;
	vl::Vec2d clickedPos_;
	double theta_;
};

class ScaleMousingAction
	: public MousingAction
{
public:
	ScaleMousingAction( const std::vector<TransformedElementPtr>& objects );
	~ScaleMousingAction();

	void render(const ScreenSpaceConverter& converter) const;

protected:
	bool select( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter );
	void move( const vl::Vec2d& screenPos, const ScreenSpaceConverter& converter );

	ActionPtr getAction() const;

private:
	vl::Vec3d projectPosition( const vl::Vec2d& pos, const ScreenSpaceConverter& converter ) const;

	unsigned int clickedAxis_;
	vl::Vec2d clickedPos_;
	std::vector< vl::Vec3f > backupScales_;
	vl::Vec3f offsetScale_;

	bool allowNonUniformScale_;
	bool allowUniformScale_;
};

class CameraPropertyDialog
	: public wxDialog
{
public:
	CameraPropertyDialog(wxFrame *parent, const wxString& title, const wxPoint& pos, 
	   long style = wxDEFAULT_DIALOG_STYLE);

	void OnFocalLengthChange(wxCommandEvent& event);
	void OnOK(wxCommandEvent& event);

	float nearZ() const      { return nearZ_; }
	float farZ() const       { return farZ_; }
	double fov() const        { return fov_; }
	std::string name() const;

	void setNearZ( float value );
	void setFarZ( float value );
	void setFOV( double value );
	void setName( const std::string& value );
	void enableFOV( bool enable );

private:
	double fovDeg( double focalLength );

	wxTextCtrl* nameBox_;
	wxTextCtrl* focalLengthBox_;
	wxTextCtrl* fovBox_;
	wxTextCtrl* nearZBox_;
	wxTextCtrl* farZBox_;

	enum IdEnum
	{
		NAME_BOX,
		FOCAL_LENGTH_BOX,
		FOV_BOX,
		NEAR_Z_BOX,
		FAR_Z_BOX
	};

	float nearZ_;
	float farZ_;
	float fov_;

    DECLARE_EVENT_TABLE()
};

class MaterialOperator
{
public:
	virtual float get( const Material& material ) const = 0;
	virtual void set( Material& material, float value ) const = 0;
	virtual ~MaterialOperator() {}
};

class Validator
{
public:
	virtual ~Validator() {}
	virtual std::string validate( float value ) const = 0;
};

class RangeValidator
	: public Validator
{
public:
	RangeValidator( float min, float max )
		: min_(min), max_(max) {}

	std::string validate( float value ) const;

private:
	float min_;
	float max_;
};

class PositiveValidator
	: public Validator
{
public:
	std::string validate( float value ) const;
};

class RandomizedPropertyDialog
	: public wxDialog
{
public:
	RandomizedPropertyDialog(wxWindow* parent, const wxString& title, const wxPoint& pos, 
	   long style = wxDEFAULT_DIALOG_STYLE);

	void OnOK(wxCommandEvent& event);

	void setMinValue( float value );
	void setMaxValue( float value );
	float getMinValue() const { return minValue_; }
	float getMaxValue() const { return maxValue_; }

private:
	wxTextCtrl* minValueBox_;
	wxTextCtrl* maxValueBox_;

	float minValue_;
	float maxValue_;

	enum IdEnum
	{
		MIN_VALUE_BOX,
		MAX_VALUE_BOX,
	};

    DECLARE_EVENT_TABLE()
};

#ifdef WINDOWS_MEDIA
class SettingsDialog
	: public wxPropertySheetDialog
{
public:
	SettingsDialog(wxFrame* parent);

	void setVideoCodec( size_t iCodec );
	void setAudioCodec( size_t iCodec );

	size_t getVideoCodec() const;
	size_t getAudioCodec() const; 

private:
	wxPanel* videoOptionsPanel_;
	wxChoice* videoCodecChoice_;
	wxChoice* audioCodecChoice_;
};
#endif

class PropertyDialog
	: public wxNotebook
{
public:
	PropertyDialog(wxFrame *parent);
	void updateValues();
	void updateSelection();

	void OnSetCheckBoxProperty( wxCommandEvent& event );
	void OnSetObjectProperty( wxCommandEvent& event );
	void OnLeaveObjectProperty( wxFocusEvent& event );
	void OnChangeObjectProperty( wxCommandEvent& event );
	void OnSetMaterialProperty( wxCommandEvent& event );
	void OnChangeObjectName( wxCommandEvent& event );
	void OnChangeMaterialName( wxCommandEvent& event );

	void OnCreateMaterial( wxCommandEvent& event );
	void OnDeleteMaterial( wxCommandEvent& event );
	void OnDuplicateMaterial( wxCommandEvent& event );

	void OnRandomizeProperty( wxCommandEvent& event );
	void OnUnrandomizeProperty( wxCommandEvent& event );

	void OnChooseMaterial( wxCommandEvent& event );
	void OnContextMenu( wxContextMenuEvent& event );

private:
	std::deque<PhysicsObjectPtr> physicsObjects( const std::deque<SceneGraphElementPtr>& selected );

	void createMaterial( Material baseMat );

	void setEnable( bool enabled );

	struct TextBox
	{
		TextBox( const std::string& attName, wxTextCtrl* box )
			: attributeName(attName), wxBox(box), label(0), dirty(false)
		{
		}

		wxTextCtrl* wxBox;
		std::string attributeName;
		wxStaticText* label;
		bool dirty;

		bool operator==( const TextBox& other )
		{
			return (this->attributeName == other.attributeName);
		}
	};

	struct TextBoxEquals
	{
		TextBoxEquals( const wxTextCtrl* box )
			: box_(box) {}

		bool operator()( const TextBox& box )
		{
			return box.wxBox == this->box_;
		}

	private:
		const wxTextCtrl* box_;
	};

	struct CheckBox
	{
		CheckBox( const std::string& attName, wxCheckBox* box )
			: attributeName(attName), wxBox(box)
		{
		}

		wxCheckBox* wxBox;
		std::string attributeName;

		bool operator==( const CheckBox& other )
		{
			return (this->attributeName == other.attributeName);
		}
	};

	TextBox& boxForId( wxObject* object );
	CheckBox& checkBoxForId( wxObject* object );

	struct MaterialBox
	{
		boost::shared_ptr<Validator> validator;
		boost::shared_ptr<MaterialOperator> materialOperator;
		wxTextCtrl* wxBox;
	};

	void update( const TextBox& box );

	std::deque< wxTextCtrl* > allTextBoxes_;

	wxPanel* generalProperties_;


	wxChoice* materialChoice_;
	wxTextCtrl* objectNameBox_;
	wxTextCtrl* materialNameBox_;
	wxStaticText* objectNameLabel_;
	wxStaticText* materialNameLabel_;

	wxButton* createMaterialButton_;
	wxButton* deleteMaterialButton_;
	wxButton* duplicateMaterialButton_;

	wxPanel* colorPanel_;

	boost::scoped_ptr<wxMenu> propertyPopup_;
	// we will need to keep track of the current property
	// so that we can apply the randomization to the correct one:
	std::string currentProperty_;
	RandomizedPropertyDialog* randomizedPropDialog_;

	MaterialPtr currentMaterial_;
	SceneGraphElementPtr currentObject_;

	enum IdEnum
	{
		MATERIAL_CHOICE,
		OBJECT_NAME_BOX, 
		MATERIAL_NAME_BOX,
		CREATE_MATERIAL_BUTTON,
		DELETE_MATERIAL_BUTTON,
		DUPLICATE_MATERIAL_BUTTON,
		COLOR_PANEL,
		RANDOMIZE_PROPERTY,
		CLEAR_RANDOMIZED_PROPERTY,
		END,
	};

	template <typename ValueWrapperType>
	void addTextBoxes( const std::string& name, int& top, wxWindow* parent, const std::string& plugName );
	void addBox( Validator* validator, MaterialOperator* matOp, int& top, std::string name, wxWindow* parent );

	unsigned int textHeight_;
	unsigned int textWidth_;
	unsigned int textSpacing_;
	unsigned int textLeft_;

	typedef STLEXT hash_map< int, MaterialBox > MaterialBoxMap;
	MaterialBoxMap materialBoxWrappers_;

	std::vector<TextBox> attributeBoxes_;
	std::vector<CheckBox> boolAttributeBoxes_;

	int lastId_;

	wxColour defaultTextBoxBackground_;
	wxColour randomizedTextBoxBackground_;

	wxColourDialog* colourDialog_;

    DECLARE_EVENT_TABLE()
};

class BitmapControl
	: public wxControl
{
public:
	BitmapControl( wxWindow* parent, wxWindowID id, 
		const wxPoint& pos = wxDefaultPosition, 
		const wxSize& size = wxDefaultSize,
		long style = wxSTATIC_BORDER );

	void addBitmap( wxBitmap* bitmap );
	void setBitmap( size_t iBitmap );

	void OnPaint( wxPaintEvent& event );
	void OnSize( wxSizeEvent& event );

private:
	std::vector<wxBitmap*> bitmaps_;
	wxBitmap* currentBitmap_;

	DECLARE_EVENT_TABLE()
};


class AboutDialog
	: public wxDialog
{
public:
	AboutDialog(wxWindow* parent);

private:
	wxTextCtrl* text_;

	enum
	{
		TEXT
	};

    DECLARE_EVENT_TABLE()
};


std::deque<SceneGraphElementPtr> cleanSelected( const std::deque<SceneGraphElementPtr>& selected );

} // namespace planning


