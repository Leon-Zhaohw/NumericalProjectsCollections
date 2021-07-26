// skinning.cpp : Defines the entry point for the application.
//

#include "stdafx.h"
#include "stdarg.h"
#include "odeWrappers.h"
#include "novodexWrappers.h"
#include "OBB.h"
#include "voxels.h"

#define MAX_LOADSTRING 100

#ifdef _WIN32
#include <windows.h>
#endif


#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>


#ifdef GUI
#include <wx/log.h>
#include <wx/settings.h>

#include "planningUI.h"

#include "wx/defs.h"
#include "wx/app.h"
#include "wx/image.h"
#else
#include "scene.h"
#include "simulationTree.h"
#include "parser.h"
#endif

#include <cstdio>
#include <cstdarg>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

void dLoadLCP( const char* filename );


extern "C" 
{

#ifdef GUI
void reportError( int errnum, const char *msg, va_list ap )
{
	std::string message( "ODE: error %d: ", errnum );
	message.append( msg );
	wxLogError(message.c_str(), ap);

	// error function is not supposed to return, according to docs
	assert( false );
	abort();
}

void reportDebug( int errnum, const char *msg, va_list ap )
{
	std::string message( "ODE: debug error %d: ", errnum );
	message.append( msg );
	wxLogError(message.c_str(), ap);

	// error function is not supposed to return, according to docs
	assert( false );
	abort();
}

void reportMessage( int errnum, const char *msg, va_list ap )
{
	std::vector<char> buf( 1000 );
#ifdef WIN32
	_vsnprintf( &buf[0], buf.size()-1, msg, ap );
#else
	snprintf( &buf[0], buf.size()-1, msg, ap );
#endif
	std::string message( "ODE: " );
	message.append( &buf[0] );

	wxLogVerbose(message.c_str(), ap);
}
#endif

} // extern "C"

namespace planning {

void init()
{
#ifdef GUI
	// set ODE error handlers
	// if we don't do this, the default handlers seem to crash
	//   the program
	dSetErrorHandler( &reportError );
	dSetDebugHandler( &reportDebug );
	dSetMessageHandler( &reportMessage );

	::wxInitAllImageHandlers();
#endif

	boost::filesystem::path::default_name_check( &boost::filesystem::native );

	// sometimes boost filesystem causes crashes.  We'll test for this here
	// so we find out early.
	boost::filesystem::path test( "c:\\test" );
}

} // namespace planning


#ifdef GUI

namespace planning {

// Define a new application type
class PlanningApp : public wxApp
{
public:
    bool OnInit();
};

// `Main program' equivalent, creating windows and returning main app frame
bool PlanningApp::OnInit()
{
#ifdef _DEBUG
    /*
    // look for "lcpBad.*"
    {
        const boost::filesystem::path dir_path( "." );
        boost::filesystem::directory_iterator end_itr;
        for ( boost::filesystem::directory_iterator itr( dir_path );
            itr != end_itr; ++itr )
        {
            if( !is_regular( itr->status() ) )
                continue;

            std::string filename = itr->path().leaf();
            const char* prefix = "lcpBad";
            if( !std::equal(prefix, prefix+sizeof(prefix), filename.begin()) )
                continue;

            dLoadLCP( filename.c_str() );
        }
    }
    */
#endif

	init();

	/*
	int width = wxSystemSettings::GetMetric(wxSYS_SCREEN_X);
	int height = wxSystemSettings::GetMetric(wxSYS_SCREEN_Y);
#ifdef _WIN32
	// not sure why this has to be 4:
	width = ::GetSystemMetrics(SM_CXMAXIMIZED) - 4*wxSystemSettings::GetMetric(wxSYS_EDGE_X);
	height = ::GetSystemMetrics(SM_CYMAXIMIZED) - 4*wxSystemSettings::GetMetric(wxSYS_EDGE_Y);
#endif
	*/
	int width = 720;
	int height = 480  + wxSystemSettings::GetMetric(wxSYS_CAPTION_Y) + wxSystemSettings::GetMetric(wxSYS_MENU_Y);

  // Create the main frame window
  MainFrame *frame = new MainFrame(NULL, wxT("Many-Worlds Browsing"), 
	  wxPoint(0, 0), 
	  wxSize(width, height));

  for( unsigned int i = 1; i < argc; ++i )
  {
	frame->loadAny( argv[i] );
	//frame->setMode( MainFrame::MODE_SOLVE );
  }

//	testOBB();
  testVoxels();
  testQuat();

  /*
#ifdef _DEBUG
  frame->test();
#endif
  */

  // Show the frame
  frame->Show(TRUE);
 
  return TRUE;

  _CrtDumpMemoryLeaks();
}


} // namespace planning

IMPLEMENT_APP(planning::PlanningApp)

#elif defined(CONDOR_WORKER)

#include "mwSample.h"
int main( int argc, char** argv )
{
	planning::Worker_Sample worker;

/*
	// Set up debug level here...
	set_MWprintf_level ( 75 );

	MWprintf ( 10, "A worker is starting.\n" );
	for ( int i = 1; i < argc; i++ )
		MWprintf(10, " Args %d was %s\n", i, argv[i] );
*/
	worker.go( argc, argv );

	return 0;
}

#elif defined(SAMPLE_SERVER)
#else
namespace po = boost::program_options;

int main( int argc, char** argv )
{
	using namespace planning;

	boost::filesystem::path::default_name_check( &boost::filesystem::native );

	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("input-file", po::value<std::string>(), "input .plan file")
		("output-file", po::value<std::string>(), "output .tree file")
		("num-samples", po::value<size_t>(), "number of samples")
		("frame-rate", po::value<size_t>(), "frame rate (default: 120)")
		("max-impulse-time", po::value<size_t>(), "max impulse spanning time (specified in number of frames)")
		("stop-time", po::value<float>(), "Stopping time (default: never)")
		("test", "Run tests");

	

	po::positional_options_description p;
	p.add("input-file", 1);
	p.add("output-file", 1);

	po::variables_map vm;
	try
	{
		po::store(po::command_line_parser(argc, argv).
			options(desc).positional(p).run(), vm);
		po::notify(vm);
	}
	catch( po::error& err )
	{
		std::cerr << err.what() << std::endl;
		return -1;
	}

	if( vm.count("help"))
	{
		std::cout << desc << "\n";
		return 0;
	}

	if( vm.count("input-file") == 0 )
	{
		std::cout << "Error: no input file specified.\n";
		std::cout << desc << "\n";
		return 1;
	}

	if( vm.count("output-file") == 0 )
	{
		std::cout << "Error: no output file specified.\n";
		std::cout << desc << "\n";
		return 1;
	}

	std::string inFilename = vm["input-file"].as< std::string >();
	std::string outFilename = vm["output-file"].as< std::string >();
	

	size_t frameRate = 120;
	if( vm.count("frame-rate") )
		frameRate = vm["frame-rate"].as< size_t >();

	size_t maxImpulseTime = 1;
	if( vm.count("max-impulse-time") )
		maxImpulseTime = vm["max-impulse-time"].as< size_t >();

	size_t numSamples = 500;
	if( vm.count("num-samples") )
		numSamples = vm["num-samples"].as< size_t >();

	float stopTime = boost::numeric::bounds<float>::highest();
	if( vm.count( "stop-time") )
		stopTime = vm["stop-time"].as< float >();

	try
	{
#ifdef USE_NOVODEX
		NxPhysicsSDKPtr sdk = 
			instantiateSDK(0, 0);
		sdk->setParameter(NX_MIN_SEPARATION_FOR_PENALTY, -0.001f);
#endif

		RandomStream random;
		boost::shared_ptr<Scene> scene = 
			parseSimulation( inFilename );

		SimulationTree tree(
			outFilename, 
			scene, 
			frameRate
#ifdef USE_NOVODEX
			, sdk
#endif
			);

		tree.setEndTime( stopTime );

		for( size_t i = 1; i <= numSamples; ++i )
		{
			std::cout << "Computing sample " << i << "... " << std::flush;
			tree.sample(random, false);
			std::cout << "done." << std::endl;
		}
	}
	catch( ParserException& e )
	{
		std::cerr << "Error reading '" << inFilename << "': " << std::endl;
		std::cerr << e;
		return 1;
	}
	catch( IOException& e )
	{
		std::cerr << "Error reading '" << inFilename << "': " << std::endl;
		std::cerr << e.message();
		return 1;
	}
	catch( Exception& e )
	{
		std::cerr << "Unknown exception: " << e.message() << std::endl;
		return 1;
	}
}

#endif

