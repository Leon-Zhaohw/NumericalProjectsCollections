
#include "stdafx.h"
#include "defs.h"
#include "util.h"
#include "defs_app.h"
#include "event_handlers.h"
#include "render.h"
#include "initialize.h"

int main(int argc, char **argv)
{

    // Parse arguments
    string inFileName;
    if (argc >= 3 && !strcmp(argv[1], "-b")) {
        bFrog = true;
        --argc; ++argv;
    }
    if (argc >= 2) {
        inFileName = argv[1];
    }
    else {
        FATAL_ERROR("\nusage: " << argv[0] << " [-b folder]\nOR\n" << argv[0] << " inputfile.off [-forcevsync:0]");
    }

    // 
    // Initialize an instance of the AppState struct that we'll use
    // to hand data between different callbacks.
    //
    AppState state = {0};

    // Load and setup the entire thing
    HRESULT hr = initialize(state, inFileName);

    // Set DXUT callbacks
    if( SUCCEEDED( hr ) )
    {
        OUT1_BEGIN("Main Loop..");
        hr = DXUTMainLoop();
        OUT1_DONE();
    }

    // Clean up
    clean_up(state, hr);

    return 0;
}

