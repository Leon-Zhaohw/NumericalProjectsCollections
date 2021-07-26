#ifdef WIN32
#pragma once
#endif

//#include "resource.h"


#include "wx/defs.h"
#include "wx/app.h"

namespace planning {

// Define a new application type
class PlanningApp : public wxApp
{
public:
    bool OnInit();
};

} // namespace planning
