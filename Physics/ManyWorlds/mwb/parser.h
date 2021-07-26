#ifdef WIN32
#pragma once
#endif

#ifndef __PARSER_H__
#define __PARSER_H__

#include "physicsFwd.h"

namespace planning {

boost::shared_ptr<Scene> parseSimulation( const std::string& filename );

} // namespace planning

#endif

