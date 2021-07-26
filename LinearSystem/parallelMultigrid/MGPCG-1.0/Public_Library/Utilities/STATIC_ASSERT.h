//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header STATIC_ASSERT
//#####################################################################
#ifndef __STATIC_ASSERT__
#define __STATIC_ASSERT__

#include <boost/static_assert.hpp>

#define STATIC_ASSERT(...) BOOST_STATIC_ASSERT((__VA_ARGS__))

#endif
