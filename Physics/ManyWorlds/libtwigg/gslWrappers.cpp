#include "stdafx.h"
#include "twigg/gslWrappers.h"

#include <boost/thread.hpp>
#include <sstream>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#ifdef USE_GSL

GSLSolverWrapper::GSLSolverWrapper( const gsl_odeiv_step_type* stepType, 
	unsigned int systemSize, 
	double *timestep,
	vl::Real error) 
	:	timestep_( timestep ), 
		systemSize_(systemSize)
#ifndef CONDOR
		, lock_( new boost::mutex() )
#endif // CONDOR
{
	s_ = gsl_odeiv_step_alloc (stepType, systemSize_);
	c_ = gsl_odeiv_control_y_new (error, 0.0);
	e_ = gsl_odeiv_evolve_alloc (systemSize);			
}

GSLSolverWrapper::~GSLSolverWrapper()
{
	gsl_odeiv_evolve_free(e_);
	gsl_odeiv_control_free(c_);
	gsl_odeiv_step_free(s_);
}

int GSLSolverWrapper::evolve( int (* function) (double t, const double y[], double dydt[], void * params), 
	double *currentTime, double endTime, double y[], void* params)
{
	gsl_odeiv_system sys = {function, 0, systemSize_, params};

#ifndef CONDOR
	boost::mutex::scoped_lock scoped_lock(*lock_);
#endif // CONDOR
	return gsl_odeiv_evolve_apply (e_, c_, s_,
							&sys,
							currentTime, endTime,
							timestep_, y);
}

gsl_odeiv_step* GSLSolverWrapper::steppingFunction()
{
	return s_;
}

#endif // USE_GSL
