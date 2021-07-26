#ifndef __GSLWRAPPERS_H__
#define __GSLWRAPPERS_H__

#ifdef USE_GSL

#include "util.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_math.h>
#include <boost/scoped_ptr.hpp>

namespace boost
{
	class mutex;
}

class GSLSolverWrapper
{
public:
	GSLSolverWrapper( const gsl_odeiv_step_type* stepType, 
		unsigned int systemSize, 
		double *timestep,
		double error);

	~GSLSolverWrapper();

	int evolve( int (* function) (double t, const double y[], double dydt[], void * params), 
		double *currentTime, double endTime, double y[], void* params);

	gsl_odeiv_step* steppingFunction();

private:
	gsl_odeiv_step* s_;
	gsl_odeiv_control* c_;
	gsl_odeiv_evolve* e_;

	double* timestep_;
	unsigned int systemSize_;

	mutable boost::scoped_ptr<boost::mutex> lock_;
};

#endif // USE_GSL

#endif

