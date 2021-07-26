// InterferenceVolume.h
//
// This class is where most of the iagm magic happens
//

#ifndef _INTERFERENCEVOLUME_H_
#define _INTERFERENCEVOLUME_H_

#include <Intersection.h>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <ext/hash_map>
#include <Eigen/Core>



namespace IAGM
{

class Curve;



class InterferenceVolume
{
public:
	InterferenceVolume(Curve *C, Intersections &is, __gnu_cxx::hash_map<size_t,size_t> &vertIndices);
	
	void computeHandleGradient(Curve *C, std::set<size_t> &movable);

	double getPressure()
	{
		double den = m_dVdp.dot(m_dVdp);
		if (den == 0.0)
			return 0.0;

		return (-m_volume / den);
	}
	
	Eigen::VectorXd& getHandleGradients()
	{ return m_dVdp; }

protected:

	void computeInterferenceVolume(Curve *C, Intersections &is, __gnu_cxx::hash_map<size_t,size_t> &vertIndices);
	void computeInterferenceVolumeGradient(Curve *C, Intersections &is, __gnu_cxx::hash_map<size_t,size_t> &vertIndices);

/*
	double getVolume(Eigen::Vector2d &r,
					 Eigen::Vector2d x0, Eigen::Vector2d x1, Eigen::Vector2d x2,
				     Eigen::Vector2d v0, Eigen::Vector2d v1, Eigen::Vector2d v2,
					 double A, double t);
*/

	// Determines the correct orientation for each summand based on relative velocity.
	// This is most useful during the analytic gradient computation
	//
	double getSign(size_t idx, double t, double a, Eigen::Vector2d *q, Eigen::Vector2d *v);

	void getDnDt(Eigen::Vector2d *v, Eigen::Vector2d &dndt);

	void getTimeGradient(Intersection &i, size_t idx, double *x0, double *x1, Eigen::Vector2d &dtdq);
	
	double getBarycentricArea(size_t idx, size_t nbrVertices, double *x);

    double m_volume;

	// Gradients with respect to fine DOFs
	//
	Eigen::VectorXd m_dVdq;

	// Gradient with respect to coarse DOFs (handles, control mesh, etc.)
	//
	Eigen::VectorXd m_dVdp;

};

} // end namespace IAGM

#endif

