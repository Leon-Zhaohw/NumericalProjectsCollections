// Curve.h
//

#ifndef _CURVE_H_
#define _CURVE_H_

#include <Intersection.h>
#include <vector>
#include <set>



namespace IAGM 
{
class InterferenceVolume;

typedef std::vector<InterferenceVolume> InterferenceVolumes;
typedef std::vector<InterferenceVolume>::iterator InterferenceVolumesIterator;


class Curve
{
public:
	Curve()
	 : m_maxNbrIters(50), m_minSep(0.0)
	{ }

	// The public interface to interference-awareness. Calling will
	// iterate until all interference is removed or maximum number of
	// iterations is exceeded. Supplied interface (see below) is used
	// to update geometry
	//
	void findAndRemoveInterference();


	// These purely virtual methods must be implemented by the deriving
	// class. They are the interface between interference awareness and
	// the modeling code that modifies the curves
	//


	// Project the vector in to the modeling subspace
	//
	virtual void getSubspaceVector(double *in, double *out) = 0;

	virtual size_t getNbrVertices() = 0;
	virtual size_t getNbrCtrlVertices() = 0;

	// These two vectors specify the desired trajectory introduced
	// by an edit. This trajectory will be modified to remove interference
	//
	virtual void getStartConfiguration(double *&x)   = 0;
	virtual void getEndConfiguration(double *&x)     = 0;

	virtual void getControlConfiguration(double *&x) = 0;

	// Notifies the subclass that the control mesh has been updated and
	// that the slave / limit mesh should be re-evaluated
	//
	virtual void copyControlConfiguration() = 0;

	// The iagm library needs to know what vertices it can move to remove
	// interference
	//
	virtual void getMovableControlVertices(std::set<size_t> &movable) = 0;

	virtual void getEdgeIndices(size_t &nbr, unsigned int *&e) = 0;


protected:

	// Interface methods to be implemented by deriving class
	//
	virtual void updateSlaveMesh() = 0;

	// Default collision detection, can be overridden
	//
	virtual bool getCurveIntersections(Intersections &is);


	// Helper functions
	//
	void modelInterference(Intersections &is, InterferenceVolumes &ivs);
	void removeInterference(InterferenceVolumes &ivs);

	bool getIntersection(double *x00, double *x01,
						 double *x10, double *x11, double *x20, double *x21,
						 double &s, double &t);


	// Getter / setter methods for control parameters
	//
	unsigned int getMaxNbrIterations() const
	{ return m_maxNbrIters; }
	void setMaxNbrIterations(unsigned int nbr)
	{ m_maxNbrIters = nbr; }

	double getMinimumSeparation() const
	{ return m_minSep; }
	void setMinimumSeparation(double h)
	{ m_minSep = h; }

private:

	unsigned int m_maxNbrIters;

	double m_minSep;

};

}

#endif

