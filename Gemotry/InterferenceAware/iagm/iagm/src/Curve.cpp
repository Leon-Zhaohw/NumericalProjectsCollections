// Curve.cpp
//

#include "Curve.h"
#include "InterferenceVolume.h"

#include <iostream>
#include <ext/hash_map>

#include <Eigen/Core>
#include "rpoly.h"
#include "Hash.h"

using namespace std;
using namespace __gnu_cxx;
using namespace Eigen;



namespace IAGM
{

void Curve::findAndRemoveInterference()
{
	size_t nbrIters = 0;

	Intersections is;
	while (getCurveIntersections(is) && nbrIters < m_maxNbrIters)
	{
		InterferenceVolumes ivs;
		modelInterference(is, ivs);
		removeInterference(ivs);

		copyControlConfiguration();
		updateSlaveMesh();

		nbrIters++;
		is.clear();
	}

	if (nbrIters == m_maxNbrIters)
	{
		// TODO: Practically speaking, it might be reasonable to cap
		// the number of iterations to something reasonable. If this
		// threshold is hit, then rollback to the start configuration,
		// which is guaranteed to be safe (since the start configuration
		// contains to intersections, by hypothesis).
		//
		cerr << "Warning: maximum number of iterations exceeded! "
		     << "Mesh may contain intersections." << endl;
	}
}

void Curve::modelInterference(Intersections &is, InterferenceVolumes &ivs)
{
	// We tag each vertex in the mesh with the earliest intersection in its
	// barycentric region that it is involved in.
	//
	size_t idx=0;
	hash_map<size_t, size_t> h;
	for (IntersectionsIterator iItr=is.begin(); iItr!=is.end(); ++iItr, ++idx)
	{
		if (h.find(iItr->getVertex()) == h.end() ||
				   iItr->getTime() < is[h[iItr->getVertex()]].getTime())
			h[iItr->getVertex()] = idx;

		if (iItr->getAlpha() <= 0.5 &&
		   (h.find(iItr->getEdgeVertex1()) == h.end() ||
				   iItr->getTime() < is[h[iItr->getEdgeVertex1()]].getTime()))
			h[iItr->getEdgeVertex1()] = idx;

		if (iItr->getAlpha() >= 0.5 &&
		   (h.find(iItr->getEdgeVertex2()) == h.end() ||
				   iItr->getTime() < is[h[iItr->getEdgeVertex2()]].getTime()))
			h[iItr->getEdgeVertex2()] = idx;
	}

	// TODO: At this point one could partition the intersections into disjoint
	// regions, in order to create one interference volume per region.
	// See the paper for details
	//

	// Create interference volume
	//
	ivs.push_back(InterferenceVolume(this, is, h));
}

void Curve::removeInterference(InterferenceVolumes &ivs)
{
	double *ctrl=0; // control mesh vertices
	getControlConfiguration(ctrl);

	set<size_t> movable;
	getMovableControlVertices(movable);

	for (InterferenceVolumesIterator ivItr=ivs.begin(); ivItr!=ivs.end(); ++ivItr)
	{
		ivItr->computeHandleGradient(this, movable);

		double p = ivItr->getPressure();
		for (set<size_t>::iterator sItr=movable.begin(); sItr!=movable.end(); ++sItr)
		{
			for (size_t i=0; i<2; ++i)
				ctrl[2*(*sItr)+i] += p * ivItr->getHandleGradients()[2*(*sItr)+i];
		}
	}
}

bool Curve::getCurveIntersections(Intersections &is)
{
	double *x0=0; // Start configuration
	double *x1=0; // End configuration
	size_t nbrVerts = getNbrVertices();
	getStartConfiguration(x0);
	getEndConfiguration(x1);

	unsigned int *e=0; // Edge i consists of vertices e[2*i] and e[2*i+1]
	size_t nbrEdges;
	getEdgeIndices(nbrEdges, e);

	static Hash grid;

	// Resize grid
	//
	Vector3d mn = Vector3d( 1e100,  1e100, 0.0);
	Vector3d mx = Vector3d(-1e100, -1e100, 1.0);
	for (size_t i=0; i<nbrVerts; ++i)
	{
		mn[0] = min(x0[2*i+0], min(x1[2*i+0], mn[0]));
		mn[1] = min(x0[2*i+1], min(x1[2*i+1], mn[1]));
		mx[0] = max(x0[2*i+0], max(x1[2*i+0], mx[0]));
		mx[1] = max(x0[2*i+1], max(x1[2*i+1], mx[1]));
	}

	// TODO: For increased performance cells should probably be sized based
	// on average edge length (or something)
	//
	double cellSize = 0.1;
	grid.resize(mn, mx, cellSize);

	// Insert points into grid
	//
	double *xStart = &x0[0];
	double *xEnd   = &x1[0];
	for (int i=0; i<(int)nbrVerts; ++i)
	{
		Vector3d xMin, xMax;
		for (size_t j=0; j<2; ++j)
		{
			xMin[j] = min(xStart[j], xEnd[j]) - getMinimumSeparation() - 1.0e-6;
			xMax[j] = max(xStart[j], xEnd[j]) + getMinimumSeparation() + 1.0e-6;
		}
		xMin[2] = -getMinimumSeparation();
		xMax[2] =  getMinimumSeparation();

		// To differentiate between vertices and edges, vertex indices
		// will always be negative (and offset by 1, to prevent ambiguity
		// for the 0-th element)
		//
		grid.addElement(xMin, xMax, -i-1);

		xStart += 2;
		xEnd   += 2;
	}

	// Insert edges into grid
	//
	for (size_t i=0; i<nbrEdges; ++i)
	{
		double *x10 = &x0[2*e[2*i+0]];
		double *x11 = &x1[2*e[2*i+0]];
		double *x20 = &x0[2*e[2*i+1]];
		double *x21 = &x1[2*e[2*i+1]];
		
		Vector3d xMin, xMax;
		for (size_t j=0; j<2; ++j)
		{
			xMin[j] = min(x10[j], min(x11[j], min(x20[j], x21[j]))) - getMinimumSeparation() - 1.0e-6;
			xMax[j] = max(x10[j], min(x11[j], min(x20[j], x21[j]))) + getMinimumSeparation() - 1.0e-6;
		}
		xMin[2] = -getMinimumSeparation();
		xMax[2] =  getMinimumSeparation();

		grid.addElement(xMin, xMax, i+1);
	}

	// Query potential intersections
	//
	Candidates hits;
	grid.getVertexEdgePairs(e, hits);


	// Low-level intersection test
	//
	is.clear();
	for (CandidatesIterator cItr=hits.begin(); cItr!=hits.end(); ++cItr)
	{
		double s, t;
		if (getIntersection(&x0[2*cItr->first],         &x1[2*cItr->first],
							&x0[2*e[2*cItr->second+0]], &x1[2*e[2*cItr->second+0]],
							&x0[2*e[2*cItr->second+1]], &x1[2*e[2*cItr->second+1]], s, t))
			is.push_back(Intersection(cItr->first, e[2*cItr->second], e[2*cItr->second+1], s, t));
	}

	return !is.empty();
}

bool Curve::getIntersection(double *x00, double *x01,
							double *x10, double *x11, double *x20, double *x21, double &s, double &t)
{
	Vector2d x1(x00[0], x00[1]);
	Vector2d x2(x10[0], x10[1]);
	Vector2d x3(x20[0], x20[1]);
	
	Vector2d v1(x01[0]-x00[0], x01[1]-x00[1]);
	Vector2d v2(x11[0]-x10[0], x11[1]-x10[1]);
	Vector2d v3(x21[0]-x20[0], x21[1]-x20[1]);

	Vector2d x1x2 = x1 - x2;
	Vector2d v1v2 = v1 - v2;
	Vector2d x3x2 = x3 - x2;
	Vector2d v3v2 = v3 - v2;

	swap(x3x2[0], x3x2[1]); x3x2[1] *= -1.0;
	swap(v3v2[0], v3v2[1]); v3v2[1] *= -1.0;

	double coeffs[3];
	coeffs[0] = v1v2.dot(v3v2);
	coeffs[1] = x1x2.dot(v3v2) + v1v2.dot(x3x2);
	coeffs[2] = x1x2.dot(x3x2);

	if (coeffs[0] == 0.0 && coeffs[1] == 0.0)
		return false;

	// TODO: Coefficient pruning tests could speed this up
	//

	int degree = 2;
	if (coeffs[0] == 0.0)
		degree = 1;

	// Find roots of at^2 + bt + c = 0
	//
	RootFinder rf;
	double rl[2], im[2];
	int cnt = rf.rpoly(&coeffs[2-degree], degree, rl, im);
	for (int i=0; i<cnt; ++i)
	{
		if (im[i] == 0.0 && rl[i] >= 0.0 && rl[i] <= 1.0)
		{
			t = rl[i];

			// Check rl[i] for intersection
			//
			Vector2d p = x1 + v1 * t;
			Vector2d a = x2 + v2 * t;
			Vector2d b = x3 + v3 * t;

			Vector2d ab = b - a;
			Vector2d ap = p - a;

			s = ab.dot(ap) / ab.dot(ab);

			// Clamp to edge endpoints
			//
			if (s < 0.0) s = 0.0;
			if (s > 1.0) s = 1.0;

			// NOTE: If the code misses intersections, it's probably because of
			// this tolerance. There is a Eurographics 2012 paper which attempts to
			// address this, but for now it can be tricky to choose a reliable epsilon.
			//
			double d = (p - (a + ab * s)).squaredNorm();
			if (d < 1.0e-4) 
			{
				// To ensure a conservative interference volume, slightly enlarge
				// it by returning an earlier time
				//
				t = max(0.0, t-0.01);
				return true;
			}
		}
	}

	return false;
}


} // end namespace IAGM

