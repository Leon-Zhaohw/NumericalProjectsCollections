// InterferenceVolume.cpp
//

#include "InterferenceVolume.h"
#include "Curve.h"
#include "rpoly.h"

#include <Eigen/Geometry>



namespace IAGM
{

using namespace std;
using namespace __gnu_cxx;
using namespace Eigen;



InterferenceVolume::InterferenceVolume(Curve *C, Intersections &is, hash_map<size_t,size_t> &vertIndices)
{
	computeInterferenceVolume(C, is, vertIndices);
	computeInterferenceVolumeGradient(C, is, vertIndices);
}

void InterferenceVolume::computeInterferenceVolume(Curve *C, Intersections &is, hash_map<size_t,size_t> &vertIndices)
{
	double *x0;
	double *x1;
	C->getStartConfiguration(x0);
	C->getEndConfiguration(x1);

	m_volume = 0.0;
	for (hash_map<size_t, size_t>::iterator iItr=vertIndices.begin(); iItr!=vertIndices.end(); ++iItr)
	{
		Intersection& i = is[iItr->second];

		Vector2d q[3] = { Vector2d(&x0[2*i.getVertex()]),
						  Vector2d(&x0[2*i.getEdgeVertex1()]),
						  Vector2d(&x0[2*i.getEdgeVertex2()]) };

		Vector2d p[3] = { Vector2d(&x1[2*i.getVertex()]),
						  Vector2d(&x1[2*i.getEdgeVertex1()]),
						  Vector2d(&x1[2*i.getEdgeVertex2()]) };

		Vector2d v[3] = { p[0] - q[0], p[1] - q[1], p[2] - q[2] };

		// Which vertex does this summand belong to
		//
		size_t idx=0;
		if (iItr->first == i.getEdgeVertex1())
			idx = 1;
		else if (iItr->first == i.getEdgeVertex2())
			idx = 2;

		// Intersection normal is the edge perpendicular
		//
		Vector2d e((q[2] + v[2] * i.getTime()) - (q[1] + v[1] * i.getTime()));
		Vector2d n(e[1], -e[0]);

		double tmp = (1.0 - i.getTime()) * v[idx].dot(n);

		// Always orient so that the summand is negative
		//
		if (tmp > 0.0)
			m_volume -= tmp;
		else
			m_volume += tmp;
	}
}

/*
// Finite difference gradient computation is extremely useful for debugging
//
void InterferenceVolume::computeInterferenceVolumeGradient(Curve *C, Intersections &is, hash_map<size_t, size_t> &vertIndices)
{
	double *x0;
	double *x1;
	C->getStartConfiguration(x0);
	C->getEndConfiguration(x1);

	m_dVdq = VectorXd::Zero(2 * C->getNbrVertices());
	for (hash_map<size_t, size_t>::iterator iItr=vertIndices.begin(); iItr!=vertIndices.end(); ++iItr)
	{
		Intersection& i = is[iItr->second];

		Vector2d q[3] = { Vector2d(&x0[2*i.getVertex()]),
						  Vector2d(&x0[2*i.getEdgeVertex1()]),
						  Vector2d(&x0[2*i.getEdgeVertex2()]) };

		Vector2d p[3] = { Vector2d(&x1[2*i.getVertex()]),
						  Vector2d(&x1[2*i.getEdgeVertex1()]),
						  Vector2d(&x1[2*i.getEdgeVertex2()]) };

		size_t idx=0;
		if (iItr->first == i.getEdgeVertex1())
			idx = 1;
		else if (iItr->first == i.getEdgeVertex2())
			idx = 2;

		Vector2d vTmp[3] = { p[0]-q[0], p[1]-q[1], p[2]-q[2] };
		
		double sign = getSign(idx, i.getTime(), i.getAlpha(), q, vTmp);
		
		Vector2d tmp = p[idx]-q[idx];
		double V = sign*getVolume(tmp, q[0], q[1], q[2], p[0]-q[0], p[1]-q[1], p[2]-q[2], 1.0, i.getTime());

		Vector2d dv[3];
		const double eps = 0.00001;
		for (size_t j=0; j<3; ++j)
		{
			for (size_t k=0; k<2; ++k)
			{
				Vector2d v[3] = { p[0] - q[0], p[1] - q[1], p[2] - q[2] };
				v[j][k] += eps;

				// Recompute time of intersection
				//
				double t=1.0;

				Vector2d x1x2 = q[0] - q[1];
				Vector2d v1v2 = v[0] - v[1];
				Vector2d x3x2 = q[2] - q[1];
				Vector2d v3v2 = v[2] - v[1];

				swap(x3x2[0], x3x2[1]); x3x2[1] *= -1.0;
				swap(v3v2[0], v3v2[1]); v3v2[1] *= -1.0;

				double coeffs[3];
				coeffs[0] = v1v2.dot(v3v2);
				coeffs[1] = x1x2.dot(v3v2) + v1v2.dot(x3x2);
				coeffs[2] = x1x2.dot(x3x2);

				int degree = 2;
				if (coeffs[0] == 0.0)
					degree = 1;

				double s=0.0;

				// Find roots of at^2 + bt + c = 0
				//
				RootFinder rf;
				double rl[2], im[2];
				int cnt = rf.rpoly(&coeffs[2-degree], degree, rl, im);
				for (int l=0; l<cnt; ++l)
				{
					if (im[l] == 0.0 && rl[l] >= 0.0 && rl[l] <= 1.0)
					{
						t = rl[l];

						// Check rl[i] for intersection
						//
						Vector2d pt = q[0] + v[0] * t;
						Vector2d a  = q[1] + v[1] * t;
						Vector2d b  = q[2] + v[2] * t;

						Vector2d ab = b - a;
						Vector2d ap = pt - a;

						s = ab.dot(ap) / ab.dot(ab);
						if (s < 0.0) s = 0.0;
						if (s > 1.0) s = 1.0;

						double d  = (pt - (a + ab * s)).squaredNorm();
						if (d < 1.0e-4)
						{
							t = max(0.0, t-0.01);
							break;
						}
					}
				}

				sign = getSign(idx, t, s, q, v);
				double Vnew = sign * getVolume(v[idx], q[0], q[1], q[2], v[0], v[1], v[2], 1.0, t);

				dv[j][k] = (Vnew - V) / eps;
			}
		}
		cout << dv[0].transpose() << "; " << dv[1].transpose() << "; " << dv[2].transpose() << endl;

		m_dVdq.segment<2>(2*i.getVertex())      += dv[0];
		m_dVdq.segment<2>(2*i.getEdgeVertex1()) += dv[1];
		m_dVdq.segment<2>(2*i.getEdgeVertex2()) += dv[2];

		{
			Intersection& i = is[iItr->second];

		Vector2d q[3] = { Vector2d(&x0[2*i.getVertex()]),
						  Vector2d(&x0[2*i.getEdgeVertex1()]),
						  Vector2d(&x0[2*i.getEdgeVertex2()]) };

		Vector2d p[3] = { Vector2d(&x1[2*i.getVertex()]),
						  Vector2d(&x1[2*i.getEdgeVertex1()]),
						  Vector2d(&x1[2*i.getEdgeVertex2()]) };

		Vector2d v[3] = { p[0] - q[0], p[1] - q[1], p[2] - q[2] };

		size_t idx=0;
		if (iItr->first == i.getEdgeVertex1())
			idx = 1;
		else if (iItr->first == i.getEdgeVertex2())
			idx = 2;

		// Intersection normal is the edge perpendicular
		//
		Vector2d e((q[2] + v[2] * i.getTime()) - (q[1] + v[1] * i.getTime()));
		Vector2d n(e[1], -e[0]);

		Vector2d dndt;
		getDnDt(i.getTime(), q, v, dndt);

		// Compute dV_i/dt = (1-t) \delta r^T dn/dt - \delta r^T n
		//
		double dVdt = (1.0 - i.getTime()) * v[idx].dot(dndt) - v[idx].dot(n);

		Vector2d dtdq1, dtdq2, dtdq3;
		getTimeGradient(i, 0, x0, x1, dtdq1);
		getTimeGradient(i, 1, x0, x1, dtdq2);
		getTimeGradient(i, 2, x0, x1, dtdq3);
		
		double sign = getSign(idx, i.getTime(), i.getAlpha(), q, v);

		Vector2d dvdq[3] = { sign * dVdt * dtdq1, sign*dVdt*dtdq2, sign*dVdt*dtdq3 };

		cout << dvdq[0].transpose() << "; " << dvdq[1].transpose() << "; " << dvdq[2].transpose() << endl << endl;

		}
	}
}
*/

void InterferenceVolume::computeInterferenceVolumeGradient(Curve *C, Intersections &is, hash_map<size_t, size_t> &vertIndices)
{
	double *x0;
	double *x1;
	C->getStartConfiguration(x0);
	C->getEndConfiguration(x1);

	// These gradients are with respect to low-level DOFs
	//
	m_dVdq = VectorXd::Zero(2 * C->getNbrVertices());
	for (hash_map<size_t, size_t>::iterator iItr=vertIndices.begin(); iItr!=vertIndices.end(); ++iItr)
	{
		Intersection& i = is[iItr->second];

		Vector2d q[3] = { Vector2d(&x0[2*i.getVertex()]),
						  Vector2d(&x0[2*i.getEdgeVertex1()]),
						  Vector2d(&x0[2*i.getEdgeVertex2()]) };

		Vector2d p[3] = { Vector2d(&x1[2*i.getVertex()]),
						  Vector2d(&x1[2*i.getEdgeVertex1()]),
						  Vector2d(&x1[2*i.getEdgeVertex2()]) };

		Vector2d v[3] = { p[0] - q[0], p[1] - q[1], p[2] - q[2] };

		size_t idx=0;
		if (iItr->first == i.getEdgeVertex1())
			idx = 1;
		else if (iItr->first == i.getEdgeVertex2())
			idx = 2;

		// Intersection normal is the edge perpendicular
		//
		Vector2d e((q[2] + v[2] * i.getTime()) - (q[1] + v[1] * i.getTime()));
		Vector2d n(e[1], -e[0]);

		Vector2d dndt;
		getDnDt(v, dndt);

		// To compute the derivative a volume summand V_i with respect to a vertex
		// q_j we do something different from what's in the paper. Using the chain
		// rule we can write it as
		// dV_i / dq_j = dV_i/dt * dt/dq_j,
		// which is much easier than the gradient in the paper.
		//

		// Compute dV_i/dt = (1-t) \delta r^T dn/dt - \delta r^T n
		//
		double dVdt = (1.0 - i.getTime()) * v[idx].dot(dndt) - v[idx].dot(n);

		// Derivative of each vertex involved in the intersection with respect to t
		//
		Vector2d dtdq1, dtdq2, dtdq3;
		getTimeGradient(i, 0, x0, x1, dtdq1);
		getTimeGradient(i, 1, x0, x1, dtdq2);
		getTimeGradient(i, 2, x0, x1, dtdq3);
		
		double sign = getSign(idx, i.getTime(), i.getAlpha(), q, v);

		m_dVdq.segment<2>(2*i.getVertex())      += sign * dVdt * dtdq1;
		m_dVdq.segment<2>(2*i.getEdgeVertex1()) += sign * dVdt * dtdq2;
		m_dVdq.segment<2>(2*i.getEdgeVertex2()) += sign * dVdt * dtdq3;
	}
}

double InterferenceVolume::getSign(size_t idx, double t, double a, Vector2d *q, Vector2d *v)
{
	Vector2d z1 = q[1] + v[1] * t;
	Vector2d z2 = q[2] + v[2] * t;
	Vector2d e(z2 - z1);
	Vector2d n(e[1], -e[0]);

	double sign = (idx) ? -1.0 : 1.0;
	double rv = (v[0] - (1.0 - a) * v[1] - a * v[2]).dot(n);
	if (rv > 0.0)
		return -sign;
	
	// Should never have an intersection between two elements
	// with zero relative velocity
	//
	assert(rv != 0.0);
	
	return sign;
}

void InterferenceVolume::computeHandleGradient(Curve *C, std::set<size_t> &movable)
{
	// Project our gradient into the modeling subspace
	//
	m_dVdp = VectorXd::Zero(2*C->getNbrCtrlVertices());
	C->getSubspaceVector(m_dVdq.data(), m_dVdp.data());

	set<size_t> all;
	for (size_t i=0; i<C->getNbrCtrlVertices(); ++i)
		all.insert(i);
	
	// Get the immovable vertices
	//
	set<size_t> immovable;
	set_difference(all.begin(), all.end(), movable.begin(), movable.end(),
				   inserter(immovable, immovable.end()));

	// Zero out gradient for unselected vertices (they don't move)
	//
	for (set<size_t>::iterator sItr=immovable.begin(); sItr!=immovable.end(); ++sItr)
		for (size_t i=0; i<2; ++i)
			m_dVdp[2*(*sItr)+i] = 0.0;
}

double InterferenceVolume::getBarycentricArea(size_t idx, size_t nbrVertices, double *x)
{
	// TODO: This assumes one spline, essentially. Ideally it should check if the vertex
	// is at the beginning or end of a new spline, not just at 0 or nbrVertices
	//
	double A = 0.0;

	if (idx > 0)
		A += 0.5 * Vector2d(x[2*(idx-1)] - x[2*idx], x[2*(idx-1)+1] - x[2*idx+1]).norm();
	
	if (idx < (nbrVertices-1))
		A += 0.5 * Vector2d(x[2*idx] - x[2*(idx+1)], x[2*idx+1] - x[2*(idx+1)+1]).norm();

	return A;
}

/*
double InterferenceVolume::getVolume(Vector2d &r,
									 Vector2d x0, Vector2d x1, Vector2d x2,
				     				 Vector2d v0, Vector2d v1, Vector2d v2,
									 double A, double t)
{
	Vector2d z1 = x1 + v1 * t;
	Vector2d z2 = x2 + v2 * t;
	Vector2d e(z2 - z1);
	Vector2d n(e[1], -e[0]);

	return (1.0 - t) * r.dot(n) * A;
}
*/

void InterferenceVolume::getDnDt(Vector2d *v, Vector2d &dndt)
{
	Vector2d e = v[2] - v[1];
	dndt = Vector2d(e[1], -e[0]);
}

void InterferenceVolume::getTimeGradient(Intersection &i, size_t idx, double *x0, double *x1, Vector2d &dtdq)
{
	Vector2d q1(x0[2*i.getVertex()],      x0[2*i.getVertex()+1]);
	Vector2d q2(x0[2*i.getEdgeVertex1()], x0[2*i.getEdgeVertex1()+1]);
	Vector2d q3(x0[2*i.getEdgeVertex2()], x0[2*i.getEdgeVertex2()+1]);
	Vector2d p1(x1[2*i.getVertex()],      x1[2*i.getVertex()+1]);
	Vector2d p2(x1[2*i.getEdgeVertex1()], x1[2*i.getEdgeVertex1()+1]);
	Vector2d p3(x1[2*i.getEdgeVertex2()], x1[2*i.getEdgeVertex2()+1]);

	Vector2d v1(p1 - q1);
	Vector2d v2(p2 - q2);
	Vector2d v3(p3 - q3);

	Vector2d x1x2 = q1 - q2;
	Vector2d v1v2 = v1 - v2;
	Vector2d x3x2 = q3 - q2;
	Vector2d v3v2 = v3 - v2;
				
	swap(x3x2[0], x3x2[1]); x3x2[1] *= -1.0;
	swap(v3v2[0], v3v2[1]); v3v2[1] *= -1.0;

	double a = v1v2.dot(v3v2);
	double b = x1x2.dot(v3v2) + v1v2.dot(x3x2);

	double scale = -1.0 / (2.0 * a * i.getTime() + b);
	
	Vector2d z1 = q2 + v2 * i.getTime();
	Vector2d z2 = q3 + v3 * i.getTime();
	Vector2d e(z2 - z1);
	Vector2d n(e[1], -e[0]);

	Vector2d num;
	if (idx==0)
		num = n;
	else if (idx==1)
		num = -(1.0 - i.getAlpha()) * n;
	else
		num = -i.getAlpha() * n;

/*
	switch (idx)
	{
	case 0:
		num = v3v2 * i.getTime() * i.getTime() + x3x2 * i.getTime();
		break;
	case 1:
		num = (-v1v2 - v3v2) * i.getTime() * i.getTime() - (x1x2 + x3x2) * i.getTime();
		break;
	case 2:
		num = v1v2 * i.getTime() * i.getTime() + x1x2 * i.getTime();
		break;
	}
*/

	dtdq = scale * num;

//	cout << dtdq.transpose() << "; ";

/*
{
	Vector2d q[3] = { Vector2d(x0[2*i.getVertex()],      x0[2*i.getVertex()+1]),
					  Vector2d(x0[2*i.getEdgeVertex1()], x0[2*i.getEdgeVertex1()+1]),
					  Vector2d(x0[2*i.getEdgeVertex2()], x0[2*i.getEdgeVertex2()+1]) };

	Vector2d p[3] = { Vector2d(x1[2*i.getVertex()],      x1[2*i.getVertex()+1]),
					  Vector2d(x1[2*i.getEdgeVertex1()], x1[2*i.getEdgeVertex1()+1]),
					  Vector2d(x1[2*i.getEdgeVertex2()], x1[2*i.getEdgeVertex2()+1]) };

	Vector2d dt;
	const double eps = 0.0001;
	for (size_t j=0; j<2; ++j)
	{
		// Perturb end position by epsilon
		//
		Vector2d v[3] = { p[0] - q[0], p[1] - q[1], p[2] - q[2] };
		v[idx][j] += eps;

		Vector2d x1x2 = q[0] - q[1];
		Vector2d v1v2 = v[0] - v[1];
		Vector2d x3x2 = q[2] - q[1];
		Vector2d v3v2 = v[2] - v[1];
	
		swap(x3x2[0], x3x2[1]); x3x2[1] *= -1.0;
		swap(v3v2[0], v3v2[1]); v3v2[1] *= -1.0;

		double coeffs[3];
		coeffs[0] = v1v2.dot(v3v2);
		coeffs[1] = x1x2.dot(v3v2) + v1v2.dot(x3x2);
		coeffs[2] = x1x2.dot(x3x2);

		if (coeffs[0] == 0.0 && coeffs[1] == 0.0)
		{
			//dt[j] = 0.0;
			assert(0);
			dtdq[j] = (1.0 - i.getTime()) / eps;
			continue;
		}

		int degree = 2;
		if (coeffs[0] == 0.0)
			degree = 1;

		// Find roots of at^2 + bt + c = 0
		//
		RootFinder rf;
		double rl[2], im[2];
		int cnt = rf.rpoly(&coeffs[2-degree], degree, rl, im);
		for (int k=0; k<cnt; ++k)
		{
			if (im[k] == 0.0 && rl[k] >= 0.0 && rl[k] <= 1.0)
			{
				double t = rl[k];

				// Check rl[i] for intersection
				//
				Vector2d pt = q[0] + v[0] * t;
				Vector2d a = q[1] + v[1] * t;
				Vector2d b = q[2] + v[2] * t;

				Vector2d ab = b - a;
				Vector2d ap = pt - a;

				double s = ab.dot(ap) / ab.dot(ab);

				if (s < 0.0) s = 0.0;
				if (s > 1.0) s = 1.0;

				double d = (pt - (a + ab * s)).squaredNorm();
				if (d < 1.0e-4)
				{
					t = max(0.0, t-0.01);
					dtdq[j] = (t - i.getTime()) / eps;
					break;
				}
			}
		}
	}
	}
*/
}

} // end namespace IAGM

