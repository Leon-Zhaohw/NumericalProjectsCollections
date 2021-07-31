// Copyright (c) 2011, Regents of the University of Utah
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the <organization> nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef SMOOTHINGGRID_is_defined
#define SMOOTHINGGRID_is_defined

#include "slMatrix.H"
#include "kdTree.H"
#include "slVector.H"
#include "slUtil.H"
#include <vector>

class SmoothingGrid {
public:
    static const unsigned int nneighbors = 32;
    static const unsigned int VARIABLE_RADIUS = 2;
    static const unsigned int NEIGHBOR_ANISOTROPY = 4;
    static const unsigned int VELOCITY_ANISOTROPY = 8;
    static const unsigned int VERBOSE = 16;
	  double h;
    int nx, ny, nz, count;
    SlVector3 bbMin, bbMax;
    SlArray3D<double> phi, phiB;
    ~SmoothingGrid();
    SmoothingGrid(double h, double r_min, double r_max, double r_init, double aniGain, double maxStretch, unsigned int flags,
                  const std::vector<SlVector3> &particles, const std::vector<double> &radii, const std::vector<SlVector3> &velocities);
	  bool doMeanCurvatureSmoothing(double dt, int ntimesteps, int redistanceFrequency);
    bool doLaplacianSmoothing(double dt, int ntimesteps, int redistanceFrequency);
	  bool doBiharmonicSmoothing(double dt, int ntimesteps, int redistanceFrequency, double cgThreshold, int cgMaxIters);
    bool doBiharmonicSmoothing(double dt, int ntimesteps, int redistanceFrequency);

private:
    double rmin, rmax, rinit;
    unsigned int flags;
    KDTree *kdtree;
    SlArray3D<double> biharmonic, laplacian, meanCurvature, phi_min, phi_max, tempPhi, prevPhi;
    SlArray3D<char> accepted, snapStatus;

    bool computeG(const std::vector<SlVector3> &particles, double rmax, std::vector<SlMatrix3x3> &G, std::vector<double> &factors, double maxStretch);
    bool computeG(const std::vector<SlVector3> &particles, const std::vector<SlVector3> &velocities, double gain, std::vector<SlMatrix3x3> &G, std::vector<double> &factors, double maxStretch);
    bool rasterize(const std::vector<SlVector3> &particles);
    bool rasterize(const std::vector<SlVector3> &particles, const std::vector<double> &radii);
    bool rasterize(const std::vector<SlVector3> &particles, const std::vector<SlMatrix3x3> &Gs, const std::vector<double> &factors);
    bool computeLaplacian();
    bool computeLaplacianCG(const SlArray3D<double> &x, SlArray3D<unsigned char> &marked);
    bool computeBiharmonic();
    bool computeMeanCurvature();
    bool stepMeanCurvature(double dt);
    bool stepLaplacian(double dt);
    double stepBiharmonic(double dt);

	  bool CG(const SlArray3D<double> &b, SlArray3D<double> &x, SlArray3D<unsigned char> &markedB, SlArray3D<unsigned char> &markedL,double &tol, double lambdadt, double threshold, int cgMaxIters);
    bool applyBiharmonic(const SlArray3D<double> &x, SlArray3D<double> &y, SlArray3D<unsigned char> &markedB, double dt) const;
    bool sameActiveSet(SlArray3D<unsigned char> &x, SlArray3D<unsigned char> &y);
    double innerProduct(SlArray3D<double> &x, SlArray3D<double> &y, SlArray3D<unsigned char> &marked);
    double sqrTwoNorm(SlArray3D<double> &x, SlArray3D<unsigned char> &marked);
};

#endif
