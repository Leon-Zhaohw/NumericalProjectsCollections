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


#include "smoothingGrid.H"
#include "slUtil.H"
#include <sys/time.h>
#include <float.h>
#include <fstream>
#include <iostream>
#include <cmath>

double SmoothingGrid::innerProduct(SlArray3D<double> &x, SlArray3D<double> &y, SlArray3D<unsigned char> &marked) {
    double sum = 0.0;
    int n = nx*ny*nz;
    double *dptrx = &(x(0,0,0)), *dptry = &(y(0,0,0));
    unsigned char *dptrm = &(marked(0,0,0));
    for(int i=0; i<n; i++, dptrx++, dptry++, dptrm++){
        if (*dptrm) {
            sum += (*dptrx) * (*dptry);
        }
    }
    return sum;
}

bool SmoothingGrid::sameActiveSet(SlArray3D<unsigned char> &x, SlArray3D<unsigned char> &y) {
    int sum;
    int n = nx*ny*nz;
    unsigned char *dptrx = &(x(0,0,0)), *dptry = &(y(0,0,0));
    for(int i=0; i<n; i++, dptrx++, dptry++){
        sum = (*dptrx)^(*dptry);
        if (sum) return false;
    }
    return true;
}

inline double sqr(double x) { return (x * x); };

double SmoothingGrid::sqrTwoNorm(SlArray3D<double> &x, SlArray3D<unsigned char> &marked) {
    double twoNorm = 0.0;
    int n = nx*ny*nz;
    double *dptrx = &(x(0,0,0));
    unsigned char *dptrm = &(marked(0,0,0));
    for(int i=0; i<n; i++, dptrx++, dptrm++){
        if (*dptrm) {
					  twoNorm += sqr(*dptrx);
        }
    }
    return twoNorm;
}

inline double cdX(int i, int j, int k, SlArray3D<double> &phi, int nx, double h) {
    return (phi(i+1,j,k) - phi(i-1,j,k)) / (2*h);
    return (-phi(i+2, j, k) + 8 * (phi(i+1, j, k) - 8 * phi(i-1, j, k)) + phi(i-2, j, k)) / (12 * h);
}

inline double cdY(int i, int j, int k, SlArray3D<double> &phi, int ny, double h) {
    return (phi(i,j+1,k) - phi(i,j-1,k)) / (2*h);
    return (-phi(i, j+2, k) + 8 * (phi(i, j+1, k) - 8 * phi(i, j-1, k)) + phi(i, j-2, k)) / (12 * h);
}

inline double cdZ(int i, int j, int k, SlArray3D<double> &phi, int nz, double h) {
    return (phi(i,j,k+1) - phi(i,j,k-1)) / (2*h);
    return (-phi(i, j, k+2) + 8 * (phi(i, j, k+1) - 8 * phi(i, j, k-1)) + phi(i, j, k-2)) / (12 * h);
}

bool SmoothingGrid::computeG(const std::vector<SlVector3> &particles, const std::vector<SlVector3> &velocities,
                             double gain, std::vector<SlMatrix3x3> &G, std::vector<double> &factors, double maxStretch) {
    double maxS = cbrt(maxStretch) - 1;

    G.resize(particles.size());
    factors.resize(particles.size());

    for (unsigned int i = 0; i<particles.size(); i++) {
        SlVector3 v(velocities[i]);
        double s = mag(velocities[i]);
        if (s > 0.01) {
            v /= s;
        } else {
            G[i].setIdentity();
            factors[i] = 1;
            continue;
        }

        s = std::min(maxS, s*gain);

        double e12 = 1.0+s;
        double e0 = 1.0/sqr(e12);

        // find two vectors orthogonal to v
        SlVector3 v1(cross(v,SlVector3(1,0,0)));
        double m = sqrMag(v1);
        if (m < 0.01) {
            v1.set(cross(v,SlVector3(0,1,0)));
            m = sqrMag(v1);
        }
        v1 /= sqrt(m);
        SlVector3 v2 = cross(v, v1);
        normalize(v2);


        G[i].set(e0*v[0], e0*v[1], e0*v[2],
                 e12*v1[0], e12*v1[1], e12*v1[2],
                 e12*v2[0], e12*v2[1], e12*v2[2]);
        // don't really need the left rotation, it doesn't change lengths

        factors[i] = 1/e0;
    }
    return true;
}

bool SmoothingGrid::computeG(const std::vector<SlVector3> &particles, double rmax, std::vector<SlMatrix3x3> &G, std::vector<double> &factors, double maxStretch) {
    // The following three parameters could be exposed to the user
    double minneighbors = 20.1; // minimum number of neighbors before G = I
    double r = 4.0*rmax; // This is the search radius, default is 4.

    G.resize(particles.size());
    factors.resize(particles.size());
    double *weights = new double[nneighbors];
    std::vector<int> neighbors;
    neighbors.reserve(nneighbors);

    for (unsigned int i = 0; i<particles.size(); i++) {
        SlVector3 xiw(0.0);
        double weightsum = 0.0;
        const SlVector3 &pos = particles[i];

        kdtree->neighbors(particles, pos, nneighbors, r, neighbors);

        for (unsigned int j=0; j<neighbors.size(); j++) {
            int &k = neighbors[j];
            const SlVector3 &npos = particles[k];
            double d = mag(pos - npos);
            double ratio = d/r;
            weights[j] = 1 - (ratio*ratio*ratio);
            weightsum += weights[j];
            xiw += weights[j]*npos;
        }

        if (weightsum < 0.01) {
            G[i].setIdentity();
            factors[i] = 1.0;
            continue;
        }

        xiw /= weightsum;
        G[i] = 0.0;
        for (unsigned int j=0; j<neighbors.size(); j++) {
            int k = neighbors[j];
            const SlVector3 &npos = particles[k];
            SlVector3 d = npos - xiw;
            G[i] += weights[j]*SlMatrix3x3(d[0]*d[0], d[0]*d[1], d[0]*d[2],
                                           d[1]*d[0], d[1]*d[1], d[1]*d[2],
                                           d[2]*d[0], d[2]*d[1], d[2]*d[2]);

        }

        G[i] /= weightsum;
        SlMatrix3x3 vecs;
        SlVector3 vals;
        SlSymetricEigenDecomp(G[i], vals, vecs);

        int maxEval = 0;
        if (vals[1] > vals[0]) maxEval = 1;
        if (vals[2] > vals[maxEval]) maxEval = 2;
        double temp = vals[0];  vals[0] = vals[maxEval]; vals[maxEval] = temp;
        SlVector3 tempV(vecs(0, 0), vecs(1, 0), vecs(2, 0));
        vecs(0,0) = vecs(0, maxEval);vecs(1,0) = vecs(1, maxEval);vecs(2,0) = vecs(2,maxEval);
        vecs(0, maxEval) = tempV[0]; vecs(1, maxEval) = tempV[1]; vecs(2, maxEval) = tempV[2];

        double maxT = vals[0] / maxStretch;
        vals[0] = 1.0 / (vals[0]);
        vals[1] = 1.0 / (std::max<double>(vals[1], maxT));
        vals[2] = 1.0 / (std::max<double>(vals[2], maxT));
				
        vals /= cbrt(vals[0]*vals[1]*vals[2]);  // make sure det(G) = 1
        // smooth falloff with decreasing # of neighbors
        double alpha = std::max(0.0, (neighbors.size()-minneighbors)/(nneighbors - minneighbors));
        vals[0] = pow(vals[0],alpha);
        vals[1] = pow(vals[1],alpha);
        vals[2] = pow(vals[2],alpha);
        vals /= cbrt(vals[0]*vals[1]*vals[2]);
        // don't really need the left rotation, it doesn't change lengths
        G[i] = vecs * diagonal(vals) * transpose(vecs);
        factors[i] = 1.0/vals[0];
    }
    return true;
}

bool SmoothingGrid::rasterize(const std::vector<SlVector3> &particles) {
    timeval startTime, endTime;
    gettimeofday(&startTime, NULL);
    int width = (int) ceil(rmax / h) + 1;
    phi = DBL_MAX;
    //phiB = DBL_MAX;
    for(unsigned int p = 0; p < particles.size(); p++) {

        SlVector3 pos(particles[p]);
        SlInt3 bin;
        bin[0] = (int) ((pos[0] - bbMin[0]) / h);
        bin[1] = (int) ((pos[1] - bbMin[1]) / h);
        bin[2] = (int) ((pos[2] - bbMin[2]) / h);

        unsigned int imax = std::min(bin[0] + width + 2, nx);
        unsigned int jmax = std::min(bin[1] + width + 2, ny);
        unsigned int kmax = std::min(bin[2] + width + 2, nz);

        for(unsigned int i = std::max(bin[0] - width, 0); i < imax; i++) {
            for(unsigned int j = std::max(bin[1] - width, 0); j < jmax; j++) {

                // for the inner most loop, the i and j coordinates do not change
                // as we are looking for sqrmag of the distance, we can precomupte
                // the first two terms and only worry about the last term during the loop
                // this reduces the computation to 4 adds, one multiply, and a min in each
                // iteration of the inner most loop.
                unsigned int k = std::max(bin[2]-width,0);
                double d = k*h+bbMin[2]-pos[2];
                double psum = sqr(i*h+bbMin[0]-pos[0]) + sqr(j*h+bbMin[1]-pos[1]);
                double *dptr = &(phi(i,j,k));

                for(; k < kmax; k++, d+=h, dptr++) {
                    (*dptr) = std::min((*dptr), psum + sqr(d));
                }
            }
        }
    }

    for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < nz; k++) {
                if (phi(i,j,k) == DBL_MAX) {
                    phi_min(i,j,k) = DBL_MAX;
                    phi_max(i,j,k) = DBL_MAX;
                } else {
                    double dist = sqrt(phi(i, j, k));
                    phi(i, j, k) = dist - rinit;
                    phi_min(i, j, k) = dist - rmin;
                    phi_max(i, j, k) = dist - rmax;
                }
            }
        }
    }
    gettimeofday(&endTime, NULL);
    return true;
}

bool SmoothingGrid::rasterize(const std::vector<SlVector3> &particles,
                              const std::vector<double> &radii) {
    phi = DBL_MAX;
    phi_min = DBL_MAX;
    phi_max = DBL_MAX;
    for(unsigned int p = 0; p < particles.size(); p++) {
        SlVector3 pos(particles[p]);
        double r = radii[p];
        double pinit = r*rinit;
        double pmax = r*rmax;
        double pmin = r*rmin;
        int width = (int) ceil(pmax / h) + 1;

        SlInt3 bin;
        bin[0] = (int) ((pos[0] - bbMin[0]) / h);
        bin[1] = (int) ((pos[1] - bbMin[1]) / h);
        bin[2] = (int) ((pos[2] - bbMin[2]) / h);

        unsigned int imax = std::min(bin[0] + width + 2, nx);
        unsigned int jmax = std::min(bin[1] + width + 2, ny);
        unsigned int kmax = std::min(bin[2] + width + 2, nz);

        for(unsigned int i = std::max(bin[0] - width, 0); i < imax; i++) {
            for(unsigned int j = std::max(bin[1] - width, 0); j < jmax; j++) {

                // for the inner most loop, the i and j coordinates do not change
                // as we are looking for sqrmag of the distance, we can precomupte
                // the first two terms and only worry about the last term during the loop
                // this reduces the computation to 4 adds, one multiply, and a min in each
                // iteration of the inner most loop.
                unsigned int k = std::max(bin[2]-width,0);
                double d = k*h+bbMin[2]-pos[2];
                double psum = sqr(i*h+bbMin[0]-pos[0]) + sqr(j*h+bbMin[1]-pos[1]);
                double *pptr = &(phi(i,j,k));
                double *mxptr = &(phi_max(i,j,k));
                double *mnptr = &(phi_min(i,j,k));

                for(; k < kmax; k++, d+=h, pptr++, mxptr++, mnptr++) {
                    double dist = sqrt(psum + sqr(d));
                    double val = dist - pinit;
                    if (val < (*pptr)) {
                        (*pptr) = val;
                        (*mxptr) = dist - pmax;
                        (*mnptr) = dist - pmin;
                    }
                }
            }
        }
    }

    return true;
}

bool SmoothingGrid::rasterize(const std::vector<SlVector3> &particles,
                              const std::vector<SlMatrix3x3> &Gs,
                              const std::vector<double> &factors) {
    phi = DBL_MAX;

    for(unsigned int p = 0; p < particles.size(); p++) {
        SlVector3 pos(particles[p]);
        const SlMatrix3x3 &G = Gs[p];

        SlVector3 Gc0(G(0,0), G(1,0), G(2,0));
        SlVector3 Gc1(G(0,1), G(1,1), G(2,1));
        SlVector3 Gc2(G(0,2), G(1,2), G(2,2));

        int width = (int)ceil(factors[p] * rmax / h) + 1 + 3;

        SlInt3 bin;
        bin[0] = (int) ((pos[0] - bbMin[0]) / h);
        bin[1] = (int) ((pos[1] - bbMin[1]) / h);
        bin[2] = (int) ((pos[2] - bbMin[2]) / h);

        unsigned int imax = std::min(bin[0] + width + 2, nx);
        unsigned int jmax = std::min(bin[1] + width + 2, ny);
        unsigned int kmax = std::min(bin[2] + width + 2, nz);

        for(unsigned int i = std::max(bin[0] - width, 0); i < imax; i++) {
            for(unsigned int j = std::max(bin[1] - width, 0); j < jmax; j++) {

                // for the inner most loop, the i and j coordinates do not change
                // as we are looking for sqrmag of the distance, we can precomupte
                // the first two terms and only worry about the last term during the loop
                // Here the terms are columns from the G matrix, we add a h* the third col
                // every time we increas k
                unsigned int k = std::max(bin[2]-width,0);
                SlVector3 psum = (i*h+bbMin[0]-pos[0])*Gc0 + (j*h+bbMin[1]-pos[1])*Gc1 + (k*h+bbMin[2]-pos[2])*Gc2;
                SlVector3 hGc2 = h*Gc2;
                double *dptr = &(phi(i,j,k));

                for(; k < kmax; k++, psum+=hGc2, dptr++) {
                    (*dptr) = std::min((*dptr), sqrMag(psum));
                }
            }
        }
    }


    for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < nz; k++) {
                if (phi(i,j,k) == DBL_MAX) {
                    phi_min(i, j, k) = DBL_MAX;
                    phi_max(i, j, k) = DBL_MAX;
                } else {
                    double dist = sqrt(phi(i, j, k));
                    phi(i, j, k) = dist - rinit;
                    phi_min(i, j, k) = dist - rmin;
                    phi_max(i, j, k) = dist - rmax;
                }
            }
        }
    }
    return true;
}

void redistance(SlArray3D<double>&phi, SlArray3D<double>&newPhi, SlArray3D<char>&accepted, double h);

SmoothingGrid::SmoothingGrid(double h, double rmin, double rmax, double rinit, double gain, double maxStretch, unsigned int flags,
                             const std::vector<SlVector3> &particles, const std::vector<double> &radii, const std::vector<SlVector3> &velocities) {
    this->h = h;
    this->rmin = rmin;
    this->rmax = rmax;
    this->rinit = rinit;
    this->flags = flags;

    // compute Gs if necessary
    std::vector<SlMatrix3x3> G;
    std::vector<double> factors;
    if (flags & NEIGHBOR_ANISOTROPY) {
        kdtree = new KDTree(particles);
        computeG(particles, rmax, G, factors, maxStretch);
    } else if (flags & VELOCITY_ANISOTROPY) {
        computeG(particles, velocities, gain, G, factors, maxStretch);
    }

    // compute bounding box and grid dimensions
    bbMin[0] = bbMin[1] = bbMin[2] = DBL_MAX;
    bbMax[0] = bbMax[1] = bbMax[2] = -DBL_MAX;

    for (std::vector<SlVector3>::const_iterator i=particles.begin(); i!=particles.end(); i++) {
        bbMin[0] = std::min(bbMin[0], (*i)[0]);
        bbMin[1] = std::min(bbMin[1], (*i)[1]);
        bbMin[2] = std::min(bbMin[2], (*i)[2]);
        bbMax[0] = std::max(bbMax[0], (*i)[0]);
        bbMax[1] = std::max(bbMax[1], (*i)[1]);
        bbMax[2] = std::max(bbMax[2], (*i)[2]);
    }

    // increase the bounding box by rmax + something a little bigger than the stencil size
    double maxFactor = 1;
    if (flags & SmoothingGrid::NEIGHBOR_ANISOTROPY) {
        maxFactor = sqrt(maxStretch);
    } else if (flags & SmoothingGrid::VELOCITY_ANISOTROPY) {
        maxFactor *= cbrt(maxStretch);
    }

    bbMin -= 5*h + maxFactor*rmax + 10*h;
    bbMax += 5*h + maxFactor*rmax + 10*h;

    bbMin[0]=floor((bbMin[0]/h))*h;
    bbMin[1]=floor((bbMin[1]/h))*h;
    bbMin[2]=floor((bbMin[2]/h))*h;
    bbMax[0]=ceil((bbMax[0]/h))*h;
    bbMax[1]=ceil((bbMax[1]/h))*h;
    bbMax[2]=ceil((bbMax[2]/h))*h;

    nx =((bbMax[0]-bbMin[0]) / h);
    ny =((bbMax[1]-bbMin[1]) / h);
    nz =((bbMax[2]-bbMin[2]) / h);



    if (flags & VERBOSE) std::cout<<"h = "<<h<<" nx = "<<nx<<" ny = "<<ny<<" nz = "<<nz<<std::endl;
    if (flags & VERBOSE) std::cout<<"Bounding box is "<<bbMin<<" X "<<bbMax<<std::endl;
    phi.allocate(nx,ny,nz);
    biharmonic.allocate(nx,ny,nz);
    laplacian.allocate(nx, ny, nz);
    tempPhi.allocate(nx, ny, nz);
    accepted.allocate(nx, ny, nz);
    phi_max.allocate(nx, ny, nz);
    phi_min.allocate(nx, ny, nz);

    if (flags & VARIABLE_RADIUS) {
        rasterize(particles, radii);
    } else if (flags & (VELOCITY_ANISOTROPY | NEIGHBOR_ANISOTROPY)) {
        rasterize(particles, G, factors);
    } else {
        rasterize(particles);
    }

    redistance(phi_min, tempPhi, accepted, h);
    redistance(phi_max, tempPhi, accepted, h);
    redistance(phi, tempPhi, accepted, h);
}

SmoothingGrid::~SmoothingGrid() {
}

bool SmoothingGrid::computeLaplacian() {
    double divisor = sqr(h), updateBand = 4 * h;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++){
                if (fabs(phi(i, j, k)) <= updateBand) {
                    laplacian(i, j, k) = (phi(i+1, j, k) + phi(i-1, j, k) + phi(i, j+1, k)
                                          + phi(i, j-1, k) + phi(i, j, k+1) + phi(i, j, k-1)
                                          - 6 * phi(i, j, k)) / divisor;
                }
            }
        }
    }
    return true;
}

bool SmoothingGrid::computeLaplacianCG(const SlArray3D<double> &x, SlArray3D<unsigned char> &marked) {
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int k = 1; k < nz - 1; k++){
                if (marked(i,j,k) ) {
                    laplacian(i, j, k) = (x(i+1, j, k) + x(i-1, j, k) + x(i, j+1, k)
                                          + x(i, j-1, k) + x(i, j, k+1) + x(i, j, k-1)
                                          - 6 * x(i, j, k));
                }
            }
        }
    }
    return true;
}

bool SmoothingGrid::computeBiharmonic() {
    double divider = sqr(h), updateBand = 3 * h;
    for (int i = 2; i < nx - 2; i++) {
        for (int j = 2; j < ny - 2; j++) {
            for (int k = 2; k < nz - 2; k++) {
                if (fabs(phi(i, j, k)) <= updateBand){
                    biharmonic(i, j, k) = (laplacian(i+1, j, k) + laplacian(i-1, j, k) + laplacian(i, j+1, k)
                                           + laplacian(i, j-1, k) + laplacian(i, j, k+1) + laplacian(i, j, k-1)
                                           - 6 * laplacian(i, j, k)) / divider;
                }
            }
        }
    }
    return true;
}

double SmoothingGrid::stepBiharmonic(double dt) {
    double change = 0.0, updateBand = 3 * h;

    for (int i = 2; i < nx - 2; i++) {
        for (int j = 2; j < ny - 2; j++) {
            for (int k =  2; k < nz - 2; k++) {
                if (fabs(phi(i, j, k)) <= updateBand) {
                    double phix = cdX(i, j, k, phi, nx, h), phiy = cdY(i, j, k, phi, ny, h), phiz = cdZ(i, j, k, phi, nz, h);
                    double gradMag = sqrt(sqr(phix) + sqr(phiy) + sqr(phiz));
                    double val = biharmonic(i,j,k);
                    double updatedPhi = 0.0;
                    updatedPhi = phi(i, j, k) - val * dt * gradMag;
                    updatedPhi = std::min(updatedPhi, phi_min(i, j, k));
                    updatedPhi = std::max(updatedPhi, phi_max(i, j, k));
                    phi(i, j, k) = updatedPhi;
                    change += fabs(val);
                }
            }
        }
    }
    if (flags & VERBOSE) std::cout<<"Change in this iteration "<<change<<std::endl;
    return change;
}

bool SmoothingGrid::doBiharmonicSmoothing(double dt, int ntimesteps, int redistanceFrequency) {
    redistance(phi, tempPhi, accepted, h);
		for(int i = 0; i < ntimesteps; i++) {
			  computeLaplacian();
				computeBiharmonic();
				stepBiharmonic(dt);
				if (i % redistanceFrequency == 0 && redistanceFrequency != 0 && i != 0) {
					  redistance(phi, tempPhi, accepted, h);
				}
    }
    return true;
}

bool SmoothingGrid::computeMeanCurvature(){
    double d1 = 12 * sqr(h), d2 = 48 * sqr(h), updateBand = 3 * h;
    for (int i = 2; i < nx - 2; i++) {
        for (int j = 2; j < ny - 2; j++) {
            for (int k = 2; k < nz - 2; k++){
                if (fabs(phi(i,j,k)) <= updateBand) {
                    double phix = cdX(i, j, k, phi, nx, h), phiy = cdY(i, j, k, phi, ny, h), phiz = cdZ(i, j, k, phi, nz, h);
                    double gradMag = sqrt(sqr(phix) + sqr(phiy) + sqr(phiz));
                    double phixx, phiyy, phixy, phizz, phixz, phiyz;
                    int ip1 = i + 1, im1 = i - 1, im2 = i - 2, ip2 = i + 2, jp1 = j + 1, jm1 = j - 1, jm2 = j - 2, jp2 = j + 2,
                            kp1 = k + 1, km1 = k - 1, km2 = k - 2, kp2 = k + 2;

                    phixx = (-phi(ip2, j, k) + 16 * phi(ip1, j, k) - 30 * phi(i,j, k) + 16 * phi(im1, j, k) - phi(im2, j, k)) / d1;
                    phiyy = (-phi(i, jp2, k) + 16 * phi(i, jp1, k) - 30 * phi(i,j, k) + 16 * phi(i, jm1, k) - phi(i, jm2, k)) / d1;
                    phizz = (-phi(i, j, kp2) + 16 * phi(i, j, kp1) - 30 * phi(i,j, k) + 16 * phi(i, j, km1) - phi(i, j, km2)) / d1;
                    phixy = (-phi(ip2, jp2, k) + 16 * phi(ip1, jp1, k) + phi(im2, jp2, k) - 16 * phi(im1, jp1, k) + phi(ip2, jm2, k)
                             - 16 * phi(ip1, jm1, k) - phi(im2, jm2, k) + 16 * phi(im1, jm1, k)) / d2;
                    phiyz = (-phi(i, jp2, kp2) + 16 * phi(i, jp1, kp1) + phi(i, jm2, kp2) - 16 * phi(i, jm1, kp1) + phi(i, jp2, km2)
                             - 16 * phi(i, jp1, km1) - phi(i, jm2, km2) + 16 * phi(i, jm1, km1)) / d2;
                    phixz = (-phi(ip2, j, kp2) + 16 * phi(ip1, j, kp1) + phi(im2, j, kp2) - 16 * phi(im1, j, kp1) + phi(ip2, j, km2)
                             - 16 * phi(ip1, j, km1) - phi(im2, j, km2) + 16 * phi(im1, j, km1)) / d2;

                    meanCurvature(i, j, k) = (sqr(phix) * phiyy - 2 * phix * phiy * phixy + sqr(phiy) * phixx + sqr(phix) * phizz
                                              - 2 * phix * phiz * phixz + sqr(phiz) * phixx + sqr(phiy) * phizz
                                              - 2 * phiy * phiz * phiyz + sqr(phiz) * phiyy) / ( sqr(gradMag));
                }
            }
        }
    }
    return true;
}

bool SmoothingGrid::stepMeanCurvature(double dt) {
    double change = 0.0, updateBand = 3 * h;
    for (int i = 2; i < nx - 2; i++) {
        for (int j = 2; j < ny - 2; j++) {
            for (int k =  2; k < nz - 2; k++) {
                if (fabs(phi(i, j, k)) <= updateBand) {
                    double val = meanCurvature(i, j, k);
                    double updatedPhi = phi(i, j, k) + val * dt;
                    updatedPhi = std::min(updatedPhi, phi_min(i, j, k));
                    updatedPhi = std::max(updatedPhi, phi_max(i, j, k));
                    phi(i, j, k) = updatedPhi;
                    change += fabs(val);
                }
            }
        }
    }
    if (flags & VERBOSE) std::cout<<"Change in this iteration "<<change<<std::endl;
    return true;
}

bool SmoothingGrid::stepLaplacian(double dt) {
    double change = 0.0, updateBand = 4 * h;
    for (int i = 2; i < nx - 2; i++) {
        for (int j = 2; j < ny - 2; j++) {
            for (int k =  2; k < nz - 2; k++) {
                if (fabs(phi(i, j, k)) <= updateBand) {
                    double phix = cdX(i, j, k, phi, nx, h), phiy = cdY(i, j, k, phi, ny, h), phiz = cdZ(i, j, k, phi, nz, h);
                    double val = laplacian(i,j,k);
                    double gradMag = sqrt(sqr(phix) + sqr(phiy) + sqr(phiz));
                    double updatedPhi;
                    updatedPhi = phi(i, j, k) + val * dt * gradMag;
                    updatedPhi = std::min(updatedPhi, phi_min(i,j,k));
                    updatedPhi = std::max(updatedPhi, phi_max(i,j,k));
                    phi(i, j, k) = updatedPhi;
                    change += fabs(val);
                }
            }
        }
    }
    if (flags & VERBOSE) std::cout<<"Change in this iteration "<<change<<std::endl;
    return true;
}

bool SmoothingGrid::doMeanCurvatureSmoothing(double dt, int ntimesteps, int redistanceFrequency) {
    for(int i = 0; i < ntimesteps; i++) {
        computeMeanCurvature();
        stepMeanCurvature(dt);
        if (i > 0 && i % redistanceFrequency == 0) redistance(phi, tempPhi, accepted, h);
    }
    return true;
}

bool SmoothingGrid::doLaplacianSmoothing(double dt, int ntimesteps, int redistanceFrequency) {
		for(int i = 0; i < ntimesteps; i++) {
			  computeLaplacian();
			  stepLaplacian(dt);
			  if (i > 0 && i % redistanceFrequency == 0) redistance(phi, tempPhi, accepted, h);
    }
    return true;
}

void sort3val(double &a, double &b, double &c){
    double temp;
    if (a>b){ temp=a; a=b; b=temp; }
    if (a>c){ temp=a; a=c; c=temp; }
    if (b>c){ temp=b; b=c; c=temp; }
}

void sweepPoint(SlArray3D<double> &newPhi, SlArray3D<char> &accepted, int i, int j, int k, double h) {
    int s = accepted(i,j,k) + accepted(i-1,j,k) + accepted(i+1,j,k) +
            accepted(i,j-1,k) + accepted(i,j+1,k) + accepted(i,j,k-1) + accepted(i,j,k+1);
    if (!s) return;

    double a = std::min<double>(fabs(newPhi(i-1,j,k)), fabs(newPhi(i+1,j,k)));
    double b = std::min<double>(fabs(newPhi(i,j-1,k)), fabs(newPhi(i,j+1,k)));
    double c = std::min<double>(fabs(newPhi(i,j,k-1)), fabs(newPhi(i,j,k+1)));
    sort3val(a,b,c);
    double x = a + h;
    if (x > b) {
        x = 0.5*(a+b+sqrt(2*sqr(h)-sqr(a-b)));
        if (x > c) {
            x = (a+b+c+sqrt(3*sqr(h)-sqr(a-b)-sqr(b-c)-sqr(c-a)))/3.0;
        }
    }
    newPhi(i,j,k) = sign(s)*x;
    accepted(i,j,k) = sign(s)*2;
}

void redistance(SlArray3D<double>&phi, SlArray3D<double> &newPhi, SlArray3D<char> &accepted, double h) {
    int nx = phi.nx();
    int ny = phi.ny();
    int nz = phi.nz();
    newPhi = DBL_MAX;
    accepted = 0;
    // determine the inital distance of neighbor grid point
    for (int i=1; i<nx - 1; i++) {
        for (int j=1; j<ny - 1; j++) {
            for (int k=1; k<nz - 1; k++) {
                double x = phi(i,j,k);
                double y = phi(i+1,j,k);
                if (sign(x) != sign(y)) {
                    if (fabs(x)+fabs(y)>0.9*h) {
                        newPhi(i,j,k) = absmin2(newPhi(i,j,k), h*x/fabs(x-y));
                        newPhi(i+1,j,k) = absmin2(newPhi(i+1,j,k), h*y/fabs(x-y));
                        accepted(i,j,k) = sign(newPhi(i,j,k));
                        accepted(i+1,j,k) = sign(newPhi(i+1,j,k));
                    } else {
                        if (accepted(i,j,k)==0) {
                            newPhi(i,j,k) = x;
                            accepted(i,j,k) = sign(newPhi(i,j,k))*3;
                        }
                        if (accepted(i+1,j,k)==0) {
                            newPhi(i+1,j,k) = y;
                            accepted(i+1,j,k) = sign(newPhi(i+1,j,k))*3;
                        }
                    }
                }
                y = phi(i,j+1,k);
                if (sign(x) != sign(y)) {
                    if (fabs(x)+fabs(y)>0.9*h) {
                        newPhi(i,j,k) = absmin2(newPhi(i,j,k), h*x/fabs(x-y));
                        newPhi(i,j+1,k) = absmin2(newPhi(i,j+1,k), h*y/fabs(x-y));
                        accepted(i,j,k) = sign(newPhi(i,j,k));
                        accepted(i,j+1,k) = sign(newPhi(i,j+1,k));
                    } else {
                        if (fabs(accepted(i,j,k))==0) {
                            newPhi(i,j,k) = x;
                            accepted(i,j,k) = sign(newPhi(i,j,k))*3;
                        }
                        if (fabs(accepted(i,j+1,k))==0) {
                            newPhi(i,j+1,k) = y;
                            accepted(i,j+1,k) = sign(newPhi(i,j+1,k))*3;
                        }
                    }
                }
                y = phi(i,j,k+1);
                if (sign(x) != sign(y)) {
                    if (fabs(x)+fabs(y)>0.9*h) {
                        newPhi(i,j,k) = absmin2(newPhi(i,j,k), h*x/fabs(x-y));
                        newPhi(i,j,k+1) = absmin2(newPhi(i,j,k+1), h*y/fabs(x-y));
                        accepted(i,j,k) = sign(newPhi(i,j,k));
                        accepted(i,j,k+1) = sign(newPhi(i,j,k+1));
                    } else {
                        if (fabs(accepted(i,j,k))==0) {
                            newPhi(i,j,k) = x;
                            accepted(i,j,k) = sign(newPhi(i,j,k))*3;
                        }
                        if (fabs(accepted(i,j,k+1))==0) {
                            newPhi(i,j,k+1) = y;
                            accepted(i,j,k+1) = sign(newPhi(i,j,k+1))*3;
                        }
                    }
                }
            }
        }
    }

    // sweeping
    for (int i=1; i<nx - 1; i++) for (int j=1; j<ny - 1; j++) for (int k=1; k<nz - 1; k++)
        if (abs(accepted(i,j,k))%2!=1) sweepPoint(newPhi,accepted,i,j,k,h);
    for (int i=1; i<nx - 1; i++) for (int j=1; j<ny - 1; j++) for (int k=nz-2; k>0; k--)
        if (abs(accepted(i,j,k))%2!=1) sweepPoint(newPhi,accepted,i,j,k,h);
    for (int i=1; i<nx - 1; i++) for (int j=ny-2; j>0; j--) for (int k=1; k<nz - 1; k++)
        if (abs(accepted(i,j,k))%2!=1) sweepPoint(newPhi,accepted,i,j,k,h);
    for (int i=1; i<nx - 1; i++) for (int j=ny-2; j>0; j--) for (int k=nz-2; k>0; k--)
        if (abs(accepted(i,j,k))%2!=1) sweepPoint(newPhi,accepted,i,j,k,h);
    for (int i=nx-2; i>0; i--) for (int j=1; j<ny - 1; j++) for (int k=1; k<nz - 1; k++)
        if (abs(accepted(i,j,k))%2!=1) sweepPoint(newPhi,accepted,i,j,k,h);
    for (int i=nx-2; i>0; i--) for (int j=1; j<ny - 1; j++) for (int k=nz-2; k>0; k--)
        if (abs(accepted(i,j,k))%2!=1) sweepPoint(newPhi,accepted,i,j,k,h);
    for (int i=nx-2; i>0; i--) for (int j=ny-2; j>0; j--) for (int k=1; k<nz - 1; k++)
        if (abs(accepted(i,j,k))%2!=1) sweepPoint(newPhi,accepted,i,j,k,h);
    for (int i=nx-2; i>0; i--) for (int j=ny-2; j>0; j--) for (int k=nz-2; k>0; k--)
        if (abs(accepted(i,j,k))%2!=1) sweepPoint(newPhi,accepted,i,j,k,h);

    phi = newPhi;
}

bool SmoothingGrid::applyBiharmonic(const SlArray3D<double> &x, SlArray3D<double> &y, SlArray3D<unsigned char> &markedB, double dt) const{
	  double factor = dt/sqr(h)*sqr(h);
    for (int i = 2; i < nx - 2; i++) {
        for (int j = 2; j < ny - 2; j++) {
            for (int k = 2; k < nz - 2; k++) {
                if (markedB(i,j,k)) {
                    y(i, j, k) = x(i,j,k)+factor*(laplacian(i+1, j, k) + laplacian(i-1, j, k) + laplacian(i, j+1, k)
																									+ laplacian(i, j-1, k) + laplacian(i, j, k+1) + laplacian(i, j, k-1)
																									- 6 * laplacian(i, j, k));
                }
            }
        }
    }

    return true;
}

bool SmoothingGrid::CG(const SlArray3D<double> &b, SlArray3D<double> &x, SlArray3D<unsigned char> &markedB,
											 SlArray3D<unsigned char> &markedL, double &tol, double dt, double cgThreshold, int cgMaxIters) {
    SlArray3D<double> r, d, q, maxDiff, minDiff;
    d.allocate(nx,ny,nz);
    r.allocate(nx,ny,nz);
    q.allocate(nx,ny,nz);
    maxDiff.allocate(nx,ny,nz);
    minDiff.allocate(nx,ny,nz);
    laplacian = 0.0;
    q = 0.0;
    d = 0.0;
    r = 0.0;

    computeLaplacianCG(x, markedL);
    applyBiharmonic(x, q, markedB, dt);

    for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
            for(int k=0; k<nz; k++) {
                if (markedB(i,j,k)){
                    r(i,j,k) = b(i,j,k) - q(i,j,k);
                    d(i,j,k) = r(i,j,k);
                }
                laplacian(i,j,k) = 0.0;
                q(i,j,k) = 0.0;
            }
        }
    }

    double deltaNew = sqrTwoNorm(r, markedB);
    if (tol < 0.0) tol = sqr(cgThreshold) * deltaNew;
    if (flags & VERBOSE) std::cout<<"delta: "<<deltaNew<<std::endl;
    int iter = 0;

    while(iter++ < cgMaxIters && deltaNew > tol) {
        maxDiff = 0.0;
        minDiff = 0.0;
        double sum = 0.0;
        computeLaplacianCG(d, markedL);
        applyBiharmonic(d, q, markedB, dt);
        double alpha = deltaNew / innerProduct(d, q, markedB);

        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
                    if (markedB(i,j,k)) {
                        x(i,j,k) = x(i,j,k) + alpha*d(i,j,k);
                        r(i,j,k) = r(i,j,k) - alpha*q(i,j,k);
                    }
                }
            }
        }

        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
                    if (markedB(i,j,k)) minDiff(i,j,k) = std::max(x(i,j,k)-phi_min(i,j,k),0.0);
                }
            }
        }

        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
                    if (markedB(i,j,k)) sum += minDiff(i,j,k)*minDiff(i,j,k);
                }
            }
        }
        if (sum >= 1E-10) return true;

        sum = 0.0;
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
                    if (markedB(i,j,k)) maxDiff(i,j,k) = std::min(x(i,j,k)-phi_max(i,j,k),0.0);
                }
            }
        }

        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
                    if (markedB(i,j,k)) sum += maxDiff(i,j,k)  *maxDiff(i,j,k);
                }
            }
        }
        if (sum >= 1E-10) return true;

        double deltaOld=deltaNew;
        if (flags & VERBOSE) std::cout<<"delta: "<<deltaNew<<" iter: "<<iter<<std::endl;

        deltaNew = innerProduct(r,r,markedB);
        double beta = deltaNew / deltaOld;

        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
                    if (markedB(i,j,k)) {
                        d(i,j,k) = r(i,j,k) + beta*d(i,j,k);
                    }
                }
            }
        }
    }
    return true;
}

bool SmoothingGrid::doBiharmonicSmoothing(double dt, int ntimesteps, int redistanceFrequency, double cgThreshold, int cgMaxIters) {
    redistance(phi,tempPhi,accepted,h);
    SlArray3D<double> x, phiZero;
    SlArray3D<unsigned char> activeSetP,activeSetC, markedB, markedL;
    x.allocate(nx,ny,nz);
    phiZero.allocate(nx,ny,nz);
    markedB.allocate(nx,ny,nz);
    markedL.allocate(nx,ny,nz);
    activeSetP.allocate(nx,ny,nz);
    activeSetC.allocate(nx,ny,nz);
    double updatebandB=3*h, updatebandL=4*h;

    for(int titer = 0; titer < ntimesteps; titer++) {
			  activeSetP = 1;
        activeSetC = 0;
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
                    phiZero(i,j,k) = phi(i,j,k);
                    x(i,j,k) = phi(i,j,k);
                    laplacian(i,j,k) = 0.0;
                    markedB(i,j,k) = 0;
                    markedL(i,j,k) = 0;
                }
            }
        }
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
                    double val = fabs(phi(i,j,k));
                    if (val<=updatebandB) markedB(i,j,k)=1;
                    if (val<=updatebandL) markedL(i,j,k)=1;
                }
            }
        }

        double tol = -1;
				int count = 0;
        while(!sameActiveSet(activeSetC, activeSetP)) {
            count++;
            CG(phiZero, x, markedB, markedL, tol, dt, cgThreshold, cgMaxIters);
            activeSetP = activeSetC;
            activeSetC = 0;
            for(int i=0; i<nx; i++) {
                for(int j=0; j<ny; j++) {
                    for(int k=0; k<nz; k++) {
                        if (markedB(i,j,k)) {
                            if (x(i,j,k) > phi_min(i,j,k)) {
                                activeSetC(i,j,k) = 1;
                                x(i,j,k) = phi_min(i,j,k);
                                markedB(i,j,k) = 0;
                            }
                            else if (x(i,j,k) < phi_max(i,j,k)) {
                                activeSetC(i,j,k) = 1;
                                x(i,j,k) = phi_max(i,j,k);
                                markedB(i,j,k) = 0;
                            }
                        }
                    }
                }
            }

        }
        for(int i=0; i<nx; i++) {
            for(int j=0; j<ny; j++) {
                for(int k=0; k<nz; k++) {
                    if (fabs(phiZero(i,j,k)) <= updatebandB) phi(i,j,k) = x(i,j,k);
                }
            }
        }
        if (titer%redistanceFrequency==0) redistance(phi,tempPhi,accepted,h);
		}
    return true;
}




