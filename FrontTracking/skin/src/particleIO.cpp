
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

#include <iostream>
#include <fstream>
#include "particleIO.H"

#ifdef OPENVDB
#include <openvdb/openvdb.h>
#include <openvdb/tools/ParticlesToLevelSet.h>

bool readfile(char *infname, MyParticleList &particles, double &radius) {
    std::ifstream in(infname, std::ios::in | std::ios::binary);
    double p0, p1, p2, r;
    while (!in.eof()) {
  			in>>p0>>p1>>p2;
        if (!in.eof()) particles.add(openvdb::Vec3R(p0, p1, p2), radius);
    }
    return true;
}

bool readfile(char *infname, MyParticleList &particles, MyParticleList &particlesmin, MyParticleList &particlesmax,
                                   double &rinit, double &rmin, double &rmax) {
    int numPoints;
    std::ifstream in(infname, std::ios::in | std::ios::binary);
    in.read((char*)&numPoints,sizeof(numPoints));
    for (int i = 0; i < numPoints; i++) {
        float x, y, z, r, lev, v1, v2, v3;
        in.read((char*)&x, sizeof(x));
        in.read((char*)&y, sizeof(y));
        in.read((char*)&z, sizeof(z));
        in.read((char*)&r, sizeof(r));
        in.read((char*)&lev, sizeof(lev));
        in.read((char*)&v1, sizeof(v1));
        in.read((char*)&v2, sizeof(v2));
        in.read((char*)&v3, sizeof(v3));
        particles.add(openvdb::Vec3R(x, y, z), r*rinit);
        particlesmax.add(openvdb::Vec3R(x, y, z), r*rmax);
        particlesmin.add(openvdb::Vec3R(x, y, z), r*rmin);
    }
    return true;
}


bool rasterize(openvdb::DoubleGrid::Ptr grid, double rinit, char *infname) {
    MyParticleList particles;
    readfile(infname, particles, rinit);

    openvdb::tools::ParticlesToLevelSet<openvdb::DoubleGrid> raster(*grid);
    raster.rasterizeSpheres(particles);
    return true;
}

bool rasterizeVariableRadius(openvdb::DoubleGrid::Ptr grid, openvdb::DoubleGrid::Ptr maxgrid, openvdb::DoubleGrid::Ptr mingrid,
                                      double rinit, double rmin, double rmax, char *infname) {
    MyParticleList particles, particlesmax, particlesmin;
    readfile(infname, particles, particlesmin, particlesmax, rinit, rmin, rmax);

		openvdb::tools::ParticlesToLevelSet<openvdb::DoubleGrid> raster(*grid);
    raster.rasterizeSpheres(particles);

		openvdb::tools::ParticlesToLevelSet<openvdb::DoubleGrid> rastermax(*maxgrid);
    rastermax.rasterizeSpheres(particlesmax);

    openvdb::tools::ParticlesToLevelSet<openvdb::DoubleGrid> rastermin(*mingrid);
    rastermin.rasterizeSpheres(particlesmin);
    return true;
}

#else // !OPENVDB

bool readfile(char *infname, std::vector<SlVector3> &particles, std::vector<double> &radii, std::vector<SlVector3> &velocities) {
    std::ifstream in(infname, std::ios::in);
    SlVector3 p;
    double r=1;
    while (!in.eof()) {
  			in>>p[0]>>p[1]>>p[2]>>r;
        if (!in.eof()){
            particles.push_back(p);
            radii.push_back(fabs(r));
        }
    }
    return true;
}

#endif
