
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

#ifndef PARTICLEIO_is_defined
#define PARTICLEIO_is_defined

#ifdef OPENVDB
#include <openvdb/openvdb.h>

class MyParticleList {

protected:
    struct MyParticle {
        openvdb::Vec3R p, v;
        double  r;
    };
    double mRadiusScale, mVelocityScale;
    std::vector<MyParticle> mParticleList;

public:
    MyParticleList(double rScale=1, double vScale=1)
        : mRadiusScale(rScale), mVelocityScale(vScale) {}

    void add(const openvdb::Vec3R &p, const double &r,
             const openvdb::Vec3R &v=openvdb::Vec3R(0,0,0)) {
        MyParticle pa;
        pa.p = p;
        pa.r = r;
        pa.v = v;
        mParticleList.push_back(pa);
    }
    /// @return coordinate bbox in the space of the specified transfrom

    openvdb::CoordBBox getBBox(const openvdb::GridBase& grid) {
        openvdb::CoordBBox bbox;
        openvdb::Coord &min= bbox.min(), &max = bbox.max();
        openvdb::Vec3R pos;
        double rad, invDx = 1/grid.voxelSize()[0];
        for (size_t n = 0, e = this->size(); n < e; ++n) {
            this->getPosRad(n, pos, rad);
            const openvdb::Vec3d xyz = grid.worldToIndex(pos);
            const double   r  = rad * invDx;
            for (int i = 0; i < 3; ++i) {
                min[i] = openvdb::math::Min(min[i], openvdb::math::Floor(xyz[i] - r));
                max[i] = openvdb::math::Max(max[i], openvdb::math::Ceil( xyz[i] + r));
            }
        }
        return bbox;
    }

    openvdb::Vec3R pos(int n)   const {return mParticleList[n].p;}

    openvdb::Vec3R vel(int n)   const {return mVelocityScale*mParticleList[n].v;}

    double radius(int n) const {return mRadiusScale*mParticleList[n].r;}

    //////////////////////////////////////////////////////////////////////////////
    /// The methods below are the only ones required by tools::ParticleToLevelSet
    /// @note We return by value since the radius and velocities are modified
    /// by the scaling factors! Also these methods are all assumed to
    /// be thread-safe.

    /// Return the total number of particles in list.
    ///  Always required!
    size_t size() const { return mParticleList.size(); }

    /// Get the world space position of n'th particle.
    /// Required by ParticledToLevelSet::rasterizeSphere(*this,radius).
    void getPos(size_t n,  openvdb::Vec3R&pos) const { pos = mParticleList[n].p; }

    void getPosRad(size_t n,  openvdb::Vec3R& pos, double& rad) const {
        pos = mParticleList[n].p;
        rad = mRadiusScale*mParticleList[n].r;
    }

    void getPosRadVel(size_t n,  openvdb::Vec3R& pos, double& rad, openvdb::Vec3R& vel) const {
        pos = mParticleList[n].p;
        rad = mRadiusScale*mParticleList[n].r;
        vel = mVelocityScale*mParticleList[n].v;
    }
    // The method below is only required for attribute transfer
    void getAtt(size_t n, openvdb::Index32& att) const { att = n; }
};

bool rasterize(openvdb::DoubleGrid::Ptr grid, double rinit, char *infname);
bool rasterizeVariableRadius(openvdb::DoubleGrid::Ptr grid, openvdb::DoubleGrid::Ptr maxgrid, openvdb::DoubleGrid::Ptr mingrid,
								 double rinit, double rmin, double rmax, char *infname);


#else // !OPENVDB

#include "slVector.H"
#include <vector>

bool readfile(char *infname, std::vector<SlVector3> &particles, std::vector<double> &radii, std::vector<SlVector3> &velocities);

#endif // OPENVDB

#endif // PARTICLEIO_IS_DEFINED
