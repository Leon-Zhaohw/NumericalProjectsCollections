/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Vortex Particle classes
 *
 *****************************************************************************/

#ifndef DDFVORTEXPART_H

#include <list>
#include "globals.h"
#include "fluidsolver.h"
#include "randomstream.h"

#define DEBUG_VORTEXPATH 1

namespace DDF {

#ifdef DEBUG_VORTEXPATH
// for debugging: draws a line behind the particles
struct LineElement {
	Vec3 p0, p1;
	Real r,g,b,a,s;
};
extern std::list<LineElement> gVortexPath;
#endif


class VortexParticle{
	public:
		VortexParticle(const Vec3& p, const Vec3& s, Real irad, Real orad);
		virtual ~VortexParticle() {};

		// check if particle should split
		static inline bool doDecay() { return msDecayTime != 0; }
		inline Vec3& getPos() {return mPos; };
		inline Vec3& getStrength() {return mStrength; };
		inline const Vec3& getPos() const {return mPos; };
		// radius in which particles merge		
		inline Real getMergeR2() { return mIRad * mIRad; }
		inline Real getRadius() { return mIRad; }
		inline Vec3& getForce() { return mForce; };
		inline Real& getTime() { return mTime; };
		int& getFlags() { return mFlags; };
		bool isDeleted() { return (mFlags & FDelete) == 1; };
		bool canMerge() { return (mFlags & (FDelete | FInertial)) == 0; }
		// mark for deletation in the next timestep
		void markDelete() { mFlags |= 1; };
		// use this to advect the particle
		inline void advance(const Vec3& v);
		void setDecayTimer();
		
		virtual void scale(Real fac) = 0;
		// merge this particle with particle v
		virtual bool merge(VortexParticle* v) = 0;
		// split into two particles, radius controlled by msRadiusCascade
		virtual VortexParticle* split() = 0;
		// add particles vorticity to desired vorticity grid (w_D)
		virtual void buildReference(Grid<Vec3>* pTemp, FlagGrid* flags, Real mult, bool doForces) = 0;
		// apply vortex stretching
		virtual void vortexStretch(Grid<Vec3>* pVel, Real dt, Real mult);
		// apply particles velocity kernel to simulation
		virtual void applyForce(Grid<Vec3>* pVel,Grid<Vec3>* pVort,Grid<Vec3>* pTemp, FlagGrid* flags, 
				Real mult, Real fadein, bool doForces, Grid<Vec3>* velHelper = NULL) = 0;
		virtual std::string toString();

		static const Real msRepulsion, msAttenuation; // factor for repulsion, attenuation forces when particles is near a wall
		static const bool msStretchingPreserve; // preserve vortex strength in vortex stretching (only rotate axis) ?
		static Real msDecayTime; // Factor to the decay timescale
		static Real msInitialTime; // Factor to the transition time to inertial subrange
		static Real msRadiusCascade; // Factor to wavenumber in decay (i.e. 0.5 -> radius is halfed). Controls granularity of decay.
		static Real msDissipateRadius; // Particle below this radius are dissipated (usually 0.5..1 grid cells)
		static Vec3 msU0; // Global velocity scale (inflow etc.)

		enum Flag { FDelete = 1, FInertial = 2};

 	protected:
		Vec3 mPos, mStrength, mForce;
		int mFlags;
		Real mIRad, mORad, mTime, mAttenuate;
};

class VortexParticleGaussian : public VortexParticle{
	public:
		VortexParticleGaussian(const Vec3& pos, const Vec3& strength, Real irad);
		
		void scale(Real fac);
		bool merge(VortexParticle* v);
		VortexParticle* split();
		void buildReference(Grid<Vec3>* pTemp, FlagGrid* flags, Real mult, bool doForces);
		void applyForce(Grid<Vec3>* pVel,Grid<Vec3>* pVort,Grid<Vec3>* pTemp, FlagGrid* flags, 
				Real mult, Real fadein, bool doForces, Grid<Vec3>* velHelper = NULL);
		static Real getAreaMultiplier() { return msAreaMult; }

		static const Real msLogCutoff, msAreaMult;
	protected:
		Real mSigma;
};

class VorticitySystem
{
	public:
		typedef std::list<VortexParticle*> VList;

		VorticitySystem(FluidSolver *fs); 
		~VorticitySystem();

		inline VList& getParticles() { return mParts; }
		static inline Real rnd() { return randStream.getReal(); }
		static inline Real rndg() { return randStream.getRandNorm(0.,1.); }

		// apply all particles to velocity grid
		void applyForces(Grid<Vec3>* pVel, Grid<Vec3>* pVort,  Grid<Vec3>* pTemp, FlagGrid* flags, Real mult, bool doForces, Grid<Vec3>* velHelper = NULL);
		void advectParticles(Grid<Vec3>* pVel, FlagGrid* flags, Real dt, Real multiplier);
		void registerPrecomp(vector<Real>* dat) { mpPrecomp = dat; }
		// merge & split particles
		void merge (FlagGrid* flags, Grid<Real>* ndist,Grid<Vec3>* vel, Real dt, Real mergeDist);
		vector<Real>* getPrecomp() { return mpPrecomp; }

		Real getFadeIn() const { return msFadeIn; }
		void setFadeIn(Real set) { msFadeIn = set; }

	protected:
		VList mParts;
		DDF::FluidSolver * mpFsolver;
		vector<Real> * mpPrecomp;
		static RandomStream randStream;
		Real msFadeIn;
};


}

#define DDFVORTEXPART_H
#endif
