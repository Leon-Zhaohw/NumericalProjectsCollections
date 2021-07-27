/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 *  Advection plugins:
 *  semi lagrange, mac cormack
 *
 *****************************************************************************/

#include "fluidsolver.h"
#include "solverplugin.h"
#include "paramset.h"

// safety boundary for semi lagrange advection step
#define SLADVBOUND 2

namespace DDF { 

//*****************************************************************************
// velocity access helpers

class getVelCentered {
	private:
		getVelCentered() {};
		~getVelCentered() {};
	public:
		static inline Vec3 getVelocity(GridAccessor<Vec3,1> &gaVel, int i,int j, int k) {
			Vec3 v = Vec3( 
				0.5*( gaVel(i  ,j  ,k)[0] + gaVel(i+1,j  ,k  )[0] ),
				0.5*( gaVel(i  ,j  ,k)[1] + gaVel(i  ,j+1,k  )[1] ),
				0.5*( gaVel(i  ,j  ,k)[2] + gaVel(i  ,j  ,k+1)[2] )   );
			if(DDF_DIMENSION==2) { v[2] = 0.; } 
			return v;
		}
}; // getVelCentered

class getVelX {
	private:
		getVelX() {};
		~getVelX() {};
	public:
		static inline Vec3 getVelocity(GridAccessor<Vec3,1> &gaVel, int i,int j, int k) {
			const Real u = gaVel(i,j,k)[0];
			const Real v = 0.25*(
					gaVel(i  ,j  ,k)[1] + gaVel(i-1,j  ,k)[1] +
					gaVel(i  ,j+1,k)[1] + gaVel(i-1,j+1,k)[1] );
#			if DDF_DIMENSION==3
			const Real w = 0.25*(
					gaVel(i  ,j,k  )[2] + gaVel(i-1,j,k  )[2] +
					gaVel(i  ,j,k+1)[2] + gaVel(i-1,j,k+1)[2] );
#			else
			const Real w = 0.;
#			endif
			return Vec3(u,v,w);
		}
}; // getVelX

class getVelY {
	private:
		getVelY() {};
		~getVelY() {};
	public:
		static inline Vec3 getVelocity(GridAccessor<Vec3,1> &gaVel, int i,int j, int k) {
			const Real u = 0.25*(
					gaVel(i  ,j  ,k)[0] + gaVel(i  ,j-1,k)[0] +
					gaVel(i+1,j  ,k)[0] + gaVel(i+1,j-1,k)[0] );
			const Real v = gaVel(i,j,k)[1];
#			if DDF_DIMENSION==3
			const Real w = 0.25*(
					gaVel(i,j  ,k  )[2] + gaVel(i,j-1,k  )[2] +
					gaVel(i,j  ,k+1)[2] + gaVel(i,j-1,k+1)[2] );
#			else
			const Real w = 0.;
#			endif
			return Vec3(u,v,w);
		}
}; // getVelY

class getVelZ {
	private:
		getVelZ() {};
		~getVelZ() {};
	public:
		static inline Vec3 getVelocity(GridAccessor<Vec3,1> &gaVel, int i,int j, int k) {
			const Real u = 0.25*(
					gaVel(i  ,j,k  )[0] + gaVel(i  ,j,k-1)[0] +
					gaVel(i+1,j,k  )[0] + gaVel(i+1,j,k-1)[0] );
			const Real v = 0.25*(
					gaVel(i,j  ,k  )[1] + gaVel(i,j  ,k-1)[1] +
					gaVel(i,j+1,k  )[1] + gaVel(i,j+1,k-1)[1] );
			const Real w = gaVel(i,j,k)[2];
			return Vec3(u,v,w);
		}
}; // getVelZ


//*****************************************************************************
// semi lagrange advection
class fsSlAdvection : public GridOpBaseFlagbord<SLADVBOUND> {
	public:
		fsSlAdvection(FlagGrid *flags, Grid<Vec3> *vel, Grid<Vec3> *dstvel, Vec3 &minv, Vec3 &maxv, Real dt) : 
			GridOpBaseFlagbord<SLADVBOUND>(), mpVel(vel), mpDstVel(dstvel), mDt(dt) { 
			mpFlags = flags;
			mMinVel = Vec3( vec2R(mpFlags->getSize()) * 100. );
			mMaxVel = Vec3(0.);
			applyOperatorToGrids(this);
			minv = mMinVel;
			maxv = mMaxVel;
		};
		~fsSlAdvection() { }; 
		void resetVariables() { };
		void buildCallList() { 
			gaVel.gridAccInit(mpVel, AM_READ, gaCalls); 
			gaDstVel.gridAccInit(mpDstVel, AM_WRITE, gaCalls); 
			setFlags(mpFlags);
		};


#		define PINIT_INTERPOL \
				Real srcpi = (Real)i-(u*mDt); \
				Real srcpj = (Real)j-(v*mDt); \
				Real srcpk = (Real)k-(w*mDt); \
				int srci = (int)srcpi; \
			 	int	srcj = (int)srcpj; \
			 	int	srck = (int)srcpk; \
				CLAMP_TO_GRID(srci,srcj,srck,mpVel); \
				const float s1 = srcpi-(float)srci, s0 = 1.-s1; \
				const float t1 = srcpj-(float)srcj, t0 = 1.-t1; \
				const float f1 = srcpk-(float)srck, f0 = 1.-f1; \
		/* end init */

		inline Real getSladvVel(const int a, const int b, const int c, const int index, const int src_i, const int src_j, const int src_k) {
			const int flag = getFlagAcc()(a,b,c);
			if (fgIsFluid(flag) ) return gaVel(a,b,c)[index];
			if (fgIsObstacle(flag) ) return 0.; 
			// default, return own source vel
			return gaVel(src_i,src_j,src_k)[index];
		}
		inline Real interpolateSladvVel(const int i, const int j, const int k, const int index, 
				const Real f0, const Real f1, const Real s0, const Real s1, const Real t0, const Real t1) {
			
			int si = i, sip1 = i+1;
			int sj = j, sjp1 = j+1;
			int sk = k, skp1 = k+1; 

			if(DDF_DIMENSION==3) {
				return ( f0*( 
							s0*(t0*getSladvVel(si  ,sj,sk  , index, si,sj,sk)  + t1*getSladvVel(si  ,sjp1,sk  , index, si,sj,sk) ) + 
							s1*(t0*getSladvVel(sip1,sj,sk  , index, si,sj,sk)  + t1*getSladvVel(sip1,sjp1,sk  , index, si,sj,sk) ) ) 
						+ f1 * ( 
							s0*(t0*getSladvVel(si  ,sj,skp1, index, si,sj,sk)  + t1*getSladvVel(si  ,sjp1,skp1, index, si,sj,sk) ) + 
							s1*(t0*getSladvVel(sip1,sj,skp1, index, si,sj,sk)  + t1*getSladvVel(sip1,sjp1,skp1, index, si,sj,sk) ) ) 
						) ;
			} else {
				// NOTE! srck not used for 2d!
				return ( f0*( 
							s0*(t0*getSladvVel(si  ,sj,sk, index, si,sj,sk)  + t1*getSladvVel(si  ,sjp1,sk, index, si,sj,sk) ) + 
							s1*(t0*getSladvVel(sip1,sj,sk, index, si,sj,sk)  + t1*getSladvVel(sip1,sjp1,sk, index, si,sj,sk) ) ) );
			} // 2d
		}

		inline void operator() (int i, int j, int k) { 
			// always init...?
			if (fgIsObstacle(getFlagAcc()(i,j,k)) ) {
				gaDstVel.write(i,j,k)=Vec3(0.);
				return;
			}
			if (!fgIsFluid(getFlagAcc()(i,j,k)) ) {
				gaDstVel.write(i,j,k) = gaVel(i,j,k); // copy, not necessary?
				return;
			}

			const bool debugAdv = false;

			if (1) { // advect x
				const Real u = gaVel(i,j,k)[0];
				const Real v = 0.25*(
						gaVel(i  ,j  ,k)[1] + gaVel(i-1,j  ,k)[1] +
						gaVel(i  ,j+1,k)[1] + gaVel(i-1,j+1,k)[1] );
#				if DDF_DIMENSION==3
				const Real w = 0.25*(
						gaVel(i  ,j,k  )[2] + gaVel(i-1,j,k  )[2] +
						gaVel(i  ,j,k+1)[2] + gaVel(i-1,j,k+1)[2] );
#				else
				const Real w = 0.;
#				endif
				PINIT_INTERPOL;
				gaDstVel.write(i,j,k)[0] = interpolateSladvVel(srci,srcj,srck,0,f0,f1,s0,s1,t0,t1);
				if (debugAdv) debMsg("sgAdvectX","vel="<<PRINT_VEC(u,v,w)<<" mDt="<<mDt<<" weights="<<s0<<","<<s1<<","<<t0<<","<<t1<<","<<f0<<","<<f1 <<" src"<<PRINT_VEC(srci,srcj,srck)<<" srcp"<<PRINT_VEC(srcpi,srcpj,srcpk)<<" ijk"<<PRINT_VEC(i,j,k) );
				if (debugAdv) debMsg("sgAdvectX","vals: "<<gaVel(srci  ,srcj,srck)[0] <<", "<<gaVel(srci  ,srcj+1,srck)[0] <<", "<<gaVel(srci+1,srcj,srck)[0]  <<", "<<gaVel(srci+1,srcj+1,srck)[0]  );
				if (debugAdv && DDF_DIMENSION==3) debMsg("sgAdvectX","vals: "<<gaVel(srci  ,srcj,srck+1)[0] <<", "<<gaVel(srci  ,srcj+1,srck+1)[0] <<", "<<gaVel(srci+1,srcj,srck+1)[0]  <<", "<<gaVel(srci+1,srcj+1,srck+1)[0]  );
			}

			if (1) { // advect y
				const Real u = 0.25*(
						gaVel(i  ,j  ,k)[0] + gaVel(i  ,j-1,k)[0] +
						gaVel(i+1,j  ,k)[0] + gaVel(i+1,j-1,k)[0] );
				const Real v = gaVel(i,j,k)[1];
#				if DDF_DIMENSION==3
				const Real w = 0.25*(
						gaVel(i,j  ,k  )[2] + gaVel(i,j-1,k  )[2] +
						gaVel(i,j  ,k+1)[2] + gaVel(i,j-1,k+1)[2] );
#				else
				const Real w = 0.;
#				endif
				PINIT_INTERPOL;
				gaDstVel.write(i,j,k)[1] = interpolateSladvVel(srci,srcj,srck,1,f0,f1,s0,s1,t0,t1);
				if (debugAdv) debMsg("sgAdvectY","vel="<<PRINT_VEC(u,v,w)<<" mDt="<<mDt<<" weights="<<s0<<","<<s1<<","<<t0<<","<<t1<<","<<f0<<","<<f1 <<" src"<<PRINT_VEC(srci,srcj,srck) );
				if (debugAdv) debMsg("sgAdvectY","vals: "<<gaVel(srci  ,srcj,srck)[1] <<", "<<gaVel(srci  ,srcj+1,srck)[1] <<", "<<gaVel(srci+1,srcj,srck)[1]  <<", "<<gaVel(srci+1,srcj+1,srck)[1]  );
			}

			if (gDim==3){ // advect z
				const Real u = 0.25*(
						gaVel(i  ,j,k  )[0] + gaVel(i  ,j,k-1)[0] +
						gaVel(i+1,j,k  )[0] + gaVel(i+1,j,k-1)[0] );
				const Real v = 0.25*(
						gaVel(i,j  ,k  )[1] + gaVel(i,j  ,k-1)[1] +
						gaVel(i,j+1,k  )[1] + gaVel(i,j+1,k-1)[1] );
				const Real w = gaVel(i,j,k)[2];
				PINIT_INTERPOL;
				gaDstVel.write(i,j,k)[2] = interpolateSladvVel(srci,srcj,srck,2,f0,f1,s0,s1,t0,t1);
				if (debugAdv) debMsg("sgAdvectZ","vel="<<PRINT_VEC(u,v,w)<<" mDt="<<mDt<<" weights="<<s0<<","<<s1<<","<<t0<<","<<t1<<","<<f0<<","<<f1 <<" src"<<PRINT_VEC(srci,srcj,srck) );
				if (debugAdv) debMsg("sgAdvectZ","vals: "<<gaVel(srci  ,srcj,srck)[2] <<", "<<gaVel(srci  ,srcj+1,srck)[2] <<", "<<gaVel(srci+1,srcj,srck)[2]  <<", "<<gaVel(srci+1,srcj+1,srck)[2]  );
			}
			if (DDF_DEBUG) {
				Vec3 v = gaDstVel(i,j,k);
				if (normNoSqrt(v)>normNoSqrt(mMaxVel)) mMaxVel=v;
				if (normNoSqrt(v)<normNoSqrt(mMinVel)) mMinVel=v;
			}
		};
		void reduce(fsSlAdvection &op) { 
			// update min/max
			if (normNoSqrt(op.mMaxVel)>normNoSqrt(mMaxVel)) mMaxVel=op.mMaxVel;
			if (normNoSqrt(op.mMinVel)<normNoSqrt(mMinVel)) mMinVel=op.mMinVel;
		};
	protected:
		Grid<Vec3> *mpVel;
		Grid<Vec3> *mpDstVel;
		Real mDt;
		GridAccessor<Vec3,SLADVBOUND> gaVel; 
		GridAccessor<Vec3,0> gaDstVel;
		Vec3 mMinVel, mMaxVel;
}; // fsSlAdvection */

#undef PINIT_INTERPOL 

template < class VEL >
class fsSlAdvectionReal : public GridOpBaseFlagbord<SLADVBOUND> {
	public:
		fsSlAdvectionReal(FlagGrid *flags, Grid<Vec3> *vel, Grid<Real> *srcreal, 
				Grid<Real> *dstreal, Real dt, bool advectLs=false) : 
			GridOpBaseFlagbord<SLADVBOUND>(), mpVel(vel), mpSrcReal(srcreal), 
			mpDstReal(dstreal), mDt(dt), mAdvectLs(advectLs) 
		{ 
			mpFlags = flags;
			//debMsg("fsSlAdvectionReal","Grid "<<srcreal->getName()<<", dt="<<mDt ); // DEBUGold
			applyOperatorToGrids(this);
		};
		~fsSlAdvectionReal() { }; 
		void resetVariables() { };
		void buildCallList() { 
			gaVel.gridAccInit(mpVel, AM_READ, gaCalls); 
			gaSrcReal.gridAccInit(mpSrcReal, AM_READ, gaCalls); 
			gaDstReal.gridAccInit(mpDstReal, AM_WRITE, gaCalls); 
			setFlags(mpFlags);
		}; 

		// apply changes for vec sladv later on...  
		inline Real getSladvReal(
				const int a, const int b, const int c, 
				const int index, const int src_i, const int src_j, const int src_k) {
			// default, return value
			return gaSrcReal(a,b,c);
		}
		inline Real interpolateSladvReal(const int i, const int j, const int k, const int index, 
				const Real f0, const Real f1, const Real s0, const Real s1, const Real t0, const Real t1) {

			int si = i, sip1 = si+1;
			int sj = j, sjp1 = sj+1;
			int sk = k, skp1 = sk+1; 

			if(DDF_DIMENSION==3) {
				return ( f0*( 
							s0*(t0*getSladvReal(si  ,sj,sk  , index, si,sj,sk)  + t1*getSladvReal(si  ,sjp1,sk  , index, si,sj,sk) ) + 
							s1*(t0*getSladvReal(sip1,sj,sk  , index, si,sj,sk)  + t1*getSladvReal(sip1,sjp1,sk  , index, si,sj,sk) ) ) 
					 	 + f1 * ( 
							s0*(t0*getSladvReal(si  ,sj,skp1, index, si,j,k)  + t1*getSladvReal(si  ,sjp1,skp1, index, si,sj,sk) ) + 
							s1*(t0*getSladvReal(sip1,sj,skp1, index, si,j,k)  + t1*getSladvReal(sip1,sjp1,skp1, index, si,sj,sk) ) ) 
						) ;
			} else {
				// NOTE! srck not used for 2d!
				return ( f0*( 
							s0*(t0*getSladvReal(si  ,sj,sk, index, si,j,k)  + t1*getSladvReal(si  ,sjp1,sk, index, si,sj,sk) ) + 
							s1*(t0*getSladvReal(sip1,sj,sk, index, si,j,k)  + t1*getSladvReal(sip1,sjp1,sk, index, si,sj,sk) ) ) );
			} // 2d
		}

		inline void operator() (int i, int j, int k) { 
			// always init...?
			if (fgIsObstacle(getFlagAcc()(i,j,k)) ) {
				gaDstReal.write(i,j,k)=(0.);
				return;
			}
			if ((!mAdvectLs) && (!fgIsFluid(getFlagAcc()(i,j,k))) ) {
				gaDstReal.write(i,j,k) = gaSrcReal(i,j,k); // copy, not necessary?
				return;
			} 
			const bool debugAdv = false;
			const Vec3 vel = VEL::getVelocity(gaVel, i,j,k);

			Real srcpi = (Real)i-(vel[0]*mDt); 
			Real srcpj = (Real)j-(vel[1]*mDt); 
			Real srcpk = (Real)k-(vel[2]*mDt); 
			int srci = (int)srcpi; 
			int	srcj = (int)srcpj; 
			int	srck = (int)srcpk; 
			CLAMP_TO_GRID(srci,srcj,srck,mpVel); 
			const float s1 = srcpi-(float)srci, s0 = 1.-s1; 
			const float t1 = srcpj-(float)srcj, t0 = 1.-t1; 
			const float f1 = srcpk-(float)srck, f0 = 1.-f1; 

			gaDstReal.write(i,j,k) = interpolateSladvReal(srci,srcj,srck,0,f0,f1,s0,s1,t0,t1);
		};
		void reduce(fsSlAdvectionReal &op) { 
			// ...
		};
	protected:
		Grid<Vec3> *mpVel;
		Grid<Real> *mpSrcReal;
		Grid<Real> *mpDstReal;
		Real mDt;
		GridAccessor<Vec3,1> gaVel; 
		GridAccessor<Real,SLADVBOUND> gaSrcReal; 
		GridAccessor<Real,0> gaDstReal;
		// levelset special, also advect non-fluid cells
		bool mAdvectLs;

}; // fsSlAdvectionReal */

// interface for advection
class spSemiLagrangeAdvVec3 : public SolverPlugin {
	public:
		spSemiLagrangeAdvVec3() : SolverPlugin(), mNameFlags("flags"),mNameCurrVel("vel-curr"),mNameOldVel("vec-temp"),mNameVel("vel-curr"),
   				 mTemp1("temp1"),mTemp2("temp2"),mMAC(true) {
		};
		~spSemiLagrangeAdvVec3() { };

		// init plugin, return failure
		virtual bool parseParams(const ParamSet& params) {
			debMsg("spSemiLagrangeAdvVec3","parse");
			mNameFlags = params.FindOneString("flags",mNameFlags);
			mNameVel = params.FindOneString("vel",mNameVel);
			mTemp1 = params.FindOneString("temp1",mTemp1);
			mTemp2 = params.FindOneString("temp2",mTemp1);
			mNameCurrVel = params.FindOneString("vel-src",mNameCurrVel);
			mNameOldVel = params.FindOneString("vel-dst",mNameOldVel);			
			mMAC = params.FindOneInt("mac",1) != 0;
			return true;
		};
		virtual bool initPlugin() {
			debMsg("spSemiLagrangeAdvVec3","init");
			return true;
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("spSemiLagrangeAdvVec3","step "<<dt);
			Vec3 mMinv,mMaxv; 
			FlagGrid*   flags   = mpPlParams->getGridInt(mNameFlags);
			Grid<Vec3>* Currvel = mpPlParams->getGridVec3(mNameCurrVel);
			Grid<Vec3>* Oldvel  = mpPlParams->getGridVec3(mNameOldVel);
			Grid<Real>* tmp1 = mpPlParams->getGridReal(mTemp1);
			Grid<Real>* tmp2 = mpPlParams->getGridReal(mTemp2);
			Grid<Vec3>* vel      = mpPlParams->getGridVec3(mNameVel);

			if (mMAC)
				fsSlAdvection(flags, Currvel, Oldvel, mMinv,mMaxv, dt);
			else
			{
				// X
				goCopyVec3ToScalar<Real,0>(tmp1,Currvel);
				fsSlAdvectionReal<getVelCentered>(flags, vel, tmp1, tmp2, dt, false);
				goCopyScalarToVec3<Real,0>(Oldvel,tmp2);
				// Y
				goCopyVec3ToScalar<Real,1>(tmp1,Currvel);
				fsSlAdvectionReal<getVelCentered>(flags, vel, tmp1, tmp2, dt, false);
				goCopyScalarToVec3<Real,1>(Oldvel,tmp2);
				// Z
				goCopyVec3ToScalar<Real,2>(tmp1,Currvel);
				fsSlAdvectionReal<getVelCentered>(flags, vel, tmp1, tmp2, dt, false);
				goCopyScalarToVec3<Real,2>(Oldvel,tmp2);
			}
			swapGrids(mpPlParams, mNameCurrVel, mNameOldVel);
			return true;
		};

	protected:
		// grid access
		std::string mNameFlags, mNameCurrVel, mNameOldVel, mNameVel, mTemp1, mTemp2;
		bool mMAC;
};

// interface for real advection
class spSemiLagrangeAdvReal : public SolverPlugin {
	public:
		spSemiLagrangeAdvReal() : 
					SolverPlugin(), mNameFlags("flags"),
   				mNameVel("vel-curr"), mNameCurrReal("-unnamed1-"), 
					mNameOldReal("real-temp") ,
					mAdvectLs(false) { }; 
		~spSemiLagrangeAdvReal() { };

		// init plugin, return failure
		virtual bool parseParams(const ParamSet& params) {
			debMsg("spSemiLagrangeAdvReal","parse");
			mNameFlags    = params.FindOneString("flags",mNameFlags);
			mNameVel      = params.FindOneString("vel",mNameVel);
			mNameCurrReal = params.FindOneString("real-src",mNameCurrReal);
			mNameOldReal  = params.FindOneString("real-dst",mNameOldReal);
			mAdvectLs     = params.FindOneInt("advect-ls",mAdvectLs);
			return true;
		};
		virtual bool initPlugin() {
			debMsg("spSemiLagrangeAdvReal","init");
			return true;
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("spSemiLagrangeAdvReal","step "<<dt);
			FlagGrid*   flags    = mpPlParams->getGridInt(mNameFlags);
			Grid<Vec3>* vel      = mpPlParams->getGridVec3(mNameVel);
			Grid<Real>* CurrReal = mpPlParams->getGridReal(mNameCurrReal);
			Grid<Real>* newReal  = mpPlParams->getGridReal(mNameOldReal);

			fsSlAdvectionReal<getVelCentered>(flags, vel, CurrReal, newReal, dt, mAdvectLs);
			swapGrids(mpPlParams, mNameCurrReal, mNameOldReal);
			return true;
		};

	protected:
		// grid access
		std::string mNameFlags, mNameVel, mNameCurrReal, mNameOldReal;
		bool mAdvectLs;
};

// mac cormack adevection - correct using the back and forth sl steps
template < class VEL >
class fsMacCormackCorrect : public GridOpBaseFlagbord<1> {
	public:
		fsMacCormackCorrect(FlagGrid *flags, Grid<Real> *dst, Grid<Real> *old, Grid<Real> *temp1, 
				Grid<Real> *temp2, Real strength, bool advLs) : 
			GridOpBaseFlagbord<1>(), mpDst(dst), mpOld(old), mpTemp1(temp1), mpTemp2(temp2) 
			, mStrength(strength), mAdvectLs(advLs)
		{ 
			mpFlags = flags;
			applyOperatorToGrids(this);
		};

		~fsMacCormackCorrect() { }; 

		void resetVariables() { };

		void buildCallList() {
			gaDst.gridAccInit(mpDst, AM_WRITE, gaCalls); 
			gaOld.gridAccInit(mpOld, AM_READ, gaCalls); 
			gaTemp1.gridAccInit(mpTemp1, AM_READ, gaCalls); 
			gaTemp2.gridAccInit(mpTemp2, AM_READ, gaCalls); 
			setFlags(mpFlags);
		};

		inline void operator() (int i, int j, int k) { 
			if (!fgIsFluid(getFlagAcc()(i,j,k)) ) {
				if(mAdvectLs)  gaDst.write(i,j,k) = gaTemp1(i,j,k); // levelset -> copy value
				else           gaDst.write(i,j,k) = 0.; // default, reset
				return;
			}
			gaDst.write(i,j,k) = gaTemp1(i,j,k) + 
				(gaOld(i,j,k)
				 - gaTemp2(i,j,k)) * 0.5;

			// interpolate between SL and MC
			if(mStrength<1.) {
				gaDst.write(i,j,k) = (1.-mStrength)*gaTemp1(i,j,k) + mStrength*gaDst(i,j,k);
			}
		};

		void reduce(fsMacCormackCorrect &op) { };

	protected:
		Grid<Real> *mpDst, *mpOld, *mpTemp1, *mpTemp2;
		GridAccessor<Real,0> gaDst;
		GridAccessor<Real,1> gaOld, gaTemp1, gaTemp2;
		Real mStrength;
		bool mAdvectLs;
}; 

// mac cormack adevection - clamp sl advection 
template < class VEL >
class fsMacCormackClamp : public GridOpBaseFlagbord<SLADVBOUND> {
	public:
		fsMacCormackClamp(FlagGrid *flags, Grid<Vec3> *vel, Grid<Real> *srcreal, 
				Grid<Real> *dstreal, Real dt) : 
			GridOpBaseFlagbord<SLADVBOUND>(), mpVel(vel), mpSrcReal(srcreal), 
			mpDstReal(dstreal), mDt(dt) { 

			mpFlags = flags;
			applyOperatorToGrids(this);
		};
		~fsMacCormackClamp() { }; 
		void resetVariables() { };
		void buildCallList() { 
			gaVel.gridAccInit(mpVel, AM_READ, gaCalls); 
			gaSrcReal.gridAccInit(mpSrcReal, AM_READ, gaCalls); 
			gaDstReal.gridAccInit(mpDstReal, AM_WRITE, gaCalls); 
			setFlags(mpFlags);
		}; 

		//! hack to fix negative rounding, this depends on always rounding down...	define P2IOFF 200

		// apply changes for vec sladv later on...  
		inline Real getSladvReal(
				const int a, const int b, const int c, 
				const int src_i, const int src_j, const int src_k) {
			// default, return value
			return gaSrcReal(a,b,c);
		}
		inline Real clampValue(const int i, const int j, const int k, 
				Real val) { 
			int si = i, sip1 = i+1;
			int sj = j, sjp1 = j+1;
			int sk = k, skp1 = k+1; 

			const Real v000 = getSladvReal(si  ,sj  ,sk , si,sj,sk);
			const Real v010 = getSladvReal(si  ,sjp1,sk , si,sj,sk);
			const Real v100 = getSladvReal(sip1,sj  ,sk , si,sj,sk);
			const Real v110 = getSladvReal(sip1,sjp1,sk , si,sj,sk);

			// compute minima/maxima
			Real max = v000, min=v000;
			if(v010>max) max = v010;
			if(v100>max) max = v100;
			if(v110>max) max = v110;

			if(v010<min) min = v010;
			if(v100<min) min = v100;
			if(v110<min) min = v110;

			if(DDF_DIMENSION==3) {
				const Real v001 = getSladvReal(si  ,sj  ,skp1, si,j,k);
				const Real v011 = getSladvReal(si  ,sjp1,skp1, si,sj,sk);
				const Real v101 = getSladvReal(sip1,sj  ,skp1, si,j,k);
				const Real v111 = getSladvReal(sip1,sjp1,skp1, si,sj,sk);

				if(v001>max) max = v001;
				if(v011>max) max = v011;
				if(v101>max) max = v101;
				if(v111>max) max = v111;
				if(v001<min) min = v001;
				if(v011<min) min = v011;
				if(v101<min) min = v101;
				if(v111<min) min = v111;
			}
			if(val<min) return min;
			if(val>max) return max;
			return val;
		}

		// same as semi-lagrange advection, only clamp instead of interpolation
		inline void operator() (int i, int j, int k) { 
			// always init...?
			if (fgIsObstacle(getFlagAcc()(i,j,k)) ) {
				return;
			}
			if (!fgIsFluid(getFlagAcc()(i,j,k)) ) {
				return;
			} 
			const Vec3 vel = VEL::getVelocity(gaVel, i,j,k);

			Real srcpi = (Real)i-(vel[0]*mDt); 
			Real srcpj = (Real)j-(vel[1]*mDt); 
			Real srcpk = (Real)k-(vel[2]*mDt); 

			int srci = (int)srcpi; 
			int	srcj = (int)srcpj; 
			int	srck = (int)srcpk; 
			CLAMP_TO_GRID(srci,srcj,srck,mpVel); 
		
			gaDstReal.write(i,j,k) = clampValue(srci,srcj,srck,  gaDstReal(i,j,k) );
		};
		void reduce(fsMacCormackClamp &op) { 
			// ...
		};
	protected:
		Grid<Vec3> *mpVel;
		Grid<Real> *mpSrcReal;
		Grid<Real> *mpDstReal;
		Real mDt;
		GridAccessor<Vec3,1> gaVel; 
		GridAccessor<Real,SLADVBOUND> gaSrcReal; 
		GridAccessor<Real,0> gaDstReal;
}; // fsMacCormackClamp */


// mac cormack adevection - clamp boundaries
template < class VEL >
class fsMacCormackBoundaryClamp : public GridOpBaseFlagbord<SLADVBOUND> {
	public:
		fsMacCormackBoundaryClamp(FlagGrid *flags, Grid<Vec3> *vel, 
				Grid<Real> *dstreal, Grid<Real> *srcorg, Real dt) : 
			GridOpBaseFlagbord<SLADVBOUND>(), 
			mpVel(vel), 
			mpSrcOrg(srcorg), 
			mpDstReal(dstreal), 
			mDt(dt) 
		{	  	
			mpFlags = flags;
			applyOperatorToGrids(this);
		};

		~fsMacCormackBoundaryClamp() { }; 

		void resetVariables() { };

		void buildCallList() { 
			gaVel.gridAccInit(mpVel, AM_READ, gaCalls); 
			gaDstReal.gridAccInit(mpDstReal, AM_WRITE, gaCalls); 
			gaSrcOrg.gridAccInit(mpSrcOrg, AM_READ, gaCalls); 
			setFlags(mpFlags);
		}; 

		// same as semi-lagrange advection, only clamp instead of interpolation
		inline void operator() (int i, int j, int k) { 
			// always init...?
			if (fgIsObstacle(getFlagAcc()(i,j,k)) ) {
				return;
			}
			if (!fgIsFluid(getFlagAcc()(i,j,k)) ) {
				// NT fix, preserve surrounding values from original grid everywhere in empty space
				gaDstReal.write(i,j,k) = gaSrcOrg(i,j,k); // fix , SURROUND_RESET
				return;
			} 
			const Vec3 vel = VEL::getVelocity(gaVel, i,j,k);
			bool doReset = false;

			// forward lookup
			Real srcpiForw = (Real)i-(vel[0]*mDt); 
			Real srcpjForw = (Real)j-(vel[1]*mDt); 
			Real srcpkForw = (Real)k-(vel[2]*mDt); 
			// backward lookup
			Real srcpiBack = (Real)i+(vel[0]*mDt); 
			Real srcpjBack = (Real)j+(vel[1]*mDt); 
			Real srcpkBack = (Real)k+(vel[2]*mDt); 

			int srciForw = (int)srcpiForw; 
			int srcjForw = (int)srcpjForw; 
			int srckForw = (int)srcpkForw; 
			if( CLAMP_TO_GRID_BOOL(srciForw,srcjForw,srckForw,mpVel) ) doReset = true;

			int srciBack = (int)srcpiBack; 
			int srcjBack = (int)srcpjBack; 
			int srckBack = (int)srcpkBack; 
			if(!doReset) {
				if( CLAMP_TO_GRID_BOOL(srciBack,srcjBack,srckBack,mpVel) ) doReset = true;
			} // not yet reset

			if(!doReset) {
				if ( fgIsObstacle(getFlagAcc()(srciForw,srcjForw,srckForw)) ||
					 fgIsObstacle(getFlagAcc()(srciBack,srcjBack,srckBack)) 
						) {
					doReset = true;
				}
			} // not yet reset

			if(doReset) {
				gaDstReal.write(i,j,k) = gaSrcOrg(i,j,k);
			} 
		};
		void reduce(fsMacCormackBoundaryClamp &op) { 
			// ...
		};
	protected:
		Grid<Vec3> *mpVel;
		Grid<Real> *mpSrcOrg;
		Grid<Real> *mpDstReal;
		Real mDt;
		GridAccessor<Vec3,1> gaVel; 
		GridAccessor<Real,0> gaDstReal;
		GridAccessor<Real,0> gaSrcOrg;
}; // fsMacCormackBoundaryClamp */


// interface for MacCormack advection of a scalar field
class spMacCormackReal : public SolverPlugin {
	public:
		spMacCormackReal() : 
					SolverPlugin(), mNameFlags("flags"),
   				mNameVel("vel-curr"), 
					mNameCurrReal("-unnamed1-"), mNameOldReal("real-temp") ,
					mNameTemp1("temp1"), mNameTemp2("temp2") ,
					mStrength(0.5), mAdvectLs(false)
				{ }; 
		~spMacCormackReal() { };

		// init plugin, return failure
		virtual bool parseParams(const ParamSet& params) {
			debMsg("spMacCormackReal","parse");
			mNameFlags = params.FindOneString("flags",mNameFlags);
			mNameVel   = params.FindOneString("vel",mNameVel);
			mNameCurrReal = params.FindOneString("real-src",mNameCurrReal);
			mNameOldReal  = params.FindOneString("real-dst",mNameOldReal);
			mNameTemp1    = params.FindOneString("temp1",mNameTemp1);
			mNameTemp2    = params.FindOneString("temp2",mNameTemp2);
			mAdvectLs     = params.FindOneInt("advect-ls",mAdvectLs);
			mStrength     = params.FindOneFloat("strength",mStrength);
			return true;
		};
		virtual bool initPlugin() {
			debMsg("spMacCormackReal","init");
			return true;
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("spMacCormackReal","step "<<dt);

			FlagGrid*   flags    = mpPlParams->getGridInt(mNameFlags);
			Grid<Vec3>* vel      = mpPlParams->getGridVec3(mNameVel);
			Grid<Real>* CurrReal = mpPlParams->getGridReal(mNameCurrReal);
			Grid<Real>* newReal  = mpPlParams->getGridReal(mNameOldReal);
			Grid<Real>* temp1    = mpPlParams->getGridReal(mNameTemp1);
			Grid<Real>* temp2    = mpPlParams->getGridReal(mNameTemp2);

			// normal forward step
			fsSlAdvectionReal<getVelCentered>(flags, vel, CurrReal, temp1, dt, mAdvectLs);

			// backwards step
			fsSlAdvectionReal<getVelCentered>(flags, vel, temp1, temp2, -1. * dt, mAdvectLs);

			// compute correction
			fsMacCormackCorrect<getVelCentered>(flags, newReal, CurrReal, temp1, temp2, mStrength, mAdvectLs);

			// clamp values
			fsMacCormackClamp<getVelCentered>(flags, vel, CurrReal, newReal, dt);
			//fsMacCormackBoundaryClamp(flags, vel, CurrReal, temp1, dt);
			fsMacCormackBoundaryClamp<getVelCentered>(flags, vel, newReal, temp1, dt);

			swapGrids(mpPlParams, mNameCurrReal, mNameOldReal);
			return true;
		};

	protected:
		// grid access
		std::string mNameFlags, mNameVel, mNameCurrReal, mNameOldReal, mNameTemp1, mNameTemp2;
		Real mStrength;
		bool mAdvectLs;
};

// interface for maccormack advection of the velocities
class spMacCormackVec3 : public SolverPlugin {
	public:
		spMacCormackVec3() : 
					SolverPlugin(), mNameFlags("flags"),
   				mNameVel("vel-curr"), 
					mNameCurrVec3("vel-curr"), mNameOldVec3("vec-temp") ,
					mNameTemp1("temp1"), mNameTemp2("temp2"), mMAC(true)
			{ }; 
		~spMacCormackVec3() { };

		// init plugin, return failure
		virtual bool parseParams(const ParamSet& params) {
			debMsg("spMacCormackVec3","parse");
			mNameFlags = params.FindOneString("flags",mNameFlags);
			mNameCurrVec3 = params.FindOneString("vel-src",mNameCurrVec3);
			mNameOldVec3  = params.FindOneString("vel-dst",mNameOldVec3);
			mMAC = params.FindOneInt("mac",1) != 0;
			
			// reuse currvec3, if given
			mNameVel      = params.FindOneString("vel",mNameCurrVec3);

			mNameTemp1    = params.FindOneString("temp1",mNameTemp1);
			mNameTemp2    = params.FindOneString("temp2",mNameTemp2);
			return true;
		};
		virtual bool initPlugin() {
			debMsg("spMacCormackVec3","init");
			return true;
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("spMacCormackVec3","step "<<dt);
			const Real mStrength = 0.5;

			FlagGrid*   flags    = mpPlParams->getGridInt(mNameFlags);
			Grid<Vec3>* advVel   = mpPlParams->getGridVec3(mNameVel);
			Grid<Vec3>* currVec3  = mpPlParams->getGridVec3(mNameCurrVec3);
			Grid<Vec3>* oldVec3   = mpPlParams->getGridVec3(mNameOldVec3);
			Grid<Real>* temp1    = mpPlParams->getGridReal(mNameTemp1);
			Grid<Real>* temp2    = mpPlParams->getGridReal(mNameTemp2);

			Grid<Real>* realTmp = mpPlParams->getGridReal("residual");
			Grid<Real>* realNew  = mpPlParams->getGridReal("tmp");

			// normal forward step
			bool advectLs = false;
			//fsSlAdvectionReal<getVelCentered>(flags, currVec3, realTmp, temp1, dt, advectLs);

			if (mMAC)
			{
				// X
				goCopyVec3ToScalar<Real,0>(realTmp,currVec3); 
				// normal forward step
				fsSlAdvectionReal<getVelX>(flags, advVel, realTmp, temp1, dt, advectLs);
				// backwards step
				fsSlAdvectionReal<getVelX>(flags, advVel, temp1, temp2, -1. * dt, advectLs); 
				// compute correction
				fsMacCormackCorrect<getVelX>(flags, realNew, realTmp, temp1, temp2, mStrength, advectLs); 
				// clamp values
				fsMacCormackClamp<getVelX>(flags, advVel, realTmp, realNew, dt);
				fsMacCormackBoundaryClamp<getVelX>(flags, advVel, realNew, temp1, dt);
				goCopyScalarToVec3<Real,0>(oldVec3,realNew); 

				// Y
				goCopyVec3ToScalar<Real,1>(realTmp,currVec3); 
				// normal forward step
				fsSlAdvectionReal<getVelY>(flags, advVel, realTmp, temp1, dt, advectLs);
				// backwards step
				fsSlAdvectionReal<getVelY>(flags, advVel, temp1, temp2, -1. * dt, advectLs); 
				// compute correction
				fsMacCormackCorrect<getVelY>(flags, realNew, realTmp, temp1, temp2, mStrength, advectLs); 
				// clamp values
				fsMacCormackClamp<getVelY>(flags, advVel, realTmp, realNew, dt);
				fsMacCormackBoundaryClamp<getVelY>(flags, advVel, realNew, temp1, dt);		
				goCopyScalarToVec3<Real,1>(oldVec3,realNew); 
		
				// Z
				goCopyVec3ToScalar<Real,2>(realTmp,currVec3); 
				// normal forward step
				fsSlAdvectionReal<getVelZ>(flags, advVel, realTmp, temp1, dt, advectLs);
				// backwards step
				fsSlAdvectionReal<getVelZ>(flags, advVel, temp1, temp2, -1. * dt, advectLs); 
				// compute correction
				fsMacCormackCorrect<getVelZ>(flags, realNew, realTmp, temp1, temp2, mStrength, advectLs); 
				// clamp values
				fsMacCormackClamp<getVelZ>(flags, advVel, realTmp, realNew, dt);
				fsMacCormackBoundaryClamp<getVelZ>(flags, advVel, realNew, temp1, dt);
				goCopyScalarToVec3<Real,2>(oldVec3,realNew); 		
			} 
			else 
			{
				for (int i=0;i<3;i++)
				{
					if (i==0) goCopyVec3ToScalar<Real,0>(realTmp,currVec3);
					if (i==1) goCopyVec3ToScalar<Real,1>(realTmp,currVec3); 
					if (i==2) goCopyVec3ToScalar<Real,2>(realTmp,currVec3); 
					// normal forward step
					fsSlAdvectionReal<getVelCentered>(flags, advVel, realTmp, temp1, dt, advectLs);
					// backwards step
					fsSlAdvectionReal<getVelCentered>(flags, advVel, temp1, temp2, -1. * dt, advectLs); 
					// compute correction
					fsMacCormackCorrect<getVelCentered>(flags, realNew, realTmp, temp1, temp2, mStrength, advectLs); 
					// clamp values
					fsMacCormackClamp<getVelCentered>(flags, advVel, realTmp, realNew, dt);
					fsMacCormackBoundaryClamp<getVelCentered>(flags, advVel, realNew, temp1, dt);
					if (i==0) goCopyScalarToVec3<Real,0>(oldVec3,realNew); 		
					if (i==1) goCopyScalarToVec3<Real,1>(oldVec3,realNew); 		
					if (i==2) goCopyScalarToVec3<Real,2>(oldVec3,realNew); 		
				}
			}
			swapGrids(mpPlParams, mNameCurrVec3, mNameOldVec3);
			return true;
		};

	protected:
		// grid access
		std::string mNameFlags, mNameVel, mNameCurrVec3, mNameOldVec3, mNameTemp1, mNameTemp2;
		bool mMAC;
		//bool mAdvectLs;
};

//*****************************************************************************
// external accessor helper function

void doSemiLagrangeReal(FlagGrid* flags, Grid<Vec3>* vel, Grid<Real>* src, Grid<Real>* dst, Real dt)
{
	fsSlAdvectionReal<getVelCentered>(flags, vel, src, dst, dt, false);
}

//*****************************************************************************

// create one of the standard hard coded plugins
SolverPlugin* MakeAdvectionPlugin(std::string name) {

	if(name.compare("semi-lagr-advect-vec3")==0) {
		return new spSemiLagrangeAdvVec3;
	} else if(name.compare( "semi-lagr-advect-real")==0) {
		return new spSemiLagrangeAdvReal;

	} else if(name.compare( string("maccormack-advect-real") )==0) {
		return new spMacCormackReal;
	} else if(name.compare( string("maccormack-advect-vec3") )==0) {
		return new spMacCormackVec3;
	}

	return NULL;
}


} // end namespace DDF 

