/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Plugins to setup fluid grids
 *
 *****************************************************************************/

// lsdebug
#include "fluidsolver.h"
#include "solverplugin.h"
#include "paramset.h"


namespace DDF { 

//*****************************************************************************
// simple, set constant value

class sinitSetConstant : public SolverPlugin {
	public:
		sinitSetConstant() : SolverPlugin(),
   				mGrid("-unnamed-"), mVal(0.) { };
		~sinitSetConstant() { };

		virtual bool parseParams(const ParamSet& params) {
			mGrid = params.FindOneString("gridname", mGrid );
			mVal = params.FindOneFloat("value", mVal );
			return true;
		};
		virtual bool initPlugin() { return true; };

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("sinitSetConstant","step "<<dt);
			Grid<int>* flags = mpPlParams->getGridInt("flags");

			// distinguish data types
			if(mpPlParams->haveGridInt(mGrid) && 1) {
				Grid<int>* grid = mpPlParams->getGridInt(mGrid);
				fsSetConstant<int>(flags, grid, (int)mVal);
			} else if(mpPlParams->haveGridReal(mGrid) && 1) {
				Grid<Real>* grid = mpPlParams->getGridReal(mGrid);
				fsSetConstant<Real>(flags, grid, (Real)mVal);
			} else if(mpPlParams->haveGridVec3(mGrid) && 1) {
				Grid<Vec3>* grid = mpPlParams->getGridVec3(mGrid);
				fsSetConstant<Vec3>(flags, grid, Vec3(mVal));				
			} else { 
				errFatal("sinitSetConstant::performStep","Grid not found "<<mGrid , SIMWORLD_GRIDERROR);
			}
			return true;
		};

	protected:
		// grid names to swap
		std::string mGrid;
		// value
		float mVal;
};

template<class T>
class fsSetConditional : public GridOpBase {
	public:
		fsSetConditional(FlagGrid *flags, Grid<T> *dst , T target, int flag) :
				GridOpBase(), mpDst(dst), mTargetVal(target), mCFlag(flag) { 
			mpFlags = flags;
			applyOperatorToGrids(this);
		};
		~fsSetConditional() { }; 
		void resetVariables() { };
		void buildCallList() {
			setFlags(mpFlags);
			gaDst.gridAccInit(mpDst, AM_WRITE, gaCalls); 
		};
		// add forces and update empty cells free surface boundaries
		inline void operator() (int i, int j, int k) { 
			if ((getFlagAcc()(i,j,k) & mCFlag)==mCFlag)
				gaDst.write(i,j,k) = mTargetVal;
		};
		void reduce(fsSetConditional &op) { };

	protected:
		Grid<T> *mpDst;
		GridAccessor<T,0> gaDst;
		T mTargetVal;
		int mCFlag;
}; // fsSetConditional */

class spluginSetConditional : public SolverPlugin {
	public:
		spluginSetConditional() : SolverPlugin(),
   			mGrid("-unnamed-"), mFlag(4),  mTargetVec(0.) {
			debMsg("spluginSetConditional","cons");
		};
		~spluginSetConditional() {
			debMsg("spluginSetConditional","des");
		};

		virtual bool parseParams(const ParamSet& params) {
			debMsg("spluginSetConditional","parse");
			mGrid = params.FindOneString("gridname", mGrid);
			mFlags = params.FindOneString("flaggrid", "flags");
			mTargetVec = params.FindOneVector("target-vec", mTargetVec );
			mTargetReal = params.FindOneFloat("target-real", 0. );
			mFlag = params.FindOneInt("flag", mFlag );
			return true;
		};
		virtual bool initPlugin() {
			debMsg("spluginSetConditional","init");
			return true;
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("spluginSetConditional"," dt="<<dt<<" dest:"<<mGrid );
			FlagGrid* flagg = mpPlParams->getGridInt(mFlags);

			if (mpPlParams->haveGridVec3(mGrid)) {
				fsSetConditional<Vec3>(flagg, mpPlParams->getGridVec3(mGrid), mTargetVec, mFlag);
			} else if (mpPlParams->haveGridReal(mGrid)) {
				fsSetConditional<Real>(flagg, mpPlParams->getGridReal(mGrid), mTargetReal, mFlag);						
			} else { 
				errFatal("spluginSetConditional::performStep","Grid not found "<<mGrid , SIMWORLD_GRIDERROR);
			}
			return true;
		};

	protected:
		std::string mGrid,mFlags;
		int mFlag;
		Vec3 mTargetVec;
		Real mTargetReal;
};

//*****************************************************************************
// initializes the domain with empty cells and a single cell border (or
// more depending on mBorder param) of obstacles

class goInitBoxDomain : public GridOpBase { 
	public:
		goInitBoxDomain(Grid<int>* dst, int border, int setempty, int flagInside = FEMPTY, int flagOutside = FOBSTACLE, int flagFloor = 0) :
				GridOpBase(), mpFlagMod(dst), mBorder(border), mDoSetEmpty(setempty), mFlagInside(flagInside), mFlagOutside(flagOutside), mFlagFloor(flagFloor) {
			mpFlags = NULL; // unused
			mS = mpFlagMod->getSize()-nVec3i(1,1,1);
			applyOperatorToGridsWithoutFlags(this, mpFlagMod);
		};

		~goInitBoxDomain() { }; 
		void resetVariables() { };
		void buildCallList() {
			gaFlagMod.gridAccInit(mpFlagMod, AM_WRITE, gaCalls); 
		};

		// add forces and update empty cells free surface boundaries
		inline void operator() (int i, int j, int k) { 


			gaFlagMod.write(i,j,k);
			if ((i<1+mBorder) || (i>=mS.x-mBorder)
				|| (j<1+mBorder) || (j>=mS.y-mBorder)
#				if DDF_DIMENSION==3
				|| (k<1+mBorder) || (k>=mS.z-mBorder)
#				endif
				) {
				gaFlagMod.write(i,j,k) = mFlagOutside;
			} else {
				if(mDoSetEmpty) gaFlagMod.write(i,j,k) = mFlagInside;
			}
			if (j==0 && mFlagFloor > 0) gaFlagMod.write(i,j,k) = mFlagFloor;
		};
		void reduce(goInitBoxDomain &op) { };

	protected:
		Grid<int> *mpFlagMod;
		GridAccessor<int,0> gaFlagMod;
		nVec3i mS;
		int mBorder, mDoSetEmpty, mFlagInside, mFlagOutside, mFlagFloor;
}; // goInitBoxDomain 

class sinitBoxDomain : public SolverPlugin {
	public:
		sinitBoxDomain()
			: SolverPlugin(), mGrid("flags"), mBorder(0) , mDoSetEmpty(1)
			{ };
		~sinitBoxDomain() { };

		virtual bool parseParams(const ParamSet& params) {
			mGrid = params.FindOneString("gridname", mGrid );
			mBorder = params.FindOneInt("border", mBorder );
			mDoSetEmpty = params.FindOneInt("set-empty", mDoSetEmpty );
			mFlag1 = params.FindOneInt("flag-inside", FEMPTY );
			mFlag2 = params.FindOneInt("flag-border", FOBSTACLE );
			mFloor = params.FindOneInt("flag-floor", 0);
			return true;
		};
		virtual bool initPlugin() { return true; };

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("sinitBoxDomain","step "<<dt);
			Grid<int>* grflags = mpPlParams->getGridInt(mGrid);

			goInitBoxDomain(grflags, mBorder,mDoSetEmpty,mFlag1,mFlag2,mFloor);
			return true;
		};

	protected:
		std::string mGrid;
		int mBorder, mDoSetEmpty, mFlag1, mFlag2, mFloor;
};

//*****************************************************************************
// creates a fluid box spanned by lowCorner and highCorner (both specified
// in range [0,1])

class goSetBoxFluid : public GridOpBase { 
	public:
		goSetBoxFluid(Grid<int> *dst, const Vec3& lowCorner, const Vec3& highCorner, int type, bool brd,
				Grid<Vec3>* vel = NULL, Vec3 set = Vec3(0.), Grid<Real>* dist = NULL, Grid<Vec3>* norm = NULL ) :
				GridOpBase(), mpFlagMod(dst), mType(type), mBorder(brd),
				mpVel(vel), mVelSet(set), mpDist(dist), mpNorm(norm) {
			mpFlags = NULL;
			mILen = 1. / dst->getMaxSize();
			mLowCorner = Vec3(lowCorner.x * float(mpFlagMod->getSizeX()), lowCorner.y * float(mpFlagMod->getSizeY()), lowCorner.z * float(mpFlagMod->getSizeZ()));
			mHighCorner = Vec3(highCorner.x * float(mpFlagMod->getSizeX()), highCorner.y * float(mpFlagMod->getSizeY()), highCorner.z * float(mpFlagMod->getSizeZ()));
			for (int i=0;i<3;i++) 
				if(mLowCorner[i] > mHighCorner[i]) {
					Real a=mLowCorner[i]; mLowCorner[i]=mHighCorner[i]; mHighCorner[i]=a; }
			applyOperatorToGridsWithoutFlags(this, mpFlagMod);
		};

		~goSetBoxFluid() { }; 
		void resetVariables() { };
		void buildCallList() {
			gaFlagMod.gridAccInit(mpFlagMod, AM_WRITE, gaCalls); 
			if(mpVel) gaVel.gridAccInit(mpVel, AM_WRITE, gaCalls); 
			if(mpDist) gaDist.gridAccInit(mpDist, AM_READWRITE, gaCalls); 
			if(mpNorm) gaNorm.gridAccInit(mpNorm, AM_READWRITE, gaCalls); 
		};

		// add forces and update empty cells free surface boundaries
		inline void operator() (int i, int j, int k) { 
			if (mBorder && (i==0 || j==0 || k==0 || i==mpFlagMod->getSizeX()-1 || j==mpFlagMod->getSizeY()-1 || k==mpFlagMod->getSizeZ()-1))
				gaFlagMod.write(i,j,k) = FOBSTACLE;
			if (gaFlagMod.write(i,j,k) == FOBSTACLE)
				return;
			if ((i>=mLowCorner.x) && (i<=mHighCorner.x)
				&& (j>=mLowCorner.y) && (j<=mHighCorner.y)
#if DDF_DIMENSION==3
				&& (k>=mLowCorner.z) && (k<=mHighCorner.z)
#endif
				)
			{
				if(mpVel) {
					Vec3 &v = gaVel.write(i,j,k);
					if(gaFlagMod.write(i,j,k) == mType) {
						//v += mVelSet;
						//v *= 0.5; // warning, depends on init order...
						// NT CJW test, dont overwrite!?
					} else {
						v = mVelSet;
					}
				}

				if(fgIsObstacle(mType)) 
					gaFlagMod.write(i,j,k) = mType; 
				else
					gaFlagMod.write(i,j,k) |= mType; 
			}
			// build distance field
			if (mpDist && mpNorm) {
				int numIn=0;
				Vec3 pos(i,j,k), ref(0.);
				for (int e=0;e<3;e++) { 
					if ((pos[e]>=mLowCorner[e]) && (pos[e]<=mHighCorner[e])) {
						ref[e] = pos[e];
						numIn++;
					} else if (pos[e]>mLowCorner[e])
						ref[e] = mHighCorner[e];
					else
						ref[e] = mLowCorner[e];
				}
				Vec3 nrm = pos-ref;
				Real newDist = normalize(nrm)*mILen, oldDist = gaDist(i,j,k);
				if (oldDist == 0 || newDist < oldDist) {
					gaDist.write(i,j,k) = newDist;
					gaNorm.write(i,j,k) = nrm;
				}
			}
		
		};
		void reduce(goSetBoxFluid &op) { };

	protected:
		Grid<int> *mpFlagMod;
		GridAccessor<int,0> gaFlagMod;
		Vec3 mLowCorner, mHighCorner;
		int mType;
		bool mBorder;

		// optional velocity
		Grid<Vec3> *mpVel;
		GridAccessor<Vec3,0> gaVel;
		Vec3 mVelSet;
		Real mILen;
		Grid<Real> *mpDist;
		Grid<Vec3> *mpNorm;
		GridAccessor<Vec3,0> gaNorm;
		GridAccessor<Real,0> gaDist;		
}; // goSetBoxFluid 

//*****************************************************************************
// creates a fluid ball with given radius at given position (both specified
// in range [0,1])
// optionally, also initializes velocity

class goSetBallFluid : public GridOpBase { 
	public:
		goSetBallFluid(Grid<int> *dst, const Vec3& center, const Vec3& scale, Real radius, int type, int cyl, 
				Grid<Vec3>* vel = NULL, Vec3 set = Vec3(0.), Grid<Real>* dist = NULL, Grid<Vec3>* norm = NULL ) :
					GridOpBase(), mpFlagMod(dst), mpDist(dist), mpNorm(norm), mType(type), mCylinderAxis(cyl), 
					mpVel(vel), mVelSet(set) {
			mpFlags = NULL;
			mScale = Vec3(1.) / scale;
			mCenter = Vec3(
					center.x * float(mpFlagMod->getSizeX()), 
					center.y * float(mpFlagMod->getSizeY()), 
					center.z * float(mpFlagMod->getSizeZ()));
			mRadiusSqr = radius * float(mpFlagMod->getSizeX());
			mRadiusSqr = mRadiusSqr*mRadiusSqr;
			mILen = 1. / (Real)dst->getMaxSize();			
			applyOperatorToGridsWithoutFlags(this, mpFlagMod);
		};

		~goSetBallFluid() { }; 
		void resetVariables() { };
		void buildCallList() {
			gaFlagMod.gridAccInit(mpFlagMod, AM_WRITE, gaCalls); 
			if(mpVel) gaVel.gridAccInit(mpVel, AM_WRITE, gaCalls); 
			if(mpDist) gaDist.gridAccInit(mpDist, AM_READWRITE, gaCalls); 
			if(mpNorm) gaNorm.gridAccInit(mpNorm, AM_READWRITE, gaCalls); 
		};

		inline void operator() (int i, int j, int k) { 
			if (gaFlagMod.write(i,j,k) == FOBSTACLE)
				return;
			//Real d[3];
			Vec3 d;
			d[0] = mCenter.x - i;
			d[1] = mCenter.y - j;
#if DDF_DIMENSION==2
			//const Real d = dx*dx + dy*dy;
			d[2] = 0.;
#else
			d[2] = mCenter.z - k;
#endif
			d *= mScale;
			
			//const Real d = dx*dx + dy*dy + dz*dz;
			if(mCylinderAxis>=0) d[mCylinderAxis] = 0.;
			const Real dist = normNoSqrt(d);

			if (dist < mRadiusSqr)
			{
				if(mpVel) {
					Vec3 &v = gaVel.write(i,j,k);
					if(gaFlagMod.write(i,j,k) == mType) {
						//v += mVelSet;
						//v *= 0.5; // warning, depends on init order...
						// NT CJW test, dont overwrite!?
					} else {
						v = mVelSet;
					}
				}

				// set flag
				//gaFlagMod.write(i,j,k) = mType; 
				if(fgIsObstacle(mType)) 
					gaFlagMod.write(i,j,k) = mType; 
				else
					gaFlagMod.write(i,j,k) |= mType; 

			} else if (mpDist && mpNorm) {
				Real ds = sqrt(dist);
				Real oldDist = gaDist(i,j,k), newDist = (ds - sqrt(mRadiusSqr)) * mILen;
				if (oldDist==0 || newDist<oldDist) {
					gaDist.write(i,j,k) = newDist;
					gaNorm.write(i,j,k) = -d / ds;
				}
			}

		};
		void reduce(goSetBallFluid &op) { };

	protected:
		Grid<int> *mpFlagMod;
		Grid<Real> *mpDist;
		Grid<Vec3> *mpNorm;
		GridAccessor<int,0> gaFlagMod;
		GridAccessor<Real,0> gaDist;
		GridAccessor<Vec3,0> gaNorm;
		Vec3 mCenter;
		Vec3 mScale;
		Real mRadiusSqr,mILen;
		int mType, mCylinderAxis;

		// optional velocity
		Grid<Vec3> *mpVel;
		GridAccessor<Vec3,0> gaVel;
		Vec3 mVelSet;
}; // goSetBallFluid 

class sinitSphere : public SolverPlugin {
	public:
		sinitSphere() : SolverPlugin() { };
		~sinitSphere() { };

		virtual bool parseParams(const ParamSet& params) {
			mGrid = params.FindOneString("gridname", "flags");
			mCenter = params.FindOneVector("center", Vec3(0.5,0.5,0.5));
			mRadius = params.FindOneFloat("radius", 0.2);
			mType = params.FindOneInt("type", FOBSTACLE);
			mDist = params.FindOneString("dist", "");
			mNorm = params.FindOneString("norm", "");
			return true;
		};
		virtual bool initPlugin() { return true; };

		virtual bool performStep(Real dt) {
			debMsg("sinitSphere","step "<<dt);
			Grid<int>* grflags = mpPlParams->getGridInt(mGrid);
			Grid<Real>* dist = mDist.empty() ? NULL : mpPlParams->getGridReal(mDist);
			Grid<Vec3>* norm = mNorm.empty() ? NULL : mpPlParams->getGridVec3(mNorm);

			goSetBallFluid(grflags, mCenter, 1., mRadius, mType, -1, NULL, Vec3(0.), dist, norm);
			return true;
		};
	protected:
		std::string mGrid, mDist, mNorm;
		Vec3 mCenter;
		Real mRadius;
		int mType;
};

class sinitBox : public SolverPlugin {
	public:
		sinitBox() : SolverPlugin() { };
		~sinitBox() { };

		virtual bool parseParams(const ParamSet& params) {
			mGrid = params.FindOneString("gridname", "flags");
			mPos1 = params.FindOneVector("pos1", Vec3(0.));
			mPos2 = params.FindOneVector("pos2", Vec3(0.));
			mType = params.FindOneInt("type", FOBSTACLE);
			mDist = params.FindOneString("dist", "");
			mNorm = params.FindOneString("norm", "");
			return true;
		};
		virtual bool initPlugin() { return true; };

		virtual bool performStep(Real dt) {
			debMsg("sinitBox","step "<<dt);
			Grid<int>* grflags = mpPlParams->getGridInt(mGrid);
			Grid<Real>* dist = mDist.empty() ? NULL : mpPlParams->getGridReal(mDist);
			Grid<Vec3>* norm = mNorm.empty() ? NULL : mpPlParams->getGridVec3(mNorm);

			goSetBoxFluid(grflags, mPos1, mPos2, mType, false, NULL, Vec3(0.), dist, norm);
			return true;
		};
	protected:
		std::string mGrid, mDist, mNorm;
		Vec3 mPos1, mPos2;
		int mType;
};

//*****************************************************************************

SolverPlugin* MakeInitPlugin(std::string name) {

	if(name.compare( string("set-constant") )==0) {
		return new sinitSetConstant;
	} else if(name.compare( string("set-conditional") )==0) {
		return new spluginSetConditional;
	} else if (name.compare( string("init-box-domain") )==0) {
		return new sinitBoxDomain;
	} else if (name.compare( string("init-sphere") )==0) {
		return new sinitSphere;
	} else if (name.compare( string("init-box") )==0) {
		return new sinitBox;
	} 
	return NULL;
}

}; // DDF

