/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Plugins to access vortex particles
 *
 *****************************************************************************/

#include "fluidsolver.h"
#include "solverplugin.h"
#include "paramset.h"

#include "randomstream.h"
#include "vortexpart.h"
#include <zlib.h>

using namespace std;

namespace DDF { 

// declare external accesors
extern void closeDatabaseWriter();
extern std::vector<Vec3> initDatabaseWriter(const string&, int, int);

#if DDF_DIMENSION==3
#define FOR_KERNEL(rad) for (int z=-rad;z<=rad;z++) for (int y=-rad;y<=rad;y++) for (int x=-rad;x<=rad;x++)
#else	
#define FOR_KERNEL(rad) for (int y=-rad,z=0;y<=rad;y++) for (int x=-rad;x<=rad;x++)
#endif

//*****************************************************************************
// Invoke vortex particle advection, vortex stretching
class spluginAdvectVPart : public SolverPlugin {
	public:
		spluginAdvectVPart() : SolverPlugin(), mGrid("vel-curr")
		{ };
		~spluginAdvectVPart() { };

		virtual bool parseParams(const ParamSet& params) {
			mGrid = params.FindOneString("gridname",mGrid);
			mSolver = params.FindOneString("from-solver","");	// use particles hostsed on a different solver
			return true;
		};
		virtual bool initPlugin() {
			debMsg("spluginAdvectVPart","init");
			return true;
		};

		// perform step with given dt
		virtual bool performStep(Real dt) {
			debMsg("spluginAdvectVPart","step "<<dt);
			FlagGrid* flags = mpPlParams->getGridInt("flags");
			Grid<Vec3>* vel = mpPlParams->getGridVec3(mGrid);
	
			FluidSolver* solver = mSolver.empty() ? mpPlParams->getFluidSolver() : ddfWorldFindSolver(mSolver);
			Real multiplier = mpPlParams->getMultiplier();
			if (multiplier <=0) multiplier = 1.;
			solver->getVorticitySys()->advectParticles(vel, flags, dt/multiplier, multiplier);
			return true;
		};

	protected:
	string mGrid, mSolver;
};


//*****************************************************************************
// Gridop : compute vorticity of velocity grid
// output vorticity is centered on ijk
class fsComputeVort : public GridOpBaseFlagbord<1> {
	public:
		fsComputeVort(FlagGrid *flags, Grid<Vec3> *vel, Grid<Vec3> *vort, Real dx) : 
				GridOpBaseFlagbord<1>(), mpVel(vel), mpVort(vort) { 
			mDxInv = 1.0;
			mpFlags = flags;
			applyOperatorToGrids(this);
		};
		~fsComputeVort() { }; 
		void resetVariables() { };
		void buildCallList() {
			gaVel.gridAccInit(mpVel, AM_READ, gaCalls); 
			gaVort.gridAccInit(mpVort, AM_WRITE, gaCalls); 
			setFlags(mpFlags);
		};

		inline void operator() (int i, int j, int k)
		{ 
			int flag = getFlagAcc()(i,j,k);
		   if(!fgIsFluid(flag)) {
				gaVort.write(i,j,k) = Vec3(0.);
				return;
			}
			#define getv(i,j,k) (!fgIsFluid(getFlagAcc()(i,j,k))?v0:gaVel(i,j,k))
			Vec3 w(0.),v0 (gaVel(i,j,k));
			for (int l=0;l<4;l++) // averaging kernel to get center value
			{
				if (gDim==3)
				{
					w[0] += (Real)_b[l] * (getv(i,j+_b[l],k+_a[l])[2] - getv(i,j+_a[l],k+_b[l])[1]);
					w[1] += (Real)_b[l] * (getv(i+_a[l],j,k+_b[l])[0] - getv(i+_b[l],j,k+_a[l])[2]);
				}
				w[2] += (Real)_b[l] * (getv(i+_b[l],j+_a[l],k)[1] - getv(i+_a[l],j+_b[l],k)[0]);
			}
			#undef getv

			gaVort.write(i,j,k) = w*(mDxInv*0.25);
		};
		void reduce(fsComputeVort &op) { };

	protected:
		Grid<Vec3> *mpVel;
		Grid<Vec3> *mpVort;
		GridAccessor<Vec3,1> gaVel;
		GridAccessor<Vec3,0> gaVort;
		
		static const int _a[4],_b[4];
		Real mDxInv;
}; 
const int fsComputeVort::_a[4]={ 1, 0, 1, 0};
const int fsComputeVort::_b[4]={ 1, 1,-1,-1};


//*****************************************************************************
// Plugin to invoke vorticity computation of velocity grid
class spluginComputeVorticity : public SolverPlugin {
	public:
		spluginComputeVorticity() : SolverPlugin(),
   				mGrid("vel-curr"),mVort("-unnamed-") {
		};
		~spluginComputeVorticity() { };

		virtual bool parseParams(const ParamSet& params) {
			mGrid = params.FindOneString("velocity", mGrid );
			mVort = params.FindOneString("vorticity", mVort );
			return true;
		};
		virtual bool initPlugin() {
			debMsg("spluginComputeVorticity","init");
			return true;
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("spluginComputeVorticity"," dt="<<dt<<" vel:"<<mGrid );
			Grid<Vec3>* vel = mpPlParams->getGridVec3(mGrid); 
			Grid<Vec3>* vort = mpPlParams->getGridVec3(mVort);
			fsComputeVort(mpPlParams->getFluidSolver()->getGridFlags(),vel,vort,mpPlParams->getDeltaX());
			return true;
		};

	protected:
		std::string mGrid, mVort;
};

//*****************************************************************************
// Apply velocity kernel of the vortex particles to velocity grid
class spluginApplyVParts : public SolverPlugin {
	public:
		spluginApplyVParts() : SolverPlugin(),
   				mGrid("vel-curr"), mVort("unnamed"),mTemp("vec-temp") {
		};
		~spluginApplyVParts() { };

		virtual bool parseParams(const ParamSet& params) {
			mGrid = params.FindOneString("velocity", mGrid );
			mVort = params.FindOneString("vorticity", mVort );	// current measured vorticity
			mTemp = params.FindOneString("temp", mTemp );		// temporary grid for storing reference vorticity field
			mSolver = params.FindOneString("from-solver", "");	// use vortex particles hosted on a different solver
			mDoForces = params.FindOneInt("doforces", 1) != 0;	// activate repulsion forces
			return true; 
		};
		virtual bool initPlugin() {return true;};

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("spluginApplyVParts"," dt="<<dt<<" vel:"<<mGrid );
			FlagGrid* flags = mpPlParams->getGridInt("flags");
			Grid<Vec3>* vel = mpPlParams->getGridVec3(mGrid);
			Grid<Vec3>* vort = mpPlParams->getGridVec3(mVort);
			Grid<Vec3>* temp = mpPlParams->getGridVec3(mTemp);
			Grid<Vec3>* velhelp = NULL;
			if (mpPlParams->haveGridVec3("vel-helper")) mpPlParams->getGridVec3("vel-helper");		// helper grid shows only induced velocity
			FluidSolver* solver = mSolver.empty() ? mpPlParams->getFluidSolver() : ddfWorldFindSolver(mSolver);
			Real multiplier = mpPlParams->getMultiplier();
			if (multiplier <0) multiplier = 1.;	// multiplier due to grid rescaling
			temp->clearGrid();		
	
			solver->getVorticitySys()->applyForces(vel,vort,temp, flags,multiplier,mDoForces, velhelp);
			
			return true;
		};

	protected:
		std::string mGrid, mVort,mTemp,mSolver;
		bool mDoForces;
};

//*****************************************************************************
// Merge, split, dissipate particles
class spluginMergeVParts : public SolverPlugin {
	public:
		spluginMergeVParts() : SolverPlugin() {};
		~spluginMergeVParts() { };

		virtual bool parseParams(const ParamSet& params) {
			mDist = params.FindOneString("ndist", "unnamed" );  // Distance field to the obstacles
			mInitTime = params.FindOneFloat("init-time", 50.);	// Factor to the lifespan in model-dependant range
			mDecayTime = params.FindOneFloat("decay-time", 250. );	// Factor to the time between particle splitting
			mMergeDist = params.FindOneFloat("merge-dist", 1.); ;	// Factor to the merge radius
			mDissipateRadius = params.FindOneFloat("dissipate-radius", 1. ); // Radius below which particles are deleted
			mRadiusCascade = params.FindOneFloat("radius-cascade", 2.); // On split, particle decay into 2 particles of with radius old_radius/radius_cascade
																		// a low value leads to frequent splitting
			return true;
		};
		virtual bool initPlugin() {
			debMsg("spluginMergeVParts","init");
			return true;
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("spluginMergeVParts"," dt="<<dt );
			FlagGrid* flags = mpPlParams->getGridInt("flags");
			Grid<Vec3>* vel = mpPlParams->getGridVec3("vel-curr");
			VortexParticle::msDecayTime = mDecayTime;
			VortexParticle::msInitialTime = mInitTime;
			VortexParticle::msRadiusCascade = mRadiusCascade;
			VortexParticle::msDissipateRadius = mDissipateRadius;
			mpPlParams->getFluidSolver()->getVorticitySys()->merge(flags, mpPlParams->getGridReal(mDist), vel, dt,mMergeDist);
			
			return true;
		};

	protected:
	std::string mDist;
	Real mInitTime,mDecayTime,mMergeDist,mRadiusCascade,mDissipateRadius;
};

//*****************************************************************************
// initialize a grid from a grid of another solver (possibly with different resolution)
class spluginInterpolateGridFrom : public SolverPlugin {
	public:
		spluginInterpolateGridFrom() : SolverPlugin(),
   				mName("-unnamed1-"), mGridSrc("-unnamed2-"), mGridDst("-unnamed3-"),
  					mDoFlagsTreatment(0), mCorrectScale(1), mCorrectAdd(0)	{
		};
		~spluginInterpolateGridFrom() { };

		virtual bool parseParams(const ParamSet& params) {
			mName = params.FindOneString("name", mName );      // solver name
			mGridSrc = params.FindOneString("src", mGridSrc ); // source grid
			mGridDst = params.FindOneString("dst", mGridDst ); // dest. grid
			mDoFlagsTreatment = params.FindOneInt("do-flags", mDoFlagsTreatment );			// adjust flags ?
			mCorrectScale = params.FindOneFloat("correctScale", mCorrectScale);
			mCorrectAdd = params.FindOneFloat("correctAdd", mCorrectAdd);
			return true;
		};
		virtual bool initPlugin() { return true; };

		virtual bool performStep(Real org_dt);
		bool performStepFlags(FluidSolver* src, FluidSolver* dst);
	protected:
		// name of solver to init from
		std::string mName;
		// name of src/dst grid
		std::string mGridSrc,mGridDst;

		int mDoFlagsTreatment;

		void initSizeScale(GridBase* src, GridBase* dst);
		nVec3i m_srcSize;
		nVec3i m_dstSize;
		Vec3 m_sizeScale;

		float mCorrectScale;
		float mCorrectAdd;
};

// Helper function to adjust scale size
void spluginInterpolateGridFrom::initSizeScale(GridBase* src, GridBase* dst) {
	m_srcSize = src->getSize();
	m_dstSize = dst->getSize();
	m_sizeScale = Vec3(
			(Real)m_srcSize[0] / (Real)m_dstSize[0],
			(Real)m_srcSize[1] / (Real)m_dstSize[1],
			(Real)m_srcSize[2] / (Real)m_dstSize[2] );
	debMsg("spluginInterpolateGridFrom ","scale size = "<<m_sizeScale<<", src="<<m_srcSize<<" to dst="<<m_dstSize);
}

// NYI, TODO, loop over range, use priority of flags: obstacle>fluid>empty
bool spluginInterpolateGridFrom::performStepFlags(FluidSolver* srcSolver, FluidSolver* dstSolver)
{
	// similar to real grid treatment
	if( ( srcSolver->getParams()->haveGridInt(mGridSrc)) &&
	    ( mpPlParams->haveGridInt(mGridDst) ) ) {
		FlagGrid* src = srcSolver->getParams()->getGridInt(mGridSrc);
		FlagGrid* dst = mpPlParams->getGridInt(mGridDst);

		initSizeScale(src,dst);

		FOR_IJK_GRID_BND(dst,1) {
			// linearly interpolate values
			Vec3 p = Vec3(i,j,k) * m_sizeScale; 

			// for safety reasons, don't overwrite obstacle flags
			if( fgIsObstacle(dst->getGlobal(i,j,k)) )
				continue;
			// downsizing, no interpol
			dst->getGlobal(i,j,k) = src->getGlobal( (int)p[0],  (int)p[1],  (int)p[2]); 
		} 

	} else {
		errFatal("spluginInterpolateGridFrom","Do-flags Grids not found! "<<mGridSrc<<" "<<mGridDst, SIMWORLD_PLUGINERROR);
		return false;
	}
	return true;
}

// interpolate grid values from another grid of another solver
bool spluginInterpolateGridFrom::performStep(Real org_dt) {

	//Real dt = org_dt * mpPlParams->getDeltaX(); 
	debMsg("spluginInterpolateGridFrom","name:"<<mName<<"; from="<<mName<<" "<<mGridSrc<<"->"<<mGridDst );

	FluidSolver* srcSolver = ddfWorldFindSolver(mName);
	FluidSolver* dstSolver = mpPlParams->getFluidSolver();
	if(!srcSolver) {
		errMsg("spluginInterpolateGridFrom","Source solver '"<<mName<<"' not found!");
		return false;
	}

	bool doneSth = false;

	if(mDoFlagsTreatment) {
		return performStepFlags(srcSolver, dstSolver);
	}

	if( ( srcSolver->getParams()->haveGridVec3(mGridSrc)) &&
	    ( mpPlParams->haveGridVec3(mGridDst) ) ) {

		Grid<Vec3>* src = srcSolver->getParams()->getGridVec3(mGridSrc);
		Grid<Vec3>* dst = mpPlParams->getGridVec3(mGridDst);

		//const nVec3i grs = mpPlParams->getFluidSolver()->get GridInitDim();
		const nVec3i grs = src->getSize();

		int ks = 1;
		int ke = grs[2]-1;
		if(gDim==2) {
			ks = mpPlParams->getFluidSolver()->get2dKstart();
			ke = ks+1;
		}

		nVec3i srcSize = src->getSize();
		nVec3i dstSize = dst->getSize();
		Vec3 sizeScale = Vec3(
				(Real)srcSize[0] / (Real)dstSize[0],
				(Real)srcSize[1] / (Real)dstSize[1],
				(Real)srcSize[2] / (Real)dstSize[2] );
		debMsg("scale ","size = "<<sizeScale);

		const Real dtFac = 1./1.; 

		FOR_IJK_GRID(dst) { 
			// linearly interpolate values
			Real pi = (Real)(i)*sizeScale[0]; 
			Real pj = (Real)(j)*sizeScale[1]; 
			Real pk = (Real)(k)*sizeScale[2]; 
			dst->getGlobal(i,j,k) = src->getInterpolated(pi, pj, pk) * dtFac * mCorrectScale; // */
		} 
		doneSth = true;
	} // VEC3

	if( ( srcSolver->getParams()->haveGridReal(mGridSrc)) &&
	    ( mpPlParams->haveGridReal(mGridDst) ) ) {
		Grid<Real>* src = srcSolver->getParams()->getGridReal(mGridSrc);
		Grid<Real>* dst = mpPlParams->getGridReal(mGridDst);

		nVec3i srcSize = src->getSize();
		nVec3i dstSize = dst->getSize();
		Vec3 sizeScale = Vec3(
				(Real)srcSize[0] / (Real)dstSize[0],
				(Real)srcSize[1] / (Real)dstSize[1],
				(Real)srcSize[2] / (Real)dstSize[2] );
		debMsg("scale ","size = "<<sizeScale);

		FOR_IJK_GRID(dst) {
			// linearly interpolate values
			Real pi = (Real)(i)*sizeScale[0]; 
			Real pj = (Real)(j)*sizeScale[1]; 
			Real pk = (Real)(k)*sizeScale[2]; 
			dst->getGlobal(i,j,k) = src->getInterpolated(pi, pj, pk) * mCorrectScale + mCorrectAdd; 
		} 

		doneSth = true;
	} // REAL

	if(!doneSth) {
		errFatal("spluginInterpolateGridFrom","Grids not found! "<<mGridSrc<<" "<<mGridDst, SIMWORLD_PLUGINERROR);
		return false;
	}

	return true;
};


//*****************************************************************************
// GridOp: update artificial boundary layer
class fsUpdatePdf : public GridOpBase {
	public:
		fsUpdatePdf(FlagGrid *flags, Grid<Vec3> *ref, Grid<Vec3> *flow, Real scaleFlow, Grid<Vec3> *rans, bool israndom,Real premult) : 
				GridOpBase(), mScaleFlow(scaleFlow), mVortPreMult(premult), mpRef(ref), mpFlow(flow),mpRans(rans),mIsRandom(israndom) { 
			mpFlags = flags;
			applyOperatorToGrids(this);
		};
		~fsUpdatePdf() { }; 
		void resetVariables() { };
		void buildCallList() {
			gaRef.gridAccInit(mpRef, AM_READ, gaCalls); 
			gaFlow.gridAccInit(mpFlow, AM_READWRITE, gaCalls);
			if (mpRans != NULL) gaRans.gridAccInit(mpRans, AM_READ, gaCalls);
			setFlags(mpFlags);
		};

		inline void operator() (int i, int j, int k)
		{ 
			if(!fgIsFluid(getFlagAcc()(i,j,k)))
				return;
			
			if (mIsRandom) // hard-coded test-case
				gaFlow.write(i,j,k) = normNoSqrt(gaRef(i,j,k))==0 ? Vec3(0,0,0) : Vec3(1e3,1e3,1e3);
			else
			{
				Vec3 ref = gaRef(i,j,k);
				Real ref2 = normNoSqrt(ref);
				if (ref2 > 1e-10) // only set if src is not empty
				{
					// use source from rans solver, ref contains surface normal
					if (mpRans != NULL)
					{
						ref = cross(ref, gaRans(i,j,k));
						ref2 = normNoSqrt(ref);
					}
					ref *= mVortPreMult;
					ref2 *= mVortPreMult*mVortPreMult;
					Real flow = dot(gaFlow(i,j,k),ref), flow2 = flow*flow;
					
					// min bounding on reference grid
					if (flow2 < ref2 || flow<0)
							gaFlow.write(i,j,k) = ref;
				}
				gaFlow.write(i,j,k) *= mScaleFlow;
			}	
		}
		void reduce(fsUpdatePdf &op) { };

	protected:
		Real mScaleFlow,mVortPreMult;
		Grid<Vec3> *mpRef, *mpFlow, *mpRans;
		GridAccessor<Vec3,0> gaRef, gaFlow, gaRans;
		bool mIsRandom;
};

//*****************************************************************************
// GridOp : Evaluate pdf, generate vortex particles
class fsGenVPart : public GridOpBase {
	public:
		fsGenVPart(FlagGrid *flags, Grid<Real> *pdf, Grid<Vec3> *flow, Grid<Real> *dist, Real thresVort, VorticitySystem* vsys, Real minRad, Real maxRad, Real mdst, Real threspdf, Real multpdf, Real u02, RandomStream& rnd, 
				bool israndom, Real randPdf, Real maxbl) : 
				GridOpBase(), randStream(rnd), mThresVort(thresVort), mThresPdf(threspdf), mMultPdf(multpdf), mMaxRad(maxRad), mMinRad(minRad), mU02(u02), mMinDist(mdst), mMaxBL(maxbl), mpVSys(vsys),mpPdf(pdf), mpDist(dist),mpFlow(flow),
				mRandom(israndom), mRandomPdf(randPdf) { 
			mpFlags = flags;
			applyOperatorToGridsSimple(this);
		};
		~fsGenVPart() { }; 
		void resetVariables() { };
		void buildCallList() {
			gaFlow.gridAccInit(mpFlow, AM_READWRITE, gaCalls); 
			gaPdf.gridAccInit(mpPdf, AM_WRITE, gaCalls); 
			gaDist.gridAccInit(mpDist, AM_READ, gaCalls); 
			setFlags(mpFlags);
		};

		inline void operator() (int i, int j, int k)
		{
			if (!fgIsFluid(getFlagAcc()(i,j,k))) return;
			
			// obtain squared vorticity
			Vec3 vvort = gaFlow(i,j,k);
			Real vort2 = normNoSqrt(vvort), dist = gaDist(i,j,k);
			if (dist <= 0 || vort2 < mThresVort*mThresVort) 
			{
				gaPdf.write(i,j,k) = 0;
				return;
			}
			// obtain reynolds stress field, pdf
			Real vort=sqrt(vort2);
			Real reyf = vort2 * dist*dist* mU02;
			Real pdf = 2. * (fabs(reyf) - mThresPdf) * mMultPdf;
			
			if (dist < 0 || dist > mMaxBL) pdf=0; // outside allowed range
			if (mRandom) pdf = mRandomPdf; // bypass pdf
			gaPdf.write(i,j,k) = pdf;
			if (pdf<0) pdf=0;
			
			// random test
			Real rnd = randStream.getReal();
			if (rnd >= pdf) return;
			
			// create particle !
			
			// suggested radius
			Real radius = randStream.getReal()*(mMaxRad-mMinRad) + mMinRad, radius2=radius*radius;
			
			// bound kernel radius to wall
			int ceilRad = (int)ceil(radius);
			FOR_KERNEL(ceilRad)
			{					
				Real cr2=x*x+y*y+z*z;
				if (!mRandom && cr2<radius2 && (!mpFlags->checkIndexValid(x+i,y+j,z+k) || fgIsObstacle(mpFlags->getGlobal(x+i,y+j,z+k)) || fgIsInflow(mpFlags->getGlobal(x+i,y+j,z+k))))
					radius2 = cr2;	
			}
			Real sugRad = radius;
			radius = sqrt(radius2);			

			// reduce radius if particle is too big and weak
			Real strength=0;
			Vec3 ndir = vvort * (1./vort);			
			for (;;) // loop to obtain next-best radius
			{
				if (radius < mMinDist) return;
				
				// sum vorticity within particle radius
				Real wsum = 0, rsum=0; 
				radius2=radius*radius;
				ceilRad = (int)ceil(radius);
				FOR_KERNEL(ceilRad)
				{
					int cr2 = x*x+y*y+z*z;
					if (cr2 <= radius2)
					{
						Real w = dot(gaFlow(i+x,j+y,z+k),ndir);
						if (w < 1e-8) continue;	
						wsum += w;
					}
				}
#if DDF_DIMENSION == 3			
				Real mint = VortexParticleGaussian::getAreaMultiplier()*radius*radius2; // gauss inner kernel integral
#else
				Real mint = M_PI*radius2;
#endif
				strength = wsum/mint;				
				
				// test if particle is too big for its strength
				if (strength > vort || mRandom)
					break;
				else
					radius-= 1.0;
			}
			
			// remove vorticity from ABL
			FOR_KERNEL(ceilRad)
			{
				int cr2 = x*x+y*y+z*z;
				if (cr2 <= radius2)
				{
					Real w = dot(gaFlow(i+x,j+y,z+k),ndir);
					if (w < 1e-8) continue;	
					gaFlow.write(i+x,j+y,z+k) -= ndir*w;
				}
			}
			
			// seed particle
			Vec3 str = ndir * (strength);
			if (mRandom)
			{
				// hard coded values for test case
				const Real vmin=0.01,vmax=0.40;
				Real nstr = randStream.getReal()* (vmax-vmin)+vmin;
				str = Vec3(randStream.getReal(),randStream.getReal(),randStream.getReal())*2.-1.;
				str *= nstr/norm(str);
			}
			mpVSys->getParticles().push_back(new VortexParticleGaussian(Vec3(i,j,k), str, radius));
		}
		void reduce(fsGenVPart &op) { };

	protected:
		RandomStream& randStream;

		Real mThresVort,mThresPdf,mMultPdf;
		Real mMaxRad, mMinRad, mU02, mMinDist, mMaxBL;
		VorticitySystem* mpVSys;
		Grid<Real> *mpPdf, *mpDist;
		Grid<Vec3> *mpFlow;
		GridAccessor<Real,0> gaPdf, gaDist;
		GridAccessor<Vec3,0> gaFlow;
		bool mRandom;
		Real mRandomPdf;
};

//*****************************************************************************
// maintain artificial boundary layer, generate vortex particles
class spluginGenVPart : public SolverPlugin {
	public:
		spluginGenVPart() : SolverPlugin(), mRandomStream(13322223) {};
		~spluginGenVPart() { };

		virtual bool parseParams(const ParamSet& params) {
			mGridRef = params.FindOneString("source", "" );				// precomputed ABL source field or surface normal, when using RANS
			mGridFlow = params.FindOneString("flow", "" );				// ABL grid
			mRans = params.FindOneString("rans", "" );					// RANS solver name, if used
			mGridDist = params.FindOneString("dist", "" );				// distance field
			mGridPdf = params.FindOneString("pdf", "" );				// pdf output grid
			mMultPdf = params.FindOneFloat("mult-pdf", 0 );				// c_P constant in the paper
			mThresPdf = params.FindOneFloat("thres-pdf", 0 );			// threshold for pdf below which no seeding will occur
			mMaxBL = params.FindOneFloat("max-bl", 0.1 );				// maximal distance to object in which boundary layer is evaluated
			mThresVort = params.FindOneFloat("thres-vort", 0 );			// threshold for confined vorticity below which no seeding will occur
			mScaleFlow = params.FindOneFloat("scale-flow", 0.99 );		// fade-out for ABL plume
			mMinRad = params.FindOneFloat("min-rad", 0 );				// minimum radius for vortex particles
			mMinDist = params.FindOneFloat("min-dist", mMinRad );		// minimal distance to wall
			mMaxRad = params.FindOneFloat("max-rad", 0 );				// maximum radius for vortex particles
			mIsRandom = params.FindOneInt("random", 0 ) != 0;			// bypass pdf, use random particles
			mRandPdf = params.FindOneFloat("random-pdf", 1.5e-4 );		// probability for random pdf
			mVortPreMult = params.FindOneFloat("vortex-gain", 1. );		// vortex gain, beta constant in the paper
			mFadein = params.FindOneFloat("fade-in", -1 ); 		// vortex particle fade-in constant
			
			if(mIsRandom) debMsg("spluginGenVPart","Randomzied vparts on! Pdf="<<mRandPdf);

			return true;
		};
		virtual bool initPlugin() {
			if (mFadein >= 0) {
				VorticitySystem* vsys = mpPlParams->getFluidSolver()->getVorticitySys();				
				vsys->setFadeIn(mFadein);
			}
			return true; 
		};

		bool performStep(Real org_dt)
		{
			FlagGrid *flags = mpPlParams->getGridInt("flags");
			Grid<Vec3>* ref = mpPlParams->getGridVec3(mGridRef);
			Grid<Real>* pdf = mpPlParams->getGridReal(mGridPdf);
			Grid<Real>* dist = mpPlParams->getGridReal(mGridDist);
			Grid<Vec3>* flow = mpPlParams->getGridVec3(mGridFlow);
			
			Grid<Vec3>* rans = mRans.empty() ? NULL : ddfWorldFindSolver(mRans)->getParams()->getGridVec3("vel-curr");
			VorticitySystem* vsys = mpPlParams->getFluidSolver()->getVorticitySys();
			Real u0 = norm(mpPlParams->mU0);
			if (u0==0) { errMsg("spluginGenVPart","u0 not set !"); exit(1); }

			fsUpdatePdf(flags,ref,flow,mScaleFlow,rans,mIsRandom,mVortPreMult);
			fsGenVPart (flags,pdf,flow,dist,mThresVort * u0, vsys, mMinRad, mMaxRad, mMinDist, mThresPdf,mMultPdf,1./(u0*u0),
					mRandomStream, mIsRandom, mRandPdf, mMaxBL );

			return true;	
		}
	protected:
		// name of src/dst grid
		RandomStream mRandomStream;
		std::string mGridRef,mGridFlow,mGridPdf,mGridDist,mRans;
		Real mThresPdf, mMultPdf, mScaleFlow, mThresVort,mMinDist;
		Real mMinRad, mMaxRad, mVortPreMult, mFadein, mMaxBL;
		bool mIsRandom; //,test
		Real mRandPdf;
};

//*****************************************************************************
// GridOp : calculate boundary layer
class fsCalcABL : public GridOpBaseFlagbord<1> {
	public:
		fsCalcABL(FlagGrid *flags, Grid<Real> *ndist, Grid<Vec3> *norm, Grid<Vec3> *abl, Grid<Vec3> *vel, Real lABL) : 
				GridOpBaseFlagbord<1>(), mpDist(ndist), mpNorm(norm), mpABL(abl), mpVel(vel) { 
			mpFlags = flags;
			mMin = (lABL - .5) / (Real)flags->getMaxSize();
			mMax = (lABL + .5) / (Real)flags->getMaxSize();
			applyOperatorToGridsSimple(this);			
		};
		~fsCalcABL() { }; 
		void resetVariables() { };
		void buildCallList() {
			gaDist.gridAccInit(mpDist, AM_READ, gaCalls); 
			gaVel.gridAccInit(mpVel, AM_READ, gaCalls); 
			gaNorm.gridAccInit(mpNorm, AM_READ, gaCalls); 
			gaABL.gridAccInit(mpABL, AM_WRITE, gaCalls); 
			setFlags(mpFlags);
		};

		inline void operator() (int i, int j, int k)
		{
			if (!fgIsFluid(getFlagAcc()(i,j,k))) return;
			// inside BL
			if (gaDist(i,j,k) > mMin && gaDist(i,j,k) < mMax)
			{
				// get centered vel
				Vec3 v = gaVel(i,j,k);
				if (mpFlags->checkIndexValid(i+1,j,k) && fgIsFluid(getFlagAcc()(i+1,j,k))) v.x = 0.5*(v.x + gaVel(i+1,j,k).x);
				if (mpFlags->checkIndexValid(i,j+1,k) && fgIsFluid(getFlagAcc()(i,j+1,k))) v.y = 0.5*(v.y + gaVel(i,j+1,k).y);
				if (mpFlags->checkIndexValid(i,j,k+1) && fgIsFluid(getFlagAcc()(i,j,k+1))) v.z = 0.5*(v.z + gaVel(i,j,k+1).z);
			
				// compute ABL
				gaABL.write(i,j,k) = -cross(v, gaNorm(i,j,k));
			}
		}
		void reduce(fsCalcABL &op) { };

	protected:
		Real mMin, mMax;
		Grid<Real> *mpDist;
		Grid<Vec3>  *mpNorm, *mpABL, *mpVel;
		GridAccessor<Real,0> gaDist;
		GridAccessor<Vec3,0> gaABL,gaNorm;
		GridAccessor<Vec3,1> gaVel;
};

//*****************************************************************************
// Plugin to calculate boundary layer
class spluginCalcABL : public SolverPlugin {
	public:
		spluginCalcABL() : SolverPlugin() {};
		~spluginCalcABL() { };

		virtual bool parseParams(const ParamSet& params) {
			mGridDist = params.FindOneString("dist", "" );				// distance field
			mGridNorm = params.FindOneString("normal", "" );			// surface normals
			mGridABL = params.FindOneString("abl", "" );				// ABL
			mGridVel = params.FindOneString("mean-vel", "" );			// averaged velocity grid
			mD = params.FindOneFloat("d", 1.5 );						// boundary layer thickness			
			return true;
		};
		virtual bool initPlugin() {return true; };

		bool performStep(Real org_dt)
		{
			FlagGrid *flags = mpPlParams->getGridInt("flags");
			Grid<Vec3>* vel = mpPlParams->getGridVec3(mGridVel);
			Grid<Real>* dist = mpPlParams->getGridReal(mGridDist);
			Grid<Vec3>* norm = mpPlParams->getGridVec3(mGridNorm);
			Grid<Vec3>* abl = mpPlParams->getGridVec3(mGridABL);
			
			fsCalcABL (flags,dist,norm,abl,vel,mD);

			return true;	
		}
	protected:
		// name of src/dst grid
		std::string mGridDist, mGridNorm, mGridVel, mGridABL;
		Real mD;
};

//*****************************************************************************
// GridOp : write database
class fsWriteDatabase : public GridOpBase {
	public:
		fsWriteDatabase(FlagGrid *flags, Grid<Vec3> *norm, Grid<Vec3> *abl, vector<Vec3>& pos, vector<Vec3>& n, vector<Vec3>& val) : 
				GridOpBase(), mpABL(abl), mpNorm(norm), mPos(pos), mNormal(n), mValue(val) { 
			mpFlags = flags;
			applyOperatorToGridsSimple(this);			
		};
		~fsWriteDatabase() { }; 
		void resetVariables() { };
		void buildCallList() {
			gaNormal.gridAccInit(mpNorm, AM_READ, gaCalls); 
			gaABL.gridAccInit(mpABL, AM_READ, gaCalls); 
			setFlags(mpFlags);
		};

		inline void operator() (int i, int j, int k)
		{
			if (!fgIsFluid(getFlagAcc()(i,j,k))) return;
			Vec3 v = gaABL(i,j,k);
			if (v.x == 0 && v.y == 0 && v.z == 0) return;
			
			// store point in pointset
			mPos.push_back(Vec3(i,j,k));
			mNormal.push_back(gaNormal(i,j,k));
			mValue.push_back(v);
		}
		void reduce(fsWriteDatabase &op) { };

	protected:
		Grid<Vec3> *mpABL, *mpNorm;
		GridAccessor<Vec3,0> gaABL,gaNormal;
		vector<Vec3>& mPos;
		vector<Vec3>& mNormal;
		vector<Vec3>& mValue;
};


//*****************************************************************************
// Plugin to write database
class spluginWriteDatabase : public SolverPlugin {
	public:
		spluginWriteDatabase() : SolverPlugin() {};
		~spluginWriteDatabase() { };

		virtual bool parseParams(const ParamSet& params) {
			mGrid = params.FindOneString("grid", "" );				// ABL Grid
			mNormal = params.FindOneString("normal", "" );			// Surface normals
			mU0 = params.FindOneVector("u0", Vec3(0.) );				// inflow/rotation axis
			return true;
		};
		virtual bool initPlugin() {return true; };

		bool performStep(Real org_dt)
		{
			FlagGrid *flags = mpPlParams->getGridInt("flags");
			Grid<Vec3>* abl = mpPlParams->getGridVec3(mGrid);
			Grid<Vec3>* norm = mpPlParams->getGridVec3(mNormal);

			// obtain sparse list
			vector<Vec3> pos, normal, val;
			fsWriteDatabase(flags, norm, abl, pos, normal, val);
			
			if (msNeedHeader) {
				// separate header for translation and rotation. Here only identical resolutions are supported.
				float grid_divisor = 1;
				Vec3 size = vec2R(flags->getSize());
				Vec3 scaling (1.,1.,1.), offset(0,0,0); // in this version these parameters are ignored
				int numPoints = pos.size();
				for (int t=0;t<2;t++) {
					gzwrite(msGzf, &numPoints, sizeof(int));
					gzwrite(msGzf, &grid_divisor, sizeof(float));
					gzwrite(msGzf, &(size[0]), sizeof(Vec3));
					gzwrite(msGzf, &(scaling[0]), sizeof(Vec3));
					gzwrite(msGzf, &(offset[0]), sizeof(Vec3));
					for (int i=0; i<numPoints; i++)	gzwrite(msGzf, &(pos[i][0]), sizeof(Vec3));
					for (int i=0; i<numPoints; i++)	gzwrite(msGzf, &(normal[i][0]), sizeof(Vec3));
					for (int i=0; i<(int)msIndexMat.size(); i++)	{ int idx=msIndexMat[i]+t; gzwrite(msGzf, &(idx), sizeof(int)); }					
				}
				msNeedHeader=false;
			}
			// write ABL data
			gzwrite(msGzf, &(mU0[0]), sizeof(Vec3));
			for (int i=0; i<(int)val.size(); i++)	gzwrite(msGzf, &(val[i][0]), sizeof(Vec3));

			return true;	
		}
		
		// open file, write main header
		static vector<Vec3> open(const string& name, int nTheta, int nPhi) {
			// open database file
			msGzf = gzopen((name+".gz").c_str(), "wb9" );
			if (!msGzf) 
				errFatal("spluginDatabaseWriter", "can't open database file", SIMWORLD_INITERROR);
			
			// write header
			int elements = 2*nTheta*nPhi;
			gzwrite(msGzf, &elements,sizeof(int));
			gzwrite(msGzf, &nTheta,sizeof(int));
			gzwrite(msGzf, &nPhi,sizeof(int));
			
			// generate angle vectors
			float *theta = new float[nTheta];
			for (int i=0;i<nTheta;++i) theta[i] = M_PI * (float)i/(float)(nTheta-1);
			gzwrite(msGzf, theta, sizeof(float)*nTheta);
			float *phi = new float[nPhi];
			for (int i=0;i<nPhi;++i) phi[i] = 2.*M_PI * (float)i/(float)nPhi;
			gzwrite(msGzf, phi, sizeof(float)*nPhi);
			
			// generate index matrices+axes
			vector<Vec3> axes;
			msIndexMat.resize(nTheta*nPhi);
			for (int i=0, idx=0;i<nTheta;i++)
				for (int j=0;j<nPhi;j++) {
					axes.push_back( Vec3( sin(theta[i])*cos(phi[j]), cos(theta[i]), sin(theta[i])*sin(phi[j])) );
					msIndexMat[nPhi*i+j] = idx;	
					idx+=2;
			}
			
			delete[] phi;
			delete[] theta;
			msNeedHeader = true;
			msIndex = 0;
			return axes;
		}
		static void close() { 
			msNeedHeader = false;
			gzclose(msGzf); 
		}
	protected:
		std::string mGrid, mNormal;		
		static gzFile msGzf;
		static bool msNeedHeader;
		static std::vector<int> msIndexMat;
		static int msIndex;
		Vec3 mU0;
};
gzFile spluginWriteDatabase::msGzf = NULL;
bool spluginWriteDatabase::msNeedHeader = false;
vector<int> spluginWriteDatabase::msIndexMat;
int spluginWriteDatabase::msIndex = 0;

// static callers
vector<Vec3> initDatabaseWriter(const string& file, int nTheta, int nPhi) {
	return spluginWriteDatabase::open(file, nTheta, nPhi);
}
void closeDatabaseWriter() {
	spluginWriteDatabase::close();
}


//*****************************************************************************

SolverPlugin* MakeVortexPlugin(std::string name)
{
	if(name.compare( string("advect-vpart") )==0)
		return new spluginAdvectVPart;
	if(name.compare( string("apply-vpart") )==0)
		return new spluginApplyVParts;
	if(name.compare( string("merge-vpart") )==0)
		return new spluginMergeVParts;
	if(name.compare( string("compute-vorticity") )==0)
		return new spluginComputeVorticity;
	if(name.compare( string("interpolate-grid-from") )==0)
		return new spluginInterpolateGridFrom;
	if(name.compare( string("gen-vpart") )==0)
		return new spluginGenVPart;
	if(name.compare( string("calc-abl") )==0)
		return new spluginCalcABL;
	if(name.compare( string("add-database") )==0)
		return new spluginWriteDatabase;
	return NULL;
}

}; // DDF

