/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 *  Plugins to animate boundaries and ABL fields.
 *
 *****************************************************************************/

#include <fstream>
#include "fluidsolver.h"
#include "solverplugin.h"
#include "paramset.h"
#include "quaternion.h"
#include <zlib.h>

#include "matrixbase.h"

using namespace std;

// safety boundary for semi lagrange advection step
#define SLADVBOUND 2

namespace DDF {

class ObjectTransformation {
	public:
		ObjectTransformation(Real v) :
			center(v), pos(v), posVel(v),
			rotMat(v), rotVel(v), 
			scale(v), scaleVel(v) {};

		Vec3 center;
		Vec3 pos, posVel;
		Mat4 rotMat;
		Vec3 rotVel;
		Vec3 scale, scaleVel;
};

//*****************************************************************************
// Gridop : update moving obstacles
class opSetMovingObstacleBC : public GridOpBaseFlagbord<SLADVBOUND> {
	public:
		opSetMovingObstacleBC(FlagGrid *flags, Grid<Vec3> *vel, FlagGrid *obs,
				ObjectTransformation trafo, Real dt, bool deleteOld, Grid<Real> *distsrc, Grid<Real> *distdst, const vector<Vec3>& if1, const vector<Vec3>& if2, bool doforces) : 
			GridOpBaseFlagbord<SLADVBOUND>(), 
			mpVel(vel), mTrafo(trafo), mDt(dt) ,
			mpObs(obs), mDeleteOld(deleteOld), mIF1(), mIF2(), mpDistSrc(distsrc), mpDistDst(distdst), mDoForces(doforces)
		{ 
			mpFlags = flags;

			// don't scale to infinity
			if (normNoSqrt(mTrafo.scale) < 1e-16) mTrafo.scale = 1.;

			// transform positions to src coordinates
			Vec3 srcSize = vec2R(mpFlags->getSize());
			Vec3 destSize = vec2R(mpObs->getSize());
			mCenterFlag = (mTrafo.center + mTrafo.pos) * srcSize;
			mTrafo.pos *= srcSize;
			mTrafo.center *= destSize;
			mTrafo.scale = Vec3(destSize.max()/srcSize.max()) / mTrafo.scale;
			mTrafo.posVel *= srcSize;
			
			for (unsigned i=0; i<if1.size(); i++) mIF1.push_back( if1[i] * destSize);
			for (unsigned i=0; i<if2.size(); i++) mIF2.push_back( if2[i] * destSize);

			cout << "pos " << mTrafo.pos << " scale " << mTrafo.scale << " center " << mTrafo.center << endl; 
			debMsg("opSetMovingObstacleBC","pos="<<mTrafo.pos << " ");

			// assume an additional layer of obstacle cells, so BND=2
			applyOperatorWithBoundToGrid(this, mpFlags, 1);			
		};
		

		~opSetMovingObstacleBC() { }; 
		void resetVariables() { };
		void reduce(opSetMovingObstacleBC &op) { };

		void buildCallList() { 
			this->gaFlags.gridAccInit(mpFlags, AM_WRITE, gaCalls); 
			gaVel.gridAccInit(mpVel, AM_WRITE, gaCalls); 
			gaObs.gridAccInit(mpObs, AM_READ, gaCalls); 

			if(mpDistSrc) gaDistSrc.gridAccInit(mpDistSrc, AM_READ, gaCalls);
			if(mpDistDst) gaDistDst.gridAccInit(mpDistDst, AM_READWRITE, gaCalls);
			//setFlags(mpFlags);
		}; 

		// point transformation
		inline void trafoPoint(Vec3 &pos, Mat4 &mat) {
			pos -= mCenterFlag;
			pos = (mat * pos) * mTrafo.scale;
			pos += mTrafo.center; // to dest coordinates
		}

		inline void operator() (int i, int j, int k) { 
			int srcFlag = FEMPTY; 
			Vec3 currVel = mTrafo.posVel;
			Vec3 pos (i,j,k);
			nVec3i src (1,1,1);

			// with rotation, use trafo matrix
			if (mDoForces)
			{
				Vec3 r = (pos-mCenterFlag);
				currVel += cross( mTrafo.rotVel, r );
			}	
			trafoPoint( pos, mTrafo.rotMat);
			src = vecround( pos);

			// check validity of lookup
			bool isvalid = mpObs->checkIndexValidWithBounds(src.x,src.y,src.z, 1);
			if(isvalid)
				srcFlag = gaObs( src.x,src.y,src.z );
			
			// write new flag grid
			if (mDeleteOld && fgIsObstacle( getFlagAcc()(i,j,k) ) )
				this->gaFlags.write(i,j,k) = FEMPTY;
			
			if (!isvalid) return;

			if (fgIsObstacle(srcFlag) ) {
				gaFlags.write(i,j,k) = FOBSTACLE;
				if (mDoForces) gaVel.write(i,j,k) = currVel; // velocity boundary condition
			}
			else if (mpDistSrc != NULL) {
				// join distance fields, if the exist
				Real ndst = gaDistDst(i,j,k), nsrc = gaDistSrc(src.x,src.y,src.z);
				if (nsrc > 0 && (ndst == 0 || nsrc < ndst)) 
					gaDistDst.write(i,j,k) = nsrc;
			}

			// initialize density source (moving smoke-box)
			for (unsigned n=0; n<mIF1.size(); n++)
			{ 
				if (src.x >= mIF1[n].x && src.y >= mIF1[n].y && src.z >= mIF1[n].z &&
					 src.x <= mIF2[n].x && src.y <= mIF2[n].y && src.z <= mIF2[n].z)
					this->gaFlags.write(i,j,k) = gaFlags(i,j,k) | FDENSITYSOURCE;
			}

			// done...
		};
	protected:
		Grid<Vec3> *mpVel;
		GridAccessor<Vec3,1> gaVel; 

		ObjectTransformation mTrafo;

		Real mDt;
		FlagGrid *mpObs;
		GridAccessor<int,0> gaObs;
		Vec3 mCenterFlag;
		bool mDeleteOld;

		Mat4 mMatRot,mMatRotNext;

		vector<Vec3> mIF1, mIF2;
		Grid<Real> *mpDistSrc, *mpDistDst;
		GridAccessor<Real,1> gaPhiDebug; 
		GridAccessor<Real,0> gaDistSrc, gaDistDst; 
		bool mDoForces;
}; 

//*****************************************************************************
// Plugin for moving obstacles with linear speed
class spSetMovingObstacleBC : public SolverPlugin {
	public:
		spSetMovingObstacleBC() : 
					SolverPlugin(), mNameFlags("flags"),
   				mNameVel("vel-curr"), mObstacle(""),mTrafo(0.)
			{ }; 
		~spSetMovingObstacleBC() { };

		// init plugin, return failure
		virtual bool parseParams(const ParamSet& params) {
			mNameFlags = params.FindOneString("flags",mNameFlags);                   // flag grid
			mNameFlagsSrc = params.FindOneString("flags-src","");                    // static background flag grid
			mNameVel = params.FindOneString("vel","vel-curr");                       // velocity grid
			mObstacle = params.FindOneString("obstacle",mObstacle);                  // flag grid of the moving object
			
			mTrafo.posVel    = params.FindOneVector("obs-vel",mTrafo.posVel);        // linear velocity
			mTrafo.pos     = params.FindOneVector("obs-pos",mTrafo.pos);             // initial position

			// center of rotation/scale
			mTrafo.center   = params.FindOneVector("obs-center",mTrafo.center);      // rotation pivot
			mRotAxis 		 = params.FindOneVector("obs-rot-axis",Vec3(0.));        // rotation axis
			mRotAngle       = params.FindOneFloat("obs-rot-angle",0.);               // initial rotation angle (degrees)
			mRotVel 			 = params.FindOneFloat("obs-rot-vel",0.);            // rotation velocity      (degrees/timestep)
			mTrafo.scaleVel = params.FindOneVector("obs-scale-vel",mTrafo.scaleVel); // scale velocity
			mTrafo.scale    = params.FindOneVector("obs-scale",mTrafo.scale);        // initial scaling
			mRotAngle *= M_PI/180.;
			mRotVel *= M_PI/180.;
			normalize(mRotAxis);
			mTrafo.rotVel = mRotAxis * mRotVel;
			return true;
		};
		virtual bool initPlugin() {
			return true;
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("spSetMovingObstacleBC","step "<<dt);
			const Real mStrength = 0.5;

			// get grids
			FlagGrid*   flags    = mpPlParams->getGridInt(mNameFlags);
			FlagGrid*   flagsrc  = mNameFlagsSrc.empty() ? NULL : mpPlParams->getGridInt(mNameFlagsSrc);
			Grid<Vec3>* vel      = mpPlParams->getGridVec3(mNameVel);
			FlagGrid*   obs      = mpPlParams->getGridInt(mObstacle);
			Grid<Real>* phiCurr  = NULL; // mpPlParams->getGridReal("phi-curr");
			
			// copy background grid
			if (flagsrc) goCopyGrid<int>(flags,flagsrc);
			
			// set flags & velocities inside cells
			Quat q(mRotAxis,mRotAngle);
			mTrafo.rotMat = q.getRotMat();
			opSetMovingObstacleBC(flags, vel, obs, mTrafo, dt, flagsrc==NULL, NULL, NULL, vector<Vec3>(),vector<Vec3>(),true);
			
			// initialize obstacle boundaries with velocities
			FOR_IJK_GRID(flags) {
				if (!fgIsObstacle( flags->getGlobal(i,j,k) ) ) {
					if (i>0 && fgIsObstacle( flags->getGlobal(i-1,j,k) ) ) {
						vel->getGlobal(i,j,k)[0] = vel->getGlobal(i-1,j,k)[0];
					}
					if (j>0 && fgIsObstacle( flags->getGlobal(i,j-1,k) ) ) {
						vel->getGlobal(i,j,k)[1] = vel->getGlobal(i,j-1,k)[1];
					}
					if(gDim==3 && (k>0 && fgIsObstacle( flags->getGlobal(i,j,k-1) ) )) {
						vel->getGlobal(i,j,k)[2] = vel->getGlobal(i,j,k-1)[2];
					}
				} 
				else {
					// reset level set inside of the obstacle?
					if (phiCurr) {
						if(phiCurr->getGlobal(i,j,k)<0.) phiCurr->getGlobal(i,j,k) = 0.5;
					}
				}
			}

			// update pos, rot, scale
			mTrafo.pos += mTrafo.posVel * dt;
			mRotAngle += mRotVel * dt;
			mTrafo.scale += mTrafo.scaleVel * dt;

			return true;
		};

	protected:
		std::string mNameFlags, mNameVel, mObstacle, mNameFlagsSrc;
		ObjectTransformation mTrafo;
		Vec3 mRotAxis;
		Real mRotAngle, mRotVel;
};

//*****************************************************************************
// A set of precalculated boundary layers
class VortexDBSet {
	public:
	// load from file
	VortexDBSet(gzFile& gzf, int numPhi, int numTheta, const Vec3& sizeFlags)
	{
		Real gd;
		if (gzread(gzf,&numPoints,sizeof(int)) != sizeof(int)) { errFatal("spAnimate::initPlugin","gzread failed", SIMWORLD_INITERROR);  }
		if (gzread(gzf,&gd,sizeof(Real)) != sizeof(Real)) { errFatal("spAnimate::initPlugin","gzread failed", SIMWORLD_INITERROR);  }		
		if (gzread(gzf,&mSize[0],sizeof(Vec3)) != sizeof(Vec3)) { errFatal("spAnimate::initPlugin","gzread failed", SIMWORLD_INITERROR);  }
		if (gzread(gzf,&mVScale[0],sizeof(Vec3)) != sizeof(Vec3)) { errFatal("spAnimate::initPlugin","gzread failed", SIMWORLD_INITERROR);  }
		if (gzread(gzf,&mVOffs[0],sizeof(Vec3)) != sizeof(Vec3)) { errFatal("spAnimate::initPlugin","gzread failed", SIMWORLD_INITERROR);  }
		mPos.resize(numPoints);
		mNoDir.resize(numPoints);
		for (int i=0;i<numPoints;i++)
			if (gzread(gzf,&mPos[i][0],sizeof(Vec3)) != sizeof(Vec3)) { errFatal("spAnimate::initPlugin","gzread failed", SIMWORLD_INITERROR);  }
		for (int i=0;i<numPoints;i++)
			if (gzread(gzf,&mNoDir[i][0],sizeof(Vec3)) != sizeof(Vec3)) { errFatal("spAnimate::initPlugin","gzread failed", SIMWORLD_INITERROR);  }
		mLookup.resize(numTheta);
		for (int i=0;i<numTheta;i++)
		{
			mLookup[i] = vector<int>(numPhi);
			if (gzread(gzf,&mLookup[i][0],numPhi*sizeof(int)) != (int)(numPhi*sizeof(int))) { errFatal("spAnimate::initPlugin","gzread failed", SIMWORLD_INITERROR);  }			
		}
		gd /= sizeFlags.max();
		mIGriddiv2 = 1./ ( gd*gd);
		mIScale = sizeFlags / (mVScale * mSize);
		//cout << "load set : N " << mSize << " vscale " << mVScale<< " offs " << mVOffs << endl << numPoints << " items" << endl;
	}

	// bilinear interpolation of the nearest datapoints in the set
	void bilinearChooseSets(const vector<Real>& mPhi, const vector<Real>& mTheta, const vector<Vec3>& mSetDir, Real& a11, Real& a12, Real& a21, Real& a22, int& idx11, int& idx12, int& idx21, int& idx22, const Vec3& vel, Real& mult) const
	{
		Vec3 dir = vel;
		mult = normalize(dir); cout << " norm dir " << mult;
		if (mult < 1e-10){
			mult=0; return;
		}
		mult /= mNorm; cout << " set norm " << mNorm << " total " << mult << endl;
		Real phi, theta;
		vecToAngle(dir, phi, theta);
		
		// find linear coefficients
		int idxp1=0, idxp2=mPhi.size()-1;
		for (unsigned i=0;i<mPhi.size();i++)
		{
			if (mPhi[i] <= phi && mPhi[i] > mPhi[idxp1]) idxp1=i;
			if (mPhi[i] >= phi && mPhi[i] < mPhi[idxp2]) idxp2=i;
		}		
		Real ap = (idxp1 == idxp2) ? 0 : (phi - mPhi[idxp1]) / (mPhi[idxp2] - mPhi[idxp1]);
		
		int idxt1=0, idxt2=mTheta.size()-1;
		for (unsigned i=0;i<mTheta.size();i++)
		{
			if (mTheta[i] <= theta && mTheta[i] > mTheta[idxt1]) idxt1=i;
			if (mTheta[i] >= theta && mTheta[i] < mTheta[idxt2]) idxt2=i;
		}
		Real at = (idxt1 == idxt2) ? 0 : (theta - mTheta[idxt1]) / (mTheta[idxt2] - mTheta[idxt1]);

		idx11 = mLookup[idxt1][idxp1]; a11 = (1.-at)*(1.-ap);
		idx12 = mLookup[idxt1][idxp2]; a12 = (1.-at)*(ap);
		idx21 = mLookup[idxt2][idxp1]; a21 = (at)*(1.-ap);
		idx22 = mLookup[idxt2][idxp2]; a22 = (at)*(ap);

		//printf("phi %g theta %g : chose phi= %g(%g)%g ,theta= %g(%g)%g\n",phi*180./M_PI,theta*180./M_PI,mPhi[idxp1]*180./M_PI,ap,mPhi[idxp2]*180./M_PI,mTheta[idxt1]*180./M_PI,at,mTheta[idxt2]*180./M_PI );
		//cout << vel << " to " << mSetDir[idx11]*mult << " " << mSetDir[idx12]*mult << " " << mSetDir[idx21]*mult << " " << mSetDir[idx22]*mult << endl;
		//cout << idx11 << "," << idx12 << "," << idx21 << "," << idx22 << "," <<endl;
		//cout << "final : " << mult*(a11*mSetDir[idx11]+a12*mSetDir[idx12]+a21*mSetDir[idx21]+a22*mSetDir[idx22]) << " to orig " << vel <<endl;
		
	}

	// map confined vorticity to simulation grid
	void mapdown(FlagGrid *flags, Grid<Vec3>* dst, const Mat4& rot, const Vec3& pos, const Vec3& nvec, const vector<vector<Vec3> >& mSetVort, const vector<Real>& mPhi, const vector<Real>& mTheta, const vector<Vec3>& mSetDir, const Vec3& cnt, const Vec3& sc)
	{
		Real mult, a11, a12, a21, a22;
		int idx11, idx12, idx22, idx21;
		bilinearChooseSets(mPhi,mTheta,mSetDir, a11,a12,a21,a22,idx11,idx12,idx21,idx22,nvec,mult);
		if (mult == 0) return;

		Vec3 srcSize = vec2R(flags->getSize());
		Vec3 centerSrc = (cnt + pos) * srcSize;
		Vec3 centerDst = cnt * mSize;
		Vec3 iscale = Vec3(srcSize.max()/mSize.max()) * sc;
		cout << "mult " << mult << " gridd " << iscale << " total " << mult*mIGriddiv2 <<endl;
		//mult /= mIGriddiv2;
		Real pt,pp;
		vecToAngle(nvec,pp,pt);
		for (unsigned n=0;n<mPos.size();n++)
		{
			Vec3 p = (mPos[n] - centerDst) * iscale;
			p = rot * p;
			p += centerSrc;
			// transfer vorticity sources
			nVec3i npos = vecround(p);
			
			if (!flags->checkIndexValid(npos.x,npos.y,npos.z)) continue;
			//printf("point %d : set %g %g %g grid %g %g %g\n",n,mPos[n].x,mPos[n].y,mPos[n].z,p.x,p.y,p.z);
			
			if (fgIsFluid(flags->getGlobal(npos.x,npos.y,npos.z)))
			{
				Vec3 vsi = a11*mSetVort[idx11][n];
				vsi += a12*mSetVort[idx12][n];
				vsi += a21*mSetVort[idx21][n];
				vsi += a22*mSetVort[idx22][n];
				Vec3 vs = rot * vsi;
				(*dst)(npos.x,npos.y,npos.z) += -vs * mult;
				//cout << dst->getName() << ":" << npos << " <<= " <<-vs*mult << " : " << (*dst)(npos.x,npos.y,npos.z) << endl;
			}
		}
	}

	vector<vector<int> > mLookup;	
	vector<Vec3> mPos, mNoDir;
	Vec3 mVScale, mVOffs, mSize, mIScale;
	Real mIGriddiv2, mNorm;
	int numPoints;
};

//*****************************************************************************
// Interface for precalculated DB
class VortexDB
{
	public:
	// load from file
	VortexDB(const string& file, const Vec3& sizeFlags)
	{
		// load database
		int numSets, numPhi, numTheta;
		gzFile gzf = gzopen(file.c_str(),"rb");
		if (gzf == NULL) { errFatal("spAnimate::initPlugin","can't open vortex source database",SIMWORLD_INITERROR);exit(1);  }
		if (gzread(gzf,&numSets,sizeof(int)) != sizeof(int)) { errFatal("spAnimate::initPlugin","gzread numSets failed",SIMWORLD_INITERROR);  }
		if (gzread(gzf,&numTheta,sizeof(int)) != sizeof(int)) { errFatal("spAnimate::initPlugin","gzread numTheta failed",SIMWORLD_INITERROR);  }
		if (gzread(gzf,&numPhi,sizeof(int)) != sizeof(int)) { errFatal("spAnimate::initPlugin","gzread numPhi failed",SIMWORLD_INITERROR);  }
		mTheta.resize(numTheta);
		mPhi.resize(numPhi);
		if (gzread(gzf,&mTheta[0],numTheta*sizeof(Real)) != (int)(numTheta*sizeof(Real))) { errFatal("spAnimate::initPlugin","gzread theta failed",SIMWORLD_INITERROR);  }
		if (gzread(gzf,&mPhi[0],numPhi*sizeof(Real)) != (int)(numPhi*sizeof(Real))) { errFatal("spAnimate::initPlugin","gzread phi failed",SIMWORLD_INITERROR);  }		 
		mTrans = new VortexDBSet(gzf, numPhi, numTheta, sizeFlags);
		mRot = new VortexDBSet(gzf, numPhi, numTheta, sizeFlags);
		
		mSetDir.resize(numSets);
		mSetVort.resize(numSets);
		for (int i=0;i<numSets;i++)
		{
			VortexDBSet *cset = ((i & 1) == 0) ? mTrans : mRot;
			if (gzread(gzf,&mSetDir[i][0],sizeof(Vec3)) != sizeof(Vec3)) { errFatal("spAnimate::initPlugin","gzread u0 failed",SIMWORLD_INITERROR);  }
			
			mSetVort[i] = vector<Vec3>(cset->numPoints);
			for (int j=0;j<cset->numPoints;j++)
				if (gzread(gzf,&mSetVort[i][j][0],sizeof(Vec3)) != sizeof(Vec3)) { errFatal("spAnimate::initPlugin","gzread data failed",SIMWORLD_INITERROR);  }			
		}
		mTrans->mNorm = norm(mSetDir[0]);
		mRot->mNorm = norm(mSetDir[1]);

		gzclose(gzf);
		debMsg("VortexDB::VortexDB()", "VDB loaded, " << numSets << " sets");
	}
	~VortexDB()
	{
		delete mRot;
		delete mTrans;
	}

	// map confined vorticity to simulation grid
	void project(FlagGrid *flags, Grid<Vec3>* dst, const ObjectTransformation& trafo, const Vec3& uvel)
	{
		Mat4 matRotBack = trafo.rotMat;
		matRotBack.transpose();
		Vec3 uv = trafo.rotMat * uvel;
		Vec3 urot = matRotBack * trafo.rotVel;
		cout << "linear " << uv << " (unrot " << uvel << ") rot " << urot << " (unrot " << trafo.rotVel << ")" <<endl;

		mTrans->mapdown(flags, dst, matRotBack, trafo.pos, uv, mSetVort, mPhi, mTheta, mSetDir, trafo.center, trafo.scale);	
		mRot->mapdown(flags, dst, matRotBack, trafo.pos, urot, mSetVort, mPhi, mTheta, mSetDir, trafo.center, trafo.scale);	
	}

	private:
	vector<Real> mPhi, mTheta;		
	vector<vector<Vec3> > mSetVort;	
	vector<Vec3>  mSetDir;
	VortexDBSet *mRot, *mTrans;
};

// Helper class for loading animation data
class AnimFrame
{
	public:
	AnimFrame(ifstream& ifs)
	{ 
		Real rx,ry,rz; Mat4 mat;
		ifs >> time >> pos.x >> pos.y >> pos.z >> rx >> ry >> rz; ifs.ignore(1000,'\n');
		rot = Quat(rx*M_PI/180.,ry*M_PI/180.,rz*M_PI/180.);
	}

	Vec3 pos, pvel, rvel;
	Quat rot;
	Real time;
};

//*****************************************************************************
// load and playback animation
class spAnimate : public SolverPlugin {
	public:
		spAnimate() : 
					SolverPlugin(),mFrames(),mTime(0.),mVDB(NULL), mIflow1(), mIflow2()
			{ }; 
		~spAnimate() { 
			if (mVDB!=NULL) delete mVDB;
		};

		// init plugin, return failure
		virtual bool parseParams(const ParamSet& params) {
			mNameSrc = params.FindOneString("flags-src","unnamed");      // static background flags
			mNameDst = params.FindOneString("flags-dst","flags");        // flag grid
			mDistSrc = params.FindOneString("ndist-src","");             // static background distance field
			mDistDst = params.FindOneString("ndist-dst","");             // distance field
			mVortSrc = params.FindOneString("vort-src","");              // static background ABL source
			mVortDst = params.FindOneString("vort-dst","");              // ABL source field
			mVortFile = params.FindOneString("vort-file","");            // vortex DB filename
			mDoForces = params.FindOneInt("do-forces",1)!=0;             // enforce velocities on the moving object boundary
			mNameVel = params.FindOneString("vel","vel-curr");           // velocity field
			mFile = params.FindOneString("file","no-filename");          // animation control file
			mAnimParams = params.FindOneString("params","unnamed");      // params where animated object is stored
			return true;
		};
		virtual bool initPlugin() {
			
			// load animation stream from text file
			ifstream ifs (mFile.c_str());
			if (!ifs) {
				errMsg("spAnimate::initPlugin","can't load animation file !"); return false;
			}
			int numFrames, numInflows;
			Vec3 lastPos(0.), lastAxis(0.);
			Real delta, lastTime(0.);
			ifs >> numFrames; ifs.ignore(1000,'\n');
			ifs >> delta; ifs.ignore(1000,'\n');
			ifs >> mScale; ifs.ignore(1000,'\n');
			ifs >> mRotCenter.x >> mRotCenter.y >> mRotCenter.z; ifs.ignore(1000,'\n');
			// obtain frame data
			for (int n=0;n<numFrames;n++)
			{
				// load frame, calculate speed by finite differences
				AnimFrame frame(ifs);
				frame.time *= delta;
				frame.pvel = (frame.pos - lastPos) / (frame.time - lastTime);

				Vec3 axis = frame.rot.getAxis();
				frame.rvel =  (axis - lastAxis) / (frame.time - lastTime);
				mFrames.push_back(frame);
				lastPos = frame.pos;
				lastAxis = axis;
				lastTime = frame.time;
			}
			mFrames[0].pvel = mFrames[1].pvel;
			mFrames[0].rvel = mFrames[1].rvel;

			// load inflow box information
			ifs >> numInflows;  ifs.ignore(1000,'\n');
			for (int n=0;n<numInflows;n++)
			{
				Vec3 b1,b2;
				ifs >> b1.x >> b1.y >> b1.z >> b2.x >> b2.y >> b2.z;  ifs.ignore(1000,'\n');
				mIflow1.push_back(b1);
				mIflow2.push_back(b2);
			}
			// load vortex database
			if (!mVortFile.empty())
				mVDB = new VortexDB(mVortFile,vec2R(mpPlParams->getFluidSolver()->getGridFlags()->getSize()));
			return true;
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("spAnimate","step "<<dt);
			
			// linear interpolatation of the animation frames
			Real multiplier = mpPlParams->getMultiplier();
			if (multiplier < 1) multiplier = 1;
			dt /= multiplier;
			mTime += dt;
			vector<AnimFrame>::const_iterator f1=mFrames.begin(), f2=mFrames.end()-1, it=mFrames.begin();
			for (;it != mFrames.end();++it)
			{
				if (it->time <= mTime && it->time > f1->time) f1 = it;
				if (it->time >= mTime && it->time < f2->time) f2 = it;
			}
			Real a2 = (f1->time == f2->time) ? 0 : (mTime - f1->time) / (f2->time - f1->time), a1=1.-a2;
			
			// setup transformation matrix
			Vec3 finalPos = f1->pos * a1 + f2->pos * a2;
			Quat finalQ = Quat::slerp(f1->rot, f2->rot, a2);
			
			ObjectTransformation trafo(0.);
			Vec3 ax = finalQ.getAxis() * 180./M_PI; 
			//printf("QUAT %g %g %g %g\nAXIS %g %g %g\n",finalQ.x,finalQ.y,finalQ.z,finalQ.w,ax.x,ax.y,ax.z);
			trafo.rotMat = finalQ.getRotMat();
			trafo.pos = finalPos;	
			trafo.posVel = f2->pvel; 
			trafo.rotVel = f2->rvel;
			trafo.scale = mScale;
			trafo.center = mRotCenter;
			cout << "pvel : " << trafo.posVel << endl;

			// log data for output in pbrt files
			/*Vec3 rpos = (finalPos * mDomainSize) + mDomain0;
			ostringstream str;
			Real angle = normalize(ax);
			if (angle < 1e-8) ax = Vec3(0.,1.,0.);
			str << rpos.x << " " << rpos.y << " "  << rpos.z << " "  << ax.x << " "  << ax.y << " "  << ax.z << " " << -angle;
			mpPlParams->getFluidSolver()->setDebugBuffer(str.str());*/

			// get grids
			FlagGrid*   src = mpPlParams->getGridInt(mNameSrc);
			FlagGrid*   dst = mpPlParams->getGridInt(mNameDst);
			FlagGrid*   obs = ddfWorldFindSolver(mAnimParams)->getGridFlags();
			Grid<Real>*   distsrc = mDistSrc.empty() ? NULL : mpPlParams->getGridReal(mDistSrc);
			Grid<Real>*   distdst = mDistDst.empty() ? NULL : mpPlParams->getGridReal(mDistDst);
			Grid<Real>*   distobs = mDistSrc.empty() ? NULL : ddfWorldFindSolver(mAnimParams)->getParams()->getGridReal(mDistSrc);
			Grid<Vec3>*   vortsrc = mVortSrc.empty() ? NULL : mpPlParams->getGridVec3(mVortSrc);
			Grid<Vec3>*   vortdst = mVortDst.empty() ? NULL : mpPlParams->getGridVec3(mVortDst);
			Grid<Vec3>* vel = mpPlParams->getGridVec3(mNameVel);
			
			cout << "time " << mTime << " pos " << trafo.pos << endl;
			//cout << "lookup bounds: " << trafo.pos * vec2R(obs->getSize())  << " to " << mScale * vec2R(obs->getSize()) + trafo.pos * vec2R(obs->getSize()) << endl;
			
			// copy static fields
			goCopyGrid<int>(dst,src);
			if (distsrc != NULL) goCopyGrid<Real>(distdst,distsrc);
			if (vortsrc != NULL) goCopyGrid<Vec3>(vortdst,vortsrc);

			// set flags & velocities inside cells
			opSetMovingObstacleBC(dst, vel, obs, trafo, dt, false, distobs, distdst, mIflow1, mIflow2, mDoForces);

			// query vortex DB
			Vec3 uvel = -trafo.posVel * vec2R(dst->getSize()) ;//+ mpPlParams->mU0 ;
			if (mVDB != NULL) 
				mVDB->project(dst, vortdst, trafo,uvel);
 
			// initialize obstacle boundaries with velocities
			if (mDoForces) {
				FOR_IJK_GRID_BND(dst,1) {
				if (!fgIsObstacle( dst->getGlobal(i,j,k) ) ) {
					if (fgIsObstacle( dst->getGlobal(i-1,j,k) ) && !fgIsObstacle( src->getGlobal(i-1,j,k) ) ) {
						vel->getGlobal(i,j,k)[0] = vel->getGlobal(i-1,j,k)[0];
					}
					if (fgIsObstacle( dst->getGlobal(i,j-1,k) ) && !fgIsObstacle( src->getGlobal(i,j-1,k) )) {
						vel->getGlobal(i,j,k)[1] = vel->getGlobal(i,j-1,k)[1];
					}
					if(gDim==3 && (fgIsObstacle( dst->getGlobal(i,j,k-1) ) && !fgIsObstacle( src->getGlobal(i,j,k-1) ) )) {
						vel->getGlobal(i,j,k)[2] = vel->getGlobal(i,j,k-1)[2];
					}
				} 
			} }

			return true;
		};

	protected:
		// grid access
		std::string mNameSrc, mNameDst, mNameVel, mAnimParams, mFile, mDistSrc, mDistDst, mVortSrc, mVortDst, mVortFile;
		vector<AnimFrame> mFrames;
		Real mTime, mScale;
		Vec3 mRotCenter;
		VortexDB* mVDB;
		vector<Vec3> mIflow1, mIflow2;
		bool mDoForces;
};

//*****************************************************************************
// create plugins
SolverPlugin* MakeAnimPlugin(std::string name) {
	if(name.compare( "set-moving-obs-bcs" )==0) {
		return new spSetMovingObstacleBC;
	} else if(name.compare( "animate" )==0) {
		return new spAnimate;
	}
	
	return NULL;
}

}
