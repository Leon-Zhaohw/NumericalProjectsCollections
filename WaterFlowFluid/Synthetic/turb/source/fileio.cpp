/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 *  Plugins for loading and saving grids 
 *
 *****************************************************************************/

// lsdebug
#include "fluidsolver.h"
#include "solverplugin.h"
#include "paramset.h"
#include "vortexpart.h"
#include <fstream>
#include <zlib.h>

// global debugging var - reuse frame nr for badtrioutput, used in spluginDumpMeshSurf
int g_meshSurfaceOutNr = 0;
using std::endl;

namespace DDF { 

// Header for universal grid file  
typedef struct {
	char id[4];
	int dimX, dimY, dimZ;
	int frames, elements, elementType, bytesPerElement, bytesPerFrame;
} UniversalHeader;


// print timing stats
class spluginPrintTimingStats : public SolverPlugin {
	public:
		spluginPrintTimingStats() : 
			SolverPlugin(), mSort(0) { };

		~spluginPrintTimingStats() { };

		virtual bool parseParams(const ParamSet& params) {
			mSort = params.FindOneInt("sort", mSort );
			return true;
		};
		virtual bool initPlugin() { return true; };


		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {

			debMsg("spluginPrintTimingStats", mpPlParams->getFluidSolver()->getTimingStats(mSort) );

			return true;
		};

	protected:
		int mSort;
};


/******************************************************************************/
// helper functions for isosurface 



/******************************************************************************/

class spluginDumperBase : public SolverPlugin {
	public:
		spluginDumperBase() : SolverPlugin(),
   				mSimTime(-1.), mMaxFrames(-1), 
				mAnimOutCounter(0), mLastAniOutTime(0.), mStart(0.), mDumpAllFrames(false), mLastDumpTime(-10.),
  				mOutputOverrideName(""), mOutname("")
					{ };

		~spluginDumperBase() { };

		bool parseParamsDumperBase(const ParamSet& params) {
			mOutputOverrideName = params.FindOneString("override-name", mOutputOverrideName );
			mDumpAllFrames = params.FindOneInt("dump-all-frames", mDumpAllFrames );
			mStart = params.FindOneFloat("start-time", mStart );
			mMaxFrames = params.FindOneInt("max-frames", mMaxFrames );
			
			// do we have a override name param? if yes use...
			// otherwise take defaults from params
			if(mOutputOverrideName.length()>0)
				mOutname = mOutputOverrideName;
			else 
				mOutname = this->mpPlParams->mAnimOutputFile;

			mLastAniOutTime = mStart;
			return true;
		}

		// check if current frame should be dumped (note, needs valid mSimTime!)
		bool doDump() {
			// dump each frame once
			if(mDumpAllFrames) {
				if(mLastDumpTime<mSimTime) {
					mLastDumpTime = mSimTime;
					return true;
				}
				return false;
			}

			// dump according to timestep
			return (mpPlParams->mTimestepAnim>0. && 
					(mSimTime-mLastAniOutTime > mpPlParams->mTimestepAnim) &&
					(!mpPlParams->getQuit()) );
		}

		void increaseDumpCounter() {
			mAnimOutCounter++;
			if (mMaxFrames>=0) {
				// signal end
				if (mAnimOutCounter>=mMaxFrames) mpPlParams->setQuit(true);
			}
		}

	protected:
		Real mSimTime; 
		// max no of frames param
		int mMaxFrames;
		// no of output frames
		int mAnimOutCounter;
		// output time of last animation frame
		Real mLastAniOutTime, mStart; 
		// dump each frame?
		int mDumpAllFrames;
		double mLastDumpTime;

		// override name from params?
		std::string mOutputOverrideName, mOutname;
};


template<class T>
class fsLoadData : public GridOpBase { 
	public:
		fsLoadData(FlagGrid *flags, Grid<T> *dst, unsigned char *data, const UniversalHeader& head, bool setonlyFluid, int interpol) :
				GridOpBase(), mpDst(dst),mpData(data), onlyFluid(setonlyFluid), mInterpolateWithBorder(interpol) {
			mpFlags = flags;
			mSizeRead = nVec3i(head.dimX, head.dimY, head.dimZ);
			mSizeDst = mpFlags->getSize() -nVec3i(1,1,1);
			mPageSize = head.bytesPerElement * head.dimX * head.dimY * head.dimZ;
			mBytesPerElement = head.bytesPerElement;
			applyOperatorToGridsSimple(this);
		};
		~fsLoadData() {}
		
		void resetVariables() { };
		void buildCallList() {
			gaDst.gridAccInit(mpDst, AM_WRITE, gaCalls); 
			setFlags(mpFlags);
		};

		inline T getRaw (unsigned chunk);
		
		// load & safety check access of file data at i,j,k
		inline T getFileData (int i, int j, int k) { 
			CLAMP(i, 0, mSizeRead[0]-1);
			CLAMP(j, 0, mSizeRead[1]-1);
			CLAMP(k, 0, mSizeRead[2]-1);
			
			if (gDim==2) k=0;
			unsigned chunk = i + mSizeRead[0]*(j + mSizeRead[1]*k);			
			return getRaw(chunk);
		}
		
		inline void operator() (int i, int j, int k) { 
			nVec3i src = nVec3i(i,j,k);
			const nVec3i dst = nVec3i(i,j,k);

			if(mInterpolateWithBorder>-1) {
				for(int l=0; l<gDim; l++) {
					if(dst[l]<mInterpolateWithBorder) {
						// keep
					} else if(dst[l]>=mSizeDst[l]-mInterpolateWithBorder) {
						// keep upper grid boundary as well
						src[l] = mSizeRead[l]-1 + (dst[l]-mSizeDst[l]);
					} else {
						//src[l] = dst[l]* mSizeRead[l] / mSizeDst[l];
						const int b = mInterpolateWithBorder;
						src[l] = ((dst[l]-b) * (mSizeRead[l]-2*b) / (mSizeDst[l]-2*b)) + b;
					}
				}
				//std::cerr <<"Init debug" << "at dst="<<dst<<"|"<<mSizeDst<<"    file src="<<src<<", filesize"<<mSizeRead <<"\n";
			} else {
				// direct read (standard way of loading)
				// nothing to do, load from (i,j,k)
			}

			if (!onlyFluid || fgIsFluid(getFlagAcc()(i,j,k)) ) 
				gaDst.write(i,j,k) = getFileData(src[0],src[1],src[2]); 			
		};
		void reduce(fsLoadData &op) { };

	protected:
		Grid<T> *mpDst;
		GridAccessor<T,0> gaDst;
	   
		unsigned char * mpData;
		bool onlyFluid;
		nVec3i mSizeRead; // size of grid in file
		nVec3i mSizeDst;   // size of destination grid
		int mInterpolateWithBorder;
		unsigned mPageSize, mBytesPerElement;
};

template<class T>
inline T fsLoadData<T>::getRaw(unsigned chunk)
{
	return *((T*) &mpData[ mBytesPerElement * chunk]);
}
template<>
inline Vec3 fsLoadData<Vec3>::getRaw(unsigned chunk)
{
	Vec3 a;
	for (int e=0; e < 3; e++)
		a[e] = *((Real*) &mpData[ mBytesPerElement * chunk + e * mPageSize]);
	return a;
}

// query dimensions of universal grid file
nVec3i getUniversalGridSize(const string& file)
{
	UniversalHeader head;
	
	gzFile gzf = gzopen(file.c_str(),"r");
	if (gzf == NULL) {
		errFatal("spluginLoadUniversal","can't open field file " << file, SIMWORLD_INITERROR);
		return nVec3i(0);
	}
	if (gzread(gzf, &head,sizeof(UniversalHeader)) != sizeof(UniversalHeader) || strncmp(head.id, "DDF1", 4)) {
		errFatal("spluginLoadUniversal","wrong filetype in " << file,SIMWORLD_INITERROR);
		return nVec3i(0);
	}
	gzclose(gzf);
		
	return nVec3i(head.dimX, head.dimY, head.dimZ);
}

// dump any grid to universal file format
class spluginDumpUniversal : public spluginDumperBase {
	public:
		spluginDumpUniversal() : spluginDumperBase(), mFirstTime(true) { };
		~spluginDumpUniversal() { };

		virtual bool parseParams(const ParamSet& params) {
			mGrid = params.FindOneString("grid", "");
			mMAC = params.FindOneInt("mac", !mGrid.compare("vel-curr")) != 0; // MAC grid is default for 'vel-curr' only
			mFilePerFrame = params.FindOneInt("file-per-frame", 0) != 0;
			mNoExt = params.FindOneInt("no-ext", 0) != 0;
			bool singleDump = params.FindOneInt("single-dump", 0) != 0;
			parseParamsDumperBase(params);
			if (singleDump) { mNoExt = true; mFilePerFrame = true; mDumpAllFrames = true; }
			return true;
		};
		virtual bool initPlugin() {
			mSimTime = 0.;
			return true;
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real org_dt) {
			Real dt = org_dt * mpPlParams->getDeltaX(); 
			debMsg("spluginDumpUniversal","dt="<<dt<<" grid:"<<mGrid<<"; simTime="<<
					mSimTime<<", last="<<mLastAniOutTime<<", cnt="<<mAnimOutCounter<<" dump-all-frame="<<mDumpAllFrames ); 
			mSimTime += dt;

			while( this->doDump() ) 
			{
				UniversalHeader head;
				const nVec3i gridSize = mpPlParams->getGridInt("flags")->getSize();
					
				// get Grid
				Grid<int> *gI=NULL;
				Grid<Real> *gR=NULL;
				Grid<Vec3> *gV=NULL;
				if (mpPlParams->haveGridInt(mGrid)) {
					gI = mpPlParams->getGridInt(mGrid);
					head.bytesPerElement = sizeof(int);
					head.elementType = 0;
					head.elements = 1;					
				}
				else if (mpPlParams->haveGridReal(mGrid)) {
					gR = mpPlParams->getGridReal(mGrid);
					head.bytesPerElement = sizeof(Real);
					head.elementType = 1;
					head.elements = 1;
				}
				else if (mpPlParams->haveGridVec3(mGrid)) {
					gV = mpPlParams->getGridVec3(mGrid);
					head.bytesPerElement = sizeof(Real);
					head.elementType = 2;
					head.elements = 3;
				}
				else {
					errMsg("spluginDumpUniversal","Grid '"<<mGrid<<"' not found or gridtype unsupported");
					return false;
				}
				
				if (mFirstTime || mFilePerFrame) {
					// Open a file					
					std::string filename = mOutname;
					if (mFilePerFrame && !mNoExt) {
						char nrStr[5]; // nr conversion 
						snprintf(nrStr, 5, "%04d", mAnimOutCounter ); 
						filename += nrStr;
					}
					mGzf = gzopen( (filename+".gz").c_str(), "wb1");
					if (!mGzf) {
						errMsg("spluginDumpUniversal","Unable to open output '"<<filename<<"' ");
						return false; 
					}
					debMsg("spluginDumpUniversal","writing to file '" << filename << ".gz'" ); 
								
					// Write header				
					strncpy(head.id,"DDF1",4);
					head.dimX = gridSize.x;
					head.dimY = gridSize.y;
					head.dimZ = (gDim==2) ? 1:gridSize.z;
					head.frames = mMaxFrames;
					head.bytesPerFrame = head.dimX * head.dimY * head.dimZ * head.elements * head.bytesPerElement + sizeof(float);
				
					gzwrite(mGzf, &head, sizeof(UniversalHeader));
				}
				
				// write frame
				const int z0 = mpPlParams->getFluidSolver()->get2dKstart();
				const int z1 = mpPlParams->getFluidSolver()->get2dKend();
				for (int e=0; e< head.elements; e++)
				{
					// write grid of current element
					FOR_IJK(z0,z1, 0, head.dimY, 0, head.dimX) 
					{
						if (head.elementType == 0)
							gzwrite(mGzf, &(gI->getGlobal(i,j,k)), sizeof(int));
						else if (head.elementType == 1)
							gzwrite(mGzf, &(gR->getGlobal(i,j,k)), sizeof(Real));
						else if (head.elementType == 2) {
							Real data = gV->getGlobal(i,j,k)[e];
							if (mMAC) // get centered value on MAC grid
							{
								nVec3i pos(i,j,k);
								pos[e] += 1;
								if (pos[e] < gridSize[e]) {
									data += gV->getGlobal(pos[0],pos[1],pos[2])[e];
									data *= 0.5;
								}
							}
							gzwrite(mGzf, &data, sizeof(Real));
						}
					}
				}
				// write timestamp
				float tm = mSimTime;
				gzwrite(mGzf, &tm, sizeof(float));

				if (mFilePerFrame)
					gzclose(mGzf);

				this->increaseDumpCounter();
			}

			return true;
		};

		virtual void finish()
		{
			if (!mFilePerFrame && !mFirstTime) {
				gzclose(mGzf);
				debMsg("spluginDumpUniversal","file archive closed.");
			}
		};

	protected:
		// grid names to swap
		std::string mGrid;
		bool mFirstTime, mFilePerFrame, mMAC, mNoExt;
		gzFile mGzf;
};


class spluginLoadUniversal : public SolverPlugin {
	public:
		spluginLoadUniversal() : SolverPlugin() { };
		~spluginLoadUniversal() { };

		virtual bool parseParams(const ParamSet& params) {
			onlyFluid = params.FindOneInt("onlyfluid", 0 ) != 0;
			mGrid = params.FindOneString("grid", "");
			mSolver = params.FindOneString("solver", "" );
			mFile = params.FindOneString("file", "" );

			// interpolate, but keep #mInterpolBorder cells from each side unscaled
			// -1 means off, 0 means scale everything
			mInterpolBorder = params.FindOneInt("interpolate-with-border", 1);

			return true;
		};
		virtual bool initPlugin() { return true; };

		// perform step with given dt, return failure
		virtual bool performStep(Real dt) {
			debMsg("spluginLoadUniversal","file " + mFile);
			SolverParams* param = mpPlParams;
			if (!mSolver.empty()) param = ddfWorldFindSolver(mSolver)->getParams();
			Grid<int>* flags = param->getGridInt("flags");
			nVec3i gridSize = flags->getSize();
			if (gDim==2) gridSize[2] = 1;

			// load file into byte array
			UniversalHeader head;
			gzFile gzf = gzopen(mFile.c_str(),"r");
			if (gzf == NULL)
				errFatal("spluginLoadUniversal","can't open field file " << mFile, SIMWORLD_INITERROR);

			if (gzread(gzf, &head,sizeof(UniversalHeader)) != sizeof(UniversalHeader) || strncmp(head.id, "DDF1", 4)) 
				errFatal("spluginLoadUniversal","wrong filetype in " << mFile,SIMWORLD_INITERROR);

			if (gridSize[0] != head.dimX || gridSize[1] != head.dimY || gridSize[2] != head.dimZ) {
				debMsg("spluginLoadUniversal","file's gridsize (" << head.dimX << "," << head.dimY << "," << head.dimZ << ") != gridsize (" << gridSize[0] << "," << gridSize[1] << "," << gridSize[2] << ")");
				if(mInterpolBorder<=-1) 
					errFatal("spluginLoadUniversal","Gridsizes don't match, and interpolation is turned off" , SIMWORLD_GRIDERROR);
			}

			unsigned char* rawFrame = new unsigned char[head.bytesPerFrame];
			if (gzread(gzf, rawFrame, head.bytesPerFrame) != head.bytesPerFrame) 
				errFatal("spluginLoadUniversal","error loading frame data in " << mFile,SIMWORLD_INITERROR);
			
			gzclose(gzf);
		
			// distinguish data types, invoke grid loader
			if(param->haveGridInt(mGrid)) {
				Grid<int>* grid = param->getGridInt(mGrid);
				if (head.elementType != 0 || head.elements != 1) 
					errFatal("spluginLoadUniversal","Trying to load non-integer data into integer grid" , SIMWORLD_GRIDERROR);
				fsLoadData<int>(flags, grid, rawFrame, head, onlyFluid, mInterpolBorder);
			} else if(param->haveGridReal(mGrid)) {
				Grid<Real>* grid = param->getGridReal(mGrid);
				if (head.elementType != 1 || head.elements != 1) 
					errFatal("spluginLoadUniversal","Trying to load non-real data into real grid" , SIMWORLD_GRIDERROR);
				fsLoadData<Real>(flags, grid, rawFrame, head, onlyFluid, mInterpolBorder);
			} else if(param->haveGridVec3(mGrid)) {
				Grid<Vec3>* grid = param->getGridVec3(mGrid);
				if (head.elementType != 2 || head.elements != 3) 
					errFatal("spluginLoadUniversal","Trying to load non-vec3 data into vec3 grid" , SIMWORLD_GRIDERROR);				
				fsLoadData<Vec3>(flags, grid, rawFrame, head, onlyFluid, mInterpolBorder);
			} else { 
				errFatal("spluginLoadUniversal","Grid not found "<< mGrid , SIMWORLD_GRIDERROR);
			}
			return true;
		};

	protected:		
		std::string mGrid, mFile, mSolver;
		bool onlyFluid;
		int mInterpolBorder;
};

// dump density to disk as gzipped df3
class spluginDumpDf3 : public spluginDumperBase {
	public:
		spluginDumpDf3() : spluginDumperBase(),
   				mGrid("-unnamed1-"), prefix("") { };
		~spluginDumpDf3() { };

		virtual bool parseParams(const ParamSet& params) {
			//debMsg("spluginDumpDf3","parse");
			mGrid = params.FindOneString("gridname", mGrid );
			prefix = params.FindOneString("prefix", mGrid );
			mPbrt = params.FindOneInt("pbrt", 0 ) != 0;
			parseParamsDumperBase(params);
			return true;
		};
		virtual bool initPlugin() {
			mSimTime = 0.;
			return true;
		};

		// perform step with given dt, return failure
		virtual bool performStep(Real org_dt) {
			Real dt = org_dt * mpPlParams->getDeltaX(); 
			debMsg("spluginDumpDf3","dt="<<dt<<" grid:"<<mGrid<<"; simTime="<<
					mSimTime<<", last="<<mLastAniOutTime<<", cnt="<<mAnimOutCounter<<" " ); 

			mSimTime += dt;
			// crude approx if dt is larger than mTimestepAnim, 
			// might be better to reduce dt?
			while ( this->doDump() ) { 
				myTime_t dumptstart = getTime(); 
				mLastAniOutTime += mpPlParams->mTimestepAnim;
				FlagGrid*  pFlags  = mpPlParams->getGridInt( "flags" );
				nVec3i gridSize = pFlags->getSize();
				Grid<Real>* srcVal = mpPlParams->getGridReal(mGrid);

				int dimMax = VMAX(gridSize);

				const int twodKs = mpPlParams->getFluidSolver()->get2dKstart();
				const int twodKe = mpPlParams->getFluidSolver()->get2dKend();
				// TODO why are x,y swapped?
				//nVec3i min(twodKs, 0, 0);
				//nVec3i max(twodKe, gridSize[0], gridSize[1]);
				//nVec3i res = max-min;

				std::ostringstream voutfilename; 
				char nrStr[5]; // nr conversion 
				snprintf(nrStr, 5, "%04d", mAnimOutCounter ); 

				voutfilename << this->mOutname;
				if(prefix.length()>0) voutfilename << "_" << prefix;
			  	voutfilename << "_" << nrStr;

				std::cout << "FRAME DUMP " << voutfilename.str() << std::endl;

				debMsg("spluginDumpDf3","Writing '"<<voutfilename.str()<<"', resolution="<<gridSize<<" " );
				string filenameDf3 = voutfilename.str() + (mPbrt ? ".pbrt.gz" : ".df3.gz");
			
				gzFile file;
				file = gzopen(filenameDf3.c_str(), "wb1"); 
				if (file != NULL) {
					// dimensions
					const int byteSize = 2;
					const unsigned short int onx=srcVal->getSizeX(),ony=srcVal->getSizeY(),onz=srcVal->getSizeZ();
					/*unsigned short int nx,ny,nz; 
					nx = onx >> 8;
					ny = ony >> 8;
					nz = onz >> 8;
					nx += (onx << 8);
					ny += (ony << 8);
					nz += (onz << 8);
					nx = onx; ny = ony; nz = onz;
					gzwrite(file, (void*)&nx, sizeof(short));
					gzwrite(file, (void*)&ny, sizeof(short));
					gzwrite(file, (void*)&nz, sizeof(short));
					*/
					// get rid of swizzling
					nVec3i write_res = nVec3i(onx,ony,onz);					
					const int nitems = onx*ony*onz;
					const float mul = (float)( (1<<(8*byteSize))-1); 

					Real *buf = new Real[nitems];
					for (int k = 0; k < onz; k++) 
						for (int j = 0; j < ony; j++) 
							for (int i = 0; i < onx; i++) {
								float val = srcVal->getGlobal(i,j,k);
								CLAMP(val, 0.f,1.f);
								buf[k*(onx*ony)+j*onx+i] = val;
					}					
					if (mPbrt) {
						// text pbrt output
						std::ostringstream sbuf;
						Vec3 prop(onx,ony,onz);
						prop /= prop.max();
						sbuf << "Volume \"volumegrid\"" << endl;
						sbuf << " \"integer nx\" " << onx << endl;
						sbuf << " \"integer ny\" " << ony << endl;
						sbuf << " \"integer nz\" " << onz << endl;
						sbuf << " \"point p0\" [ 0 0 0 ]" << endl;
						sbuf << " \"point p1\" [ " << prop.x << " " << prop.y << " " << prop.z << " ]" << endl;
						sbuf << " \"float density\" [ " << endl;
						for (int i=0;i<nitems;++i)
							sbuf << buf[i] << " ";
						sbuf << " ] " << endl;
						string bs = sbuf.str();
						gzwrite(file, (void*)(bs.c_str()), bs.length());						
					} else {
						// binary DF3 output
						gzwrite(file, &write_res[0] , 3*sizeof(int));
						unsigned short int *dbuf = new unsigned short int[nitems];
						for (int i=0;i<nitems;++i)	
							dbuf[i] = (unsigned short int)(buf[i]*mul);
						gzwrite(file, (void*)dbuf, sizeof(unsigned short int)* nitems);
						delete[] dbuf;
					}
					gzclose(file);
					delete[] buf;
					// df3 gz written
				} else {
					errMsg("spluginDumpDf3","Unable to write to "<<filenameDf3 );
				}

				this->increaseDumpCounter();

				myTime_t dumptend = getTime(); 
				debMsg("spluginDumpDf3","took "<< getTimeString(dumptend-dumptstart)<<" for frame "<<mAnimOutCounter );
			}

			return true;
		};

	protected:
		// grid names to swap
		std::string mGrid, prefix;
		bool mPbrt;
};



//*****************************************************************************

SolverPlugin* MakeIoPlugin(std::string name) {
	if (name.compare( string("dump-df3") )==0) {
		return new spluginDumpDf3;
	} else if (name.compare( string("dump-universal") )==0) {
		return new spluginDumpUniversal;
	} else if (name.compare( string("load-universal") )==0) {
		return new spluginLoadUniversal;
	} else if (name.compare( string("print-timing-stats") )==0) {
		return new spluginPrintTimingStats;
	}
	return NULL;
}

} // end namespace DDF 

