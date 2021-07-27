/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * GUI
 *
 *****************************************************************************/

#include <list>

#include "globals.h"
#include "vectorbase.h"
#include "solverparams.h"
#include "paramset.h"

#include "fluidsolver.h"
#include "boundbox.h"
#include "vortexpart.h"

#if DDF_GLUTGUI==1
#if defined(__APPLE__)
// mac
#include <GLUT/glut.h>
#else 
#include <GL/glut.h>
#endif
#endif // DDF_GLUTGUI==1

// solver object from main
DDF::FluidSolver *gpGlutFsolver = NULL;

// velocity grid to display "string showvel" param
std::string gShowVelgridName = std::string("vel-curr");
// name of real grid to show
std::string gShowRealgridName = std::string("pressure");

// scale display of velocities with "float vel-scale" parameter for glut gui
static float gScaleVelocityDisplay = 1.;  // = g Vec3 Divider
// scale real values in grids for display, store one value per grid
std::map<std::string,float> gScaleRealDisplay;

// draw velocities at center instead of faces?
static bool gDrawCenteredVels = false;
static bool gVelocityMAC = true;
static bool gFinalizationDone = false;
// draw lines for grid (only for dims <= 128)
static bool gDrawLines = true;

// toggle flag bit string display for probing on/off
static bool gProbeShowFlags = false;

// show real values in empty cells?
const bool gHideEmptyValues = true;


/******************************************************************************/
namespace DDF { 



// set by parser externally, read here...
ParamSet gGlutGuiParams;

#if DDF_GLUTGUI==1

#if DDF_DIMENSION==2
Vec3 gViewTrafo(0., 0., -0.01);
#else
Vec3 gViewTrafo(0., 0., -0.5);
#endif
Vec3 gViewRot(0., 0., 0.);

bool gDrawVortexParticles = true;
//bool gDrawPhiRef = false;
bool gDrawVel = false;

// time stepping
bool gStep = false;
bool gPause = false;
// movement speeds
const float gMoveSpeed = 0.002;
const float gRotSpeed  = 0.2;

// region of interest
BboxVeci gGlutRoi;
int gRoiPlane = 2;
// minimal and maximal values for grid info in GUI
Real gRealgridMin = 0., gRealgridmMax = 0.;
Vec3 gVecgridMin = Vec3(0.), gVecgridMax = Vec3(0.);

bool gGlobalVPart = true;

// helper functions in guihelpers.cpp
void glutgSelectNextSolver();
void glutgSelectNextFloatArray();
void glutgSelectNextVec3Array();
void glutgSelectNextFloatArrayO();
void initFloatArrayDisplay();

void handle_display();

// mouse motion
static int g_mouseBut = -1;
static int g_specialKey = -1;
static int g_lastMouseX=-1, g_lastMouseY=-1;

/******************************************************************************/

void initGlutRoi(Real planeFac) { 
	nVec3i grsize = gpGlutFsolver->getGridFlags()->getSize();

	// roi checks are "inclusive" of boundaries, so reduce by one...
	grsize += nVec3i(-1);
	grsize[gRoiPlane] = (grsize[gRoiPlane] * planeFac);

	nVec3i grsizeMin(0);
	grsizeMin[gRoiPlane] = grsize[gRoiPlane];
	gGlutRoi = BboxVeci( grsizeMin, grsize ); 

}

inline void initGlutRoiRect(nVec3f &p0,nVec3f &p1,nVec3f &p2,nVec3f &p3, nVec3i off, Real mDX,nVec3f pOffset) {
	switch( gRoiPlane ) {
		case 2:
			p0 = nVec3f((off[0]+0)*mDX, (off[1]+0)*mDX, (off[2]+0)*mDX)+pOffset;
			p1 = nVec3f((off[0]+1)*mDX, (off[1]+0)*mDX, (off[2]+0)*mDX)+pOffset;
			p2 = nVec3f((off[0]+1)*mDX, (off[1]+1)*mDX, (off[2]+0)*mDX)+pOffset;
			p3 = nVec3f((off[0]+0)*mDX, (off[1]+1)*mDX, (off[2]+0)*mDX)+pOffset;
			break;
		case 1:
			p0 = nVec3f((off[0]+0)*mDX, (off[1])*mDX, (off[2]+0)*mDX)+pOffset;
			p1 = nVec3f((off[0]+1)*mDX, (off[1])*mDX, (off[2]+0)*mDX)+pOffset;
			p2 = nVec3f((off[0]+1)*mDX, (off[1])*mDX, (off[2]+1)*mDX)+pOffset;
			p3 = nVec3f((off[0]+0)*mDX, (off[1])*mDX, (off[2]+1)*mDX)+pOffset;
			break;
		case 0:
			p0 = nVec3f((off[0])*mDX, (off[1]+0)*mDX, (off[2]+0)*mDX)+pOffset;
			p1 = nVec3f((off[0])*mDX, (off[1]+1)*mDX, (off[2]+0)*mDX)+pOffset;
			p2 = nVec3f((off[0])*mDX, (off[1]+1)*mDX, (off[2]+1)*mDX)+pOffset;
			p3 = nVec3f((off[0])*mDX, (off[1]+0)*mDX, (off[2]+1)*mDX)+pOffset;
			break;
	}
}

void drawVortexParticles(VorticitySystem::VList& parts)
{
	if (gpGlutFsolver->getVorticitySys()->getParticles().size() == 0) return;
	
	Real dx = gpGlutFsolver->getParams()->getDeltaX();
	int a = gGlutRoi.getStart()[gRoiPlane];
	int b = gGlutRoi.getEnd()[gRoiPlane];
	
	for (VorticitySystem::VList::const_iterator it = parts.begin(); it != parts.end(); ++it)
	{
		if (((*it)->getFlags() & VortexParticle::FInertial) == 0)
			glColor3f(1., 1., 0.);
		else
			glColor3f(0.,1.,1.);
	
		Vec3 pos = (*it)->getPos() +.5;
		glPointSize((*it)->getRadius()*1.5);
			
		if (gDim == 2)
			pos[2] = 0.0;
		else {
			//pos[2] -= 0.5; // z offset to account for centered phi value
			if (!gGlobalVPart && (pos[2]<a-0.5	||pos[2]>b+0.5))
				continue;
		}
		pos *= dx;
		glBegin(GL_POINTS);glVertex3f(pos[0],pos[1],pos[2]);glEnd();	
	}
	
#if DEBUG_VORTEXPATH==1
	glEnable(GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glBegin(GL_LINES);
	for (std::list<LineElement>::const_iterator it = gVortexPath.begin(); it != gVortexPath.end();++it)
	{
		Vec3 p1 = it->p0 + .5, p2 = it->p1 + .5;
		if (gDim == 2)
			p2[2] = p1[2] = 0.0;
		else {
			int a2 = gGlutRoi.getStart()[gRoiPlane];
			int b2 = gGlutRoi.getEnd()[gRoiPlane];
			if (!gGlobalVPart && (p1[2]<a2-0.5	||p1[2]>b2+0.5 || p2[2]<a2-0.5	||p2[2]>b2+0.5))
				continue;
		}
		const Real intensity = 0.7;
		glColor4f(intensity*it->r,intensity*it->g,intensity*it->b,it->a);
		glVertex3f(dx*p1[0],dx*p1[1],dx*p1[2]);
		glVertex3f(dx*p2[0],dx*p2[1],dx*p2[2]);		
	}
	glEnd();	
	glDisable(GL_BLEND);	
#endif
	
}

// draw grid obstacles and compute min/max pressure
class fsDrawGridLines : public GridOpBase {
	public:
		fsDrawGridLines(FlagGrid *flags, BboxVeci roi) : 
			GridOpBase() { 
			mpFlags = flags;
			mRoi = roi;
			mDX = 1. / (Real)(mpFlags->getMaxSize());
			pOffset = nVec3f(0., 0., -.00001);
			glBegin(GL_LINES);
			applyOperatorToGridsSimple(this);
			glEnd();
		};
		~fsDrawGridLines() { }; 
		void resetVariables() { };
		void reduce(fsDrawGridLines &op) { };

		void buildCallList() {
			setFlags(mpFlags);
		};
		inline void operator() (int i, int j, int k) { 
			nVec3i off = nVec3i(i,j,k);
			if (!mRoi.contains(off)) return; // if (gDim==3 && off[2]!=mRoi) return;
			if (gDim==2) off[2] = 0;
			nVec3f p0,p1,p2,p3; initGlutRoiRect(p0,p1,p2,p3,off,mDX,pOffset);
			glColor3f(0.2,0.2,0.2);
			glVertex3fv(&p0[0]); glVertex3fv(&p1[0]);
			glVertex3fv(&p1[0]); glVertex3fv(&p2[0]);
			glVertex3fv(&p2[0]); glVertex3fv(&p3[0]);
			glVertex3fv(&p3[0]); glVertex3fv(&p0[0]);
		};
	protected:
		nVec3f pOffset; // drawing offset
		BboxVeci mRoi; // display a single slice in 3d
		Real mDX;
}; // fsDrawGridLines */

// draw grid obstacles and compute min/max pressure
// which scalar field to show in fsDrawObstacles drawTypes:
// 0 = eg pressure (w min/max)
// 1 = levelset (if active)
// 2 = error function
class fsDrawObstacles : public GridOpBase {
	public:
		fsDrawObstacles(FlagGrid *flags, 
				Grid<Real> *srcGrid, Grid<Real> *err, 
				Real minp, Real maxp, BboxVeci roi, int drawtype) : 
			GridOpBase(), mpSrc(srcGrid),mpErr(err), mDrawType(0) { 

			mpFlags = flags;
			mRoi = roi;
			mDX = 1. / (Real)(mpSrc->getMaxSize());
			mMinP = minp; 
			mMaxP = maxp;
			mFirstPInit = false;
			pOffset = nVec3f(0., 0., -.00001);
			mDrawType = drawtype;

			//mProjMult = gGlutProject / gScaleRealDisplay  * mDX;
			mValScale = gScaleRealDisplay[gShowRealgridName];
			mProjMult = 0;//gGlutProject / mValScale * mDX;

			SolverParams* prms = gpGlutFsolver->getParams();

			mHideEmpty = gHideEmptyValues;
			if(mpSrc->getDisplayFlags() & 2) mHideEmpty = false;

			glBegin(GL_QUADS);
			applyOperatorToGridsSimple(this);
			glEnd();
			//minp = mMinP; maxp = mMaxP;
		};
		~fsDrawObstacles() { }; 
		void resetVariables() { };
		void reduce(fsDrawObstacles &op) { };

		void buildCallList() {
			gaSrc.gridAccInit(mpSrc, AM_READ, gaCalls); 
			if (mpErr) {
				gaErr.gridAccInit(mpErr, AM_READ, gaCalls); 
			}
			setFlags(mpFlags);
		};
		inline void operator() (int i, int j, int k) { 
			const int currFlag = getFlagAcc()(i,j,k);
			nVec3i gridPos = nVec3i(i,j,k);

			if (!mRoi.contains(gridPos)) return; 
			if (gDim==2) gridPos[2] = 0;
			nVec3f p0,p1,p2,p3; 
			initGlutRoiRect(p0,p1,p2,p3,gridPos,mDX,pOffset);
						
			if (mDrawType==0) { 
				if (mHideEmpty && fgIsObstacle(currFlag)) {
					// obstacle
					glColor3f(0.15,0.15,0.15);
				} else if (1 && fgIsOutflow(currFlag)) {
					glColor3f(0.3,0.0,0.0);
				} else if (mHideEmpty && fgIsEmpty(currFlag)) {
					// empty cells
					glColor3f(0.,0.2,0.); // */
				} else {
					Real ps = 0.;

					// override with projection
					{
						// default, show field value
						const Real val = gaSrc(i,j,k);
						ps = val / mValScale ;
					}

					if (fgIsInflow(currFlag)) {
						glColor3f(ps, ps, 0.3);
					} else {
						if (ps>0)
							glColor3f(ps,0.,0.);
						else
							glColor3f(0.,0.,fabs(ps));
					}
				}
			} else if (mDrawType==1) { 
				// level set display
				if (gDim==2) gridPos[2] = gPatchSize/2;

				Real v = gaSrc(i,j,k); 

				const Real vScale = 1./5. * mValScale;
				v *= vScale;
				CLAMP(v, (Real)-1., (Real)1.); 

				if (v>=0.) {
					glColor3f(v,0.,0.5);
				} else {
					if (v<=-1000.) {
						glColor3f(0., 0., 0.);  // invalid values
					} else {
						glColor3f(0.5, 1.+v, 0.);
					}
				}
			} else if (mDrawType == 2) {
				/* a to b fades from blue to green, b to c from green to red */
				/*const Real a = 0.01;
				const Real b = 1.0;
				const Real c = 20.0;*/
				const Real a = 0.0001;
				const Real b = 0.01;
				const Real c = 1.0;
				if (gDim==2) gridPos[2] = gPatchSize/2;
				Real err = fabs(gaErr(i,j,k));
				
				if (err > b) {
					Real val = (err-b)/(c-b);
					glColor3f(val,1.0f-val,0.0f);
				} else if (err > a) {
					Real val = (err-a)/(b-a);
					glColor3f(0.0f,val,1-val);
				} else
					glColor3f(0.0f, 0.0f, 0.0f);

			} else {
				return;
			}

			glVertex3fv(&p0[0]);
			glVertex3fv(&p1[0]);
			glVertex3fv(&p2[0]);
			glVertex3fv(&p3[0]); 
		};
	protected:
		nVec3f pOffset; // drawing offset
		BboxVeci mRoi; 
		Real mDX;
		Real mMinP, mMaxP;
		bool mFirstPInit;
		Grid<Real> *mpSrc;
		Grid<Real> *mpErr;
		GridAccessor<Real,0> gaSrc;
		GridAccessor<Real,0> gaErr;

		Grid<Real> *mpPhi;
		Grid<Real> *mpPhiRef;

		//! what to display 0=grid, 1=levelset
		int mDrawType;
		int mHideEmpty;
		// multiplier for projected displays, adapted to gridsize
		Real mProjMult, mValScale;
}; // fsDrawObstacles */

// draw grid obstacles and compute min/max pressure
class fsDrawCellIndex : public GridOpBase {
	public:
		fsDrawCellIndex(FlagGrid *flags, BboxVeci roi) : 
				GridOpBase()
		{ 
			mpFlags = flags;
			mDX = 1. / (Real)(flags->getMaxSize());
			pOffset = nVec3f(0., 0., -.00001);
			mRoi = roi;

			glBegin(GL_QUADS);
			applyOperatorToGridsSimple(this);
			glEnd();
		};
		~fsDrawCellIndex() { }; 
		void resetVariables() { };
		void reduce(fsDrawCellIndex &op) { };

		void buildCallList() {
			setFlags(mpFlags);
		};
		inline void operator() (int i, int j, int k) { 
			const int currFlag = getFlagAcc()(i,j,k);
			nVec3i gridPos = /*getFlagAcc().getAccPatch()->mOffset + */nVec3i(i,j,k);
			if (!mRoi.contains(gridPos)) return; 

			if (gDim==2) gridPos[2] = 0;
			nVec3f p0,p1,p2,p3; 
			initGlutRoiRect(p0,p1,p2,p3,gridPos,mDX,pOffset);
			
			if (gDim==2) gridPos[2] = gPatchSize/2;

			Vec3 col = Vec3(
					gridPos[0]/(Real)mpFlags->getSizeX(), 
					gridPos[1]/(Real)mpFlags->getSizeY(), 
					gridPos[2]/(Real)mpFlags->getSizeZ() );
			glColor3f(col[0],col[1],col[2]);
			//debMsg("at "," "<<PRINT_IJK<<" "<<col); // DEBUG

			glVertex3fv(&p0[0]);
			glVertex3fv(&p1[0]);
			glVertex3fv(&p2[0]);
			glVertex3fv(&p3[0]); 
		};
	protected:
		nVec3f pOffset; // drawing offset
		BboxVeci mRoi; 
		Real mDX;
}; // fsDrawLevelset */



// draw centered grid vels
class fsDrawVels : public GridOpBase {
	public:
		fsDrawVels(FlagGrid *flags, Grid<Vec3> *vel, BboxVeci roi) : 
			GridOpBase(), mpVel(vel) { 
			mpFlags = flags;
			mDX = 1. / (Real)(mpVel->getMaxSize());
			mRoi = roi;

			// switch off mac display for this grid, if flags|4 is set
			//if( mpVel->getDisplayFlags() & 4) mMacValues = false;

			glBegin(GL_LINES);
			applyOperatorToGridsSimple(this);
			glEnd();
		};
		~fsDrawVels() { }; 
		void resetVariables() { };
		void reduce(fsDrawVels &op) { };

		void buildCallList() {
			gaVel.gridAccInit(mpVel, AM_READ, gaCalls); 
			setFlags(mpFlags);
		};
		inline void operator() (int i, int j, int k) { 
			const int currFlag = getFlagAcc()(i,j,k);
			const bool dxScale = true;
			//if (!fgIsFluid(currFlag)) return;
			//if (fgIsObstacle(currFlag)) return;

			nVec3i off = nVec3i(i,j,k);
			if (!mRoi.contains(off)) return;
			if (gDim==2) off[2] = 0;

			nVec3f p = nVec3f(off[0]*mDX, off[1]*mDX, off[2]*mDX);

			if (gDrawCenteredVels) { // centered
				nVec3f vOrg = vec2F( gaVel(i,j,k) );
				nVec3f v = vOrg;
				if(mpVel->checkIndexValid(i+1,j+1,k+1) && gVelocityMAC) {
					v[0] += gaVel(i+1,j,k)[0];
					v[1] += gaVel(i,j+1,k)[1];
					v[2] += gaVel(i,j,k+1)[2];
					v *= 0.5;
				}
				if(dxScale) v *= mDX;
				v *= gScaleVelocityDisplay;

				p[0] += mDX*0.5;
				p[1] += mDX*0.5;
				if (gDim==3) p[2] += mDX*0.5; // offset in z dir for display only

				glColor3f(0.,1.,0.);
				glVertex3fv(&p[0]);
				p += v;
				glColor3f(1.,1.,1.);
				glVertex3fv(&p[0]);
			} else { // side
				// dont use averaged vels!
				nVec3f v = vec2F( gaVel(i,j,k) ); // * 0.33;
				if(dxScale) v *= mDX;
				v *= gScaleVelocityDisplay;

				p = nVec3f(off[0]*mDX, off[1]*mDX, off[2]*mDX);
				p[0] += 0.5*mDX;
				glColor3f(0.,0.,1.);
				glVertex3fv(&p[0]);
				p[1] += v[1];
				glColor3f(1.,1.,1.);
				glVertex3fv(&p[0]);

				p = nVec3f(off[0]*mDX, off[1]*mDX, off[2]*mDX);
				p[1] += 0.5*mDX;
				glColor3f(1.,0.,0.);
				glVertex3fv(&p[0]);
				p[0] += v[0];
				glColor3f(1.,1.,1.);
				glVertex3fv(&p[0]);

				if (gDim==3) {
				p = nVec3f(off[0]*mDX, off[1]*mDX, off[2]*mDX);
				p[2] += 0.5*mDX;
				glColor3f(1.,0.,0.);
				glVertex3fv(&p[0]);
				p[2] += v[2];
				glColor3f(1.,1.,1.);
				glVertex3fv(&p[0]); } // 3dim
			}
		};
	protected:
		Real mDX;
		BboxVeci mRoi; 
		Grid<Vec3> *mpVel;
		GridAccessor<Vec3,1> gaVel;
		//bool mMacValues;
}; // fsDrawVels */


/******************************************************************************/

extern bool advanceAllSolvers(); // from main
extern void finalizeAllSolvers(); // from main
extern void deleteAllSolvers(); // from main

void runSimulation () {
	if (gPause) return;
	if (gFinalizationDone) return;

	bool quit = advanceAllSolvers();
	if (quit) {
		gPause=true;

		if (!gFinalizationDone) {
			finalizeAllSolvers();
			gFinalizationDone = true;
		}

		// immediatly quit for qutorun
		if(getenv("DDF_GUIAUTORUN")) {
			if(atoi(getenv("DDF_GUIAUTORUN"))>=1) {
				exit(0);
			}
		}
	} // quit

	// debug info per timestep
	gRealgridMin = gRealgridmMax = 0.;

	// compute only once each time step
	if (1) {
		Grid<Real> *realGrid = NULL; 
		realGrid = gpGlutFsolver->getParams()->getGridReal(gShowRealgridName); // DEBrealGrid
		goFindMinMax<Real> gomm = goFindMinMax<Real>( gpGlutFsolver->getGridFlags(), realGrid );
		gRealgridMin = gomm.mMinVal;
		gRealgridmMax = gomm.mMaxVal;
		debMsg("MINMAX-REAL","min="<<gRealgridMin<<" at "<<gomm.mMinPos<<" max="<<gRealgridmMax<<" at "<<gomm.mMaxPos);

		goFindMinMax<Vec3> gomv = goFindMinMax<Vec3>( gpGlutFsolver->getGridFlags(), 
				gpGlutFsolver->getParams()->getGridVec3(gShowVelgridName) );
		gVecgridMin = gomv.mMinVal;
		gVecgridMax = gomv.mMaxVal;
		debMsg("MINMAXV","min="<<gomv.mMinVal<<" at "<<gomv.mMinPos<<" max="<<gomv.mMaxVal<<" at "<<gomv.mMaxPos);

	} // end minmax debug info

	glutPostRedisplay();
	if (gStep) {
		gPause = true;
		gStep = false;
	}

	// pause sim, if an error occurred
	if(!isSimworldOk()) {
		gPause = true;
	}
}

// debug output for patch solver
void drawPatchGrid()
{
	Grid<Real> *perr = NULL; 

	int drawType = 0;
	
	Grid<Real> *realGrid = NULL; 
	realGrid = gpGlutFsolver->getParams()->getGridReal(gShowRealgridName); // DEBUG
	
	// hide everything else for some iso  disp settings
	bool drawGridInfo = true;
	
	if(drawGridInfo) {
		fsDrawObstacles( 
				gpGlutFsolver->getGridFlags(), 
				realGrid, perr,
				gRealgridMin,gRealgridmMax, gGlutRoi, drawType); 
	}

	glLineWidth(1.0);
	if(gDrawLines && drawGridInfo) {
		fsDrawGridLines( gpGlutFsolver->getGridFlags(), gGlutRoi ); 
	}

	if (gDrawVel && drawGridInfo) {
		Grid<Vec3> *vels = gpGlutFsolver->getParams()->getGridVec3(gShowVelgridName);
		fsDrawVels( gpGlutFsolver->getParams()->getGridInt("flags"), vels, gGlutRoi );
	}

	if (gDrawVortexParticles && drawGridInfo && gpGlutFsolver->getVorticitySys() != NULL) {
		drawVortexParticles(gpGlutFsolver->getVorticitySys()->getParticles());
	}
}


/******************************************************************************/
void handle_display (void) {
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE); 

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	//gluLookAt(0.5, 0.5, 1.0, 0.5, 0.5, 0.0, 0.0, 1.0, 0.0);
	//        from           to             up
	gluLookAt(0., 0., 0.5,   0., 0., -0.5,  0., 1., 0.);
	glPushMatrix();

	glTranslatef( gViewTrafo[0], gViewTrafo[1], gViewTrafo[2] );
	glRotatef(gViewRot[0],  1.,0.,0.);
	glRotatef(gViewRot[1],  0.,1.,0.);
	// center on .5 .5 .5
	glTranslatef( -0.5,-0.5,-0.5 );

	drawPatchGrid();

	glPopMatrix();
	glFlush();
	glutSwapBuffers();
}

void print_stats() { 
	std::ostringstream outstr;
	int frame = gpGlutFsolver->getFrameNum();
	int aniframe = gpGlutFsolver->getAniFrameNum();
	outstr << "Solver: "<<gpGlutFsolver->getName()<<", t="<< gpGlutFsolver->getSimTime() << " [frame " << frame << ", ani " << aniframe  <<  "]\n";

	// real grid stats
	//Real absmax=gRealgridmMax;
	//if (fabs(gRealgridMin)>fabs(absmax)) absmax=gRealgridMin; 
	outstr << "Real-grid: "<<gShowRealgridName<< " [range " << gScaleRealDisplay[gShowRealgridName]   ;
	if(gRealgridMin != 0.) outstr << " min " << gRealgridMin <<  "," ;
	outstr << " max " << gRealgridmMax <<  "]\n";
	int numP = (gpGlutFsolver->getVorticitySys()==NULL) ? 0 : gpGlutFsolver->getVorticitySys()->getParticles().size();
	if (numP > 0) outstr << numP << " vortex parts active\n";

	// vec grid stats
	if(gDrawVel) {
		outstr<<"Vec3-grid: "<<gShowVelgridName<<"";
		if(gDrawCenteredVels) outstr<<", centered ";
		else outstr<<", faces ";

		outstr << " range=" << gScaleVelocityDisplay ;
		if(norm( gVecgridMin ) > 0.) outstr << " min=" << gVecgridMin <<  ",";
		outstr << " max=" << gVecgridMax ;
		if (!gVelocityMAC) outstr << " noMAC" ;
		outstr << "\n";
	}

	debMsg("GRID", outstr.str());
}

void handle_reshape (int w, int h) {
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, (GLfloat) w/(GLfloat) h, 0.05, 20.0);

	glMatrixMode(GL_MODELVIEW);
}

void handle_keyboard(unsigned char key, int x, int y) {
	bool roiMsgs = false;
	switch (key) {
		case 'p': 
			// run simulation
			gPause = !gPause;
			break;
		case 'l': 
			// single step
			gPause = false; 
			gStep = true;
			break;

			// roi z plane
		case '-': 
		case '_': 
			gGlutRoi.getStart()[gRoiPlane]--; gGlutRoi.getEnd()[gRoiPlane]--; 
			if (roiMsgs) debMsg("key","- plane = " << gGlutRoi.getStart()[gRoiPlane]);
			if (roiMsgs) debMsg("changeGlutRoi","- roi = "<<gGlutRoi.toString()<<" plane="<<gRoiPlane ); // debug
			break;
		case '+': 
		case '=': 
			gGlutRoi.getStart()[gRoiPlane]++; gGlutRoi.getEnd()[gRoiPlane]++; 
			if (roiMsgs) debMsg("key","+ plane = " << gGlutRoi.getStart()[gRoiPlane]);
			if (roiMsgs) debMsg("changeGlutRoi","+ roi = "<<gGlutRoi.toString()<<" plane="<<gRoiPlane ); // debug
			break;
		case '*': 
			gRoiPlane++; gRoiPlane=gRoiPlane%3; 
			initGlutRoi( 0.5 ); 
			if (roiMsgs) debMsg("changeGlutRoi","* roi = "<<gGlutRoi.toString()<<" plane="<<gRoiPlane ); // debug
			break;
		case 'o':
			gDrawVortexParticles = !gDrawVortexParticles;
			break;
		case 'v':
			// toggle velocity display
			gDrawVel = !gDrawVel;
			gVelocityMAC = !(gpGlutFsolver->getParams()->getGridVec3(gShowVelgridName)->getDisplayFlags() & 4);
			break;
		case 'V':
			// toggle centered velocity display
			gDrawCenteredVels = !gDrawCenteredVels;
			break;
		case 'g':
			//gGlutProject *=-1;
			break;
		case 'G':
			gGlobalVPart =!gGlobalVPart;
			break;
		case '`': {
				// remember grid selection factor
				nVec3i grsize = gpGlutFsolver->getGridFlags()->getSize();
				Real planeFac = ( ((Real)gGlutRoi.getStart()[gRoiPlane]+1.0) / (Real)grsize[gRoiPlane]) ;
				//grsize[gRoiPlane] = (grsize[gRoiPlane] * planeFac);

				// show next solver, cycle through all available
				glutgSelectNextSolver();

				// reset struct that have to be initialized anew for solver
				initGlutRoi(planeFac); 
				debMsg("selectNextSolver","Roi factor="<<planeFac<<", plane = "<< gGlutRoi.getStart()[gRoiPlane] );
				print_stats();
			 } break;

		case '[':
			gScaleRealDisplay[gShowRealgridName]  *= 2.;
			print_stats();
			break;
		case ']':
			gScaleRealDisplay[gShowRealgridName]  /= 2.;
			print_stats();
			break;
		case '{':
			gScaleVelocityDisplay /= 2.;
			print_stats();
			break;
		case '}':
			gScaleVelocityDisplay *= 2.;
			print_stats();
			break;

		case 'z': 
			// show next float grid
			glutgSelectNextFloatArray();
			initFloatArrayDisplay();
			print_stats();
			break;
		case 'x': 
			// show next vec3 grid
			glutgSelectNextVec3Array();
			gVelocityMAC = !(gpGlutFsolver->getParams()->getGridVec3(gShowVelgridName)->getDisplayFlags() & 4);
			print_stats();
			break;
			// cam movement
		case 'a': gViewTrafo[0] += gMoveSpeed* 15; break;
		case 'd': gViewTrafo[0] += gMoveSpeed*-15; break;
		case 's': gViewTrafo[1] += gMoveSpeed* 15; break;
		case 'w': gViewTrafo[1] += gMoveSpeed*-15; break;
		case 'q': gViewTrafo[2] += gMoveSpeed* 15; break;
		case 'e': gViewTrafo[2] += gMoveSpeed*-15; break;

		case 27: // ESC
			// quit
			//cleanSimulation();
			exit(0);
			break;
	}

	glutPostRedisplay();
}


void handle_mousePress(int button, int state, int x, int y) {
	g_specialKey = glutGetModifiers();
	// if both a mouse button, and the ALT key, are pressed  then
	if ( (state == GLUT_DOWN) ) {

	      if (button == GLUT_LEFT_BUTTON) {
		      g_mouseBut = GLUT_LEFT_BUTTON;
	      } else if (button == GLUT_MIDDLE_BUTTON) {
		      g_mouseBut = GLUT_MIDDLE_BUTTON;
	      } else {
		      g_mouseBut = GLUT_RIGHT_BUTTON;
	      }
	} else {
		g_mouseBut = -1;
		g_lastMouseX=g_lastMouseY=-1;
		// gProbeGridString = std::string(""); // reset display
	}
}
void handle_mouseMotion(int x,int y) {

	if (g_lastMouseX==-1 && g_lastMouseY==-1) {
		g_lastMouseX = x;
		g_lastMouseY = y;
		return;
	}

	if (g_mouseBut==GLUT_RIGHT_BUTTON) {
		gViewTrafo[0] += (g_lastMouseX-x) * -gMoveSpeed;
		gViewTrafo[1] += (g_lastMouseY-y) *  gMoveSpeed;
		glutPostRedisplay();
	}
	if (g_mouseBut==GLUT_MIDDLE_BUTTON) {
		gViewTrafo[2] += (g_lastMouseX-x) * -gMoveSpeed;
		gViewTrafo[2] += (g_lastMouseY-y) *  gMoveSpeed;
		glutPostRedisplay();
	}
	if (g_mouseBut==GLUT_LEFT_BUTTON) {
		gViewRot[1] += (g_lastMouseX-x) * -gRotSpeed;
		gViewRot[0] += (g_lastMouseY-y) * -gRotSpeed;
		glutPostRedisplay();
	}
	g_lastMouseX = x;
	g_lastMouseY = y;
	//debMsg("mouseMotion"," at "<<PRINT_VEC(x,y,0.)<<", last "<<PRINT_VEC(g_lastMouseX,g_lastMouseY,0.)<<" but="<<g_mouseBut);
}



/******************************************************************************/
// glut gui test
int showGlutGui(const char *argv_title)
{
	int argc=1;
	char argvb[64]; strcpy(argvb,"ddfCmd");
	char* argv = argvb;
	int sx = 800, sy = 800;

	debMsg("\nglutGui","Starting...");

	// init parser parameters
	string solverName = gGlutGuiParams.FindOneString("solvername","default");
	gScaleVelocityDisplay = gGlutGuiParams.FindOneFloat("vel-scale", gScaleVelocityDisplay );
	gpGlutFsolver = ddfWorldFindSolver(solverName);
	if (!gpGlutFsolver) { return 1; }

	//if(gpGlutFsolver->getGridFlags()->getMaxSize() > 128) gDrawLines = false;
	gDrawLines = gGlutGuiParams.FindOneInt("draw-lines", gDrawLines ) > 0;

	gShowRealgridName = gGlutGuiParams.FindOneString("showgrid",gShowRealgridName);
	gShowVelgridName = gGlutGuiParams.FindOneString("showvel",gShowVelgridName);
	
	gGlutGuiParams.ReportUnused();

	if(gShowRealgridName.length()<=0 || (!gpGlutFsolver->getParams()->haveGridReal(gShowRealgridName)) ) {
		debMsg("glutGui","Invalid gShowRealgridName init '"<<gShowRealgridName<<"', re-init...");
		glutgSelectNextFloatArray();
	}
	initFloatArrayDisplay();

	// center view to grid 
	Grid<Real> *realGrid = gpGlutFsolver->getParams()->getGridReal(gShowRealgridName);
	nVec3i size = realGrid->getSize();
	double dx = 1./(double)realGrid->getMaxSize();
	gViewTrafo[0] = 0.5 - 0.5* (double)size[0] *dx;
	gViewTrafo[1] = 0.5 - 0.5* (double)size[1] *dx;

	// init GLUT
	glutInit(&argc, &argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
	glutInitWindowSize(sx, sy);
	glutInitWindowPosition(10, 10);

	std::ostringstream title;
	title << string(argv_title);
	glutCreateWindow( title.str().c_str() );

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	
	// initialize OpenGL
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClearDepth(1.0);
	//glShadeModel(GL_SMOOTH);
	glDepthFunc(GL_LESS);

	// initialize Simulation drawing
	initGlutRoi( 0.5 );
	
	// setup callbacks
	glutDisplayFunc(handle_display);
	glutReshapeFunc(handle_reshape);
	glutKeyboardFunc(handle_keyboard);
	glutIdleFunc(runSimulation);
	glutMouseFunc(handle_mousePress);
	glutMotionFunc(handle_mouseMotion);
	glutMainLoop();

	return 0;
}


#else // DDF_GLUTGUI==1
/******************************************************************************/
int showGlutGui(const char* str)
{
	// dummy
	return 0;
}
#endif // DDF_GLUTGUI==1

}; // namespace DDF

