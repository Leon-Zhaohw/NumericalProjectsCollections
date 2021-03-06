/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Grid representation
 *
 ******************************************************************************/

#include "grid.h"
#include "levelset.h"
#include "kernel.h"
#include "commonkernels.h"
#include <limits>
#include <sstream>
#include <cstring>
#include "fileio.h"

using namespace std;
namespace Manta {

//******************************************************************************
// GridBase members

GridBase::GridBase (FluidSolver* parent) 
	: PbClass(parent), mType(TypeNone)
{
	checkParent();
	m3D = getParent()->is3D();
}

//******************************************************************************
// Grid<T> members

// helpers to set type
template<class T> inline GridBase::GridType typeList() { return GridBase::TypeNone; }
template<> inline GridBase::GridType typeList<Real>()  { return GridBase::TypeReal; }
template<> inline GridBase::GridType typeList<int>()   { return GridBase::TypeInt;  }
template<> inline GridBase::GridType typeList<Vec3>()  { return GridBase::TypeVec3; }

template<class T>
Grid<T>::Grid(FluidSolver* parent, bool show)
	: GridBase(parent)
{
	mType = typeList<T>();
	mSize = parent->getGridSize();
	mData = parent->getGridPointer<T>();
	
	mStrideZ = parent->is2D() ? 0 : (mSize.x * mSize.y);
	mDx = 1.0 / mSize.max();
	clear();
	setHidden(!show);
}

template<class T>
Grid<T>::Grid(const Grid<T>& a) : GridBase(a.getParent()) {
	mSize = a.mSize;
	mType = a.mType;
	mStrideZ = a.mStrideZ;
	mDx = a.mDx;
	FluidSolver *gp = a.getParent();
	mData = gp->getGridPointer<T>();
	memcpy(mData, a.mData, sizeof(T) * a.mSize.x * a.mSize.y * a.mSize.z);
}

template<class T>
Grid<T>::~Grid() {
	mParent->freeGridPointer<T>(mData);
}

template<class T>
void Grid<T>::clear() {
	memset(mData, 0, sizeof(T) * mSize.x * mSize.y * mSize.z);    
}

template<class T>
void Grid<T>::swap(Grid<T>& other) {
	if (other.getSizeX() != getSizeX() || other.getSizeY() != getSizeY() || other.getSizeZ() != getSizeZ())
		errMsg("Grid::swap(): Grid dimensions mismatch.");
	
	T* dswap = other.mData;
	other.mData = mData;
	mData = dswap;
}

template<class T>
void Grid<T>::load(string name) {
	if (name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if (ext == ".raw")
		readGridRaw(name, this);
	else if (ext == ".uni")
		readGridUni(name, this);
	else if (ext == ".vol")
		readGridVol(name, this);
	else
		errMsg("file '" + name +"' filetype not supported");
}

template<class T>
void Grid<T>::save(string name) {
	if (name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if (ext == ".raw")
		writeGridRaw(name, this);
	else if (ext == ".uni")
		writeGridUni(name, this);
	else if (ext == ".vol")
		writeGridVol(name, this);
	else if (ext == ".txt")
		writeGridTxt(name, this);
	else
		errMsg("file '" + name +"' filetype not supported");
}

//******************************************************************************
// Grid<T> operators

//! Kernel: Compute min value of Real grid
KERNEL(idx, reduce=min) returns(Real minVal=std::numeric_limits<Real>::max())
Real CompMinReal(const Grid<Real>& val) {
	if (val[idx] < minVal)
		minVal = val[idx];
}

//! Kernel: Compute max value of Real grid
KERNEL(idx, reduce=max) returns(Real maxVal=-std::numeric_limits<Real>::max())
Real CompMaxReal(const Grid<Real>& val) {
	if (val[idx] > maxVal)
		maxVal = val[idx];
}

//! Kernel: Compute min value of int grid
KERNEL(idx, reduce=min) returns(int minVal=std::numeric_limits<int>::max())
int CompMinInt(const Grid<int>& val) {
	if (val[idx] < minVal)
		minVal = val[idx];
}

//! Kernel: Compute max value of int grid
KERNEL(idx, reduce=max) returns(int maxVal=-std::numeric_limits<int>::max())
int CompMaxInt(const Grid<int>& val) {
	if (val[idx] > maxVal)
		maxVal = val[idx];
}

//! Kernel: Compute min norm of vec grid
KERNEL(idx, reduce=min) returns(Real minVal=std::numeric_limits<Real>::max())
Real CompMinVec(const Grid<Vec3>& val) {
	const Real s = normSquare(val[idx]);
	if (s < minVal)
		minVal = s;
}

//! Kernel: Compute max norm of vec grid
KERNEL(idx, reduce=max) returns(Real maxVal=-std::numeric_limits<Real>::max())
Real CompMaxVec(const Grid<Vec3>& val) {
	const Real s = normSquare(val[idx]);
	if (s > maxVal)
		maxVal = s;
}


template<class T> Grid<T>& Grid<T>::safeDivide (const Grid<T>& a) {
	gridSafeDiv<T> (*this, a);
	return *this;
}
template<class T> Grid<T>& Grid<T>::copyFrom (const Grid<T>& a, bool copyType ) {
	assertMsg (a.mSize.x == mSize.x && a.mSize.y == mSize.y && a.mSize.z == mSize.z, "different grid resolutions "<<a.mSize<<" vs "<<this->mSize );
	memcpy(mData, a.mData, sizeof(T) * mSize.x * mSize.y * mSize.z);
	if(copyType) mType = a.mType; // copy type marker
	return *this;
}
/*template<class T> Grid<T>& Grid<T>::operator= (const Grid<T>& a) {
	note: do not use , use copyFrom instead
}*/

template<class T> void Grid<T>::add(const Grid<T>& a)  { (*this) += a; }
template<class T> void Grid<T>::sub(const Grid<T>& a)  { (*this) -= a; }
template<class T> void Grid<T>::mult(const Grid<T>& a) { (*this) *= a; }
template<class T> void Grid<T>::addConst(const T& a)   { (*this) += a; }
template<class T> void Grid<T>::multConst(const T& a)  { (*this) *= a; }
template<class T> void Grid<T>::addScaled(const Grid<T>& a, const T& factor) { 
	gridScaledAdd<T,T>(*this, a, factor);
}
template<class T> void Grid<T>::setConst(const T& a) {
	gridSetConst<T>(*this, a);
}
template<class T> void Grid<T>::clamp(const T& min, const T& max) {
	gridClamp<T>(*this, min, max);
}
template<class T> void Grid<T>::stomp(const T& threshold) {
	gridStomp<T>(*this, threshold);
}

template<> Real Grid<Real>::getMax() const {
	return CompMaxReal(*this);
}
template<> Real Grid<Real>::getMin() const {
	return CompMinReal(*this);
}
template<> Real Grid<Real>::getMaxAbs() const {
	Real amin = CompMinReal (*this);
	Real amax = CompMaxReal (*this);
	return max( fabs(amin), fabs(amax));
}
template<> Real Grid<Vec3>::getMax() const {
	return sqrt(CompMaxVec (*this));
}
template<> Real Grid<Vec3>::getMin() const {
	return sqrt(CompMinVec (*this));
}
template<> Real Grid<Vec3>::getMaxAbs() const {
	return sqrt(CompMaxVec (*this));
}
template<> Real Grid<int>::getMax() const {
	return (Real) CompMaxInt (*this);
}
template<> Real Grid<int>::getMin() const {
	return (Real) CompMinInt (*this);
}
template<> Real Grid<int>::getMaxAbs() const {
	int amin = CompMinInt (*this);
	int amax = CompMaxInt (*this);
	return max( fabs((Real)amin), fabs((Real)amax));
}
template<class T> std::string Grid<T>::getDataPointer() {
	std::ostringstream out;
	out << mData ;
	return out.str();
}

// compute maximal diference of two cells in the grid
// used for testing system
PYTHON() Real gridMaxDiff(Grid<Real>& g1, Grid<Real>& g2)
{
	double maxVal = 0.;
	FOR_IJK(g1) {
		maxVal = std::max(maxVal, (double)fabs(g1(i, j, k) - g2(i, j, k)));
	}
	return maxVal;
}
PYTHON() Real gridMaxDiffInt(Grid<int>& g1, Grid<int>& g2)
{
	double maxVal = 0.;
	FOR_IJK(g1) {
		maxVal = std::max(maxVal, (double)fabs((double)g1(i, j, k) - g2(i, j, k)));
	}
	return maxVal;
}
PYTHON() Real gridMaxDiffVec3(Grid<Vec3>& g1, Grid<Vec3>& g2)
{
	double maxVal = 0.;
	FOR_IJK(g1) {
		// accumulate differences with double precision
		// note - don't use norm here! should be as precise as possible...
		double d = 0.;
		for (int c = 0; c<3; ++c) {
			d += fabs((double)g1(i, j, k)[c] - (double)g2(i, j, k)[c]);
		}
		maxVal = std::max(maxVal, d);
		//maxVal = std::max(maxVal, (double)fabs( norm(g1(i,j,k)-g2(i,j,k)) ));
	}
	return maxVal;
}

PYTHON()
Vec3 calcGridSizeFactorWithRange(const Vec3 &from_old, const Vec3 &to_old, const Vec3 &from_new, const Vec3 &to_new) {
	return (to_new - from_new)/(to_old - from_old);
}

// simple helper functions to copy (convert) mac to vec3 , and levelset to real grids
// (are assumed to be the same for running the test cases - in general they're not!)
PYTHON() void copyMacToVec3 (MACGrid &source, Grid<Vec3>& target)
{
	FOR_IJK(target) {
		target(i,j,k) = source(i,j,k);
	}
}

PYTHON() void convertMacToVec3 (MACGrid &source , Grid<Vec3> &target) { debMsg("Deprecated - do not use convertMacToVec3... use copyMacToVec3 instead",1); copyMacToVec3(source,target); }
PYTHON() void moveMacToCen(const MACGrid &source, Grid<Vec3>& target) { GetCentered(target, source); }

//! vec3->mac grid conversion , but with full resampling 
PYTHON() void resampleVec3ToMac (Grid<Vec3>& source, MACGrid &target ) {
	FOR_IJK_BND(target,1) {
		target(i,j,k)[0] = 0.5*(source(i-1,j,k)[0]+source(i,j,k))[0];
		target(i,j,k)[1] = 0.5*(source(i,j-1,k)[1]+source(i,j,k))[1];
		if(target.is3D()) {
		target(i,j,k)[2] = 0.5*(source(i,j,k-1)[2]+source(i,j,k))[2]; }
	}
}
//! mac->vec3 grid conversion , with full resampling 
PYTHON() void resampleMacToVec3 (MACGrid &source, Grid<Vec3>& target ) {
	FOR_IJK_BND(target,1) {
		target(i,j,k) = source.getCentered(i,j,k);
	}
}

PYTHON() void copyLevelsetToReal (LevelsetGrid &source , Grid<Real> &target)
{
	FOR_IJK(target) {
		target(i,j,k) = source(i,j,k);
	}
}
PYTHON() void copyVec3ToReal (Grid<Vec3> &source, Grid<Real> &targetX, Grid<Real> &targetY, Grid<Real> &targetZ)
{
	FOR_IJK(source) {
		targetX(i,j,k) = source(i,j,k).x;
		targetY(i,j,k) = source(i,j,k).y;
		targetZ(i,j,k) = source(i,j,k).z;
	}
}

PYTHON() void copyRealToVec3 (Grid<Real> &sourceX, Grid<Real> &sourceY, Grid<Real> &sourceZ, Grid<Vec3> &target)
{
	FOR_IJK(target) {
		target(i,j,k).x = sourceX(i,j,k);
		target(i,j,k).y = sourceY(i,j,k);
		target(i,j,k).z = sourceZ(i,j,k);
	}
}
PYTHON() void convertLevelsetToReal (LevelsetGrid &source , Grid<Real> &target) { debMsg("Deprecated - do not use convertLevelsetToReal... use copyLevelsetToReal instead",1); copyLevelsetToReal(source,target); }

template<class T> void Grid<T>::printGrid(int zSlice, bool printIndex) {
	std::ostringstream out;
	out << std::endl;
	const int bnd = 1;
	FOR_IJK_BND(*this,bnd) {
		IndexInt idx = (*this).index(i,j,k);
		if(zSlice>=0 && k==zSlice) { 
			out << " ";
			if(printIndex &&  this->is3D()) out << "  "<<i<<","<<j<<","<<k <<":";
			if(printIndex && !this->is3D()) out << "  "<<i<<","<<j<<":";
			out << (*this)[idx]; 
			if(i==(*this).getSizeX()-1 -bnd) out << std::endl; 
		}
	}
	out << endl; debMsg("Printing " << this->getName() << out.str().c_str() , 1);
}

//! helper to swap components of a grid (eg for data import)
PYTHON() void swapComponents(Grid<Vec3>& vel, int c1=0, int c2=1, int c3=2) {
	FOR_IJK(vel) {
		Vec3 v = vel(i,j,k);
		vel(i,j,k)[0] = v[c1];
		vel(i,j,k)[1] = v[c2];
		vel(i,j,k)[2] = v[c3];
	}
}

// helper functions for UV grid data (stored grid coordinates as Vec3 values, and uv weight in entry zero)

// make uv weight accesible in python
PYTHON() Real getUvWeight (Grid<Vec3> &uv) { return uv[0][0]; }

// note - right now the UV grids have 0 values at the border after advection... could be fixed with an extrapolation step...

// compute normalized modulo interval
static inline Real computeUvGridTime(Real t, Real resetTime) {
	return fmod( (t / resetTime), (Real)1. );
}
// create ramp function in 0..1 range with half frequency
static inline Real computeUvRamp(Real t) {
	Real uvWeight = 2. * t; 
	if (uvWeight>1.) uvWeight=2.-uvWeight;
	return uvWeight;
}

KERNEL() void knResetUvGrid (Grid<Vec3>& target) { target(i,j,k) = Vec3((Real)i,(Real)j,(Real)k); }

PYTHON() void resetUvGrid (Grid<Vec3> &target)
{
	knResetUvGrid reset(target); // note, llvm complains about anonymous declaration here... ?
}
PYTHON() void updateUvWeight(Real resetTime, int index, int numUvs, Grid<Vec3> &uv)
{
	const Real t   = uv.getParent()->getTime();
	Real  timeOff  = resetTime/(Real)numUvs;

	Real lastt = computeUvGridTime(t +(Real)index*timeOff - uv.getParent()->getDt(), resetTime);
	Real currt = computeUvGridTime(t +(Real)index*timeOff                          , resetTime);
	Real uvWeight = computeUvRamp(currt);

	// normalize the uvw weights , note: this is a bit wasteful...
	Real uvWTotal = 0.;
	for(int i=0; i<numUvs; ++i) {
		uvWTotal += computeUvRamp( computeUvGridTime(t +(Real)i*timeOff , resetTime) );
	}
	if(uvWTotal<=VECTOR_EPSILON) { uvWeight =  uvWTotal = 1.; }
	else                           uvWeight /= uvWTotal;

	// check for reset
	if( currt < lastt ) 
		knResetUvGrid reset( uv );

	// write new weight value to grid
	uv[0] = Vec3( uvWeight, 0.,0.);

	// print info about uv weights?
	debMsg("Uv grid "<<index<<"/"<<numUvs<< " t="<<currt<<" w="<<uvWeight<<", reset:"<<(int)(currt<lastt) , 2);
}

KERNEL() template<class T> void knSetBoundary (Grid<T>& grid, T value, int w) { 
	bool bnd = (i<=w || i>=grid.getSizeX()-1-w || j<=w || j>=grid.getSizeY()-1-w || (grid.is3D() && (k<=w || k>=grid.getSizeZ()-1-w)));
	if (bnd) 
		grid(i,j,k) = value;
}

template<class T> void Grid<T>::setBound(T value, int boundaryWidth) {
	knSetBoundary<T>( *this, value, boundaryWidth );
}


KERNEL() template<class T> void knSetBoundaryNeumann (Grid<T>& grid, int w) { 
	bool set = false;
	int  si=i, sj=j, sk=k;
	if( i<=w) {
		si = w+1; set=true;
	}
	if( i>=grid.getSizeX()-1-w){
		si = grid.getSizeX()-1-w-1; set=true;
	}
	if( j<=w){
		sj = w+1; set=true;
	}
	if( j>=grid.getSizeY()-1-w){
		sj = grid.getSizeY()-1-w-1; set=true;
	}
	if( grid.is3D() ){
		 if( k<=w ) {
			sk = w+1; set=true;
		 }
		 if( k>=grid.getSizeZ()-1-w ) {
			sk = grid.getSizeZ()-1-w-1; set=true;
		 }
	}
	if(set)
		grid(i,j,k) = grid(si, sj, sk);
}

template<class T> void Grid<T>::setBoundNeumann(int boundaryWidth) {
	knSetBoundaryNeumann<T>( *this, boundaryWidth );
}

//! kernel to set velocity components of mac grid to value for a boundary of w cells
KERNEL() void knSetBoundaryMAC (Grid<Vec3>& grid, Vec3 value, int w) { 
	if (i<=w   || i>=grid.getSizeX()  -w || j<=w-1 || j>=grid.getSizeY()-1-w || (grid.is3D() && (k<=w-1 || k>=grid.getSizeZ()-1-w)))
		grid(i,j,k).x = value.x;
	if (i<=w-1 || i>=grid.getSizeX()-1-w || j<=w   || j>=grid.getSizeY()  -w || (grid.is3D() && (k<=w-1 || k>=grid.getSizeZ()-1-w)))
		grid(i,j,k).y = value.y;
	if (i<=w-1 || i>=grid.getSizeX()-1-w || j<=w-1 || j>=grid.getSizeY()-1-w || (grid.is3D() && (k<=w   || k>=grid.getSizeZ()  -w)))
		grid(i,j,k).z = value.z;
} 

//! only set normal velocity components of mac grid to value for a boundary of w cells
KERNEL() void knSetBoundaryMACNorm (Grid<Vec3>& grid, Vec3 value, int w) { 
	if (i<=w   || i>=grid.getSizeX()  -w ) grid(i,j,k).x = value.x;
	if (j<=w   || j>=grid.getSizeY()  -w ) grid(i,j,k).y = value.y;
	if ( (grid.is3D() && (k<=w   || k>=grid.getSizeZ()  -w))) grid(i,j,k).z = value.z;
} 

//! set velocity components of mac grid to value for a boundary of w cells (optionally only normal values)
void MACGrid::setBoundMAC(Vec3 value, int boundaryWidth, bool normalOnly) { 
	if(!normalOnly) knSetBoundaryMAC    ( *this, value, boundaryWidth ); 
	else            knSetBoundaryMACNorm( *this, value, boundaryWidth ); 
}

//! various sums of grid data
KERNEL(idx, reduce=+) returns(double result=0.0)
double knGridTotalSum(const Grid<Real>& a, const FlagGrid* flags) {
	if(flags) { if(flags->isFluid(idx)) result += a[idx]; }
	else      { result += a[idx]; }
}
KERNEL(idx, reduce=+) returns(double result=0.0)
double knGridTotalMagSum(const Grid<Real>& a, const FlagGrid* flags) {
	if(flags) { if(flags->isFluid(idx)) result += std::fabs(a[idx]); }
	else      { result += std::fabs(a[idx]); }
}
KERNEL(idx, reduce=+) returns(double result=0.0)
double knGridTotalVecMagSum(const Grid<Vec3>& a, const FlagGrid* flags) {
	if(flags) { if(flags->isFluid(idx)) result += norm(a[idx]); }
	else      { result += norm(a[idx]); }
}
KERNEL(idx, reduce=+) returns(double result=0.0)
double knGridTotalSqrSum(const Grid<Vec3>& a, const FlagGrid* flags) {
	if(flags) { if(flags->isFluid(idx)) result += normSquare(a[idx]); }
	else      { result += normSquare(a[idx]); }
}
PYTHON() Real getGridTotalSum(const Grid<Real>& a, const FlagGrid* flags=NULL) { return knGridTotalSum(a, flags); }
PYTHON() Real getGridTotalMagSum(const Grid<Real>& a, const FlagGrid* flags=NULL) { return knGridTotalMagSum(a, flags); }
PYTHON() Real getGridTotalVecMagSum(const Grid<Vec3>& a, const FlagGrid* flags=NULL) { return knGridTotalVecMagSum(a, flags); }
PYTHON() Real getGridTotalSqrSum(const Grid<Vec3>& a, const FlagGrid* flags=NULL) { return knGridTotalSqrSum(a, flags); }

KERNEL(idx, reduce=+) returns(int numEmpty=0)
int knCountFluidCells(const FlagGrid& flags) { if (flags.isFluid(idx) ) numEmpty++; }
PYTHON() int getCountFluidCells(const FlagGrid& flags) { return knCountFluidCells(flags); }

//! averaged value for all cells (if flags are given, only for fluid cells)
PYTHON() Real getGridAvg(Grid<Real>& source, FlagGrid* flags=NULL) 
{
	double sum = knGridTotalSum(source, flags);

	double cells;
	if(flags) { cells = knCountFluidCells(*flags); }
	else      { cells = source.getSizeX()*source.getSizeY()*source.getSizeZ(); }

	if(cells>0.) sum *= 1./cells;
	else         sum = -1.;
	return sum;
}

//! transfer data between real and vec3 grids

KERNEL(idx) void knGetComponent(Grid<Vec3>& source, Grid<Real>& target, int component) { 
	target[idx] = source[idx][component]; 
}
PYTHON() void getComponent(Grid<Vec3>& source, Grid<Real>& target, int component) { knGetComponent(source, target, component); }

KERNEL(idx) void knSetComponent(Grid<Real>& source, Grid<Vec3>& target, int component) { 
	target[idx][component] = source[idx]; 
}
PYTHON() void setComponent(Grid<Real>& source, Grid<Vec3>& target, int component) { knSetComponent(source, target, component); }

//******************************************************************************
// Specialization classes

void FlagGrid::InitMinXWall(const int &boundaryWidth, Grid<Real>& phiWalls) {
	const int w = boundaryWidth;
	FOR_IJK(phiWalls) {
		phiWalls(i,j,k) = std::min(i - w - .5, (double)phiWalls(i,j,k));
	}
}

void FlagGrid::InitMaxXWall(const int &boundaryWidth, Grid<Real>& phiWalls) {
	const int w = boundaryWidth;
	FOR_IJK(phiWalls) {
		phiWalls(i,j,k) = std::min(mSize.x-i-1.5-w, (double)phiWalls(i,j,k));
	}
}

void FlagGrid::InitMinYWall(const int &boundaryWidth, Grid<Real>& phiWalls) {
	const int w = boundaryWidth;
	FOR_IJK(phiWalls) {
		phiWalls(i,j,k) = std::min(j - w - .5, (double)phiWalls(i,j,k));
	}
}

void FlagGrid::InitMaxYWall(const int &boundaryWidth, Grid<Real>& phiWalls) {
	const int w = boundaryWidth;
	FOR_IJK(phiWalls) {
		phiWalls(i,j,k) = std::min(mSize.y-j-1.5-w, (double)phiWalls(i,j,k));
	}
}

void FlagGrid::InitMinZWall(const int &boundaryWidth, Grid<Real>& phiWalls) {
	const int w = boundaryWidth;
	FOR_IJK(phiWalls) {
		phiWalls(i,j,k) = std::min(k - w - .5, (double)phiWalls(i,j,k));
	}
}

void FlagGrid::InitMaxZWall(const int &boundaryWidth, Grid<Real>& phiWalls) {
	const int w = boundaryWidth;
	FOR_IJK(phiWalls) {
		phiWalls(i,j,k) = std::min(mSize.z-k-1.5-w, (double)phiWalls(i,j,k));
	}
}

void FlagGrid::initDomain(const int &boundaryWidth,
			  const string &wallIn,
			  const string &openIn,
			  const string &inflowIn,
			  const string &outflowIn,
			  Grid<Real>* phiWalls) {
	int  types[6] = {0};
	bool set  [6] = {false};
	// make sure we have at least 6 entries
	string wall    = wallIn;    wall.append("      ");
	string open    = openIn;    open.append("      ");
	string inflow  = inflowIn;  inflow.append("      ");
	string outflow = outflowIn; outflow.append("      ");

	if(phiWalls) phiWalls->setConst(1000000000);

	for (char i = 0; i<6; ++i) {
		//min x-direction
		if(!set[0]) {
			if(open[i]=='x')         {types[0] = TypeOpen;set[0] = true;}
			else if(inflow[i]=='x')  {types[0] = TypeInflow;set[0] = true;}
			else if(outflow[i]=='x') {types[0] = TypeOutflow;set[0] = true;}
			else if(wall[i]=='x') {
				types[0]    = TypeObstacle;
				if(phiWalls) InitMinXWall(boundaryWidth, *phiWalls);
				set[0] = true;
			}			
		}
		//max x-direction
		if(!set[1]) {
			if(open[i]=='X')         {types[1] = TypeOpen;set[1] = true;}
			else if(inflow[i]=='X')  {types[1] = TypeInflow;set[1] = true;}
			else if(outflow[i]=='X') {types[1] = TypeOutflow;set[1] = true;}
			else if(wall[i]=='X')  {
				types[1]    = TypeObstacle;
				if(phiWalls) InitMaxXWall(boundaryWidth, *phiWalls);
				set[1] = true;
			}			
		}
		//min y-direction
		if(!set[2]) {
			if(open[i]=='y')         {types[2] = TypeOpen;set[2] = true;}
			else if(inflow[i]=='y')  {types[2] = TypeInflow;set[2] = true;}
			else if(outflow[i]=='y') {types[2] = TypeOutflow;set[2] = true;}
			else if(wall[i]=='y') {
				types[2]    = TypeObstacle;
				if(phiWalls) InitMinYWall(boundaryWidth, *phiWalls);
				set[2] = true;
			}			
		}
		//max y-direction
		if(!set[3]) {
			if(open[i]=='Y')         {types[3] = TypeOpen;set[3] = true;}
			else if(inflow[i]=='Y')  {types[3] = TypeInflow;set[3] = true;}
			else if(outflow[i]=='Y') {types[3] = TypeOutflow;set[3] = true;}
			else if(wall[i]=='Y') {
				types[3]    = TypeObstacle;
				if(phiWalls) InitMaxYWall(boundaryWidth, *phiWalls);
				set[3] = true;
			}			
		}
		if(this->is3D()) {
		//min z-direction
			if(!set[4]) {
				if(open[i]=='z')         {types[4] = TypeOpen;set[4] = true;}
				else if(inflow[i]=='z')  {types[4] = TypeInflow;set[4] = true;}
				else if(outflow[i]=='z') {types[4] = TypeOutflow;set[4] = true;}
				else if(wall[i]=='z') {
					types[4]    = TypeObstacle;
					if(phiWalls) InitMinZWall(boundaryWidth, *phiWalls);
					set[4] = true;
				}				
			}
			//max z-direction
			if(!set[5]) {
				if(open[i]=='Z')         {types[5] = TypeOpen;set[5] = true;}
				else if(inflow[i]=='Z')  {types[5] = TypeInflow;set[5] = true;}
				else if(outflow[i]=='Z') {types[5] = TypeOutflow;set[5] = true;}
				else if(wall[i]=='Z') {
					types[5]    = TypeObstacle;
					if(phiWalls) InitMaxZWall(boundaryWidth, *phiWalls);
					set[5] = true;
				}				
			}
		}
	}

	setConst(TypeEmpty); 
	initBoundaries(boundaryWidth, types); 
}

void FlagGrid::initBoundaries(const int &boundaryWidth, const int *types) {
	const int w = boundaryWidth;
	FOR_IJK(*this) {
		bool bnd = (i <= w);
		if (bnd) mData[index(i,j,k)] = types[0];
		bnd = (i >= mSize.x-1-w);
		if (bnd) mData[index(i,j,k)] = types[1];
		bnd = (j <= w);
		if (bnd) mData[index(i,j,k)] = types[2];
		bnd = (j >= mSize.y-1-w);
		if (bnd) mData[index(i,j,k)] = types[3];
		if(is3D()) {
			bnd = (k <= w);
			if (bnd) mData[index(i,j,k)] = types[4];
			bnd = (k >= mSize.z-1-w);
			if (bnd) mData[index(i,j,k)] = types[5];
		}
	}
}

void FlagGrid::minifyFrom(const FlagGrid &flags, const Vec3i &scale) {
	FOR_IJK(*this) {
		if(isObstacle(i,j,k)) continue;
		bool fluidP = false;
		for(int si=0; si<scale.x; ++si) {
			for(int sj=0; sj<scale.y; ++sj) {
				for(int sk=0; sk<(is3D() ? scale.z : 1); ++sk) {
					fluidP |= flags.isFluid(i*scale.x+si, j*scale.y+sj, k*scale.z+sk);
				}
			}
		}

		(*this)(i,j,k) &= ~(TypeEmpty | TypeFluid); // clear empty/fluid flags
		(*this)(i,j,k) |= (fluidP) ? TypeFluid : TypeEmpty; // set resepctive flag
	}
}

void FlagGrid::updateFromLevelset(const LevelsetGrid& levelset) {
	if(getSize()!=levelset.getSize()) return updateFromLevelsetNonMatched(levelset);
	FOR_IDX(*this) {
		if (!isObstacle(idx)) {
			const Real phi = levelset[idx];
			if (phi <= levelset.invalidTimeValue()) continue;
			
			mData[idx] &= ~(TypeEmpty | TypeFluid); // clear empty/fluid flags
			mData[idx] |= (phi <= 0) ? TypeFluid : TypeEmpty; // set resepctive flag
		}
	}
}

KERNEL()
void knUpdateFromLevelsetNonMatched(FlagGrid &flags, const LevelsetGrid &levelset,
				    const Vec3 &factor, const Vec3 &offset, const int orderSpace) {
	if(flags.isObstacle(i,j,k)) return;
	Vec3 pos = Vec3(i,j,k)*factor + offset;
	if(!flags.is3D()) pos[2] = 0; // allow 2d -> 3d
	const Real phi = levelset.getInterpolatedHi(pos, orderSpace);
	if(phi <= levelset.invalidTimeValue()) return;

	flags(i,j,k) &= ~(FlagGrid::TypeEmpty | FlagGrid::TypeFluid); // clear empty/fluid flags
	flags(i,j,k) |= (phi <= 0) ? FlagGrid::TypeFluid : FlagGrid::TypeEmpty; // set resepctive flag
}

void FlagGrid::updateFromLevelsetNonMatched(const LevelsetGrid &levelset, const int orderSpace) {
	const Vec3 sourceFactor = calcGridSizeFactor(levelset.getSize(), getSize());
	const Vec3 offset       = sourceFactor*0.5;
	knUpdateFromLevelsetNonMatched(*this, levelset, sourceFactor, offset, orderSpace);
}

void FlagGrid::fillGrid(int type) {
	FOR_IDX(*this) {
		if ((mData[idx] & TypeObstacle)==0 && (mData[idx] & TypeInflow)==0&& (mData[idx] & TypeOutflow)==0&& (mData[idx] & TypeOpen)==0)
			mData[idx] = (mData[idx] & ~(TypeEmpty | TypeFluid)) | type;
	}
}

void FlagGrid::extendRegion(const int region, const int exclude, const int depth) {
	const int I=getSizeX()-1, J=getSizeY()-1, K=getSizeZ()-1;
	for(int i_depth=0; i_depth<depth; ++i_depth) {
		std::vector<int> update;
		FOR_IJK(*this) {
			if(get(i, j, k) & exclude) continue;
			if((i>0 && (get(i-1, j, k)&region)) || (i<I && (get(i+1, j, k)&region)) ||
			   (j>0 && (get(i, j-1, k)&region)) || (j<J && (get(i, j+1, k)&region)) ||
			   (is3D() && ((k>0 && (get(i, j, k-1)&region)) || (k<K && (get(i, j, k+1)&region)))))
				update.push_back(index(i, j, k));
		}

		for(std::vector<int>::const_iterator it=update.begin(); it!=update.end(); ++it) {
			mData[*it] = region;
		}
	}
}

void dfs(Grid<int> &r, const FlagGrid &flags, const IndexInt idx, const int c, const int type) {
	r(idx) = c;
	if((flags(idx-flags.getStrideX()) & type) && !r[idx-flags.getStrideX()]) dfs(r, flags, idx-flags.getStrideX(), c, type);
	if((flags(idx+flags.getStrideX()) & type) && !r[idx+flags.getStrideX()]) dfs(r, flags, idx+flags.getStrideX(), c, type);
	if((flags(idx-flags.getStrideY()) & type) && !r[idx-flags.getStrideY()]) dfs(r, flags, idx-flags.getStrideY(), c, type);
	if((flags(idx+flags.getStrideY()) & type) && !r[idx+flags.getStrideY()]) dfs(r, flags, idx+flags.getStrideY(), c, type);
	if(!flags.is3D()) return;
	if((flags(idx-flags.getStrideZ()) & type) && !r[idx-flags.getStrideZ()]) dfs(r, flags, idx-flags.getStrideZ(), c, type);
	if((flags(idx+flags.getStrideZ()) & type) && !r[idx+flags.getStrideZ()]) dfs(r, flags, idx+flags.getStrideZ(), c, type);
}

PYTHON() int getRegions(Grid<int> &r, const FlagGrid &flags, const int ctype) {
	r.clear();
	int n_regions = 0;

	FOR_IDX(flags) {
		if((flags(idx) & ctype) && !r(idx)) dfs(r, flags, idx, ++n_regions, ctype);
	}
	return n_regions;
}

PYTHON() void getRegionalCounts(Grid<int> &r, const FlagGrid &flags, const int ctype) {
	const int n_regions = getRegions(r, flags, ctype);
	std::vector<int> cnt(n_regions+1, 0);
	FOR_IDX(flags) {
		if(r[idx]>0) ++(cnt[r[idx]]);
	}
	FOR_IDX(flags) {
		r[idx] = cnt[r[idx]];
	}
}

bool isIsolatedFluidCell(const IndexInt idx, const FlagGrid &flags) {
	if(!flags.isFluid(idx)) return false;
	if(flags.isFluid(idx-flags.getStrideX())) return false;
	if(flags.isFluid(idx+flags.getStrideX())) return false;
	if(flags.isFluid(idx-flags.getStrideY())) return false;
	if(flags.isFluid(idx+flags.getStrideY())) return false;
	if(!flags.is3D()) return true;
	if(flags.isFluid(idx-flags.getStrideZ())) return false;
	if(flags.isFluid(idx+flags.getStrideZ())) return false;
	return true;
}

KERNEL(idx)
void knMarkIsolatedFluidCell(FlagGrid &flags, const int mark) {
	if(isIsolatedFluidCell(idx, flags)) flags[idx] = mark;
}

PYTHON()
void markIsolatedFluidCell(FlagGrid &flags, const int mark) {
	knMarkIsolatedFluidCell(flags, mark);
}

KERNEL(idx)
void knMarkSmallRegions(FlagGrid &flags, const Grid<int> &rcnt, const int mark, const int exclude, const int th) {
	if(flags[idx] & exclude) return;
	if(rcnt[idx] <= th) flags[idx] = mark;
}

PYTHON()
void markSmallRegions(FlagGrid &flags, const Grid<int> &rcnt, const int mark, const int exclude, const int th=1) {
	knMarkSmallRegions(flags, rcnt, mark, exclude, th);
}

PYTHON() void getGradientGrid(Grid<Vec3> &gradient, const Grid<Real> &grid) {
	GradientOp(gradient, grid);
}

PYTHON() void getDivergenceMAC(Grid<Real> &divergence, const MACGrid &grid) {
	DivergenceOpMAC(divergence, grid);
}

PYTHON() void getLaplacian(Grid<Real> &laplacian, const Grid<Real> &grid) {
	LaplaceOp(laplacian, grid);
}

// explicit instantiation
template class Grid<int>;
template class Grid<Real>;
template class Grid<Vec3>;

} //namespace
