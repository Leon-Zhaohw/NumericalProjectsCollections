/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Generate defined noise textures
 *
 *****************************************************************************/


#include "randomstream.h"
#include "vectorbase.h"
#include "paramset.h"
using DDF::Real;
using DDF::Vec3;
using DDF::nVec3i;

namespace DDF {
// global tile data, defined in globals.cpp
extern float* stdTileData;
}; 


namespace WAVELETNOISE {

#define NOISE_TILE_SIZE 128
static const int noiseTileSize = NOISE_TILE_SIZE;
static const char* stdTileName = "noiseTile3d.dqt";

//*****************************************************************************
// helper functions

static int modSlow(int x, int n) {
  int m = x % n; 
  return (m<0) ? m+n : m;
}

// warning - noiseTileSize has to be 128^3!
#define modFast128(x)  ((x) & 127)


//*****************************************************************************
// Wavelet downsampling
static float _aCoeffs[32] = {
		0.000334f,-0.001528f, 0.000410f, 0.003545f,-0.000938f,-0.008233f, 0.002172f, 0.019120f,
		-0.005040f,-0.044412f, 0.011655f, 0.103311f,-0.025936f,-0.243780f, 0.033979f, 0.655340f,
		0.655340f, 0.033979f,-0.243780f,-0.025936f, 0.103311f, 0.011655f,-0.044412f,-0.005040f,
		0.019120f, 0.002172f,-0.008233f,-0.000938f, 0.003546f, 0.000410f,-0.001528f, 0.000334f};
static void Downsample(float *from, float *to, int n, int stride){
	const float *a = &_aCoeffs[16];
	for (int i = 0; i < n / 2; i++) {
		to[i * stride] = 0;
		for (int k = 2 * i - 16; k <= 2 * i + 16; k++)
			to[i * stride] += a[k - 2 * i] * from[modFast128(k) * stride];
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
// Wavelet upsampling
//////////////////////////////////////////////////////////////////////////////////////////
static float _pCoeffs[4] = {0.25f, 0.75f, 0.75f, 0.25f};
static void Upsample(float *from, float *to, int n, int stride) {
	const float *p = &_pCoeffs[2];

	for (int i = 0; i < n; i++) {
		to[i * stride] = 0;
		for (int k = i / 2; k <= i / 2 + 1; k++)
			to[i * stride] += p[i - 2 * k] * from[modSlow(k, n / 2) * stride];
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
// load in an existing noise tile
//////////////////////////////////////////////////////////////////////////////////////////
static bool loadTile(float* const noiseTileData, std::string filename)
{
	FILE* file;
	file = fopen(filename.c_str(), "rb");

	if (file == NULL) {
		printf("loadTile: No noise tile '%s' found.\n", filename.c_str());
		return false;
	}

	// dimensions
	int gridSize = noiseTileSize * noiseTileSize * noiseTileSize;
	// noiseTileData memory is managed by caller
	int bread = fread((void*)noiseTileData, sizeof(float), gridSize, file);
	fclose(file);
	printf("Noise tile file loaded.\n");

	if (bread != gridSize) {
		printf("loadTile: Noise tile '%s' is wrong size %d.\n", filename.c_str(), bread);
		return false;
	} 
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////
// write out an existing noise tile
//////////////////////////////////////////////////////////////////////////////////////////
static void saveTile(float* const noiseTileData, std::string filename)
{
	FILE* file;
	file = fopen(filename.c_str(), "wb");

	if (file == NULL) {
		printf("saveTile: No noise tile '%s' found.\n", filename.c_str());
		return;
	} 

	unsigned len = noiseTileSize * noiseTileSize * noiseTileSize;
	if (fwrite((void*)noiseTileData, sizeof(float), len, file) != len) errFatal("saveTile","fwrite failed", SIMWORLD_INITERROR);
	fclose(file);

	printf("saveTile: Noise tile file saved.\n");
}

//////////////////////////////////////////////////////////////////////////////////////////
// create a new noise tile if necessary
//////////////////////////////////////////////////////////////////////////////////////////
static void GenerateNoiseTile(float* const noiseTileData, std::string filename) {
	// if a tile already exists, just use that
	if (loadTile(noiseTileData, filename)) return;

	const int n = noiseTileSize;
	const int n3 = n*n*n;
	std::cout <<"Generating new 3d noise tile size="<<n<<"^3 \n";
	//MTRand twister;
	DDF::RandomStream randStream = DDF::RandomStream(13322223);
	DDF::RandomStream randStream2 = DDF::RandomStream(1238201283);

	float *temp13 = new float[n3];
	float *temp23 = new float[n3];
	float *noise3 = new float[n3];

	// initialize
	for (int i = 0; i < n3; i++) {
		temp13[i] = temp23[i] =
			noise3[i] = 0.;
	}

	// Step 1. Fill the tile with random numbers in the range -1 to 1.
	for (int i = 0; i < n3; i++) {
		//noise3[i] = twister.randNorm();
		noise3[i] = (randStream.getFloat() + randStream2.getFloat()) -1.; // produces repeated values??
		noise3[i] = (randStream.getFloat() + randStream2.getFloat()) -1.;
	}

	// Steps 2 and 3. Downsample and upsample the tile
	for (int iy = 0; iy < n; iy++) 
		for (int iz = 0; iz < n; iz++) {
			const int i = iy * n + iz*n*n;
			Downsample(&noise3[i], &temp13[i], n, 1);
			Upsample  (&temp13[i], &temp23[i], n, 1);
		}
	for (int ix = 0; ix < n; ix++) 
		for (int iz = 0; iz < n; iz++) {
			const int i = ix + iz*n*n;
			Downsample(&temp23[i], &temp13[i], n, n);
			Upsample  (&temp13[i], &temp23[i], n, n);
		}
	for (int ix = 0; ix < n; ix++) 
		for (int iy = 0; iy < n; iy++) {
			const int i = ix + iy*n;
			Downsample(&temp23[i], &temp13[i], n, n*n);
			Upsample  (&temp13[i], &temp23[i], n, n*n);
		}

	// Step 4. Subtract out the coarse-scale contribution
	for (int i = 0; i < n3; i++) { 
		noise3[i] -= temp23[i];
	}

	// Avoid even/odd variance difference by adding odd-offset version of noise to itself.
	int offset = n / 2;
	if (offset % 2 == 0) offset++;

	int icnt=0;
	for (int ix = 0; ix < n; ix++)
		for (int iy = 0; iy < n; iy++)
			for (int iz = 0; iz < n; iz++) { 
				temp13[icnt] = noise3[modFast128(ix+offset) + modFast128(iy+offset)*n + modFast128(iz+offset)*n*n];
				icnt++;
			}

	for (int i = 0; i < n3; i++) {
		noise3[i] += temp13[i];
	}

	for (int i = 0; i < n3; i++) {
			noiseTileData[i] = noise3[i];
	}
	saveTile(noise3, filename); 

	delete[] temp13;
	delete[] temp23;
	std::cout <<"Generating new 3d noise done\n";
}


#define ADD_WEIGHTED(x,y,z)\
  weight = 1.0f;\
  xC = modFast128(midX + (x));\
  weight *= w[0][(x) + 1];\
  yC = modFast128(midY + (y));\
  weight *= w[1][(y) + 1];\
  zC = modFast128(midZ + (z));\
  weight *= w[2][(z) + 1];\
  result += weight * data[(zC * NOISE_TILE_SIZE + yC) * NOISE_TILE_SIZE + xC];

//////////////////////////////////////////////////////////////////////////////////////////
// derivatives of 3D noise - unrolled for performance
//////////////////////////////////////////////////////////////////////////////////////////
static inline float WNoiseDx(Vec3 p, float *data) {
	float w[3][3], t, result = 0;

	// Evaluate quadratic B-spline basis functions
	int midX = (int)ceil(p[0] - 0.5f); 
	t        =   midX - (p[0] - 0.5f);
	w[0][0] = -t;
	w[0][2] = (1.f - t);
	w[0][1] = 2.0f * t - 1.0f;

	int midY = (int)ceil(p[1] - 0.5f); 
	t        =   midY - (p[1] - 0.5f);
	w[1][0] = t * t * 0.5f; 
	w[1][2] = (1.f - t) * (1.f - t) *0.5f; 
	w[1][1] = 1.f - w[1][0] - w[1][2];

	int midZ = (int)ceil(p[2] - 0.5f); 
	t        =   midZ - (p[2] - 0.5f);
	w[2][0] = t * t * 0.5f; 
	w[2][2] = (1.f - t) * (1.f - t) *0.5f; 
	w[2][1] = 1.f - w[2][0] - w[2][2];

	// Evaluate noise by weighting noise coefficients by basis function values
	int xC, yC, zC;
	float weight = 1;

	ADD_WEIGHTED(-1,-1, -1); ADD_WEIGHTED( 0,-1, -1); ADD_WEIGHTED( 1,-1, -1);
	ADD_WEIGHTED(-1, 0, -1); ADD_WEIGHTED( 0, 0, -1); ADD_WEIGHTED( 1, 0, -1);
	ADD_WEIGHTED(-1, 1, -1); ADD_WEIGHTED( 0, 1, -1); ADD_WEIGHTED( 1, 1, -1);

	ADD_WEIGHTED(-1,-1, 0);  ADD_WEIGHTED( 0,-1, 0);  ADD_WEIGHTED( 1,-1, 0);
	ADD_WEIGHTED(-1, 0, 0);  ADD_WEIGHTED( 0, 0, 0);  ADD_WEIGHTED( 1, 0, 0);
	ADD_WEIGHTED(-1, 1, 0);  ADD_WEIGHTED( 0, 1, 0);  ADD_WEIGHTED( 1, 1, 0);

	ADD_WEIGHTED(-1,-1, 1);  ADD_WEIGHTED( 0,-1, 1);  ADD_WEIGHTED( 1,-1, 1);
	ADD_WEIGHTED(-1, 0, 1);  ADD_WEIGHTED( 0, 0, 1);  ADD_WEIGHTED( 1, 0, 1);
	ADD_WEIGHTED(-1, 1, 1);  ADD_WEIGHTED( 0, 1, 1);  ADD_WEIGHTED( 1, 1, 1);

	return result;
}

static inline float WNoise(Vec3 p, float *data) {
	float w[3][3], t, result = 0;

	// Evaluate quadratic B-spline basis functions
	int midX = (int)ceil(p[0] - 0.5f); 
	t        =   midX - (p[0] - 0.5f);
	w[0][0] = t * t * 0.5f; 
	w[0][2] = (1.f - t) * (1.f - t) *0.5f; 
	w[0][1] = 1.f - w[0][0] - w[0][2];

	int midY = (int)ceil(p[1] - 0.5f); 
	t        =   midY - (p[1] - 0.5f);
	w[1][0] = t * t * 0.5f; 
	w[1][2] = (1.f - t) * (1.f - t) *0.5f; 
	w[1][1] = 1.f - w[1][0] - w[1][2];

	int midZ = (int)ceil(p[2] - 0.5f); 
	t        =   midZ - (p[2] - 0.5f);
	w[2][0] = t * t * 0.5f; 
	w[2][2] = (1.f - t) * (1.f - t) *0.5f; 
	w[2][1] = 1.f - w[2][0] - w[2][2];

	// Evaluate noise by weighting noise coefficients by basis function values
	int xC, yC, zC;
	float weight = 1;

	ADD_WEIGHTED(-1,-1, -1); ADD_WEIGHTED( 0,-1, -1); ADD_WEIGHTED( 1,-1, -1);
	ADD_WEIGHTED(-1, 0, -1); ADD_WEIGHTED( 0, 0, -1); ADD_WEIGHTED( 1, 0, -1);
	ADD_WEIGHTED(-1, 1, -1); ADD_WEIGHTED( 0, 1, -1); ADD_WEIGHTED( 1, 1, -1);

	ADD_WEIGHTED(-1,-1, 0);  ADD_WEIGHTED( 0,-1, 0);  ADD_WEIGHTED( 1,-1, 0);
	ADD_WEIGHTED(-1, 0, 0);  ADD_WEIGHTED( 0, 0, 0);  ADD_WEIGHTED( 1, 0, 0);
	ADD_WEIGHTED(-1, 1, 0);  ADD_WEIGHTED( 0, 1, 0);  ADD_WEIGHTED( 1, 1, 0);

	ADD_WEIGHTED(-1,-1, 1);  ADD_WEIGHTED( 0,-1, 1);  ADD_WEIGHTED( 1,-1, 1);
	ADD_WEIGHTED(-1, 0, 1);  ADD_WEIGHTED( 0, 0, 1);  ADD_WEIGHTED( 1, 0, 1);
	ADD_WEIGHTED(-1, 1, 1);  ADD_WEIGHTED( 0, 1, 1);  ADD_WEIGHTED( 1, 1, 1);

	return result;
}



#define ADD_WEIGHTEDX(x,y,z)\
  weight = dw[0][(x) + 1] * w[1][(y) + 1] * w[2][(z) + 1];\
  result += weight * neighbors[x + 1][y + 1][z + 1];

#define ADD_WEIGHTEDY(x,y,z)\
  weight = w[0][(x) + 1] * dw[1][(y) + 1] * w[2][(z) + 1];\
  result += weight * neighbors[x + 1][y + 1][z + 1];

#define ADD_WEIGHTEDZ(x,y,z)\
  weight = w[0][(x) + 1] * w[1][(y) + 1] * dw[2][(z) + 1];\
  result += weight * neighbors[x + 1][y + 1][z + 1];

//////////////////////////////////////////////////////////////////////////////////////////
// compute all derivatives in at once
//////////////////////////////////////////////////////////////////////////////////////////
static inline Vec3 WNoiseVec(const Vec3& p, float *data)
{
	Vec3 final(0.);
	float w[3][3];
	float dw[3][3];
	float result = 0;
	int xC, yC, zC;
	float weight;

	int midX = (int)ceil(p[0] - 0.5f); 
	int midY = (int)ceil(p[1] - 0.5f); 
	int midZ = (int)ceil(p[2] - 0.5f);

	float t0 =   midX - (p[0] - 0.5f);
	float t1 =   midY - (p[1] - 0.5f);
	float t2 =   midZ - (p[2] - 0.5f);

	// precache all the neighbors for fast access
	float neighbors[3][3][3];
	for (int z = -1; z <=1; z++)
		for (int y = -1; y <= 1; y++)
			for (int x = -1; x <= 1; x++)
			{
				xC = modFast128(midX + (x));
				yC = modFast128(midY + (y));
				zC = modFast128(midZ + (z));
				neighbors[x + 1][y + 1][z + 1] = data[zC * NOISE_TILE_SIZE * NOISE_TILE_SIZE + yC * NOISE_TILE_SIZE + xC];
			}

	///////////////////////////////////////////////////////////////////////////////////////
	// evaluate splines
	///////////////////////////////////////////////////////////////////////////////////////
	dw[0][0] = -t0;
	dw[0][2] = (1.f - t0);
	dw[0][1] = 2.0f * t0 - 1.0f;

	dw[1][0] = -t1;
	dw[1][2] = (1.0f - t1);
	dw[1][1] = 2.0f * t1 - 1.0f;

	dw[2][0] = -t2;
	dw[2][2] = (1.0f - t2);
	dw[2][1] = 2.0f * t2 - 1.0f;

	w[0][0] = t0 * t0 * 0.5f; 
	w[0][2] = (1.f - t0) * (1.f - t0) *0.5f; 
	w[0][1] = 1.f - w[0][0] - w[0][2];

	w[1][0] = t1 * t1 * 0.5f; 
	w[1][2] = (1.f - t1) * (1.f - t1) *0.5f; 
	w[1][1] = 1.f - w[1][0] - w[1][2];

	w[2][0] = t2 * t2 * 0.5f; 
	w[2][2] = (1.f - t2) * (1.f - t2) *0.5f;
	w[2][1] = 1.f - w[2][0] - w[2][2];

	///////////////////////////////////////////////////////////////////////////////////////
	// x derivative
	///////////////////////////////////////////////////////////////////////////////////////
	result = 0.0f;
	ADD_WEIGHTEDX(-1,-1, -1); ADD_WEIGHTEDX( 0,-1, -1); ADD_WEIGHTEDX( 1,-1, -1);
	ADD_WEIGHTEDX(-1, 0, -1); ADD_WEIGHTEDX( 0, 0, -1); ADD_WEIGHTEDX( 1, 0, -1);
	ADD_WEIGHTEDX(-1, 1, -1); ADD_WEIGHTEDX( 0, 1, -1); ADD_WEIGHTEDX( 1, 1, -1);

	ADD_WEIGHTEDX(-1,-1, 0);  ADD_WEIGHTEDX( 0,-1, 0);  ADD_WEIGHTEDX( 1,-1, 0);
	ADD_WEIGHTEDX(-1, 0, 0);  ADD_WEIGHTEDX( 0, 0, 0);  ADD_WEIGHTEDX( 1, 0, 0);
	ADD_WEIGHTEDX(-1, 1, 0);  ADD_WEIGHTEDX( 0, 1, 0);  ADD_WEIGHTEDX( 1, 1, 0);

	ADD_WEIGHTEDX(-1,-1, 1);  ADD_WEIGHTEDX( 0,-1, 1);  ADD_WEIGHTEDX( 1,-1, 1);
	ADD_WEIGHTEDX(-1, 0, 1);  ADD_WEIGHTEDX( 0, 0, 1);  ADD_WEIGHTEDX( 1, 0, 1);
	ADD_WEIGHTEDX(-1, 1, 1);  ADD_WEIGHTEDX( 0, 1, 1);  ADD_WEIGHTEDX( 1, 1, 1);
	final[0] = result;

	///////////////////////////////////////////////////////////////////////////////////////
	// y derivative
	///////////////////////////////////////////////////////////////////////////////////////
	result = 0.0f;
	ADD_WEIGHTEDY(-1,-1, -1); ADD_WEIGHTEDY( 0,-1, -1); ADD_WEIGHTEDY( 1,-1, -1);
	ADD_WEIGHTEDY(-1, 0, -1); ADD_WEIGHTEDY( 0, 0, -1); ADD_WEIGHTEDY( 1, 0, -1);
	ADD_WEIGHTEDY(-1, 1, -1); ADD_WEIGHTEDY( 0, 1, -1); ADD_WEIGHTEDY( 1, 1, -1);

	ADD_WEIGHTEDY(-1,-1, 0);  ADD_WEIGHTEDY( 0,-1, 0);  ADD_WEIGHTEDY( 1,-1, 0);
	ADD_WEIGHTEDY(-1, 0, 0);  ADD_WEIGHTEDY( 0, 0, 0);  ADD_WEIGHTEDY( 1, 0, 0);
	ADD_WEIGHTEDY(-1, 1, 0);  ADD_WEIGHTEDY( 0, 1, 0);  ADD_WEIGHTEDY( 1, 1, 0);

	ADD_WEIGHTEDY(-1,-1, 1);  ADD_WEIGHTEDY( 0,-1, 1);  ADD_WEIGHTEDY( 1,-1, 1);
	ADD_WEIGHTEDY(-1, 0, 1);  ADD_WEIGHTEDY( 0, 0, 1);  ADD_WEIGHTEDY( 1, 0, 1);
	ADD_WEIGHTEDY(-1, 1, 1);  ADD_WEIGHTEDY( 0, 1, 1);  ADD_WEIGHTEDY( 1, 1, 1);
	final[1] = result;

	///////////////////////////////////////////////////////////////////////////////////////
	// z derivative
	///////////////////////////////////////////////////////////////////////////////////////
	result = 0.0f;
	ADD_WEIGHTEDZ(-1,-1, -1); ADD_WEIGHTEDZ( 0,-1, -1); ADD_WEIGHTEDZ( 1,-1, -1);
	ADD_WEIGHTEDZ(-1, 0, -1); ADD_WEIGHTEDZ( 0, 0, -1); ADD_WEIGHTEDZ( 1, 0, -1);
	ADD_WEIGHTEDZ(-1, 1, -1); ADD_WEIGHTEDZ( 0, 1, -1); ADD_WEIGHTEDZ( 1, 1, -1);

	ADD_WEIGHTEDZ(-1,-1, 0);  ADD_WEIGHTEDZ( 0,-1, 0);  ADD_WEIGHTEDZ( 1,-1, 0);
	ADD_WEIGHTEDZ(-1, 0, 0);  ADD_WEIGHTEDZ( 0, 0, 0);  ADD_WEIGHTEDZ( 1, 0, 0);
	ADD_WEIGHTEDZ(-1, 1, 0);  ADD_WEIGHTEDZ( 0, 1, 0);  ADD_WEIGHTEDZ( 1, 1, 0);

	ADD_WEIGHTEDZ(-1,-1, 1);  ADD_WEIGHTEDZ( 0,-1, 1);  ADD_WEIGHTEDZ( 1,-1, 1);
	ADD_WEIGHTEDZ(-1, 0, 1);  ADD_WEIGHTEDZ( 0, 0, 1);  ADD_WEIGHTEDZ( 1, 0, 1);
	ADD_WEIGHTEDZ(-1, 1, 1);  ADD_WEIGHTEDZ( 0, 1, 1);  ADD_WEIGHTEDZ( 1, 1, 1);
	final[2] = result;

	//debMsg("FINAL","at "<<p<<" = "<<final); // DEBUG
	return final;
}
#undef ADD_WEIGHTEDX
#undef ADD_WEIGHTEDY
#undef ADD_WEIGHTEDZ


//*****************************************************************************
static bool loadStdTile()
{
	if(!DDF::stdTileData) {
		DDF::stdTileData = new float[NOISE_TILE_SIZE * NOISE_TILE_SIZE * NOISE_TILE_SIZE];
	}
	GenerateNoiseTile(DDF::stdTileData, stdTileName);
	return true;
}

static void freeStdTile()
{
	if(DDF::stdTileData) {
		delete [] DDF::stdTileData;
	}
}

static inline float getWNoise(Vec3 p) {
#ifdef DDF_DEBUG 
#if DDF_DEBUG==1
	if(!DDF::stdTileData) {
		fprintf(stderr,"waveletnoise.h: Fatal, noise not inited!\n");
		exit(1);
	}
#endif // DDF_DEBUG==1
#endif // DDF_DEBUG 
	return WNoise(p, DDF::stdTileData);
}

static inline Vec3 getWNoiseVec(Vec3 p) {
#ifdef DDF_DEBUG 
#if DDF_DEBUG==1
	if(!DDF::stdTileData) {
		fprintf(stderr,"waveletnoise.h: Fatal, noise not inited!\n");
		exit(1);
	}
#endif // DDF_DEBUG==1
#endif // DDF_DEBUG 
	return WNoiseVec(p, DDF::stdTileData);
}


// wrapper for a parametrized field of wavelet noise
class WaveletNoiseField {
	public:
		WaveletNoiseField(std::string name, const ParamSet &params, nVec3i gridsize) :
				mName(name),
				mPosOffset(0.), mPosScale(1.),
				mValOffset(0.), mValScale(1.),
				mClamp(false), mClampNeg(0.), mClampPos(1.),
				mTimeAnim(0.) , mTime(0.)
		{

			// init params
			mPosOffset = params.FindOneVector("pos-offset", mPosOffset );
			mPosScale  = params.FindOneVector("pos-scale", mPosScale );

			mValOffset = params.FindOneFloat("val-offset", mValOffset );
			mValScale  = params.FindOneFloat("val-scale", mValScale );

			mClamp = params.FindOneBool("clamp", mClamp );
			mClampNeg = params.FindOneFloat("clamp-neg", mClampNeg );
			mClampPos = params.FindOneFloat("clamp-pos", mClampPos );

			mTimeAnim = params.FindOneFloat("time-anim", mTimeAnim );

			// init inv grid size
			mGsInvX = 1. /(Real)gridsize[0];
			mGsInvY = 1. /(Real)gridsize[1];
			mGsInvZ = 1. /(Real)gridsize[2];

			if(!loadStdTile()) {
				errFatal("WaveletNoiseField::WaveletNoiseField","Failed to laod/generate tile! Aborting...", SIMWORLD_INITERROR);
				exit(1);
			}
			debMsg("WaveletNoiseField","Created "<<this->toString() );
		};
		~WaveletNoiseField() {};

		//! evaluate noise
		inline Real evaluate(Vec3 pos) { 
			pos[0] *= mGsInvX;
			pos[1] *= mGsInvY;
			pos[2] *= mGsInvZ;

			// time anim
			pos += Vec3(mTime);

			pos[0] *= mPosScale[0];
			pos[1] *= mPosScale[1];
			pos[2] *= mPosScale[2];
			pos += mPosOffset;

			Real v = getWNoise(pos);

			v += mValOffset;
			v *= mValScale;
			if (mClamp) {
				if (v< mClampNeg) v = mClampNeg;
				if (v> mClampPos) v = mClampPos;
			}
			return v;
		}

		//! evaluate noise as a vector
		inline Vec3 evaluateVec(Vec3 pos) { 
			pos[0] *= mGsInvX;
			pos[1] *= mGsInvY;
			pos[2] *= mGsInvZ;

			// time anim
			pos += Vec3(mTime);

			pos[0] *= mPosScale[0];
			pos[1] *= mPosScale[1];
			pos[2] *= mPosScale[2];
			pos += mPosOffset;

			Vec3 v = getWNoiseVec(pos);

			v += Vec3(mValOffset);
			v *= mValScale;
			if (mClamp) {
				for(int i=0; i<3; i++) {
					if (v[i]< mClampNeg) v[i] = mClampNeg;
					if (v[i]> mClampPos) v[i] = mClampPos;
				}
			}
			return v;
		}

		// increase time
		void increaseTime(Real dt) { mTime += dt * mTimeAnim; }

		//! helper
		std::string toString() {
			std::ostringstream out;
			out <<	"NoiseField: name '"<<mName<<"' "<<
				"  pos off="<<mPosOffset<<" scale="<<mPosScale<<
				"  val off="<<mValOffset<<" scale="<<mValScale<<
				"  clamp ="<<mClamp<<" val="<<mClampNeg<<" to "<<mClampPos<<
				"  timeAni ="<<mTimeAnim<<
				"  gridInv ="<<PRINT_VEC(mGsInvX,mGsInvY,mGsInvZ) ;
			return out.str();
		}

		Real getTimeAnim() { return mTimeAnim; }
		void setTimeAnim(Real set) { mTimeAnim = set; }

	protected:

		//! name
		std::string mName;
		//! texcoord position and scale
		Vec3 mPosOffset, mPosScale;
		//! value offset & scale
		Real mValOffset, mValScale;
		//! clamp? (default 0-1)
		bool mClamp;
		Real mClampNeg, mClampPos; 
		// grid size normalization (inverse size)
		Real mGsInvX, mGsInvY, mGsInvZ;
		// animation over time
		Real mTimeAnim, mTime;
};


} // namespace WAVELETNOISE 

