#ifndef FASTLEVELSET_H
#define FASTLEVELSET_H

// A uniform levelset on [0,1]^3, supporting fairly fast
// operations (but not optimal storage).
// I use the convention that negative values are inside
// and positive values are outside.


typedef float (*AnalyticLevelSet) (float x, float y, float z);

typedef float *pointer_to_float;

class FastLevelSet {
	int N, Ncoarse, Nsub; // N = Ncoarse*Nsub
	float **coarse_raw; // Ncoarse^3 grid of float[Nsub^3] grids or NULLs
	float dx, coarsedx;
	float overdx, overcoarsedx;
	bool *inside_raw; // Ncoarse^3 grid indicating which unpopulated (coarse(i,j,k)==0) blocks are inside
	unsigned int status_counter;
	unsigned int *status_raw;
	float *temp_block;

	public:

	FastLevelSet (void);
	FastLevelSet (int N, AnalyticLevelSet phi);
	FastLevelSet (const FastLevelSet &f);
	~FastLevelSet (void);
	FastLevelSet &operator= (const FastLevelSet &f);

	void glRender (void) const;
	bool traceRay (const float orig[3], const float dir[3], float p[3]) const;
	float maxToolSize (void) const {return coarsedx;}
	void localAdd (const float center[3], float radius, float amount);
	void localSmooth (const float center[3], float radius, float amount);
	void write (const char *filename);
	void read (const char *filename);
	void initialize (AnalyticLevelSet phi);
	void refine (void);
	void enforceSharedValues (void);
	void redistance (void);

	void smoothRefined (void);
	private:

	pointer_to_float &coarse (const int i, const int j, const int k);
	const pointer_to_float &coarse (const int i, const int j, const int k) const;

	float &fine (float *block, const int i, const int j, const int k);
	const float &fine (const float *block, const int i, const int j, const int k) const;

	void resetStatus (void);
	bool &inside (const int i, const int j, const int k);
	const bool &inside (const int i, const int j, const int k) const;

	unsigned int &status (const int i, const int j, const int k);
	const unsigned int &status (const int i, const int j, const int k) const;

	void glRenderCoarse (const int i, const int j, const int k) const;
	void glRenderFine (const int i, const int j, const int k) const;
	void glRenderFineGood (const int i, const int j, const int k) const;
	void glRenderFineGreat (const int i, const int j, const int k) const;

	bool traceRayFine (int i, int j, int k, const float dir[3], float maxt, float p[3]) const;
	float evalFine (const float *block, float fx, float fy, float fz) const;

	void localRedistance (const int i, const int j, const int k);
	void enforceSharedBlock (const int i, const int j, const int k);
	void reinitBlock (int i, int j, int k, bool &containsBoundary);
	void initBlock (int i, int j, int k);
	void localDistanceSweeps (float *block);

	bool findphi (int i, int j, int k, int fi, int fj, int fk, float &phi);
	void refineNsub (void);

	void checkShared (const int i, const int j, const int k);
	void checkInsideness (const int i, const int j, const int k);
};

#endif

