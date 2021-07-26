#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <OpenGL/gl.h>
#include "fastlevelset.h"

// PLEASE NOTE: I assume sizeof(int)==sizeof(float)==4 in the read/write code


using namespace std;


//============================= inlines etc. ==================================

#if __BYTE_ORDER == __LITTLE_ENDIAN
static inline bool fwrite_32 (const void *buf, unsigned int n, FILE *stream)
{
	return (fwrite (buf, 4, n, stream) != n);
}

static inline bool fread_32 (void *buf, unsigned int n, FILE *stream)
{
	return (fread (buf, 4, n, stream) != n);
}
#else
static bool fwrite_32 (const void *buf, unsigned int n, FILE *stream)
{
	for (int i=0; i < n; ++i) {
		if (fwrite (4*i+3+(char*)buf, 1, 1, stream) != 1) return 1;
		if (fwrite (4*i+2+(char*)buf, 1, 1, stream) != 1) return 1;
		if (fwrite (4*i+1+(char*)buf, 1, 1, stream) != 1) return 1;
		if (fwrite (4*i+0+(char*)buf, 1, 1, stream) != 1) return 1;
	}
	return 0;
}

static inline bool fread_32 (void *buf, unsigned int n, FILE *stream)
{
	for (int i=0; i < n; ++i) {
		if (fread (4*i+3+(char*)buf, 1, 1, stream) != 1) return 1;
		if (fread (4*i+2+(char*)buf, 1, 1, stream) != 1) return 1;
		if (fread (4*i+1+(char*)buf, 1, 1, stream) != 1) return 1;
		if (fread (4*i+0+(char*)buf, 1, 1, stream) != 1) return 1;
	}
	return 0;
}
#endif


static inline int cube (int x)
{ return x*x*x; }


static inline float magnitude (float x, float y, float z)
{
	return sqrt(x*x+y*y+z*z);
}


static inline void divideByNormSquared (float v[3])
{
	float overmag2 = 1./sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	v[0] *= overmag2; v[1] *= overmag2; v[2] *= overmag2;
}


static inline float min2 (float a, float b)
{ if (a<b) return a; else return b; }


static inline float max2 (float a, float b)
{ if (a>b) return a; else return b; }


static inline float min3 (float a, float b, float c)
{ if (a<b) { if (a<c) return a; else return c; } else if (b<c) return b; else return c; }


static inline float max3 (float a, float b, float c)
{ if (a>b) { if (a>c) return a; else return c; } else if (b>c) return b; else return c; }


static inline float flabs (float a)
{ return (a > 0) ? a : -a; }


static inline float minmod (float a, float b)
{
	if (a <= b) {
		if (b <= 0) return b;
		else if (a >= 0) return a;
		else return 0;
	} else {
		if (a <= 0) return a;
		else if (b >= 0) return b;
		else return 0;
	}
} 


static inline float sign (float a)
{
	if (a <= 0) return -1;
	else return 1;
}


inline pointer_to_float &FastLevelSet::
coarse (const int i, const int j, const int k)
{ return coarse_raw[i+Ncoarse*(j+Ncoarse*k)]; }


inline const pointer_to_float &FastLevelSet::
coarse (const int i, const int j, const int k) const
{ return coarse_raw[i+Ncoarse*(j+Ncoarse*k)]; }


inline float &FastLevelSet::
fine (float *block, const int i, const int j, const int k)
{ return block[i+(Nsub+1)*(j+(Nsub+1)*k)]; }


inline const float &FastLevelSet::
fine (const float *block, const int i, const int j, const int k) const
{ return block[i+(Nsub+1)*(j+(Nsub+1)*k)]; }


inline bool &FastLevelSet::
inside (const int i, const int j, const int k)
{ return inside_raw[i+Ncoarse*(j+Ncoarse*k)]; }


inline const bool &FastLevelSet::
inside (const int i, const int j, const int k) const
{ return inside_raw[i+Ncoarse*(j+Ncoarse*k)]; }


inline void FastLevelSet::
resetStatus (void)
{
	++status_counter;
	if (status_counter == 0) {
		status_counter = 1;
		memset (status_raw, 0, sizeof(unsigned int)*cube(Nsub+1));
	}
}


inline unsigned int &FastLevelSet::
status (int i, int j, int k)
{ return status_raw[i+(Nsub+1)*(j+(Nsub+1)*k)]; }


inline const unsigned int &FastLevelSet::
status (int i, int j, int k) const
{ return status_raw[i+(Nsub+1)*(j+(Nsub+1)*k)]; }


//=========================== constructors etc. ==================================

FastLevelSet::
FastLevelSet (void)
	: N(0), Ncoarse(0), Nsub(0), coarse_raw(0), dx(0), coarsedx(0), overdx(0), overcoarsedx(0),
	  inside_raw(0),status_counter(0), status_raw(0), temp_block(0)
{}


FastLevelSet::
FastLevelSet (int desiredN, AnalyticLevelSet phi)
{
	Ncoarse = (int) (0.8*pow(desiredN,.75));
	if (Ncoarse < 1) Ncoarse = 1;
	Nsub = desiredN/Ncoarse;
	if (Ncoarse*Nsub < desiredN) ++Nsub;
	if (Nsub < 2) Nsub = 2;
	Ncoarse = desiredN/Nsub;
	if (Ncoarse*Nsub < desiredN) ++Ncoarse;
	if (Ncoarse < 1) Ncoarse = 1;
	N = Ncoarse*Nsub;
	printf("N = %d, Ncoarse = %d, Nsub = %d\n", N, Ncoarse, Nsub);
	dx = 1./N; overdx = N;
	coarsedx = Nsub*dx; overcoarsedx = 1./coarsedx;
	coarse_raw = new (float*)[cube(Ncoarse)];
	inside_raw = new bool[cube(Ncoarse)];
	for (int i = cube(Ncoarse)-1; i >= 0; --i) {coarse_raw[i] = 0; inside_raw[i] = false;}
	status_counter = 0;
	status_raw = new unsigned int[cube(Nsub+1)];
	memset (status_raw, 0, sizeof(unsigned int)*cube(Nsub+1));
	temp_block = new float[cube(Nsub+1)];
	initialize (phi);
}


FastLevelSet::
FastLevelSet (const FastLevelSet &f)
	: N(f.N), Ncoarse(f.Ncoarse), Nsub(f.Nsub), coarse_raw(new (float*)[cube(f.Ncoarse)]),
	  dx(f.dx), coarsedx(f.coarsedx), overdx(f.overdx), overcoarsedx(f.overcoarsedx),
	  inside_raw(new bool[cube(f.Ncoarse)]), temp_block(new float[cube(f.Nsub+1)])
{
	for (int i = cube(Ncoarse)-1; i >= 0; --i) {
		if (f.coarse_raw[i]) {
			coarse_raw[i] = new float[cube(Nsub+1)];
			memcpy (coarse_raw[i], f.coarse_raw[i], cube(Nsub+1)*sizeof(float));
		} else
			coarse_raw[i] = 0;
		inside_raw[i] = f.inside_raw[i];
	}
	status_counter = 0;
	status_raw = new unsigned int[cube(Nsub+1)];
	memset (status_raw, 0, sizeof(unsigned int)*cube(Nsub+1));
}


FastLevelSet::
~FastLevelSet (void)
{
	for (int i = cube(Ncoarse)-1; i>= 0;  --i) delete[] coarse_raw[i];
	delete[] coarse_raw;
	coarse_raw = 0;
	delete[] inside_raw;
	inside_raw = 0;
	delete[] status_raw;
	status_raw = 0;
	delete[] temp_block;
	temp_block = 0;
}


FastLevelSet &FastLevelSet::
operator= (const FastLevelSet &f)
{
	for (int i = cube(Ncoarse)-1; i>= 0;  --i) delete[] coarse_raw[i];
	delete[] coarse_raw;
	delete[] inside_raw;
	delete[] status_raw;
	delete[] temp_block;

	N = f.N;
	Ncoarse = f.Ncoarse;
	Nsub = f.Nsub;
	coarse_raw = new (float*)[cube(Ncoarse)];
	dx = f.dx;
	coarsedx = f.coarsedx;
	overdx = f.overdx;
	overcoarsedx = f.overcoarsedx;
	inside_raw = new bool[cube(Ncoarse)];
	status_counter = 0;
	status_raw = new unsigned int[cube(Nsub+1)];
	memset (status_raw, 0, sizeof(unsigned int)*cube(Nsub+1));
	temp_block = new float[cube(Nsub+1)];

	for (int i = cube(Ncoarse)-1; i >= 0; --i) {
		if (f.coarse_raw[i]) {
			coarse_raw[i] = new float[cube(Nsub+1)];
			memcpy (coarse_raw[i], f.coarse_raw[i], cube(Nsub+1)*sizeof(float));
		} else
			coarse_raw[i] = 0;
		inside_raw[i] = f.inside_raw[i];
	}
	return *this;
}


//============================== glRender stuff ================================

// bound on two-norm of upper 3x3 submatrix
static float matrixbound (float a[16])
{
	return sqrt(max3(flabs(a[0])+flabs(a[1])+flabs(a[2]),flabs(a[4])+flabs(a[5])+flabs(a[6]),flabs(a[8])+flabs(a[9])+flabs(a[10]))
	           *max3(flabs(a[0])+flabs(a[4])+flabs(a[8]),flabs(a[1])+flabs(a[5])+flabs(a[9]),flabs(a[2])+flabs(a[6])+flabs(a[10])));
}


void FastLevelSet::
glRender (void) const
{
	glEnable(GL_NORMALIZE);

	float transform[16];
	{
		GLfloat p[16], m[16];
		glGetFloatv (GL_PROJECTION_MATRIX, p);
		glGetFloatv (GL_MODELVIEW_MATRIX, m);
		transform[0] = p[0]*m[0]+p[4]*m[1]+p[8]*m[2]+p[12]*m[3];
		transform[1] = p[1]*m[0]+p[5]*m[1]+p[9]*m[2]+p[13]*m[3];
		transform[2] = p[2]*m[0]+p[6]*m[1]+p[10]*m[2]+p[14]*m[3];
		transform[3] = p[3]*m[0]+p[7]*m[1]+p[11]*m[2]+p[15]*m[3];
		transform[4] = p[0]*m[4]+p[4]*m[5]+p[8]*m[6]+p[12]*m[7];
		transform[5] = p[1]*m[4]+p[5]*m[5]+p[9]*m[6]+p[13]*m[7];
		transform[6] = p[2]*m[4]+p[6]*m[5]+p[10]*m[6]+p[14]*m[7];
		transform[7] = p[3]*m[4]+p[7]*m[5]+p[11]*m[6]+p[15]*m[7];
		transform[8] = p[0]*m[8]+p[4]*m[9]+p[8]*m[10]+p[12]*m[11];
		transform[9] = p[1]*m[8]+p[5]*m[9]+p[9]*m[10]+p[13]*m[11];
		transform[10] = p[2]*m[8]+p[6]*m[9]+p[10]*m[10]+p[14]*m[11];
		transform[11] = p[3]*m[8]+p[7]*m[9]+p[11]*m[10]+p[15]*m[11];
		transform[12] = p[0]*m[12]+p[4]*m[13]+p[8]*m[14]+p[12]*m[15];
		transform[13] = p[1]*m[12]+p[5]*m[13]+p[9]*m[14]+p[13]*m[15];
		transform[14] = p[2]*m[12]+p[6]*m[13]+p[10]*m[14]+p[14]*m[15];
		transform[15] = p[3]*m[12]+p[7]*m[13]+p[11]*m[14]+p[15]*m[15];
	}
	float m3 = matrixbound(transform);

	float pixmag;
	{
		GLint viewport[4];
		glGetIntegerv (GL_VIEWPORT, viewport);
		pixmag = (float) ((viewport[2] > viewport[3]) ? viewport[2] : viewport[3]);
	}

	enum {DRAWING_NOTHING, DRAWING_TRIANGLES, DRAWING_POINTS} drawmode = DRAWING_NOTHING;
	int pointsize = 1;
	float xyz[3], overw, screenxyz[3], voxelsize, voxelscreensize, finescreensize, *block;
	for (int k = 0; k < Ncoarse; ++k) {
		xyz[2] = coarsedx*(k+.5);
		for (int j = 0; j < Ncoarse; ++j) {
			xyz[1] = coarsedx*(j+.5);
			for (int i = 0; i < Ncoarse; ++i) {
				if (coarse(i,j,k)) {
					xyz[0] = coarsedx*(i+.5);
					overw = 1./(transform[3]*xyz[0] + transform[7]*xyz[1] + transform[11]*xyz[2] + transform[15]);
					voxelsize = .5*coarsedx*m3*overw;
					screenxyz[0] = (transform[0]*xyz[0] + transform[4]*xyz[1] + transform[8]*xyz[2] + transform[12])*overw;
					if (screenxyz[0] < -1-voxelsize || screenxyz[0] > 1+voxelsize) continue;
					screenxyz[1] = (transform[1]*xyz[0] + transform[5]*xyz[1] + transform[9]*xyz[2] + transform[13])*overw;
					if (screenxyz[1] < -1-voxelsize || screenxyz[1] > 1+voxelsize) continue;
					screenxyz[2] = (transform[2]*xyz[0] + transform[6]*xyz[1] + transform[10]*xyz[2] + transform[14])*overw;
					if (screenxyz[2] < -voxelsize || screenxyz[2] > 1+voxelsize) continue;
					// now we have to decide how to draw this voxel
					block = coarse(i,j,k);
					finescreensize = .6*pixmag*dx*m3*overw;
					if (true || finescreensize > 16) {
						if (drawmode == DRAWING_POINTS) glEnd();
						if (drawmode != DRAWING_TRIANGLES) {glBegin (GL_TRIANGLES); drawmode = DRAWING_TRIANGLES;}
						glRenderFineGreat (i, j, k);
					} else if (finescreensize > 12) {
						if (drawmode == DRAWING_POINTS) glEnd();
						if (drawmode != DRAWING_TRIANGLES) {glBegin (GL_TRIANGLES); drawmode = DRAWING_TRIANGLES;}
						glRenderFineGood (i, j, k);
					} else {
						voxelscreensize = pixmag*voxelsize;
						if (voxelscreensize > 32) {
							if (drawmode == DRAWING_TRIANGLES) glEnd();
							if (drawmode == DRAWING_POINTS && (pointsize < finescreensize || finescreensize < 5*pointsize/4)) {
								glEnd(); drawmode = DRAWING_NOTHING;
							}
							if (drawmode != DRAWING_POINTS) {
								pointsize = (int) finescreensize;
								glPointSize (pointsize);
								glBegin (GL_POINTS);
								drawmode = DRAWING_POINTS;
							}
							glRenderFine (i, j, k);
						} else {
							if (drawmode == DRAWING_TRIANGLES) glEnd();
							if (drawmode == DRAWING_POINTS && (pointsize < voxelscreensize || voxelscreensize < 5*pointsize/4)) {
								glEnd(); drawmode = DRAWING_NOTHING;
							}
							if (drawmode != DRAWING_POINTS) {
								pointsize = (int) voxelscreensize;
								glPointSize (pointsize);
								glBegin (GL_POINTS);
								drawmode = DRAWING_POINTS;
							}
							glRenderCoarse (i, j, k);
						}
					}
				}
			}
		}
	}
	if (drawmode != DRAWING_NOTHING) glEnd();
}


void FastLevelSet::
glRenderCoarse (const int i, const int j, const int k) const
{
	float *block = coarse(i,j,k);
	int mid = Nsub/2, lo = mid-1;
	float scale = .5*overdx;
	if (lo < 0) {lo = 0; scale = overdx;}
	float normal[3];
	normal[0] = scale * (fine(block,mid+1,mid,mid) - fine(block,lo,mid,mid));
	normal[1] = scale * (fine(block,mid,mid+1,mid) - fine(block,mid,lo,mid));
	normal[2] = scale * (fine(block,mid,mid,mid+1) - fine(block,mid,mid,lo));
	glNormal3f (normal[0], normal[1], normal[2]);
	divideByNormSquared(normal);
	float phi = fine(block,mid,mid,mid);
	glVertex3f (coarsedx*(i+.5)-phi*normal[0],coarsedx*(j+.5)-phi*normal[0],coarsedx*(k+.5)-phi*normal[0]);
}


void FastLevelSet::
glRenderFine (const int i, const int j, const int k) const
{
	float *block = coarse(i,j,k);
	float normal[3];
	for (int fk = 1; fk <= Nsub; ++fk)
		for (int fj = 1; fj <= Nsub; ++fj)
			for (int fi = 1; fi <= Nsub; ++fi) {
				float p111 = fine(block,fi,fj,fk),
					  p011 = fine(block,fi-1,fj,fk),
					  p101 = fine(block,fi,fj-1,fk),
					  p110 = fine(block,fi,fj,fk-1);
				if ((p111<=0 && (p011>0||p101>0||p110>0)) || (p111>0 && (p011<=0||p101<=0||p110<=0))) {
					normal[0] = overdx*(p111-p011); normal[1] = overdx*(p111-p101); normal[2] = overdx*(p111-p110);
					glNormal3f (normal[0], normal[1], normal[2]);
					glVertex3f (coarsedx*i+dx*(fi+.5), coarsedx*j+dx*(fj+.5), coarsedx*k+dx*(fk+.5));
				}
			}
}


static void rendertetgood (float x1, float y1, float z1, float p1,
						   float x2, float y2, float z2, float p2,
						   float x3, float y3, float z3, float p3,
						   float x4, float y4, float z4, float p4)
{
	if (p4 <= 0 && p1>0 && p2>0 && p3>0) {
		glVertex3f(x1,y1,z1); glVertex3f(x3,y3,z3); glVertex3f(x2,y2,z2);
	} else if (p3 <= 0 && p1>0 && p2>0 && p4>0) {
		glVertex3f(x1,y1,z1); glVertex3f(x2,y2,z2); glVertex3f(x4,y4,z4);
	} else if (p2 <= 0 && p1>0 && p3>0 && p4>0) {
		glVertex3f(x1,y1,z1); glVertex3f(x4,y4,z4); glVertex3f(x3,y3,z3);
	} else if (p1 <= 0 && p2>0 && p3>0 && p4>0) {
		glVertex3f(x2,y2,z2); glVertex3f(x3,y3,z3); glVertex3f(x4,y4,z4);
	}
}


void FastLevelSet::
glRenderFineGood (const int i, const int j, const int k) const
{
	float *block = coarse(i,j,k);
	float normal[3];
	for (int fk = 1; fk <= Nsub; ++fk) {
		for (int fj = 1; fj <= Nsub; ++fj) {
			for (int fi = 1; fi <= Nsub; ++fi) {
				float p111 = fine(block,fi,fj,fk),
					  p011 = fine(block,fi-1,fj,fk),
					  p101 = fine(block,fi,fj-1,fk),
					  p110 = fine(block,fi,fj,fk-1),
					  p100 = fine(block,fi,fj-1,fk-1),
					  p010 = fine(block,fi-1,fj,fk-1),
					  p001 = fine(block,fi-1,fj-1,fk),
					  p000 = fine(block,fi-1,fj-1,fk-1);
				bool pos = p111>0 || p011>0 || p101>0 || p110>0 || p100>0 || p010>0 || p001>0 || p000>0;
				bool neg = p111<=0 || p011<=0 || p101<=0 || p110<=0 || p100<=0 || p010<=0 || p001<=0 || p000<=0;
				if (!pos || !neg) continue;
				normal[0] = overdx*((p100+p101+p110+p111)-(p000+p001+p010+p011));
				normal[1] = overdx*((p010+p011+p110+p111)-(p000+p001+p100+p101));
				normal[2] = overdx*((p001+p011+p101+p111)-(p000+p010+p100+p110));
				glNormal3f (normal[0], normal[1], normal[2]);
				float x = coarsedx*i+dx*fi, y = coarsedx*j+dx*fj, z = coarsedx*k+dx*fk;
				if ((fi+fj+fk+Nsub*(i+j+k))%2) {
					rendertetgood (x-dx, y-dx, z-dx, p000, x, y-dx, z-dx, p100,
								   x-dx, y, z-dx, p010, x-dx, y-dx, z, p001);
					rendertetgood (x, y-dx, z-dx, p100, x, y, z-dx, p110,
							       x-dx, y, z-dx, p010, x, y, z, p111);
					rendertetgood (x, y-dx, z-dx, p100, x, y, z, p111,
							       x-dx, y-dx, z, p001, x, y-dx, z, p101);
					rendertetgood (x-dx, y-dx, z, p001, x, y, z, p111,
							       x-dx, y, z-dx, p010, x-dx, y, z, p011);
					rendertetgood (x, y-dx, z-dx, p100, x-dx, y-dx, z, p001,
							       x-dx, y, z-dx, p010, x, y, z, p111);
				} else {
					rendertetgood (x-dx, y-dx, z-dx, p000, x, y-dx, z-dx, p100,
							       x, y, z-dx, p110, x, y-dx, z, p101);
					rendertetgood (x-dx, y-dx, z-dx, p000, x, y, z-dx, p110,
							       x-dx, y, z, p011, x, y-dx, z, p101);
					rendertetgood (x-dx, y-dx, z-dx, p000, x, y-dx, z, p101,
							       x-dx, y, z, p011, x-dx, y-dx, z, p001);
					rendertetgood (x-dx, y-dx, z-dx, p000, x, y, z-dx, p110,
							       x-dx, y, z-dx, p010, x-dx, y, z, p011);
					rendertetgood (x, y-dx, z, p101, x, y, z-dx, p110,
							       x-dx, y, z, p011, x, y, z, p111);
				}
			}
		}
	}
}


static void rendertetgreat (float x1, float y1, float z1, float p1, float n1[3],
							float x2, float y2, float z2, float p2, float n2[3],
							float x3, float y3, float z3, float p3, float n3[3],
							float x4, float y4, float z4, float p4, float n4[3])
{
	if (p4 <= 0 && p1>0 && p2>0 && p3>0) {
		glNormal3f (n1[0], n1[1], n1[2]); glVertex3f(x1-p1*n1[0],y1-p1*n1[1],z1-p1*n1[2]);
		glNormal3f (n3[0], n3[1], n3[2]); glVertex3f(x3-p3*n3[0],y3-p3*n3[1],z3-p3*n3[2]);
		glNormal3f (n2[0], n2[1], n2[2]); glVertex3f(x2-p2*n2[0],y2-p2*n2[1],z2-p2*n2[2]);
	} else if (p3 <= 0 && p1>0 && p2>0 && p4>0) {
		glNormal3f (n1[0], n1[1], n1[2]); glVertex3f(x1-p1*n1[0],y1-p1*n1[1],z1-p1*n1[2]);
		glNormal3f (n2[0], n2[1], n2[2]); glVertex3f(x2-p2*n2[0],y2-p2*n2[1],z2-p2*n2[2]);
		glNormal3f (n4[0], n4[1], n4[2]); glVertex3f(x4-p4*n4[0],y4-p4*n4[1],z4-p4*n4[2]);
	} else if (p2 <= 0 && p1>0 && p3>0 && p4>0) {
		glNormal3f (n1[0], n1[1], n1[2]); glVertex3f(x1-p1*n1[0],y1-p1*n1[1],z1-p1*n1[2]);
		glNormal3f (n4[0], n4[1], n4[2]); glVertex3f(x4-p4*n4[0],y4-p4*n4[1],z4-p4*n4[2]);
		glNormal3f (n3[0], n3[1], n3[2]); glVertex3f(x3-p3*n3[0],y3-p3*n3[1],z3-p3*n3[2]);
	} else if (p1 <= 0 && p2>0 && p3>0 && p4>0) {
		glNormal3f (n2[0], n2[1], n2[2]); glVertex3f(x2-p2*n2[0],y2-p2*n2[1],z2-p2*n2[2]);
		glNormal3f (n3[0], n3[1], n3[2]); glVertex3f(x3-p3*n3[0],y3-p3*n3[1],z3-p3*n3[2]);
		glNormal3f (n4[0], n4[1], n4[2]); glVertex3f(x4-p4*n4[0],y4-p4*n4[1],z4-p4*n4[2]);
	}
}


void FastLevelSet::
glRenderFineGreat (const int i, const int j, const int k) const
{
	float *block = coarse(i,j,k);
	float *prevk, *nextk, *prevj, *nextj, *previ, *nexti;
	float mulpk, mulnk, mulpj, mulnj, mulpi, mulni;
	float n000[3], n001[3], n010[3], n011[3], n100[3], n101[3], n110[3], n111[3];
	for (int fk = 1; fk <= Nsub; ++fk) {
		mulpk = .5*overdx;
		if (fk > 1) 						prevk = block+(fk-2)*(Nsub+1)*(Nsub+1);
		else if (k > 0 && coarse(i,j,k-1)) 	prevk = coarse(i,j,k-1)+(Nsub-1)*(Nsub+1)*(Nsub+1);
		else 								{prevk = block; mulpk = overdx;}
		mulnk = .5*overdx;
		if (fk < Nsub) 								nextk = block+(fk+1)*(Nsub+1)*(Nsub+1);
		else if (k+1 < Ncoarse && coarse(i,j,k+1)) 	nextk = coarse(i,j,k+1)+(Nsub+1)*(Nsub+1);
		else 										{nextk = block+Nsub*(Nsub+1)*(Nsub+1); mulnk = overdx;}
		for (int fj = 1; fj <= Nsub; ++fj) {
			mulpj = .5*overdx;
			if (fj > 1) 						prevj = block+(fj-2)*(Nsub+1);
			else if (j > 0 && coarse(i,j-1,k)) 	prevj = coarse(i,j-1,k)+(Nsub-1)*(Nsub+1);
			else 								{prevj = block; mulpj = overdx;}
			mulnj = .5*overdx;
			if (fj < Nsub) 								nextj = block+(fj+1)*(Nsub+1);
			else if (j+1 < Ncoarse && coarse(i,j+1,k)) 	nextj = coarse(i,j+1,k)+(Nsub+1);
			else 										{nextj = block+Nsub*(Nsub+1); mulnj = overdx;}
			for (int fi = 1; fi <= Nsub; ++fi) {
				float p111 = fine(block,fi,fj,fk),
					  p011 = fine(block,fi-1,fj,fk),
					  p101 = fine(block,fi,fj-1,fk),
					  p110 = fine(block,fi,fj,fk-1),
					  p100 = fine(block,fi,fj-1,fk-1),
					  p010 = fine(block,fi-1,fj,fk-1),
					  p001 = fine(block,fi-1,fj-1,fk),
					  p000 = fine(block,fi-1,fj-1,fk-1);
				bool pos = p111>0 || p011>0 || p101>0 || p110>0 || p100>0 || p010>0 || p001>0 || p000>0;
				bool neg = p111<=0 || p011<=0 || p101<=0 || p110<=0 || p100<=0 || p010<=0 || p001<=0 || p000<=0;
				if (!pos || !neg) continue;
				mulpi = .5*overdx;
				if (fi > 1) 						previ = block+(fi-2);
				else if (i > 0 && coarse(i-1,j,k)) 	previ = coarse(i-1,j,k)+Nsub-1;
				else 								{previ = block; mulpi = overdx;}
				mulni = .5*overdx;
				if (fi < Nsub) 								nexti = block+(fi+1);
				else if (i+1 < Ncoarse && coarse(i+1,j,k)) 	nexti = coarse(i+1,j,k)+1;
				else 										{nexti = block+Nsub; mulni = overdx;}
				n000[0] = mulpi * (p100-previ[(Nsub+1)*((fj-1)+(Nsub+1)*(fk-1))]);
				n000[1] = mulpj * (p010-prevj[(fi-1)+(Nsub+1)*(Nsub+1)*(fk-1)]);
				n000[2] = mulpk * (p001-prevk[(fi-1)+(Nsub+1)*(fj-1)]);
				divideByNormSquared (n000);
				n001[0] = mulpi * (p101-previ[(Nsub+1)*((fj-1)+(Nsub+1)*fk)]);
				n001[1] = mulpj * (p011-prevj[(fi-1)+(Nsub+1)*(Nsub+1)*fk]);
				n001[2] = mulnk * (nextk[(fi-1)+(Nsub+1)*(fj-1)]-p000);
				divideByNormSquared (n001);
				n010[0] = mulpi * (p110-previ[(Nsub+1)*(fj+(Nsub+1)*(fk-1))]);
				n010[1] = mulnj * (nextj[(fi-1)+(Nsub+1)*(Nsub+1)*(fk-1)]-p000);
				n010[2] = mulpk * (p011-prevk[(fi-1)+(Nsub+1)*fj]);
				divideByNormSquared (n010);
				n011[0] = mulpi * (p111-previ[(Nsub+1)*(fj+(Nsub+1)*fk)]);
				n011[1] = mulnj * (nextj[(fi-1)+(Nsub+1)*(Nsub+1)*fk]-p001);
				n011[2] = mulnk * (nextk[(fi-1)+(Nsub+1)*fj]-p010);
				divideByNormSquared (n011);
				n100[0] = mulni * (nexti[(Nsub+1)*((fj-1)+(Nsub+1)*(fk-1))]-p000);
				n100[1] = mulpj * (p110-prevj[fi+(Nsub+1)*(Nsub+1)*(fk-1)]);
				n100[2] = mulpk * (p101-prevk[fi+(Nsub+1)*(fj-1)]);
				divideByNormSquared (n100);
				n101[0] = mulni * (nexti[(Nsub+1)*((fj-1)+(Nsub+1)*fk)]-p001);
				n101[1] = mulpj * (p111-prevj[fi+(Nsub+1)*(Nsub+1)*fk]);
				n101[2] = mulnk * (nextk[fi+(Nsub+1)*(fj-1)]-p100);
				divideByNormSquared (n101);
				n110[0] = mulni * (nexti[(Nsub+1)*(fj+(Nsub+1)*(fk-1))]-p010);
				n110[1] = mulnj * (nextj[fi+(Nsub+1)*(Nsub+1)*(fk-1)]-p100);
				n110[2] = mulpk * (p111-prevk[fi+(Nsub+1)*fj]);
				divideByNormSquared (n110);
				n111[0] = mulni * (nexti[(Nsub+1)*(fj+(Nsub+1)*fk)]-p011);
				n111[1] = mulnj * (nextj[fi+(Nsub+1)*(Nsub+1)*fk]-p101);
				n111[2] = mulnk * (nextk[fi+(Nsub+1)*fj]-p110);
				divideByNormSquared (n111);
				float x = coarsedx*i+dx*fi, y = coarsedx*j+dx*fj, z = coarsedx*k+dx*fk;
				if ((fi+fj+fk+Nsub*(i+j+k))%2) {
					rendertetgreat (x-dx, y-dx, z-dx, p000, n000, x, y-dx, z-dx, p100, n100,
					           		x-dx, y, z-dx, p010, n010, x-dx, y-dx, z, p001, n001);
					rendertetgreat (x, y-dx, z-dx, p100, n100, x, y, z-dx, p110, n110,
							   		x-dx, y, z-dx, p010, n010, x, y, z, p111, n111);
					rendertetgreat (x, y-dx, z-dx, p100, n100, x, y, z, p111, n111,
							   		x-dx, y-dx, z, p001, n001, x, y-dx, z, p101, n101);
					rendertetgreat (x-dx, y-dx, z, p001, n001, x, y, z, p111, n111,
							   		x-dx, y, z-dx, p010, n010, x-dx, y, z, p011, n011);
					rendertetgreat (x, y-dx, z-dx, p100, n100, x-dx, y-dx, z, p001, n001,
							   		x-dx, y, z-dx, p010, n010, x, y, z, p111, n111);
				} else {
					rendertetgreat (x-dx, y-dx, z-dx, p000, n000, x, y-dx, z-dx, p100, n100,
							   		x, y, z-dx, p110, n110, x, y-dx, z, p101, n101);
					rendertetgreat (x-dx, y-dx, z-dx, p000, n000, x, y, z-dx, p110, n110,
							   		x-dx, y, z, p011, n011, x, y-dx, z, p101, n101);
					rendertetgreat (x-dx, y-dx, z-dx, p000, n000, x, y-dx, z, p101, n101,
							   		x-dx, y, z, p011, n011, x-dx, y-dx, z, p001, n001);
					rendertetgreat (x-dx, y-dx, z-dx, p000, n000, x, y, z-dx, p110, n110,
							   		x-dx, y, z-dx, p010, n010, x-dx, y, z, p011, n011);
					rendertetgreat (x, y-dx, z, p101, n101, x, y, z-dx, p110, n110,
							   		x-dx, y, z, p011, n011, x, y, z, p111, n111);
				}
			}
		}
	}
}


//================================= ray tracing =================================

bool FastLevelSet::
traceRay (const float orig[3], const float dir[3], float p[3]) const
{
	// begin by clipping ray against cube
	float t = 0, maxt = 1e38, overdirx = 1e38, overdiry = 1e38, overdirz = 1e38, near, far;
	if (dir[0] != 0) {
		overdirx = 1./dir[0];
		near = -orig[0]*overdirx;
		far = (1-orig[0])*overdirx;
		if (near > far) {float temp = far; far = near; near = temp;}
		if (near > 0) t = near; else t = 0;
		maxt = far;
		if (t >= maxt) return false;
	} else if (orig[0] <= 0 || orig[0] >= 1) return false;
	if (dir[1] != 0) {
		overdiry = 1./dir[1];
		near = -orig[1]*overdiry;
		far = (1-orig[1])*overdiry;
		if (near > far) {float temp = far; far = near; near = temp;}
		if (near > t) t = near;
		if (far < maxt) maxt = far;
		if (t >= maxt) return false;
	} else if (orig[1] <= 0 || orig[1] >= 1) return false;
	if (dir[2] != 0) {
		overdirz = 1./dir[2];
		near = -orig[2]*overdirz;
		far = (1-orig[2])*overdirz;
		if (near > far) {float temp = far; far = near; near = temp;}
		if (near > t) t = near;
		if (far < maxt) maxt = far;
		if (t >= maxt) return false;
	} else if (orig[2] <= 0 || orig[2] >= 1) return false;

	// trace through the coarse grid, doing a fine trace whenever we hit a populated block
	t += .001*dx;
	p[0] = orig[0]+t*dir[0];
	p[1] = orig[1]+t*dir[1];
	p[2] = orig[2]+t*dir[2];
	int i = (int)(p[0]*overcoarsedx), j = (int)(p[1]*overcoarsedx), k = (int)(p[2]*overcoarsedx);
	if (dir[2] >= 0) {
		if (dir[1] >= 0) {
			if (dir[0] >= 0)
				while (i < Ncoarse && j < Ncoarse && k < Ncoarse) {
					float tx = ((i+1)*coarsedx-p[0])*overdirx,
						  ty = ((j+1)*coarsedx-p[1])*overdiry,
						  tz = ((k+1)*coarsedx-p[2])*overdirz;
					if (coarse(i,j,k) && traceRayFine (i, j, k, dir, min3(tx,ty,tz), p)) return true;
					if (tx <= ty && tx <= tz) {++i; p[0] = i*coarsedx; p[1] += tx*dir[1]; p[2] += tx*dir[2];}
					else if (ty <= tx && ty <= tz) {++j; p[0] += ty*dir[0]; p[1] = j*coarsedx; p[2] += ty*dir[2];}
					else {++k; p[0] += tz*dir[0]; p[1] += tz*dir[1]; p[2] = k*coarsedx;}
				}
			else 
				while (i >= 0 && j < Ncoarse && k < Ncoarse) {
					float tx = (i*coarsedx-p[0])*overdirx,
						  ty = ((j+1)*coarsedx-p[1])*overdiry,
						  tz = ((k+1)*coarsedx-p[2])*overdirz;
					if (coarse(i,j,k) && traceRayFine (i, j, k, dir, min3(tx,ty,tz), p)) return true;
					if (tx <= ty && tx <= tz) {p[0] = i*coarsedx; p[1] += tx*dir[1]; p[2] += tx*dir[2]; --i;}
					else if (ty <= tx && ty <= tz) {++j; p[0] += ty*dir[0]; p[1] = j*coarsedx; p[2] += ty*dir[2];}
					else {++k; p[0] += tz*dir[0]; p[1] += tz*dir[1]; p[2] = k*coarsedx;}
				}
		} else {
			if (dir[0] >= 0)
				while (i < Ncoarse && j >= 0 && k < Ncoarse) {
					float tx = ((i+1)*coarsedx-p[0])*overdirx,
						  ty = (j*coarsedx-p[1])*overdiry,
						  tz = ((k+1)*coarsedx-p[2])*overdirz;
					if (coarse(i,j,k) && traceRayFine (i, j, k, dir, min3(tx,ty,tz), p)) return true;
					if (tx <= ty && tx <= tz) {++i; p[0] = i*coarsedx; p[1] += tx*dir[1]; p[2] += tx*dir[2];}
					else if (ty <= tx && ty <= tz) {p[0] += ty*dir[0]; p[1] = j*coarsedx; p[2] += ty*dir[2]; --j;}
					else {++k; p[0] += tz*dir[0]; p[1] += tz*dir[1]; p[2] = k*coarsedx;}
				}
			else
				while (i >= 0 && j >= 0 && k < Ncoarse) {
					float tx = (i*coarsedx-p[0])*overdirx,
						  ty = (j*coarsedx-p[1])*overdiry,
						  tz = ((k+1)*coarsedx-p[2])*overdirz;
					if (coarse(i,j,k) && traceRayFine (i, j, k, dir, min3(tx,ty,tz), p)) return true;
					if (tx <= ty && tx <= tz) {p[0] = i*coarsedx; p[1] += tx*dir[1]; p[2] += tx*dir[2]; --i;}
					else if (ty <= tx && ty <= tz) {p[0] += ty*dir[0]; p[1] = j*coarsedx; p[2] += ty*dir[2]; --j;}
					else {++k; p[0] += tz*dir[0]; p[1] += tz*dir[1]; p[2] = k*coarsedx;}
				}
		}
	} else {
		if (dir[1] >= 0) {
			if (dir[0] >= 0)
				while (i < Ncoarse && j < Ncoarse && k >= 0) {
					float tx = ((i+1)*coarsedx-p[0])*overdirx,
						  ty = ((j+1)*coarsedx-p[1])*overdiry,
						  tz = (k*coarsedx-p[2])*overdirz;
					if (coarse(i,j,k) && traceRayFine (i, j, k, dir, min3(tx,ty,tz), p)) return true;
					if (tx <= ty && tx <= tz) {++i; p[0] = i*coarsedx; p[1] += tx*dir[1]; p[2] += tx*dir[2];}
					else if (ty <= tx && ty <= tz) {++j; p[0] += ty*dir[0]; p[1] = j*coarsedx; p[2] += ty*dir[2];}
					else {p[0] += tz*dir[0]; p[1] += tz*dir[1]; p[2] = k*coarsedx; --k;}
				}
			else
				while (i >= 0 && j < Ncoarse && k >= 0) {
					float tx = (i*coarsedx-p[0])*overdirx,
						  ty = ((j+1)*coarsedx-p[1])*overdiry,
						  tz = (k*coarsedx-p[2])*overdirz;
					if (coarse(i,j,k) && traceRayFine (i, j, k, dir, min3(tx,ty,tz), p)) return true;
					if (tx <= ty && tx <= tz) {p[0] = i*coarsedx; p[1] += tx*dir[1]; p[2] += tx*dir[2]; --i;}
					else if (ty <= tx && ty <= tz) {++j; p[0] += ty*dir[0]; p[1] = j*coarsedx; p[2] += ty*dir[2];}
					else {p[0] += tz*dir[0]; p[1] += tz*dir[1]; p[2] = k*coarsedx; --k;}
				}
		} else {
			if (dir[0] >= 0)
				while (i < Ncoarse && j >= 0 && k >= 0) {
					float tx = ((i+1)*coarsedx-p[0])*overdirx,
						  ty = (j*coarsedx-p[1])*overdiry,
						  tz = (k*coarsedx-p[2])*overdirz;
					if (coarse(i,j,k) && traceRayFine (i, j, k, dir, min3(tx,ty,tz), p)) return true;
					if (tx <= ty && tx <= tz) {++i; p[0] = i*coarsedx; p[1] += tx*dir[1]; p[2] += tx*dir[2];}
					else if (ty <= tx && ty <= tz) {p[0] += ty*dir[0]; p[1] = j*coarsedx; p[2] += ty*dir[2]; --j;}
					else {p[0] += tz*dir[0]; p[1] += tz*dir[1]; p[2] = k*coarsedx; --k;}
				}
			else 
				while (i >= 0 && j >= 0 && k >= 0) {
					float tx = (i*coarsedx-p[0])*overdirx,
						  ty = (j*coarsedx-p[1])*overdiry,
						  tz = (k*coarsedx-p[2])*overdirz;
					if (coarse(i,j,k) && traceRayFine (i, j, k, dir, min3(tx,ty,tz), p)) return true;
					if (tx <= ty && tx <= tz) {p[0] = i*coarsedx; p[1] += tx*dir[1]; p[2] += tx*dir[2]; --i;}
					else if (ty <= tx && ty <= tz) {p[0] += ty*dir[0]; p[1] = j*coarsedx; p[2] += ty*dir[2]; --j;}
					else {p[0] += tz*dir[0]; p[1] += tz*dir[1]; p[2] = k*coarsedx; --k;}
				}
		}
	}
	return false;
}


bool FastLevelSet::
traceRayFine (int i, int j, int k, const float dir[3], float maxt, float p[3]) const
{
	const float *block = coarse(i,j,k);
	float oldval = dx;
	float t = 0;
	do {
		float val = evalFine (block, ((p[0]+t*dir[0])*overcoarsedx-i)*Nsub,
								     ((p[1]+t*dir[1])*overcoarsedx-j)*Nsub,
									 ((p[2]+t*dir[2])*overcoarsedx-k)*Nsub);
		if (val < 0.01*dx) { // if we almost hit a zero or the sign changes to negative
			val = (val*oldval)/(oldval-val); // estimate where the zero is (secant rule)
			t += val;
			if (t < 0) t = 0; else if (t > maxt) t = maxt;
			p[0] += t*dir[0]; p[1] += t*dir[1]; p[2] += t*dir[2];
			return true;
		}
		if (val > dx) t += dx; else t += val;  // don't trust that it really is signed distance
		oldval = val;
	} while (t < maxt);
	return false;
}


float FastLevelSet::
evalFine (const float *block, float fx, float fy, float fz) const
{
	int i = (int)fx, j = (int)fy, k = (int)fz;
	if (i < 0) i = 0; else if (i >= Nsub) i = Nsub-1;
	if (j < 0) j = 0; else if (j >= Nsub) j = Nsub-1;
	if (k < 0) k = 0; else if (k >= Nsub) k = Nsub-1;
	float p000 = fine(block,i,j,k), p100 = fine(block,i+1,j,k),
		  p010 = fine(block,i,j+1,k), p001 = fine(block,i,j,k+1),
		  p011 = fine(block,i,j+1,k+1), p101 = fine(block,i+1,j,k+1),
		  p110 = fine(block,i+1,j+1,k), p111 = fine(block,i+1,j+1,k+1);
	float a = fx-i, b = fy-j, c = fz-k;
	return (1-a)*((1-b)*((1-c)*p000+c*p001) + b*((1-c)*p010+c*p011))
			  +a*((1-b)*((1-c)*p100+c*p101) + b*((1-c)*p110+c*p111));
}


//================================ local modifications ==========================

void FastLevelSet::
localAdd (const float center[3], float radius, float amount)
{
	for (int k = (int)((center[2]-radius)*overcoarsedx); k < (center[2]+radius)*overcoarsedx; ++k) {
		if (k < 0 || k >= Ncoarse) continue;
		for (int j = (int)((center[1]-radius)*overcoarsedx); j < (center[1]+radius)*overcoarsedx; ++j) {
			if (j < 0 || j >= Ncoarse) continue;
			for (int i = (int)((center[0]-radius)*overcoarsedx); i < (center[0]+radius)*overcoarsedx; ++i) {
				if (i < 0 || i >= Ncoarse) continue;
				float *block = coarse(i,j,k);
				if (!block) continue;
				for (int fk = 0; fk <= Nsub; ++fk) {
					float z = coarsedx*k+dx*fk;
					for (int fj = 0; fj <= Nsub; ++fj) {
						float y = coarsedx*j+dx*fj;
						for (int fi = 0; fi <= Nsub; ++fi) {
							float x = coarsedx*i+dx*fi;
							float scale = 1 - 1.001*magnitude(center[0]-x, center[1]-y, center[2]-z)/radius;
							if (scale > 0) {
								//scale = 1+scale*scale*(-3+2*scale);
								fine(block,fi,fj,fk) += amount*scale;
							}
						}
					}
				}
			}
		}
	}
	for (int k = (int)((center[2]-radius)*overcoarsedx); k < (center[2]+radius)*overcoarsedx; ++k) {
		if (k < 0 || k >= Ncoarse) continue;
		for (int j = (int)((center[1]-radius)*overcoarsedx); j < (center[1]+radius)*overcoarsedx; ++j) {
			if (j < 0 || j >= Ncoarse) continue;
			for (int i = (int)((center[0]-radius)*overcoarsedx); i < (center[0]+radius)*overcoarsedx; ++i) {
				if (i < 0 || i >= Ncoarse) continue;
				if (!coarse(i,j,k)) continue;
				localRedistance (i, j, k);
			}
		}
	}
}


void FastLevelSet::
localSmooth (const float center[3], float radius, float amount)
{
	int i, j, k, fi, fj, fk;
	float *prevk, *nextk, *prevj, *nextj, *previ, *nexti;
	for (k = (int)((center[2]-radius)*overcoarsedx); k < (center[2]+radius)*overcoarsedx; ++k) {
		if (k < 0 || k >= Ncoarse) continue;
		for (j = (int)((center[1]-radius)*overcoarsedx); j < (center[1]+radius)*overcoarsedx; ++j) {
			if (j < 0 || j >= Ncoarse) continue;
			for (i = (int)((center[0]-radius)*overcoarsedx); i < (center[0]+radius)*overcoarsedx; ++i) {
				if (i < 0 || i >= Ncoarse) continue;
				float *block = coarse(i,j,k);
				if (!block) continue;
				for (fk = 0; fk <= Nsub; ++fk) {
					float z = coarsedx*k+dx*fk;
					if (fk > 1) 						prevk = block+(fk-2)*(Nsub+1)*(Nsub+1);
					else if (k > 0 && coarse(i,j,k-1)) 	prevk = coarse(i,j,k-1)+(Nsub-1)*(Nsub+1)*(Nsub+1);
					else 								prevk = block;
					if (fk < Nsub) 								nextk = block+(fk+1)*(Nsub+1)*(Nsub+1);
					else if (k+1 < Ncoarse && coarse(i,j,k+1)) 	nextk = coarse(i,j,k+1)+(Nsub+1)*(Nsub+1);
					else 										nextk = block+Nsub*(Nsub+1)*(Nsub+1);
					for (fj = 0; fj <= Nsub; ++fj) {
						float y = coarsedx*j+dx*fj;
						if (fj > 1) 						prevj = block+(fj-2)*(Nsub+1);
						else if (j > 0 && coarse(i,j-1,k)) 	prevj = coarse(i,j-1,k)+(Nsub-1)*(Nsub+1);
						else 								prevj = block;
						if (fj < Nsub) 								nextj = block+(fj+1)*(Nsub+1);
						else if (j+1 < Ncoarse && coarse(i,j+1,k)) 	nextj = coarse(i,j+1,k)+(Nsub+1);
						else 										nextj = block+Nsub*(Nsub+1);
						for (fi = 0; fi <= Nsub; ++fi) {
							float x = coarsedx*i+dx*fi;
							float scale = 1 - 1.001*magnitude(center[0]-x, center[1]-y, center[2]-z)/radius;
							if (scale > 0) {
								//scale = 1+scale*scale*(-3+2*scale);
								if (fi > 1) 						previ = block+(fi-2);
								else if (i > 0 && coarse(i-1,j,k)) 	previ = coarse(i-1,j,k)+Nsub-1;
								else 								previ = block;
								if (fi < Nsub) 								nexti = block+(fi+1);
								else if (i+1 < Ncoarse && coarse(i+1,j,k)) 	nexti = coarse(i+1,j,k)+1;
								else 										nexti = block+Nsub;
								float avg = ( 6*fine(block,fi,fj,fk)
										     +previ[(Nsub+1)*(fj+(Nsub+1)*fk)] + nexti[(Nsub+1)*(fj+(Nsub+1)*fk)]
										     +prevj[fi+(Nsub+1)*(Nsub+1)*fk] + nextj[fi+(Nsub+1)*(Nsub+1)*fk]
										     +prevk[fi+(Nsub+1)*fj] + nextk[fi+(Nsub+1)*fj])/12.;
								scale *= amount*overdx;
								if (scale > 1) scale = 1; else if (scale < -.5) scale = -.5;
								fine(block,fi,fj,fk) += scale*(avg-fine(block,fi,fj,fk));
							}
						}
					}
				}
			}
		}
	}
	for (int k = (int)((center[2]-radius)*overcoarsedx); k < (center[2]+radius)*overcoarsedx; ++k) {
		if (k < 0 || k >= Ncoarse) continue;
		for (int j = (int)((center[1]-radius)*overcoarsedx); j < (center[1]+radius)*overcoarsedx; ++j) {
			if (j < 0 || j >= Ncoarse) continue;
			for (int i = (int)((center[0]-radius)*overcoarsedx); i < (center[0]+radius)*overcoarsedx; ++i) {
				if (i < 0 || i >= Ncoarse) continue;
				if (!coarse(i,j,k)) continue;
				localRedistance (i, j, k);
			}
		}
	}
}


void FastLevelSet::
enforceSharedValues (void)
{
	for (int k = 0; k < Ncoarse; ++k)
		for (int j = 0; j < Ncoarse; ++j)
			for (int i = 0; i < Ncoarse; ++i)
				if (coarse(i,j,k)) enforceSharedBlock (i, j, k);
}


void FastLevelSet::
redistance (void)
{
	bool containsBoundary;
	for (int k = 0; k < Ncoarse; ++k)
		for (int j = 0; j < Ncoarse; ++j)
			for (int i = 0; i < Ncoarse; ++i)
				if (coarse(i,j,k)) reinitBlock (i, j, k, containsBoundary);
}


static inline void share (float &a, float &b)
{ a = b = minmod(a,b); }


void FastLevelSet::
enforceSharedBlock (const int i, const int j, const int k)
{
	float *block = coarse(i,j,k), *nbr;
	// share common faces
	if (i > 0 && (nbr = coarse(i-1,j,k)) != 0)
		for (int fk = 0; fk <= Nsub; ++fk) for (int fj = 0; fj <= Nsub; ++fj)
			share(fine(nbr,Nsub,fj,fk), fine(block,0,fj,fk));
	if (i+1 < Ncoarse && (nbr = coarse(i+1,j,k)) != 0)
		for (int fk = 0; fk <= Nsub; ++fk) for (int fj = 0; fj <= Nsub; ++fj)
			share(fine(nbr,0,fj,fk), fine(block,Nsub,fj,fk));
	if (j > 0 && (nbr = coarse(i,j-1,k)) != 0)
		for (int fk = 0; fk <= Nsub; ++fk) for (int fi = 0; fi <= Nsub; ++fi)
			share(fine(nbr,fi,Nsub,fk), fine(block,fi,0,fk));
	if (j+1 < Ncoarse && (nbr = coarse(i,j+1,k)) != 0)
		for (int fk = 0; fk <= Nsub; ++fk) for (int fi = 0; fi <= Nsub; ++fi)
			share(fine(nbr,fi,0,fk), fine(block,fi,Nsub,fk));
	if (k > 0 && (nbr = coarse(i,j,k-1)) != 0)
		for (int fj = 0; fj <= Nsub; ++fj) for (int fi = 0; fi <= Nsub; ++fi)
			share(fine(nbr,fi,fj,Nsub), fine(block,fi,fj,0));
	if (k+1 < Ncoarse && (nbr = coarse(i,j,k+1)) != 0)
		for (int fj = 0; fj <= Nsub; ++fj) for (int fi = 0; fi <= Nsub; ++fi)
			share(fine(nbr,fi,fj,0), fine(block,fi,fj,Nsub));
	// share common edges along k
	if (i > 0 && j > 0 && (nbr = coarse(i-1,j-1,k)) != 0)
		for (int fk = 0; fk <= Nsub; ++fk)
			fine(nbr,Nsub,Nsub,fk) = fine(block,0,0,fk);
	if (i+1 < Ncoarse && j > 0 && (nbr = coarse(i+1,j-1,k)) != 0)
		for (int fk = 0; fk <= Nsub; ++fk)
			fine(nbr,0,Nsub,fk) = fine(block,Nsub,0,fk);
	if (i > 0 && j+1 < Ncoarse && (nbr = coarse(i-1,j+1,k)) != 0)
		for (int fk = 0; fk <= Nsub; ++fk)
			fine(nbr,Nsub,0,fk) = fine(block,0,Nsub,fk);
	if (i+1 < Ncoarse && j+1 < Ncoarse && (nbr = coarse(i+1,j+1,k)) != 0)
		for (int fk = 0; fk <= Nsub; ++fk)
			fine(nbr,0,0,fk) = fine(block,Nsub,Nsub,fk);
	// share common edges along j
	if (i > 0 && k > 0 && (nbr = coarse(i-1,j,k-1)) != 0)
		for (int fj = 0; fj <= Nsub; ++fj)
			fine(nbr,Nsub,fj,Nsub) = fine(block,0,fj,0);
	if (i+1 < Ncoarse && k > 0 && (nbr = coarse(i+1,j,k-1)) != 0)
		for (int fj = 0; fj <= Nsub; ++fj)
			fine(nbr,0,fj,Nsub) = fine(block,Nsub,fj,0);
	if (i > 0 && k+1 < Ncoarse && (nbr = coarse(i-1,j,k+1)) != 0)
		for (int fj = 0; fj <= Nsub; ++fj)
			fine(nbr,Nsub,fj,0) = fine(block,0,fj,Nsub);
	if (i+1 < Ncoarse && k+1 < Ncoarse && (nbr = coarse(i+1,j,k+1)) != 0)
		for (int fj = 0; fj <= Nsub; ++fj)
			fine(nbr,0,fj,0) = fine(block,Nsub,fj,Nsub);
	// share common edges along i
	if (j > 0 && k > 0 && (nbr = coarse(i,j-1,k-1)) != 0)
		for (int fi = 0; fi <= Nsub; ++fi)
			fine(nbr,fi,Nsub,Nsub) = fine(block,fi,0,0);
	if (j+1 < Ncoarse && k > 0 && (nbr = coarse(i,j+1,k-1)) != 0)
		for (int fi = 0; fi <= Nsub; ++fi)
			fine(nbr,fi,0,Nsub) = fine(block,fi,Nsub,0);
	if (j > 0 && k+1 < Ncoarse && (nbr = coarse(i,j-1,k+1)) != 0)
		for (int fi = 0; fi <= Nsub; ++fi)
			fine(nbr,fi,Nsub,0) = fine(block,fi,0,Nsub);
	if (j+1 < Ncoarse && k+1 < Ncoarse && (nbr = coarse(i,j+1,k+1)) != 0)
		for (int fi = 0; fi <= Nsub; ++fi)
			fine(nbr,fi,0,0) = fine(block,fi,Nsub,Nsub);
	// share common corners
	if (i > 0 && j > 0 && k > 0 && (nbr = coarse(i-1,j-1,k-1)) != 0)
		fine(nbr,Nsub,Nsub,Nsub) = fine(block,0,0,0);
	if (i > 0 && j > 0 && k+1 < Ncoarse && (nbr = coarse(i-1,j-1,k+1)) != 0)
		fine(nbr,Nsub,Nsub,0) = fine(block,0,0,Nsub);
	if (i > 0 && j+1 < Ncoarse && k > 0 && (nbr = coarse(i-1,j+1,k-1)) != 0)
		fine(nbr,Nsub,0,Nsub) = fine(block,0,Nsub,0);
	if (i > 0 && j+1 < Ncoarse && k+1 < Ncoarse && (nbr = coarse(i-1,j+1,k+1)) != 0)
		fine(nbr,Nsub,0,0) = fine(block,0,Nsub,Nsub);
	if (i+1 < Ncoarse && j > 0 && k > 0 && (nbr = coarse(i+1,j-1,k-1)) != 0)
		fine(nbr,0,Nsub,Nsub) = fine(block,Nsub,0,0);
	if (i+1 < Ncoarse && j > 0 && k+1 < Ncoarse && (nbr = coarse(i+1,j-1,k+1)) != 0)
		fine(nbr,0,Nsub,0) = fine(block,Nsub,0,Nsub);
	if (i+1 < Ncoarse && j+1 < Ncoarse && k > 0 && (nbr = coarse(i+1,j+1,k-1)) != 0)
		fine(nbr,0,0,Nsub) = fine(block,Nsub,Nsub,0);
	if (i+1 < Ncoarse && j+1 < Ncoarse && k+1 < Ncoarse && (nbr = coarse(i+1,j+1,k+1)) != 0)
		fine(nbr,0,0,0) = fine(block,Nsub,Nsub,Nsub);
}


void FastLevelSet::
localRedistance (const int i, const int j, const int k)
{
	bool containsBoundary;
	reinitBlock (i, j, k, containsBoundary);
	float *block = coarse(i,j,k), *nbr;
	if (containsBoundary) {
		if (i > 0 && (nbr = coarse(i-1,j,k)) == 0) {
			bool interfacepos = false, interfaceneg = false;
			for (int fk = 0; fk <= Nsub; ++fk) for (int fj = 0; fj <= Nsub; ++fj)
				if (fine(block,0,fj,fk) <= 0) interfaceneg = true; else interfacepos = true;
			if ((interfaceneg && !inside(i-1,j,k)) || (interfacepos && inside(i-1,j,k))) initBlock (i-1,j,k); // create it if so
		}
		if (i+1 < Ncoarse && (nbr = coarse(i+1,j,k)) == 0) {
			bool interfacepos = false, interfaceneg = false;
			for (int fk = 0; fk <= Nsub; ++fk) for (int fj = 0; fj <= Nsub; ++fj)
				if (fine(block,Nsub,fj,fk) <= 0) interfaceneg = true; else interfacepos = true;
			if ((interfaceneg && !inside(i+1,j,k)) || (interfacepos && inside(i+1,j,k))) initBlock (i+1,j,k); // create it if so
		}
		if (j > 0 && (nbr = coarse(i,j-1,k)) == 0) {
			bool interfacepos = false, interfaceneg = false;
			for (int fk = 0; fk <= Nsub; ++fk) for (int fi = 0; fi <= Nsub; ++fi)
				if (fine(block,fi,0,fk) <= 0) interfaceneg = true; else interfacepos = true;
			if ((interfaceneg && !inside(i,j-1,k)) || (interfacepos && inside(i,j-1,k))) initBlock (i,j-1,k); // create it if so
		}
		if (j+1 < Ncoarse && (nbr = coarse(i,j+1,k)) == 0) {
			bool interfacepos = false, interfaceneg = false;
			for (int fk = 0; fk <= Nsub; ++fk) for (int fi = 0; fi <= Nsub; ++fi)
				if (fine(block,fi,Nsub,fk) <= 0) interfaceneg = true; else interfacepos = true;
			if ((interfaceneg && !inside(i,j+1,k)) || (interfacepos && inside(i,j+1,k))) initBlock (i,j+1,k); // create it if so
		}
		if (k > 0 && (nbr = coarse(i,j,k-1)) == 0) {
			bool interfacepos = false, interfaceneg = false;
			for (int fj = 0; fj <= Nsub; ++fj) for (int fi = 0; fi <= Nsub; ++fi)
				if (fine(block,fi,fj,0) <= 0) interfaceneg = true; else interfacepos = true;
			if ((interfaceneg && !inside(i,j,k-1)) || (interfacepos && inside(i,j,k-1))) initBlock (i,j,k-1); // create it if so
		}
		if (k+1 < Ncoarse && (nbr = coarse(i,j,k+1)) == 0) {
			bool interfacepos = false, interfaceneg = false;
			for (int fj = 0; fj <= Nsub; ++fj) for (int fi = 0; fi <= Nsub; ++fi)
				if (fine(block,fi,fj,Nsub) <= 0) interfaceneg = true; else interfacepos = true;
			if ((interfaceneg && !inside(i,j,k+1)) || (interfacepos && inside(i,j,k+1))) initBlock (i,j,k+1); // create it if so
		}

	} else if (block[0] > 0) { // boundary moved outside of block, this block is now outside
		if (i > 0) {
			if ((nbr = coarse(i-1,j,k)) != 0) { // if neighbour exists, enforce shared (but limited) values.
				for (int fk = 0; fk <= Nsub; ++fk) for (int fj = 0; fj <= Nsub; ++fj)
					fine(nbr,Nsub,fj,fk) = fine(block,0,fj,fk);
			} else if (inside(i-1,j,k)) initBlock (i-1,j,k); // otherwise, if boundary moved to neighbour then create it
		}
		if (i+1 < Ncoarse) {
			if ((nbr = coarse(i+1,j,k)) != 0) { // if neighbour exists, enforce shared (but limited) values.
				for (int fk = 0; fk <= Nsub; ++fk) for (int fj = 0; fj <= Nsub; ++fj)
					fine(nbr,0,fj,fk) = fine(block,Nsub,fj,fk);
			} else if (inside(i+1,j,k)) initBlock (i-1,j,k); // otherwise, if boundary moved to neighbour then create it
		}
		if (j > 0) {
			if ((nbr = coarse(i,j-1,k)) != 0) { // if neighbour exists, enforce shared (but limited) values.
				for (int fk = 0; fk <= Nsub; ++fk) for (int fi = 0; fi <= Nsub; ++fi)
					fine(nbr,fi,Nsub,fk) = fine(block,fi,0,fk);
			} else if (inside(i,j-1,k)) initBlock (i,j-1,k); // otherwise, if boundary moved to neighbour then create it
		}
		if (j+1 < Ncoarse) {
			if ((nbr = coarse(i,j+1,k)) != 0) { // if neighbour exists, enforce shared (but limited) values.
				for (int fk = 0; fk <= Nsub; ++fk) for (int fi = 0; fi <= Nsub; ++fi)
					fine(nbr,fi,0,fk) = fine(block,fi,Nsub,fk);
			} else if (inside(i,j+1,k)) initBlock (i,j+1,k); // otherwise, if boundary moved to neighbour then create it
		}
		if (k > 0) {
			if ((nbr = coarse(i,j,k-1)) != 0) { // if neighbour exists, enforce shared (but limited) values.
				for (int fj = 0; fj <= Nsub; ++fj) for (int fi = 0; fi <= Nsub; ++fi)
					fine(nbr,fi,fj,Nsub) = fine(block,fi,fj,0);
			} else if (inside(i,j,k-1)) initBlock (i,j,k-1); // otherwise, if boundary moved to neighbour then create it
		}
		if (k+1 < Ncoarse) {
			if ((nbr = coarse(i,j,k+1)) != 0) { // if neighbour exists, enforce shared (but limited) values.
				for (int fj = 0; fj <= Nsub; ++fj) for (int fi = 0; fi <= Nsub; ++fi)
					fine(nbr,fi,fj,0) = fine(block,fi,fj,Nsub);
			} else if (inside(i,j,k+1)) initBlock (i,j,k+1); // otherwise, if boundary moved to neighbour then create it
		}
		delete[] block;
		coarse(i,j,k) = 0;
		inside(i,j,k) = false;

	} else {
		if (i > 0) {
			if ((nbr = coarse(i-1,j,k)) != 0) { // if neighbour exists, enforce shared (but limited) values.
				for (int fk = 0; fk <= Nsub; ++fk) for (int fj = 0; fj <= Nsub; ++fj)
					fine(nbr,Nsub,fj,fk) = fine(block,0,fj,fk);
			} else if (!inside(i-1,j,k)) initBlock (i-1,j,k); // otherwise, if boundary moved to neighbour then create it
		}
		if (i+1 < Ncoarse) {
			if ((nbr = coarse(i+1,j,k)) != 0) { // if neighbour exists, enforce shared (but limited) values.
				for (int fk = 0; fk <= Nsub; ++fk) for (int fj = 0; fj <= Nsub; ++fj)
					fine(nbr,0,fj,fk) = fine(block,Nsub,fj,fk);
			} else if (!inside(i+1,j,k)) initBlock (i-1,j,k); // otherwise, if boundary moved to neighbour then create it
		}
		if (j > 0) {
			if ((nbr = coarse(i,j-1,k)) != 0) { // if neighbour exists, enforce shared (but limited) values.
				for (int fk = 0; fk <= Nsub; ++fk) for (int fi = 0; fi <= Nsub; ++fi)
					fine(nbr,fi,Nsub,fk) = fine(block,fi,0,fk);
			} else if (!inside(i,j-1,k)) initBlock (i,j-1,k); // otherwise, if boundary moved to neighbour then create it
		}
		if (j+1 < Ncoarse) {
			if ((nbr = coarse(i,j+1,k)) != 0) { // if neighbour exists, enforce shared (but limited) values.
				for (int fk = 0; fk <= Nsub; ++fk) for (int fi = 0; fi <= Nsub; ++fi)
					fine(nbr,fi,0,fk) = fine(block,fi,Nsub,fk);
			} else if (!inside(i,j+1,k)) initBlock (i,j+1,k); // otherwise, if boundary moved to neighbour then create it
		}
		if (k > 0) {
			if ((nbr = coarse(i,j,k-1)) != 0) { // if neighbour exists, enforce shared (but limited) values.
				for (int fj = 0; fj <= Nsub; ++fj) for (int fi = 0; fi <= Nsub; ++fi)
					fine(nbr,fi,fj,Nsub) = fine(block,fi,fj,0);
			} else if (!inside(i,j,k-1)) initBlock (i,j,k-1); // otherwise, if boundary moved to neighbour then create it
		}
		if (k+1 < Ncoarse) {
			if ((nbr = coarse(i,j,k+1)) != 0) { // if neighbour exists, enforce shared (but limited) values.
				for (int fj = 0; fj <= Nsub; ++fj) for (int fi = 0; fi <= Nsub; ++fi)
					fine(nbr,fi,fj,0) = fine(block,fi,fj,Nsub);
			} else if (!inside(i,j,k+1)) initBlock (i,j,k+1); // otherwise, if boundary moved to neighbour then create it
		}
		delete[] block;
		coarse(i,j,k) = 0;
		inside(i,j,k) = true;
	}
}


static void resetDistance (float &newp000, float p000, float p001, float p010, float p100,
						   float p011, float p101, float p110, float p111, float dx)
{
	if (p000 == 0) newp000 = 0;
	else if (p000 > 0) {
		if (p001 <= 0) newp000 = min2 (newp000, .9*p000 + .1*(dx*p000/(p000-p001)));
		if (p010 <= 0) newp000 = min2 (newp000, .9*p000 + .1*(dx*p000/(p000-p010)));
		if (p100 <= 0) newp000 = min2 (newp000, .9*p000 + .1*(dx*p000/(p000-p100)));
		if (p011 <= 0) newp000 = min2 (newp000, .9*p000 + .1*(1.4142136*dx*p000/(p000-p011)));
		if (p101 <= 0) newp000 = min2 (newp000, .9*p000 + .1*(1.4142136*dx*p000/(p000-p101)));
		if (p110 <= 0) newp000 = min2 (newp000, .9*p000 + .1*(1.4142136*dx*p000/(p000-p110)));
		if (p111 <= 0) newp000 = min2 (newp000, .9*p000 + .1*(1.7320508*dx*p000/(p000-p111)));
	} else { // p000 < 0
		if (p001 >= 0) newp000 = max2 (newp000, .9*p000 + .1*(dx*p000/(p001-p000)));
		if (p010 >= 0) newp000 = max2 (newp000, .9*p000 + .1*(dx*p000/(p010-p000)));
		if (p100 >= 0) newp000 = max2 (newp000, .9*p000 + .1*(dx*p000/(p100-p000)));
		if (p011 >= 0) newp000 = max2 (newp000, .9*p000 + .1*(1.4142136*dx*p000/(p011-p000)));
		if (p101 >= 0) newp000 = max2 (newp000, .9*p000 + .1*(1.4142136*dx*p000/(p101-p000)));
		if (p110 >= 0) newp000 = max2 (newp000, .9*p000 + .1*(1.4142136*dx*p000/(p110-p000)));
		if (p111 >= 0) newp000 = max2 (newp000, .9*p000 + .1*(1.7320508*dx*p000/(p111-p000)));
	}
}


void FastLevelSet::
reinitBlock (int i, int j, int k, bool &containsBoundary)
{
	float *oldblock = coarse(i,j,k), *newblock = temp_block;

	// initialize with upper bounds on signed distance
	for (int r=cube(Nsub+1)-1; r>=0; --r) newblock[r] = sign(oldblock[r])*coarsedx;

	// find (and fix) boundary
	resetStatus();
	containsBoundary = false;
	for (int fk = 1; fk <= Nsub; ++fk)
		for (int fj = 1; fj <= Nsub; ++fj)
			for (int fi = 1; fi <= Nsub; ++fi) {
				float p111 = fine(oldblock,fi,fj,fk), p011 = fine(oldblock,fi-1,fj,fk),
					  p101 = fine(oldblock,fi,fj-1,fk), p110 = fine(oldblock,fi,fj,fk-1),
					  p100 = fine(oldblock,fi,fj-1,fk-1), p010 = fine(oldblock,fi-1,fj,fk-1),
					  p001 = fine(oldblock,fi-1,fj-1,fk), p000 = fine(oldblock,fi-1,fj-1,fk-1);
				bool pos = p111>0 || p011>0 || p101>0 || p110>0 || p100>0 || p010>0 || p001>0 || p000>0;
				bool neg = p111<=0 || p011<=0 || p101<=0 || p110<=0 || p100<=0 || p010<=0 || p001<=0 || p000<=0;
				if (pos && neg) {
					status(fi,fj,fk) = status_counter;
					status(fi-1,fj,fk) = status_counter;
					status(fi,fj-1,fk) = status_counter;
					status(fi-1,fj-1,fk) = status_counter;
					status(fi,fj,fk-1) = status_counter;
					status(fi-1,fj,fk-1) = status_counter;
					status(fi,fj-1,fk-1) = status_counter;
					status(fi-1,fj-1,fk-1) = status_counter;
					containsBoundary = true;
					resetDistance (fine(newblock,fi-1,fj-1,fk-1),p000,p001,p010,p100,p011,p101,p110,p111, dx);
					resetDistance (fine(newblock,fi-1,fj-1,fk),p001,p000,p011,p101,p010,p100,p111,p110, dx);
					resetDistance (fine(newblock,fi-1,fj,fk-1),p010,p011,p000,p110,p001,p111,p100,p101, dx);
					resetDistance (fine(newblock,fi,fj-1,fk-1),p100,p101,p110,p000,p111,p001,p010,p011, dx);
					resetDistance (fine(newblock,fi-1,fj,fk),p011,p010,p001,p111,p000,p110,p101,p100, dx);
					resetDistance (fine(newblock,fi,fj-1,fk),p101,p100,p111,p001,p110,p000,p011,p010, dx);
					resetDistance (fine(newblock,fi,fj,fk-1),p110,p111,p100,p010,p101,p011,p000,p001, dx);
					resetDistance (fine(newblock,fi,fj,fk),p111,p110,p101,p011,p100,p010,p001,p000, dx);
				}
			}
	if (containsBoundary) {
		localDistanceSweeps (newblock);
		//for (int r = cube(Nsub+1)-1; r >= 0; --r) oldblock[r] = .99*oldblock[r]+.01*newblock[r];
		temp_block = oldblock;
		coarse(i,j,k) = newblock;
		enforceSharedBlock (i, j, k);
	}
}


void FastLevelSet::
initBlock (int i, int j, int k)
{
	float *block = (coarse(i,j,k) = new float[cube(Nsub+1)]), *nbr;
	if (inside(i,j,k)) {for (int r=cube(Nsub+1)-1; r>=0; --r) block[r] = -coarsedx;}
	else {for (int r=cube(Nsub+1)-1; r>=0; --r) block[r] = coarsedx;}
	resetStatus();
	if (i > 0 && (nbr = coarse(i-1,j,k)) != 0) // if neighbour exists, copy shared values
		for (int fk = 0; fk <= Nsub; ++fk) for (int fj = 0; fj <= Nsub; ++fj) {
			fine(block,0,fj,fk) = fine(nbr,Nsub,fj,fk);
			status(0,fj,fk) = status_counter;
		}
	if (i+1 < Ncoarse && (nbr = coarse(i+1,j,k)) != 0) // if neighbour exists, copy shared values
		for (int fk = 0; fk <= Nsub; ++fk) for (int fj = 0; fj <= Nsub; ++fj) {
			fine(block,Nsub,fj,fk) = fine(nbr,0,fj,fk);
			status(Nsub,fj,fk) = status_counter;
		}
	if (j > 0 && (nbr = coarse(i,j-1,k)) != 0) // if neighbour exists, copy shared values
		for (int fk = 0; fk <= Nsub; ++fk) for (int fi = 0; fi <= Nsub; ++fi) {
			fine(block,fi,0,fk) = fine(nbr,fi,Nsub,fk);
			status(fi,0,fk) = status_counter;
		}
	if (j+1 < Ncoarse && (nbr = coarse(i,j+1,k)) != 0) // if neighbour exists, copy shared values
		for (int fk = 0; fk <= Nsub; ++fk) for (int fi = 0; fi <= Nsub; ++fi) {
			fine(block,fi,Nsub,fk) = fine(nbr,fi,0,fk);
			status(fi,Nsub,fk) = status_counter;
		}
	if (k > 0 && (nbr = coarse(i,j,k-1)) != 0) // if neighbour exists, copy shared values
		for (int fj = 0; fj <= Nsub; ++fj) for (int fi = 0; fi <= Nsub; ++fi) {
			fine(block,fi,fj,0) = fine(nbr,fi,fj,Nsub);
			status(fi,fj,0) = status_counter;
		}
	if (k+1 < Ncoarse && (nbr = coarse(i,j,k+1)) != 0) // if neighbour exists, copy shared values
		for (int fj = 0; fj <= Nsub; ++fj) for (int fi = 0; fi <= Nsub; ++fi) {
			fine(block,fi,fj,Nsub) = fine(nbr,fi,fj,0);
			status(fi,fj,Nsub) = status_counter;
		}
	localDistanceSweeps (block); // propagate this interface information into all of interior
	checkShared (i, j, k);
}


static void solveEikonal (float &p, float pilo, float pihi, float pjlo, float pjhi, float pklo, float pkhi, float dx)
{
	float p1, p2, p3, swap;
	if (p <= 0) {
		p1 = max2 (pilo, pihi);
		p2 = max2 (pjlo, pjhi);
		p3 = max2 (pklo, pkhi);
		if (p1 < p2) {swap = p1; p1 = p2; p2 = swap;}
		if (p2 < p3) {swap = p2; p2 = p3; p3 = swap;}
		if (p1 < p2) {swap = p1; p1 = p2; p2 = swap;}
		if (p >= p1) return;
		p = p1-dx;
		if (p >= p2) return;
		p = .5*(p1+p2-sqrt(2*dx*dx-(p2-p1)*(p2-p1)));
		if (p >= p3) return;
		p = (p1+p2+p3)/3.-.5*sqrt(4./9.*(p1+p2+p3)*(p1+p2+p3)-4./3.*(p1*p1+p2*p2+p3*p3-dx*dx));
	} else {
		p1 = min2 (pilo, pihi);
		p2 = min2 (pjlo, pjhi);
		p3 = min2 (pklo, pkhi);
		if (p1 > p2) {swap = p1; p1 = p2; p2 = swap;}
		if (p2 > p3) {swap = p2; p2 = p3; p3 = swap;}
		if (p1 > p2) {swap = p1; p1 = p2; p2 = swap;}
		if (p <= p1) return;
		p = p1+dx;
		if (p <= p2) return;
		p = .5*(p1+p2+sqrt(2*dx*dx-(p2-p1)*(p2-p1)));
		if (p <= p3) return;
		p = (p1+p2+p3)/3.+.5*sqrt(4./9.*(p1+p2+p3)*(p1+p2+p3)-4./3.*(p1*p1+p2*p2+p3*p3-dx*dx));
	}
}


void FastLevelSet::
localDistanceSweeps (float *block)
{
	int i, j, k, ilo, jlo, klo, ihi, jhi, khi;
	for (k=0; k<=Nsub; ++k) {
		if (k > 0) klo = k-1; else klo = k+1;
		if (k < Nsub) khi = k+1; else khi = k-1;
		for (j=0; j<=Nsub; ++j) {
			if (j > 0) jlo = j-1; else jlo = j+1;
			if (j < Nsub) jhi = j+1; else jhi = j-1;
			for (i=0; i<=Nsub; ++i) {
				if (status(i,j,k) == status_counter) continue;
				if (i > 0) ilo = i-1; else ilo = i+1;
				if (i < Nsub) ihi = i+1; else ihi = i-1;
				solveEikonal (fine(block,i,j,k), fine(block,ilo,j,k), fine(block,ihi,j,k),
							  fine(block,i,jlo,k), fine(block,i,jhi,k), fine(block,i,j,klo), fine(block,i,j,khi), dx);
			}
		}
	}
	for (k=Nsub; k>=Nsub; --k) {
		if (k > 0) klo = k-1; else klo = k+1;
		if (k < Nsub) khi = k+1; else khi = k-1;
		for (j=0; j<=Nsub; ++j) {
			if (j > 0) jlo = j-1; else jlo = j+1;
			if (j < Nsub) jhi = j+1; else jhi = j-1;
			for (i=0; i<=Nsub; ++i) {
				if (status(i,j,k) == status_counter) continue;
				if (i > 0) ilo = i-1; else ilo = i+1;
				if (i < Nsub) ihi = i+1; else ihi = i-1;
				solveEikonal (fine(block,i,j,k), fine(block,ilo,j,k), fine(block,ihi,j,k),
							  fine(block,i,jlo,k), fine(block,i,jhi,k), fine(block,i,j,klo), fine(block,i,j,khi), dx);
			}
		}
	}
	for (k=0; k<=Nsub; ++k) {
		if (k > 0) klo = k-1; else klo = k+1;
		if (k < Nsub) khi = k+1; else khi = k-1;
		for (j=Nsub; j>=0; --j) {
			if (j > 0) jlo = j-1; else jlo = j+1;
			if (j < Nsub) jhi = j+1; else jhi = j-1;
			for (i=0; i<=Nsub; ++i) {
				if (status(i,j,k) == status_counter) continue;
				if (i > 0) ilo = i-1; else ilo = i+1;
				if (i < Nsub) ihi = i+1; else ihi = i-1;
				solveEikonal (fine(block,i,j,k), fine(block,ilo,j,k), fine(block,ihi,j,k),
							  fine(block,i,jlo,k), fine(block,i,jhi,k), fine(block,i,j,klo), fine(block,i,j,khi), dx);
			}
		}
	}
	for (k=Nsub; k>=Nsub; --k) {
		if (k > 0) klo = k-1; else klo = k+1;
		if (k < Nsub) khi = k+1; else khi = k-1;
		for (j=Nsub; j>=0; --j) {
			if (j > 0) jlo = j-1; else jlo = j+1;
			if (j < Nsub) jhi = j+1; else jhi = j-1;
			for (i=0; i<=Nsub; ++i) {
				if (status(i,j,k) == status_counter) continue;
				if (i > 0) ilo = i-1; else ilo = i+1;
				if (i < Nsub) ihi = i+1; else ihi = i-1;
				solveEikonal (fine(block,i,j,k), fine(block,ilo,j,k), fine(block,ihi,j,k),
							  fine(block,i,jlo,k), fine(block,i,jhi,k), fine(block,i,j,klo), fine(block,i,j,khi), dx);
			}
		}
	}
	for (k=0; k<=Nsub; ++k) {
		if (k > 0) klo = k-1; else klo = k+1;
		if (k < Nsub) khi = k+1; else khi = k-1;
		for (j=0; j<=Nsub; ++j) {
			if (j > 0) jlo = j-1; else jlo = j+1;
			if (j < Nsub) jhi = j+1; else jhi = j-1;
			for (i=Nsub; i>=0; --i) {
				if (status(i,j,k) == status_counter) continue;
				if (i > 0) ilo = i-1; else ilo = i+1;
				if (i < Nsub) ihi = i+1; else ihi = i-1;
				solveEikonal (fine(block,i,j,k), fine(block,ilo,j,k), fine(block,ihi,j,k),
							  fine(block,i,jlo,k), fine(block,i,jhi,k), fine(block,i,j,klo), fine(block,i,j,khi), dx);
			}
		}
	}
	for (k=Nsub; k>=Nsub; --k) {
		if (k > 0) klo = k-1; else klo = k+1;
		if (k < Nsub) khi = k+1; else khi = k-1;
		for (j=0; j<=Nsub; ++j) {
			if (j > 0) jlo = j-1; else jlo = j+1;
			if (j < Nsub) jhi = j+1; else jhi = j-1;
			for (i=Nsub; i>=0; --i) {
				if (status(i,j,k) == status_counter) continue;
				if (i > 0) ilo = i-1; else ilo = i+1;
				if (i < Nsub) ihi = i+1; else ihi = i-1;
				solveEikonal (fine(block,i,j,k), fine(block,ilo,j,k), fine(block,ihi,j,k),
							  fine(block,i,jlo,k), fine(block,i,jhi,k), fine(block,i,j,klo), fine(block,i,j,khi), dx);
			}
		}
	}
	for (k=0; k<=Nsub; ++k) {
		if (k > 0) klo = k-1; else klo = k+1;
		if (k < Nsub) khi = k+1; else khi = k-1;
		for (j=Nsub; j>=0; --j) {
			if (j > 0) jlo = j-1; else jlo = j+1;
			if (j < Nsub) jhi = j+1; else jhi = j-1;
			for (i=Nsub; i>=0; --i) {
				if (status(i,j,k) == status_counter) continue;
				if (i > 0) ilo = i-1; else ilo = i+1;
				if (i < Nsub) ihi = i+1; else ihi = i-1;
				solveEikonal (fine(block,i,j,k), fine(block,ilo,j,k), fine(block,ihi,j,k),
							  fine(block,i,jlo,k), fine(block,i,jhi,k), fine(block,i,j,klo), fine(block,i,j,khi), dx);
			}
		}
	}
	for (k=Nsub; k>=Nsub; --k) {
		if (k > 0) klo = k-1; else klo = k+1;
		if (k < Nsub) khi = k+1; else khi = k-1;
		for (j=Nsub; j>=0; --j) {
			if (j > 0) jlo = j-1; else jlo = j+1;
			if (j < Nsub) jhi = j+1; else jhi = j-1;
			for (i=Nsub; i>=0; --i) {
				if (status(i,j,k) == status_counter) continue;
				if (i > 0) ilo = i-1; else ilo = i+1;
				if (i < Nsub) ihi = i+1; else ihi = i-1;
				solveEikonal (fine(block,i,j,k), fine(block,ilo,j,k), fine(block,ihi,j,k),
							  fine(block,i,jlo,k), fine(block,i,jhi,k), fine(block,i,j,klo), fine(block,i,j,khi), dx);
			}
		}
	}
}


//=================================== file IO ==================================

void FastLevelSet::
checkShared (const int i, const int j, const int k)
{
	float *block = coarse(i,j,k);
	if (i > 0 && coarse(i-1,j,k)) {
		float *other = coarse(i-1,j,k);
		for (int fk = 0; fk <= Nsub; ++fk)
			for (int fj = 0; fj <= Nsub; ++fj)
				if (flabs(fine(block,0,fj,fk)-fine(other,Nsub,fj,fk)) > .001*dx) {
					fprintf(stderr,"inconsistent boundary between %d %d %d and %d %d %d at (-, %d, %d)\n",i-1,j,k,i,j,k,fj,fk);
					fprintf(stderr,"   values:%g  versus %g\n",fine(block,0,fj,fk),fine(other,Nsub,fj,fk));
				}
	}
	if (i+1 < Ncoarse && coarse(i+1,j,k)) {
		float *other = coarse(i+1,j,k);
		for (int fk = 0; fk <= Nsub; ++fk)
			for (int fj = 0; fj <= Nsub; ++fj)
				if (flabs(fine(block,Nsub,fj,fk)-fine(other,0,fj,fk)) > .001*dx) {
					fprintf(stderr,"inconsistent boundary between %d %d %d and %d %d %d at (-, %d, %d)\n",i,j,k,i+1,j,k,fj,fk);
					fprintf(stderr,"   values:%g  versus %g\n",fine(block,Nsub,fj,fk),fine(other,0,fj,fk));
				}
	}
	if (j > 0 && coarse(i,j-1,k)) {
		float *other = coarse(i,j-1,k);
		for (int fk = 0; fk <= Nsub; ++fk)
			for (int fi = 0; fi <= Nsub; ++fi)
				if (flabs(fine(block,fi,0,fk)-fine(other,fi,Nsub,fk)) > .001*dx) {
					fprintf(stderr,"inconsistent boundary between %d %d %d and %d %d %d at (%d, -, %d)\n",i,j-1,k,i,j,k,fi,fk);
					fprintf(stderr,"   values:%g  versus %g\n",fine(block,fi,0,fk),fine(other,fi,Nsub,fk));
				}
	}
	if (j+1 < Ncoarse && coarse(i,j+1,k)) {
		float *other = coarse(i,j+1,k);
		for (int fk = 0; fk <= Nsub; ++fk)
			for (int fi = 0; fi <= Nsub; ++fi)
				if (flabs(fine(block,fi,Nsub,fk)-fine(other,fi,0,fk)) > .001*dx) {
					fprintf(stderr,"inconsistent boundary between %d %d %d and %d %d %d at (-, %d, %d)\n",i,j,k,i,j+1,k,fi,fk);
					fprintf(stderr,"   values:%g  versus %g\n",fine(block,fi,Nsub,fk),fine(other,fi,0,fk));
				}
	}
}


void FastLevelSet::
write (const char *filename)
{
//	enforceSharedValues();
	FILE *out = fopen(filename,"wb");
	if (out == NULL) {
		fprintf (stderr, "Error opening file '%s' for writing\n", filename);
		return;
	}
	bool error = false;
	error |= fwrite_32 (&N, 1, out);
	error |= fwrite_32 (&Ncoarse, 1, out);
	error |= fwrite_32 (&Nsub, 1, out);
	if (error) {
		fprintf (stderr, "Error writing file '%s' (preamble)\n", filename);
		fclose (out);
		return;
	}
	char zero=0, one=1, negone=-1;
	for (int k = 0; k < Ncoarse; ++k)
		for (int j = 0; j < Ncoarse; ++j)
			for (int i = 0; i < Ncoarse; ++i) {
				if (coarse(i,j,k)) {
					error |= (fwrite (&zero, 1, 1, out) != 1);
				} else if (inside(i,j,k)) {
					error |= (fwrite (&negone, 1, 1, out) != 1);
				} else {
					error |= (fwrite (&one, 1, 1, out) != 1);
				}
			}
	if (error) {
		fprintf (stderr, "Error writing file '%s' (coarse voxels)\n", filename);
		fclose (out);
		return;
	}
	for (int k = 0; k < Ncoarse; ++k)
		for (int j = 0; j < Ncoarse; ++j)
			for (int i = 0; i < Ncoarse; ++i)
				if (coarse(i,j,k)) error |= fwrite_32 (coarse(i,j,k), cube(Nsub+1), out);
	fclose (out);
	if (error) fprintf (stderr, "Error writing file '%s' (fine samples)\n", filename);
}


void FastLevelSet::
read (const char *filename)
{
	FILE *in = fopen(filename,"rb");
	if (in == NULL) {
		fprintf (stderr, "Error opening file '%s' for reading\n", filename);
		return;
	}
	bool error = false;
	int newN, newNcoarse, newNsub;
	error |= fread_32 (&newN, 1, in);
	error |= fread_32 (&newNcoarse, 1, in);
	error |= fread_32 (&newNsub, 1, in);
	printf ("N=%d, Ncoarse=%d, Nsub=%d\n", newN, newNcoarse, newNsub);
	if (error) {
		fprintf (stderr, "Error reading file '%s' (preamble)\n", filename);
		fclose (in);
		return;
	}

	for (int i = cube(Ncoarse)-1; i>= 0; --i) delete[] coarse_raw[i];
	delete[] coarse_raw;
	delete[] inside_raw;
	delete[] status_raw;
	delete[] temp_block;

	N = newN;
	Ncoarse = newNcoarse;
	Nsub = newNsub;
	dx = 1./N; overdx = 1./dx;
	coarsedx = Nsub*dx; overcoarsedx = 1./coarsedx;
	coarse_raw = new (float*)[cube(Ncoarse)];
	inside_raw = new bool[cube(Ncoarse)];
	for (int i = cube(Ncoarse)-1; i >= 0; --i) {coarse_raw[i] = 0; inside_raw[i] = false;}
	status_counter = 0;
	status_raw = new unsigned int[cube(Nsub+1)];
	memset (status_raw, 0, sizeof(unsigned int)*cube(Nsub+1));
	temp_block = new float[cube(Nsub+1)];

	char voxel;
	for (int k = 0; k < Ncoarse; ++k)
		for (int j = 0; j < Ncoarse; ++j)
			for (int i = 0; i < Ncoarse; ++i) {
				error |= (fread (&voxel, 1, 1, in) != 1);
				if (error) {
					fprintf (stderr, "Error reading file '%s' (coarse voxels)\n", filename);
					fclose (in);
					return;
				}
				if (voxel == 0) coarse(i,j,k) = new float[cube(Nsub+1)];
				inside(i,j,k) = (voxel <= 0);
			}
	for (int k = 0; k < Ncoarse; ++k)
		for (int j = 0; j < Ncoarse; ++j)
			for (int i = 0; i < Ncoarse; ++i)
				if (coarse(i,j,k)) error |= fread_32 (coarse(i,j,k), cube(Nsub+1), in);
	fclose (in);
	if (error) fprintf (stderr, "Error reading file '%s' (fine samples)\n", filename);
}


//================================= initialize =====================================

void FastLevelSet::
initialize (AnalyticLevelSet phi)
{
	int nonempties = 0;
	// first clean up
	for (int i = cube(Ncoarse)-1; i>= 0;  --i) { delete[] coarse_raw[i]; coarse_raw[i] = 0; }
	// next construct the new levelset
	float diaglen = .5*sqrt(3.)*coarsedx;
	for (int k = 0; k < Ncoarse; ++k)
		for (int j = 0; j < Ncoarse; ++j)
			for (int i = 0; i < Ncoarse; ++i) {
				float pmiddle = phi(coarsedx*(i+.5),coarsedx*(j+.5),coarsedx*(k+.5));
				if (flabs(pmiddle) < 1.01*diaglen) { // if the surface can intersect this voxel
					coarse(i,j,k) = new float[cube(Nsub+1)];
					float *block = coarse(i,j,k);
					bool pos = false, neg = false;
					for (int fk = 0; fk <= Nsub; ++fk)
						for (int fj = 0; fj <= Nsub; ++fj)
							for (int fi = 0; fi <= Nsub; ++fi) {
								fine(block,fi,fj,fk) = phi(coarsedx*i+dx*fi,coarsedx*j+dx*fj,coarsedx*k+dx*fk);
								if (fine(block,fi,fj,fk) >= 0) pos = true;
								if (fine(block,fi,fj,fk) <= 0) neg = true;
							}
					if (neg) inside(i,j,k) = true; else inside(i,j,k) = false;
					if (!pos || !neg) { delete[] coarse(i,j,k); coarse(i,j,k) = 0; }
					else ++nonempties;
				} else {
					coarse(i,j,k) = 0;
					inside(i,j,k) = (pmiddle <= 0);
				}
			}
}


//=============================== refine =================================

// try to look up phi from whatever block might define it - return true if successful
bool FastLevelSet::
findphi (int i, int j, int k, int fi, int fj, int fk, float &phi)
{
	while (fi < 0) {fi += Nsub; --i;}
	while (fi > Nsub) {fi -= Nsub; ++i;}
	while (fj < 0) {fj += Nsub; --j;}
	while (fj > Nsub) {fj -= Nsub; ++j;}
	while (fk < 0) {fk += Nsub; --k;}
	while (fk > Nsub) {fk -= Nsub; ++k;}
	if (i < 0 || i >= Ncoarse || j < 0 || j >= Ncoarse || k < 0 || k >= Ncoarse)
		return false;
	if (coarse(i,j,k)) {
		phi = fine(coarse(i,j,k),fi,fj,fk);
		return true;
	}

	if (fi == 0 && i > 0) {
		if (coarse(i-1,j,k)) {
			phi = fine(coarse(i-1,j,k),Nsub,fj,fk);
			return true;
		}
		if (fj == 0 && j > 0) {
			if (coarse(i-1,j-1,k)) {
				phi = fine(coarse(i-1,j-1,k),Nsub,Nsub,fk);
				return true;
			}
			if (fk == 0 && k-1 >= 0 && coarse(i-1,j-1,k-1)) {
				phi = fine(coarse(i-1,j-1,k-1),Nsub,Nsub,Nsub);
				return true;
			}
			if (fk == Nsub && k+1 < Ncoarse && coarse(i-1,j-1,k+1)) {
				phi = fine(coarse(i-1,j-1,k+1),Nsub,Nsub,0);
				return true;
			}
		}
		if (fj == Nsub && j+1 < Ncoarse) {
			if (coarse(i-1,j+1,k)) {
				phi = fine(coarse(i-1,j+1,k),Nsub,0,fk);
				return true;
			}
			if (fk == 0 && k-1 >= 0 && coarse(i-1,j+1,k-1)) {
				phi = fine(coarse(i-1,j+1,k-1),Nsub,0,Nsub);
				return true;
			}
			if (fk == Nsub && k+1 < Ncoarse && coarse(i-1,j+1,k+1)) {
				phi = fine(coarse(i-1,j+1,k+1),Nsub,0,0);
				return true;
			}
		}
		if (fk == 0 && k-1 >= 0 && coarse(i-1,j,k-1)) {
			phi = fine(coarse(i-1,j,k-1),Nsub,fj,Nsub);
			return true;
		}
		if (fk == Nsub && k+1 < Ncoarse && coarse(i-1,j,k+1)) {
			phi = fine(coarse(i-1,j,k+1),Nsub,fj,0);
			return true;
		}
	}
	if (fi == Nsub && i+1 < Ncoarse) {
		if (coarse(i+1,j,k)) {
			phi = fine(coarse(i+1,j,k),0,fj,fk);
			return true;
		}
		if (fj == 0 && j > 0) {
			if (coarse(i+1,j-1,k)) {
				phi = fine(coarse(i+1,j-1,k),0,Nsub,fk);
				return true;
			}
			if (fk == 0 && k-1 >= 0 && coarse(i+1,j-1,k-1)) {
				phi = fine(coarse(i+1,j-1,k-1),0,Nsub,Nsub);
				return true;
			}
			if (fk == Nsub && k+1 < Ncoarse && coarse(i+1,j-1,k+1)) {
				phi = fine(coarse(i+1,j-1,k+1),0,Nsub,0);
				return true;
			}
		}
		if (fj == Nsub && j+1 < Ncoarse) {
			if (coarse(i+1,j+1,k)) {
				phi = fine(coarse(i+1,j+1,k),0,0,fk);
				return true;
			}
			if (fk == 0 && k-1 >= 0 && coarse(i+1,j+1,k-1)) {
				phi = fine(coarse(i+1,j+1,k-1),0,0,Nsub);
				return true;
			}
			if (fk == Nsub && k+1 < Ncoarse && coarse(i+1,j+1,k+1)) {
				phi = fine(coarse(i+1,j+1,k+1),0,0,0);
				return true;
			}
		}
		if (fk == 0 && k-1 >= 0 && coarse(i+1,j,k-1)) {
			phi = fine(coarse(i+1,j,k-1),0,fj,Nsub);
			return true;
		}
		if (fk == Nsub && k+1 < Ncoarse && coarse(i+1,j,k+1)) {
			phi = fine(coarse(i+1,j,k+1),0,fj,0);
			return true;
		}
	}

	if (fj == 0 && j > 0) {
		if (coarse(i,j-1,k)) {
			phi = fine(coarse(i,j-1,k),fi,Nsub,fk);
			return true;
		}
		if (fk == 0 && k-1 >= 0 && coarse(i,j-1,k-1)) {
			phi = fine(coarse(i,j-1,k-1),fi,Nsub,Nsub);
			return true;
		}
		if (fk == Nsub && k+1 < Ncoarse && coarse(i,j-1,k+1)) {
			phi = fine(coarse(i,j-1,k+1),fi,Nsub,0);
			return true;
		}
	}
	if (fj == Nsub && j+1 < Ncoarse) {
		if (coarse(i,j+1,k)) {
			phi = fine(coarse(i,j+1,k),fi,0,fk);
			return true;
		}
		if (fk == 0 && k-1 >= 0 && coarse(i,j+1,k-1)) {
			phi = fine(coarse(i,j+1,k-1),fi,0,Nsub);
			return true;
		}
		if (fk == Nsub && k+1 < Ncoarse && coarse(i,j+1,k+1)) {
			phi = fine(coarse(i,j+1,k+1),fi,0,0);
			return true;
		}
	}

	if (fk == 0 && k-1 >= 0) {
		if (coarse(i,j,k-1)) {
			phi = fine(coarse(i,j,k-1),fi,fj,Nsub);
			return true;
		}
	}
	if (fk == Nsub && k+1 < Ncoarse) {
	 	if (coarse(i,j,k+1)) {
			phi = fine(coarse(i,j,k+1),fi,fj,0);
			return true;
		}
	}

	return false;
}


void FastLevelSet::
refineNsub (void)
{
	int newNsub = 2*Nsub;
	// first do trilinear interpolation
	for (int i=cube(Ncoarse)-1; i >= 0; --i)
		if (coarse_raw[i]) {
			float *oldblock = coarse_raw[i];
			float *newblock = new float[cube(newNsub+1)];
			for (int fk=0; fk <= Nsub; ++fk)
				for (int fj=0; fj <= Nsub; ++fj)
					for (int fi=0; fi <= Nsub; ++fi) {
						if (fi > 0 && fj > 0 && fk > 0)
							newblock[2*fi-1+(newNsub+1)*(2*fj-1+(newNsub+1)*(2*fk-1))] =
								.125*(fine(oldblock,fi-1,fj-1,fk-1)+fine(oldblock,fi,fj-1,fk-1)+
								      fine(oldblock,fi-1,fj,fk-1)+fine(oldblock,fi,fj,fk-1)+
									  fine(oldblock,fi-1,fj-1,fk)+fine(oldblock,fi,fj-1,fk)+
									  fine(oldblock,fi-1,fj,fk)+fine(oldblock,fi,fj,fk));
						newblock[2*fi+(newNsub+1)*(2*fj+(newNsub+1)*2*fk)] = fine(oldblock,fi,fj,fk);
						if (fi > 0)
							newblock[2*fi-1+(newNsub+1)*(2*fj+(newNsub+1)*2*fk)] =
								.5*(fine(oldblock,fi-1,fj,fk)+fine(oldblock,fi,fj,fk));
						if (fj > 0)
							newblock[2*fi+(newNsub+1)*(2*fj-1+(newNsub+1)*2*fk)] =
								.5*(fine(oldblock,fi,fj-1,fk)+fine(oldblock,fi,fj,fk));
						if (fi > 0 && fj > 0)
							newblock[2*fi-1+(newNsub+1)*(2*fj-1+(newNsub+1)*2*fk)] =
								.25*(fine(oldblock,fi-1,fj-1,fk)+fine(oldblock,fi,fj-1,fk)+
									 fine(oldblock,fi-1,fj,fk)+fine(oldblock,fi,fj,fk));
						if (fk > 0)
							newblock[2*fi+(newNsub+1)*(2*fj+(newNsub+1)*(2*fk-1))] =
								.5*(fine(oldblock,fi,fj,fk-1)+fine(oldblock,fi,fj,fk));
						if (fi > 0 && fk > 0)
							newblock[2*fi-1+(newNsub+1)*(2*fj+(newNsub+1)*(2*fk-1))] =
								.25*(fine(oldblock,fi-1,fj,fk-1)+fine(oldblock,fi,fj,fk-1)+
									 fine(oldblock,fi-1,fj,fk)+fine(oldblock,fi,fj,fk));
						if (fj > 0 && fk > 0)
							newblock[2*fi+(newNsub+1)*(2*fj-1+(newNsub+1)*(2*fk-1))] =
								.25*(fine(oldblock,fi,fj-1,fk-1)+fine(oldblock,fi,fj,fk-1)+
									 fine(oldblock,fi,fj-1,fk)+fine(oldblock,fi,fj,fk));
					}
			delete[] coarse_raw[i];
			coarse_raw[i] = newblock;
		}

	delete[] status_raw;
	status_raw = new unsigned int[cube(newNsub+1)];
	memset (status_raw, 0, sizeof(unsigned int)*cube(newNsub+1));
	status_counter = 0;
	delete[] temp_block;
	temp_block = new float[cube(newNsub+1)];
	Nsub = newNsub;
	N = Ncoarse*Nsub;
	dx *= .5; overdx *= 2;
}


void FastLevelSet::
smoothRefined (void)
{
	// figure out the smoothed vertex positions
	for (int k = 0; k < Ncoarse; ++k) for (int j = 0; j < Ncoarse; ++j) for (int i = 0; i < Ncoarse; ++i) {
		float *block = coarse(i,j,k);
		if (!block) continue;
		for (int fk = 0; fk <= Nsub; ++fk) for (int fj = 0; fj <= Nsub; ++fj) for (int fi = 0; fi <= Nsub; ++fi)
			if (fi%2 + fj%2 + fk%2 == 0) {
				float center0, center1, center2, center3, center4, center5, center6, center7;
				float face0, face1, face2, face3, face4, face5, face6, face7, face8, face9, face10, face11;
				float edge0, edge1, edge2, edge3, edge4, edge5;
				center0 = center1 = center2 = center3 = center4 = center5 = center6 = center7 = face0 = face1 = face2 = face3
					= face4 = face5 = face6 = face7 = face8 = face9 = face10 = face11
					= edge0 = edge1 = edge2 = edge3 = edge4 = edge5 = fine(block,fi,fj,fk);
				findphi(i,j,k,fi-1,fj-1,fk-1,center0);
				findphi(i,j,k,fi+1,fj-1,fk-1,center1);
				findphi(i,j,k,fi-1,fj+1,fk-1,center2);
				findphi(i,j,k,fi+1,fj+1,fk-1,center3);
				findphi(i,j,k,fi-1,fj-1,fk+1,center4);
				findphi(i,j,k,fi+1,fj-1,fk+1,center5);
				findphi(i,j,k,fi-1,fj+1,fk+1,center6);
				findphi(i,j,k,fi+1,fj+1,fk+1,center7);
				findphi(i,j,k,fi-1,fj-1,fk,face0);
				findphi(i,j,k,fi+1,fj-1,fk,face1);
				findphi(i,j,k,fi-1,fj+1,fk,face2);
				findphi(i,j,k,fi+1,fj+1,fk,face3);
				findphi(i,j,k,fi-1,fj,fk-1,face4);
				findphi(i,j,k,fi+1,fj,fk-1,face5);
				findphi(i,j,k,fi-1,fj,fk+1,face6);
				findphi(i,j,k,fi+1,fj,fk+1,face7);
				findphi(i,j,k,fi,fj-1,fk-1,face8);
				findphi(i,j,k,fi,fj+1,fk-1,face9);
				findphi(i,j,k,fi,fj-1,fk+1,face10);
				findphi(i,j,k,fi,fj+1,fk+1,face11);
				findphi(i,j,k,fi-1,fj,fk,edge0);
				findphi(i,j,k,fi+1,fj,fk,edge1);
				findphi(i,j,k,fi,fj-1,fk,edge2);
				findphi(i,j,k,fi,fj+1,fk,edge3);
				findphi(i,j,k,fi,fj,fk-1,edge4);
				findphi(i,j,k,fi,fj,fk+1,edge5);
				fine(block,fi,fj,fk) = (fine(block,fi,fj,fk)+3*(edge0+edge1+edge2+edge3+edge4+edge5)/6
										+3*(face0+face1+face2+face3+face4+face5+face6+face7+face8+face9+face10+face11)/12
										+(center0+center1+center2+center3+center4+center5+center6+center7)/8)/8;
			}
		}
	// figure out the smoothed edge positions
	for (int k = 0; k < Ncoarse; ++k) for (int j = 0; j < Ncoarse; ++j) for (int i = 0; i < Ncoarse; ++i) {
		float *block = coarse(i,j,k);
		if (!block) continue;
		for (int fk = 0; fk <= Nsub; ++fk) for (int fj = 0; fj <= Nsub; ++fj) for (int fi = 0; fi <= Nsub; ++fi)
			if (fi%2 + fj%2 + fk%2 == 1) {
				float center0, center1, center2, center3;
				float face0, face1, face2, face3;
				center0 = center1 = center2 = center3 = face0 = face1 = face2 = face3 = fine(block,fi,fj,fk);
				if (fi%2) {
					findphi(i,j,k,fi,fj-1,fk-1,center0);
					findphi(i,j,k,fi,fj+1,fk-1,center1);
					findphi(i,j,k,fi,fj-1,fk+1,center2);
					findphi(i,j,k,fi,fj+1,fk+1,center3);
					findphi(i,j,k,fi,fj-1,fk,face0);
					findphi(i,j,k,fi,fj+1,fk,face1);
					findphi(i,j,k,fi,fj,fk-1,face2);
					findphi(i,j,k,fi,fj,fk+1,face3);
				} else if (fj%2) {
					findphi(i,j,k,fi-1,fj,fk-1,center0);
					findphi(i,j,k,fi+1,fj,fk-1,center1);
					findphi(i,j,k,fi-1,fj,fk+1,center2);
					findphi(i,j,k,fi+1,fj,fk+1,center3);
					findphi(i,j,k,fi-1,fj,fk,face0);
					findphi(i,j,k,fi+1,fj,fk,face1);
					findphi(i,j,k,fi,fj,fk-1,face2);
					findphi(i,j,k,fi,fj,fk+1,face3);
				} else {
					findphi(i,j,k,fi-1,fj-1,fk,center0);
					findphi(i,j,k,fi+1,fj-1,fk,center1);
					findphi(i,j,k,fi-1,fj+1,fk,center2);
					findphi(i,j,k,fi+1,fj+1,fk,center3);
					findphi(i,j,k,fi-1,fj,fk,face0);
					findphi(i,j,k,fi+1,fj,fk,face1);
					findphi(i,j,k,fi,fj-1,fk,face2);
					findphi(i,j,k,fi,fj+1,fk,face3);
				}
				fine(block,fi,fj,fk) = (4*fine(block,fi,fj,fk)+(center0+center1+center2+center3)+2*(face0+face1+face2+face3))/16;
			}
		}
	// figure out the smoothed face positions
	for (int k = 0; k < Ncoarse; ++k) for (int j = 0; j < Ncoarse; ++j) for (int i = 0; i < Ncoarse; ++i) {
		float *block = coarse(i,j,k);
		if (!block) continue;
		for (int fk = 0; fk <= Nsub; ++fk) for (int fj = 0; fj <= Nsub; ++fj) for (int fi = 0; fi <= Nsub; ++fi)
			if (fi%2 + fj%2 + fk%2 == 2) {
				float center0, center1;
				center0 = center1 = fine(block,fi,fj,fk);
				if (fi%2 && fj%2) {
					findphi(i,j,k,fi,fj,fk-1,center0);
					findphi(i,j,k,fi,fj,fk+1,center1);
				} else if (fi%2 && fk%2) {
					findphi(i,j,k,fi,fj-1,fk,center0);
					findphi(i,j,k,fi,fj+1,fk,center1);
				} else {
					findphi(i,j,k,fi-1,fj,fk,center0);
					findphi(i,j,k,fi+1,fj,fk,center1);
				}
				fine(block,fi,fj,fk) = (2*fine(block,fi,fj,fk)+center0+center1)/4;
			}
		}
	// now check for holes that moving the boundary may have caused
	for (int d = 0; d < 3; ++d) { // we loop three times to cover common faces, edges, and corners
		for (int k = 0; k < Ncoarse; ++k) for (int j = 0; j < Ncoarse; ++j) for (int i = 0; i < Ncoarse; ++i) {
			float *block = coarse(i,j,k);
			if (!block) continue;
			if (i > 0 && !coarse(i-1,j,k)) {
				bool pos = false, neg = false;
				for (int fk = 0; fk <= Nsub; ++fk) for (int fj = 0; fj <= Nsub; ++fj) {
					if (fine(block,0,fj,fk) >= 0) pos = true;
					if (fine(block,0,fj,fk) <= 0) neg = true;
				}
				if (pos && neg) initBlock (i-1,j,k);
			}
			if (i+1 < Ncoarse && !coarse(i+1,j,k)) {
				bool pos = false, neg = false;
				for (int fk = 0; fk <= Nsub; ++fk) for (int fj = 0; fj <= Nsub; ++fj) {
					if (fine(block,Nsub,fj,fk) >= 0) pos = true;
					if (fine(block,Nsub,fj,fk) <= 0) neg = true;
				}
				if (pos && neg) initBlock (i+1,j,k);
			}
			if (j > 0 && !coarse(i,j-1,k)) {
				bool pos = false, neg = false;
				for (int fk = 0; fk <= Nsub; ++fk) for (int fi = 0; fi <= Nsub; ++fi) {
					if (fine(block,fi,0,fk) >= 0) pos = true;
					if (fine(block,fi,0,fk) <= 0) neg = true;
				}
				if (pos && neg) initBlock (i,j-1,k);
			}
			if (j+1 < Ncoarse && !coarse(i,j+1,k)) {
				bool pos = false, neg = false;
				for (int fk = 0; fk <= Nsub; ++fk) for (int fi = 0; fi <= Nsub; ++fi) {
					if (fine(block,fi,Nsub,fk) >= 0) pos = true;
					if (fine(block,fi,Nsub,fk) <= 0) neg = true;
				}
				if (pos && neg) initBlock (i,j+1,k);
			}
			if (k > 0 && !coarse(i,j,k-1)) {
				bool pos = false, neg = false;
				for (int fj = 0; fj <= Nsub; ++fj) for (int fi = 0; fi <= Nsub; ++fi) {
					if (fine(block,fi,fj,0) >= 0) pos = true;
					if (fine(block,fi,fj,0) <= 0) neg = true;
				}
				if (pos && neg) initBlock (i,j,k-1);
			}
			if (k+1 < Ncoarse && !coarse(i,j,k+1)) {
				bool pos = false, neg = false;
				for (int fj = 0; fj <= Nsub; ++fj) for (int fi = 0; fi <= Nsub; ++fi) {
					if (fine(block,fi,fj,Nsub) >= 0) pos = true;
					if (fine(block,fi,fj,Nsub) <= 0) neg = true;
				}
				if (pos && neg) initBlock (i,j,k+1);
			}
		}
	}
	// and also check for unnecessary blocks
	for (int k = 0; k < Ncoarse; ++k) for (int j = 0; j < Ncoarse; ++j) for (int i = 0; i < Ncoarse; ++i) {
		float *block = coarse(i,j,k);
		if (!block) continue;
		bool pos = false, neg = false;
		for (int r = cube(Nsub+1)-1; r >= 0; --r) {
			if (block[r] >= 0) pos = true;
			if (block[r] <= 0) neg = true;
		}
		if (!neg) { inside(i,j,k) = false; delete[] block; coarse(i,j,k) = 0; }
		else if (!pos) { inside(i,j,k) = true; delete[] block; coarse(i,j,k) = 0; }
	}
}


static bool copycube (int newNsub, pointer_to_float &newblock, int oldNsub, const float *oldblock)
{
	newblock = new float[cube(newNsub+1)];
	bool pos = false, neg = false;
	for (int k = 0; k <= newNsub; ++k) {
		for (int j = 0; j <= newNsub; ++j) {
			for (int i = 0; i <= newNsub; ++i) {
				float value = oldblock[i+(oldNsub+1)*(j+(oldNsub+1)*k)];
				if (value > 0) pos = true; else neg = true;
				newblock[i+(newNsub+1)*(j+(newNsub+1)*k)] = value;
			}
		}
	}
	if (pos && neg) return true;
	else {
		delete[] newblock;
		newblock = 0;
		return neg;
	}
}


void FastLevelSet::
refine (void)
{
	enforceSharedValues();
	refineNsub();
	smoothRefined();
	if (2*Ncoarse < 0.8*pow(N,.75)) {
		int nc = 2*Ncoarse, ns = Nsub/2;
		float **newcoarse_raw = new (float*)[cube(nc)];
		bool *newinside_raw = new bool[cube(nc)];
		for (int k = 0; k < Ncoarse; ++k)
			for (int j = 0; j < Ncoarse; ++j)
				for (int i = 0; i < Ncoarse; ++i) {
#define newcoarse(a,b,c) (newcoarse_raw[(a)+nc*((b)+nc*(c))])
#define newinside(a,b,c) (newinside_raw[(a)+nc*((b)+nc*(c))])
					float *block = coarse(i,j,k);
					if (block) {
						newinside(2*i,2*j,2*k) = copycube (ns, newcoarse(2*i,2*j,2*k), Nsub, block);
						newinside(2*i+1,2*j,2*k) = copycube (ns, newcoarse(2*i+1,2*j,2*k), Nsub, block+ns);
						newinside(2*i,2*j+1,2*k) = copycube (ns, newcoarse(2*i,2*j+1,2*k), Nsub, block+(Nsub+1)*ns);
						newinside(2*i+1,2*j+1,2*k) = copycube (ns, newcoarse(2*i+1,2*j+1,2*k), Nsub, block+ns+(Nsub+1)*ns);
						newinside(2*i,2*j,2*k+1) = copycube (ns, newcoarse(2*i,2*j,2*k+1), Nsub, block+(Nsub+1)*(Nsub+1)*ns);
						newinside(2*i+1,2*j,2*k+1) = copycube (ns, newcoarse(2*i+1,2*j,2*k+1), Nsub, block+ns+(Nsub+1)*(Nsub+1)*ns);
						newinside(2*i,2*j+1,2*k+1) = copycube (ns, newcoarse(2*i,2*j+1,2*k+1), Nsub, block+(Nsub+1)*(ns+(Nsub+1)*ns));
						newinside(2*i+1,2*j+1,2*k+1) = copycube (ns, newcoarse(2*i+1,2*j+1,2*k+1), Nsub, block+ns+(Nsub+1)*(ns+(Nsub+1)*ns));
						delete[] coarse(i,j,k); coarse(i,j,k) = 0;
					} else {
						newcoarse(2*i,2*j,2*k) = 0;
						newcoarse(2*i+1,2*j,2*k) = 0;
						newcoarse(2*i,2*j+1,2*k) = 0;
						newcoarse(2*i+1,2*j+1,2*k) = 0;
						newcoarse(2*i,2*j,2*k+1) = 0;
						newcoarse(2*i+1,2*j,2*k+1) = 0;
						newcoarse(2*i,2*j+1,2*k+1) = 0;
						newcoarse(2*i+1,2*j+1,2*k+1) = 0;
						newinside(2*i,2*j,2*k) = inside(i,j,k);
						newinside(2*i+1,2*j,2*k) = inside(i,j,k);
						newinside(2*i,2*j+1,2*k) = inside(i,j,k);
						newinside(2*i+1,2*j+1,2*k) = inside(i,j,k);
						newinside(2*i,2*j,2*k+1) = inside(i,j,k);
						newinside(2*i+1,2*j,2*k+1) = inside(i,j,k);
						newinside(2*i,2*j+1,2*k+1) = inside(i,j,k);
						newinside(2*i+1,2*j+1,2*k+1) = inside(i,j,k);
					}
#undef newcoarse
#undef newinside
				}
		delete[] coarse_raw;
		delete[] inside_raw;
		coarse_raw = newcoarse_raw;
		inside_raw = newinside_raw;

		delete[] status_raw;
		status_raw = new unsigned int[cube(ns+1)];
		memset (status_raw, 0, sizeof(unsigned int)*cube(ns+1));
		status_counter = 0;
		delete[] temp_block;
		temp_block = new float[cube(ns+1)];

		Ncoarse = nc;
		Nsub = ns;
		coarsedx *= .5; overcoarsedx *= 2;
	}
	printf("N = %d, Ncoarse = %d, Nsub = %d\n", N, Ncoarse, Nsub);
}

