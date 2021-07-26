//
// Source code for the paper
// 
// D. Heck and T. Schloemer and O. Deussen, "Blue Noise Sampling with
// Controlled Aliasing", ACM Trans. Graph., 2013, in press
//
// Copyright (C) 2012,2013 Daniel Heck and Thomas Schloemer
//
#ifndef COMMON_HH_INCLUDED
#define COMMON_HH_INCLUDED

#include <string>
#include <vector>
#include <cmath>
#include <stdint.h>

// -------------------- Math --------------------

template <class T, class U>
inline T Clamp(T x, U min, U max) {
    if (x < min)
	return min;
    else if (x > max)
	return max;
    else
	return x;
}


// -------------------- Utility --------------------

bool HasSuffix(const std::string &s, const std::string &suffix);

enum Endianness {
    BigEndian, LittleEndian
};

Endianness SystemEndianness();
void SwapEndian4(uint8_t *buf, int n);
bool WriteFloats(FILE *fp, const float *p, int n, Endianness endian);
bool ReadFloats(FILE *fp, float *p, int n, Endianness endian);

void ReadPFM(int &width, int &height, float *&data, const char *fname);
void WritePFM(int width, int height, const float *data, const char *fname);

// Output floating point matrix as PGM image
void SavePGM(float *p, int rows, int cols, const char *fname);


// -------------------- Point --------------------

struct Point {
    float x, y;
    Point(float x = 0, float y = 0) : x(x), y(y) {}
};

void LoadPoints(const std::string &fname, std::vector<Point> &points);
void WritePoints(const std::string &fname, std::vector<Point> &points);

// -------------------- Curve --------------------

class Curve {
    std::vector<float> y;
public:
    float x0, x1, dx;
    Curve() { x0 = x1 = dx = 0.0f; }
    Curve(const Curve &c);
    Curve(int size, float x0, float x1);

    float &operator[](size_t i) { return y[i]; }
    const float &operator[](size_t i) const { return y[i]; }

    int size() const { return (int)y.size(); }

    int ToIndex(float x) const {
        return static_cast<int>((x-x0)/dx);
    }

    float ToXCenter(int index) const {
        return x0 + (index+0.5f)*dx;
    }
    float ToX(int index) const {
        return x0 + index*dx;
    }

    // Evaluate curve at point 'x'. Use linear interpolation if 'x' falls
    // between two bins.
    float At(float x) const;

    void Write(const std::string &fname);
    static Curve Read(const std::string &fname);

    void FilterGauss(const Curve &source, float sigma);
};


inline Curve FilterGauss(int nbins, const Curve &source, float sigma) {
    Curve c(nbins, source.x0, source.x1);
    c.FilterGauss(source, sigma);
    return c;
}

float Integrate(const Curve &c, float x0, float x1, float *E = NULL);

inline float Integrate(const Curve &c) {
    return Integrate(c, c.x0, c.x1);
}


// -------------------- Toroidal distance --------------------

inline static float Sqr(float x) {
    return x*x;
}
  
// Compute toroidal distance between points a and b. Both a and b must
// lie inside the unit square.
inline static float DistSqrToroidal(float const *a, float const *b) {
    float dx = a[0]-b[0], dy = a[1]-b[1];
    return Sqr(dx + (dx < -.5f) - (dx > .5f)) + 
        Sqr(dy + (dy < -.5f) - (dy > .5f));
}

inline static float DistToroidal(float const *a, float const *b) {
    return sqrtf(DistSqrToroidal(a, b));
}

// -------------------- RDF calculations --------------------

void CalcRDF(Curve &c, size_t npts, const float *pts, 
        float (*distfunc)(const float *, const float *) = DistToroidal);
Curve CalcRDF(int numbins, size_t npts, const float *pts, 
        float sigma = 0, float (*distfunc)(const float *, const float *) = DistToroidal);

Curve RDF2Power(int npts, const Curve &rdf, int nbins, float x0, float x1);
Curve Power2RDF(int npts, const Curve &power, int nbins, float x0, float x1, float smoothing);

#endif
