#ifndef DEFS_H
#define DEFS_H


#include <iomanip>
#include <cmath>

using namespace std;

#include "stdafx.h"
#include "SDKmisc.h"
#include "DXUTCamera.h"
#include "dxut.h"
#include "dxutmisc.h"

#include "defs_polymesh.h"

////////////////////////////////////////////////////////////
// For debugging and status output
////////////////////////////////////////////////////////////

#define DEBUG_LEVEL 2

const int STATUS_INDENT = 4;
inline int make_nonneg(int v) { return (v < 0 ? 0 : v); }

// Output functions with indenting capability
int g_statusIndentLevel = 0;
#define     INDENT ( STATUS_INDENT * g_statusIndentLevel   )
#define INC_INDENT ( STATUS_INDENT * g_statusIndentLevel++ )
#define DEC_INDENT ( STATUS_INDENT * (g_statusIndentLevel > 0 ? --g_statusIndentLevel : 0) )
#define OUT0(text)        cout << setw(make_nonneg(    INDENT)) << "" << text << endl, cout.flush()
#define OUT0_BEGIN(text)  cout << setw(make_nonneg(INC_INDENT)) << "" << text << endl, cout.flush()
#define OUT0_DONE()       cout << setw(make_nonneg(DEC_INDENT)) << "" << "DONE\n", cout.flush()
//#define RETURN_IF_FAILED(hr) if( FAILED(hr) ) return (DEC_INDENT, hr)

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define RETURN_IF_FAILED(hr) if( FAILED(hr) ) return ( OUT2("ERROR in FILE " __FILE__ " at line " TOSTRING(__LINE__) "."), hr )
#define WARN(text)        OUT0("WARNING: " << text)
#define ERR(text)         OUT0("ERROR: " << text)
#define FATAL_ERROR(text) exit( (ERR(text), 1) )

// Debug level 1
#if (DEBUG_LEVEL >= 1)
#define OUT1(text)        OUT0(text)
#define OUT1_BEGIN(text)  OUT0_BEGIN(text)
#define OUT1_DONE()       OUT0_DONE()
#else
#define OUT1(text)
#define OUT1_BEGIN(text)
#define OUT1_DONE()
#endif

// Debug level 2
#if (DEBUG_LEVEL >= 2)
#define OUT2(text)        OUT0(text)
#define OUT2_BEGIN(text)  OUT0_BEGIN(text)
#define OUT2_DONE()       OUT0_DONE()
#else
#define OUT2(text)
#define OUT2_BEGIN(text)
#define OUT2_DONE()
#endif

////////////////////////////////////////////////////////////

// Shader flags bit (MATCH WITH SHADER)
const int BIT_GROUP_COLOR = 0;
const int BIT_DARK_BG = 1;

// Make sure this matches the shader
const int MAX_FACE_SIZE = 5;
const int POLAR = 1;
const int REG = 2;
const int NGON = 3;

const int MAX_TYPE = MAX_FACE_SIZE;
const int NUM_SIDES[MAX_TYPE+1] = { 0, 3, 4, 3, 4, 5 };
const int MAX_CONTROL_POINTS[MAX_TYPE+1] = { 0, 16, 16, 19, 25, 31 };

const REAL PI = 3.141592653589793f;
const int Vs = 4;
const int MAX_VALENCE = 8;
const int MAX_ALL = 8192;
const unsigned int NUM_OF_INSTANCES = 1;
const unsigned int NUM_PASSES = 2;

// For encoding the offset
const int NBITS_OFFSET = 4; // match the value in the shader

// Number of frames in the frog animation.
const int FROG_FRAMES = 100;
const int MORPH_FRAMES = 100;

using namespace std;

////////////////////////////////////////////////////////////

// Create a string with an integer index in it.
string str_index(const char *pre, int i, const char *post = "", int wid = 1)
{
    assert(0 < wid && wid <= 100); // arbitrary upper bound

    const size_t len = size_t( strlen(pre) + strlen(post) );
    assert(len < 2048); // an artificial bound for security
    char s[4097];
    sprintf_s(s, 4096, "%s%0*d%s", pre, wid, i, post);
    return string(s);
}

////////////////////////////////////////////////////////////

// Bit manipulation
inline void set_bit(UINT32 &flag, int bit) {
    assert(0 <= bit && bit < sizeof(UINT32)*8);
    flag |= (1 << bit);
}
// Clear bit
inline void clear_bit(UINT32 &flag, int bit) {
    assert(0 <= bit && bit < sizeof(UINT32)*8);
    flag &= ~(1 << bit);
}
// Toggle
inline void toggle_bit(UINT32 &flag, int bit) {
    assert(0 <= bit && bit < sizeof(UINT32)*8);
    flag ^= (1 << bit);
}
// Get bit
inline UINT get_bit(UINT32 flag, int bit) {
    assert(0 <= bit && bit < sizeof(UINT32)*8);
    return (flag >> bit) & 1;
}

#endif // DEFS_H

