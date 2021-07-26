
#define PHONG 1
#define BLEND 0

#define MAX_VALENCE 8 
#define PI 3.141592653589793f

// Make sure this matches the C++ code
#define MAX_FACE_SIZE 5
#define DEFAULT 0
#define POLAR   1
#define REG     2
#define MAX_TYPE MAX_FACE_SIZE
const int MAX_CONTROL_POINTS[MAX_TYPE+1] = { 0, 16, 16, 19, 25, 31 };

#define NBITS_OFFSET 4          // match the value in defs.h
#define NBITS_OFFSET_MASK 0xF   // corresponding mask

typedef float4 Point;

//data textures
Texture2D<uint> gOffsetData;
Buffer<float4> gControlPoints1;
Buffer<float4> gControlPoints2;
Buffer<float4> gControlPoints3;
Buffer<float4> gControlPoints4;
Buffer<float4> gControlPoints5;
#define gControlPoints_polar gControlPoints1
#define gControlPoints_reg   gControlPoints2
#define gControlPoints_irreg gControlPoints4
Texture2D gVertexLocation;
Texture2D gRingIndex;

// The formula determining the translation based on the instance ID
inline float4 instance_translation(int iID)
{
//  return float4 (4.0*(iID%3)-4.0, 4.0*floor(iID/3)-4.0, 0, 0);
    return float4 (0, 0, 0, 0);
}

const float4 gDiffuseMtrl[MAX_TYPE+1] = {
//  float4(0.55f, 0.47f, 0.14f, 1.0f),  // 0
//  float4(0.55f, 0.47f, 0.14f, 1.0f),  // 1
//  float4(0.55f, 0.47f, 0.14f, 1.0f),  // 2
//  float4(    0, 0.47f, 0.30f, 1.0f),  // 3
//  float4(0.14f, 0.31f, 0.55f, 1.0f),  // 4
//  float4(0.45f, 0.45f, 0.45f, 1.0f)   // 5
//
    // BezierView colors
    float4(0.75164 , 0.60648 , 0.22648 , 1.0), // gold
    float4(0.714   , 0.4284  , 0.18144 , 1.0), // bronze
    float4(0.75164 , 0.60648 , 0.22648 , 1.0), // gold
//  float4(0.75164 , 0.60648 , 0.22648 , 1.0), // gold
    float4(0.07568 , 0.61424 , 0.07568 , 1.0), // emerald
    float4(0.61424 , 0.04136 , 0.04136 , 1.0), // ruby
    float4(0.50754 , 0.50754 , 0.50754 , 1.0), // silver
//  float4(0.714   , 0.4284  , 0.18144 , 1.0), // bronze
//  float4(0.5     , 0.5     , 0.0     , 1.0), // yellow plastic
//  float4(0.4     , 0.4     , 0.4     , 1.0), // chrome
//  float4(0.427451, 0.470588, 0.541176, 1.0), // Pewter
//  float4(0.0     , 0.51    , 0.51    , 1.0), // cyan plastic
//  float4(0.54    , 0.89    , 0.630   , 1.0), // Jade
//  float4(0.0     , 0.0     , 0.0     , 0.0), // random
//  float4(0.0     , 0.0     , 0.0     , 0.0), // transparent
};

// Shader flags bit (MATCH WITH SHADER)
const int BIT_GROUP_COLOR = 0;
const int BIT_DARK_BG = 1;

DepthStencilState RenderWithStencilState
{
    DepthEnable = true;
    DepthFunc = less;
};
//-----------------------------------------------------------------------------------------
// Constant Buffers (where we store variables by frequency of update)
//-----------------------------------------------------------------------------------------
cbuffer cbEveryFrame
{
    matrix matWorldView;
    matrix matWorldViewProj; 
    float3 gLightVecW1;
    float3 gEyePosW;
};


cbuffer cbOccasional
{
    uint gStart[MAX_TYPE+1];
};

cbuffer cbPerFrame
{
    uint gFlags = 0;
    uint gRowTex = 0;
};

const float WEIGHT[] = {0, 0, 0, 2, 1, -3};

//-----------------------------------------------------------------------------------------
// constant value for Cosine [valence-3][i]
//-----------------------------------------------------------------------------------------
const float cCosf[11][12]=
{
    1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0,
    1.0, -0.5, -0.5, 1.0, -0.5, -0.5, 1.0, -0.5, -0.5, 1.0, -0.5, -0.5,
    1.0, 0.0, -1.0, 0.0,  1.0, 0.0, -1.0, 0.0,  1.0, 0.0, -1.0, 0.0,
    1.0, 0.309016994374947, -0.809016994374947, -0.809016994374948, 0.309016994374947,  1.0, 0.309016994374947, -0.809016994374947, -0.809016994374948, 0.309016994374947,  1.0, 0.309016994374947,
    1.0, 0.5, -0.5, -1.0, -0.5, 0.5,  1.0, 0.5, -0.5, -1.0, -0.5, 0.5,
    1.0, 0.623489801858734, -0.222520933956314, -0.900968867902419, -0.900968867902419, -0.222520933956315, 0.623489801858733,  1.0, 0.623489801858734, -0.222520933956314, -0.900968867902419, -0.900968867902419, 
    1.0, 0.707106781186548, 0.0, -0.707106781186547, -1.0, -0.707106781186548, 0.0, 0.707106781186547, 1.0, 0.707106781186548, 0.0, -0.707106781186547, 
    1.0, 0.766044443118978, 0.17364817766693, -0.5, -0.939692620785908, -0.939692620785908, -0.5, 0.17364817766693, 0.766044443118978,     1.0, 0.766044443118978, 0.17364817766693,
    1.0, 0.809016994374947, 0.309016994374947, -0.309016994374947, -0.809016994374947, -1.0, -0.809016994374948, -0.309016994374948, 0.309016994374947, 0.809016994374947, 1.0, 0.809016994374947,
    1.0, 0.841253532831181, 0.415415013001886, -0.142314838273285, -0.654860733945285, -0.959492973614497, -0.959492973614497, -0.654860733945285, -0.142314838273285, 0.415415013001886, 0.841253532831181, 1.0,
    1.0, 0.866025403784439, 0.5, 0.0, -0.5, -0.866025403784439, -1.0, -0.866025403784439, -0.5, 0.0, 0.5, 0.866025403784438,
};

//-----------------------------------------------------------------------------------------
// constant value for Sine [valence-3][i]
//-----------------------------------------------------------------------------------------
const float cSinf[11][12]=
{
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.866025403784439, -0.866025403784438, 0.0, 0.866025403784439, -0.866025403784438, 0.0, 0.866025403784439, -0.866025403784438, 0.0, 0.866025403784439, -0.866025403784438,
    0.0, 1.0, 0.0, -1.0,  0.0, 1.0, 0.0, -1.0,  0.0, 1.0, 0.0, -1.0,
    0.0, 0.951056516295154, 0.587785252292473, -0.587785252292473, -0.951056516295154,  0.0, 0.951056516295154, 0.587785252292473, -0.587785252292473, -0.951056516295154,  0.0, 0.951056516295154, 
    0.0, 0.866025403784439, 0.866025403784439, 0.0, -0.866025403784438, -0.866025403784439,     0.0, 0.866025403784439, 0.866025403784439, 0.0, -0.866025403784438, -0.866025403784439,
    0.0, 0.78183148246803, 0.974927912181824, 0.433883739117558, -0.433883739117558, -0.974927912181824, -0.78183148246803,   0.0, 0.78183148246803, 0.974927912181824, 0.433883739117558, -0.433883739117558, 
    0.0, 0.707106781186547, 1.0, 0.707106781186548, 0.0, -0.707106781186547, -1.0, -0.707106781186548,  0.0, 0.707106781186547, 1.0, 0.707106781186548, 
    0.0, 0.642787609686539, 0.984807753012208, 0.866025403784439, 0.342020143325669, -0.342020143325669, -0.866025403784438, -0.984807753012208, -0.64278760968654,  0.0, 0.642787609686539, 0.984807753012208,
    0.0, 0.587785252292473, 0.951056516295154, 0.951056516295154, 0.587785252292473, 0.0, -0.587785252292473, -0.951056516295154, -0.951056516295154, -0.587785252292473, 0.0, 0.587785252292473, 
    0.0, 0.540640817455598, 0.909631995354518, 0.989821441880933, 0.755749574354258, 0.28173255684143, -0.281732556841429, -0.755749574354258, -0.989821441880933, -0.909631995354519, -0.540640817455597, 0.0,
    0.0, 0.5, 0.866025403784439, 1.0, 0.866025403784439, 0.5, 0.0, -0.5, -0.866025403784438, -1.0, -0.866025403784439, -0.5,
};


inline float cCos(const uint val, const uint i)
{

    return cCosf[val-2][i];
   // return cos(2*PI*i/val);

}

inline float cSin(const uint val, const uint i)
{
    return cSinf[val-2][i];  
//    return sin(2*PI*i/val);
}

// For construction
struct VS_INPUT
{
    uint n                    : BLENDINDICES0; 
};
struct VS_OUTPUT_reg
{
    Point eop                 : SV_POSITION;
    Point e0                  : POSITION1;
    Point e1                  : POSITION2;
    Point a11[4]              : TANGENT0;
};
struct VS_OUTPUT_irreg
{
    uint n                    : BLENDINDICES1;
    Point eop                 : SV_POSITION;
    Point e0                  : POSITION1;
    Point e1                  : POSITION2;
    Point a11[MAX_VALENCE]    : TANGENT0;
};
struct GS_OUTPUT
{
    float4 bc                       : POSITION;
};
#define STREAM_OUT(val) output.bc = float4(val.xyz,1), Stream.Append(output), Stream.RestartStrip()

// For evaluation
struct VS_INPUT_eval
{
    uint patchID              : BLENDINDICES0;
    float4 pos                : SV_POSITION;
};
struct VS_OUTPUT_eval
{
    float4 pos                : SV_POSITION;
    float3 normal             : TANGENT0;
    float4 color              : COLOR;
};

// For mesh display
typedef VS_INPUT VS_INPUT_mesh;
struct VS_OUTPUT_mesh
{
    Point pos                 : SV_POSITION;
    float3 normal             : TANGENT0;
    float4 color              : COLOR;
};

typedef float4 Corner[3][3];
typedef Corner Patch[4];

// To support n-gons
struct PartPatch {
    Point b112;
    Point b211, b121;
    Point b300, b210, b120; // Note: deg 3 boundary
};
typedef PartPatch Gon3Patch[3];
typedef PartPatch Gon4Patch[4];
typedef PartPatch Gon5Patch[5];

RasterizerState EnableCulling
{
    CullMode = BACK;
};


////////////////////////////////////////////////////////////////////////////////

VS_OUTPUT_reg VSConstruction_reg( VS_INPUT input, uint vID : SV_VertexID, uint iID : SV_InstanceID )
{
    uint i;
    VS_OUTPUT_reg vout;
    float4 direct[MAX_VALENCE], diag[MAX_VALENCE];
    float x, y, z, temp[2][MAX_VALENCE];
    float4 vLocation;
    uint directVal[MAX_VALENCE];

    vLocation = float4( gVertexLocation.Load(int3(vID,gRowTex,0)).xyz, 1);

    // instance 
    const float4 loc = instance_translation(iID);
    vLocation += loc;

    const uint n = 4;

    const float Af = (float)(4.0/9.0);
    const float Bf = (float)(2.0/9.0);
    const float Cf = (float)(1.0/9.0);
    const float Ae = (float)(2.0/18.0);
    const float Be = (float)(8.0/18.0);
    const float Ce = (float)(1.0/18.0);
    const float De = (float)(4.0/18.0);

    [unroll]
    for (i = 0; i < MAX_VALENCE; ++i) {
        if (i < n) {
            float ftmp = gRingIndex.Load(int3(i*3, vID, 0));
            float4 vtmp = gVertexLocation.Load(int3((uint)ftmp, gRowTex, 0));
            direct[i] = float4(vtmp.xyz, 1);
            directVal[i] = uint(vtmp.w);

            float ftmp1 = gRingIndex.Load(int3(i*3+1, vID, 0));
            float ftmp2 = gRingIndex.Load(int3(i*3+2, vID, 0));
            vtmp = 0.5 * ( gVertexLocation.Load(int3((uint)ftmp1, gRowTex, 0))
                         + gVertexLocation.Load(int3((uint)ftmp2, gRowTex, 0)) );
            diag[i] = float4(vtmp.xyz, 1);
            
            // instance
            direct[i] += loc;
            diag[i] += loc;
        }
    }

    // Compute the a11 values around the vertex.
    vout.eop = float4(0,0,0,0);
    [unroll]
    for (i = 0; i < MAX_VALENCE; ++i) {
        if (i < n) {
            uint ip = (i+1) % n;
            vout.a11[i] = Af * vLocation + Bf * (direct[i] + direct[ip]) + Cf * diag[i];

#if (BLEND)
            vout.a11[i] = (1+input.blend)*(1-input.blend) * vLocation
                          + input.blend*(1-input.blend) * 0.5 * ( vLocation + direct[i] )
                          + input.blend*(1-input.blend) * 0.5 * ( vLocation + direct[ip] )
                          + input.blend*input.blend * 0.25 * ( vLocation + direct[i] + direct[ip] + diag[i] );
#endif
            vout.eop += vout.a11[i];
        }
    }
    vout.eop = ( vout.eop*(9.0/float(n)) + (float(n)-4)*vLocation ) * (1.0/(float(n)+5));

//  vout.e0  = float4(0,0,0,0);
//  vout.e1  = float4(0,0,0,0);
//  [unroll]
//  for (i = 0; i < n; ++i) {
//      uint im = (i+n-1) % n;
//      Point e = 0.5 * ( vout.a11[i] + vout.a11[im] );
//      vout.e0 += cCos(n,i) * e;
//      vout.e1 += cSin(n,i) * e;
//  }

    // Projection
    vout.e0 = 0.5 * ( vout.a11[0] + vout.a11[3] - vout.a11[2] - vout.a11[1]);
    vout.e1 = 0.5 * ( vout.a11[1] + vout.a11[0] - vout.a11[3] - vout.a11[2]);

    const float lambda = 0.5;
    const float factor = 1.0/(lambda*n);
    vout.e0 *= factor;
    vout.e1 *= factor;

    return vout;
}

////////////////////////////////////////////////////////////////////////////////

VS_OUTPUT_irreg VSConstruction_irreg( VS_INPUT input, uint vID : SV_VertexID, uint iID : SV_InstanceID  )
{
    uint i;
    VS_OUTPUT_irreg vout;
    float4 direct[MAX_VALENCE], diag[MAX_VALENCE];
    float x, y, z, temp[2][MAX_VALENCE];
    float4 vLocation;
    uint directVal[MAX_VALENCE];
    
    vLocation = float4( gVertexLocation.Load(int3(vID,gRowTex,0)).xyz, 1);

    // instance 
    const float4 loc = instance_translation(iID);
    vLocation += loc;

    const uint n = input.n;

    const float Af = (float)(4.0/9.0);
    const float Bf = (float)(2.0/9.0);
    const float Cf = (float)(1.0/9.0);
    const float Ae = (float)(2.0/18.0);
    const float Be = (float)(8.0/18.0);
    const float Ce = (float)(1.0/18.0);
    const float De = (float)(4.0/18.0);

    vout.n   = input.n;

    [unroll]
    for (i = 0; i < MAX_VALENCE; ++i) {
        if (i < n) {
            float ftmp = gRingIndex.Load(int3(i*3, vID, 0));
            float4 vtmp = gVertexLocation.Load(int3((uint)ftmp, gRowTex, 0));
            direct[i] = float4(vtmp.xyz, 1);
            directVal[i] = uint(vtmp.w);

            float ftmp1 = gRingIndex.Load(int3(i*3+1, vID, 0));
            float ftmp2 = gRingIndex.Load(int3(i*3+2, vID, 0));
            vtmp = 0.5 * ( gVertexLocation.Load(int3((uint)ftmp1, gRowTex, 0))
                         + gVertexLocation.Load(int3((uint)ftmp2, gRowTex, 0)) );
            diag[i] = float4(vtmp.xyz, 1);

            // instance
            direct[i] += loc;
            diag[i] += loc;
        }
    }

    // Compute the a11 values around the vertex.
    vout.eop = float4(0,0,0,0);
    [unroll]
    for (i = 0; i < MAX_VALENCE; ++i) {
        if (i < n) {
            uint ip = (i+1) % n;
            vout.a11[i] = Af * vLocation + Bf * (direct[i] + direct[ip]) + Cf * diag[i];

#if (BLEND)
            vout.a11[i] = (1+input.blend)*(1-input.blend) * vLocation
                          + input.blend*(1-input.blend) * 0.5 * ( vLocation + direct[i] )
                          + input.blend*(1-input.blend) * 0.5 * ( vLocation + direct[ip] )
                          + input.blend*input.blend * 0.25 * ( vLocation + direct[i] + direct[ip] + diag[i] );
#endif
            vout.eop += vout.a11[i];
        }
    }
    vout.eop = ( vout.eop*(9.0/float(n)) + (float(n)-4)*vLocation ) * (1.0/(float(n)+5));

    const float c = cCos(n, 1);
    const float s = cSin(n, 1);
    // Check for (n==4) for exact watertight computation with regular shader
    const float lambda = (n==4) ? 0.5 : (c + 5 + sqrt((c + 9) * (c + 1))) / 16;
    const float a = (n == 3) ? 0.53 : float(0.25/lambda);

    vout.e0  = float4(0,0,0,0);
    vout.e1  = float4(0,0,0,0);
    [unroll]
    for (i = 0; i < MAX_VALENCE; ++i) {
        if (i < n) {
            uint im = (i+n-1) % n;
            Point e = 0.5 * ( vout.a11[i] + vout.a11[im] );
            vout.e0 += cCos(n,i) * e;
            vout.e1 += cSin(n,i) * e;
        }
    }
    const float factor = 1.0/lambda;
    vout.e0 *= factor / n;
    vout.e1 *= factor / n;

    return vout;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

    [maxvertexcount(16)]
void GSConstruction_polar( triangleadj VS_OUTPUT_irreg input[6], inout PointStream<GS_OUTPUT> Stream, uint p_ID : SV_PrimitiveID )
{   //can handle 1024 dwords for outstream ==> less than 128 points (pos, color) to evaluate
    uint i, j, k;

    GS_OUTPUT output; 

    VS_OUTPUT_irreg inp[3];
    inp[0] = input[0];
    inp[1] = input[1];
    inp[2] = input[2];

    const uint num_sides = 3;

    uint pID = p_ID + gStart[POLAR];
    uint rot_packed = gOffsetData.Load(int3(pID, 0, 0));
    uint rot_off[num_sides];
    [unroll]
    for (i = 0; i < num_sides; ++i)
        rot_off[i] = (rot_packed >> (NBITS_OFFSET*i)) & NBITS_OFFSET_MASK;

    Patch pat;
  
    [unroll] // unroll needed no matter what
    for (k = 0; k <= 1; ++k) {
        const uint n = 4;
        const uint off  = rot_off[k];
        const uint offp = (off + 1) % n;

        pat[k][0][0] = inp[k].eop;
        pat[k][1][0] = inp[k].eop + inp[k].e0*cCos(n,off ) + inp[k].e1*cSin(n,off );
        pat[k][0][1] = inp[k].eop + inp[k].e0*cCos(n,offp) + inp[k].e1*cSin(n,offp);
        pat[k][1][1] = inp[k].a11[off];
    }

    { // Handle creases here!
        const uint k = 2;
        const uint j  = k;
        const uint jp = (j+1) % 4;
        const uint n = inp[k].n;
        const uint off  = rot_off[k];
        const uint offp = (off + 1) % n;
        const float c = cCos(n,1);
        const float factor = float(1.0/(2.0+c));
        pat[jp][0][0] = inp[k].eop;
        pat[jp][0][1] = inp[k].eop;
        pat[j ][1][0] = inp[k].eop;
        pat[j ][0][0] = inp[k].eop;
        pat[jp][1][0] = inp[k].eop + inp[k].e0*cCos(n,off ) + inp[k].e1*cSin(n,off );
        pat[j ][0][1] = inp[k].eop + inp[k].e0*cCos(n,offp) + inp[k].e1*cSin(n,offp);
        pat[jp][1][1] = factor * (2*pat[jp][1][0] + pat[j ][0][1] + (c-1)*inp[k].eop);
        pat[j ][1][1] = factor * (2*pat[j ][0][1] + pat[jp][1][0] + (c-1)*inp[k].eop);

//      // apply blend ratios
//      float alpha0 = float(3.0/2.0)*inp[k].a11[off ].w; // should be gotten via w coordinate in GPU implementation
//      float alpha1 = float(3.0/2.0)*inp[k].a11[offp].w;
//      pat[jp][1][1] = alpha0*pat[jp][1][1] + (1-alpha0)*pat[jp][1][0];
//      pat[j ][1][1] = alpha1*pat[j ][1][1] + (1-alpha1)*pat[j ][0][1];

        // Move the singular vertex epsilon out to keep the normals
        // non-degenerate.
//      const float epsilon = 1.0e-6;
//      pat[jp][0][0] = inp[k].eop * (1 - epsilon) + pat[jp][0][0] * epsilon;
//      pat[jp][0][1] = inp[k].eop * (1 - epsilon) + pat[jp][0][1] * epsilon;
//      pat[j ][1][0] = inp[k].eop * (1 - epsilon) + pat[j ][1][0] * epsilon;
//      pat[j ][0][0] = inp[k].eop * (1 - epsilon) + pat[j ][0][0] * epsilon;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //   Stream Out Bi-3
    ////////////////////////////////////////////////////////////////////////////////

    float4 bc[4][4];

    bc[0][3]=pat[3][0][0]; bc[1][3]=pat[3][0][1]; bc[2][3]=pat[2][1][0]; bc[3][3]=pat[2][0][0];
    bc[0][2]=pat[3][1][0]; bc[1][2]=pat[3][1][1]; bc[2][2]=pat[2][1][1]; bc[3][2]=pat[2][0][1];   
    bc[0][1]=pat[0][0][1]; bc[1][1]=pat[0][1][1]; bc[2][1]=pat[1][1][1]; bc[3][1]=pat[1][1][0];
    bc[0][0]=pat[0][0][0]; bc[1][0]=pat[0][1][0]; bc[2][0]=pat[1][0][1]; bc[3][0]=pat[1][0][0];

    ///////////////////////////////////////////////////////////
    // degree lowering from 5 to 3 for b11
    ///////////////////////////////////////////////////////////

    // Stream out
    [unroll]
    for (i=0; i<4; i++) 
        for (j=0; j<4; j++)
            STREAM_OUT(bc[i][j]);
}

////////////////////////////////////////////////////////////////////////////////

    [maxvertexcount(16)]
void GSConstruction_reg( triangleadj VS_OUTPUT_reg input[6], inout PointStream<GS_OUTPUT> Stream, uint p_ID : SV_PrimitiveID )
{   //can handle 1024 dwords for outstream ==> less than 128 points (pos, color) to evaluate
    uint i, j, k;

    GS_OUTPUT output; 

    VS_OUTPUT_reg inp[4];
    inp[0] = input[0];
    inp[1] = input[1];
    inp[2] = input[2];
    inp[3] = input[3];

    const uint num_sides = 4;
    const uint n = 4;

    uint pID = p_ID + gStart[REG];
    uint rot_packed = gOffsetData.Load(int3(pID, 0, 0));
    uint rot_off[num_sides];
    [unroll]
    for (i = 0; i < num_sides; ++i)
        rot_off[i] = (rot_packed >> (NBITS_OFFSET*i)) & NBITS_OFFSET_MASK;

    Patch pat;
  
    [unroll] // unroll needed no matter what
    for (k = 0; k < num_sides; ++k) {
        const uint off  = rot_off[k];
        const uint offp = (off + 1) % n;

        pat[k][0][0] = inp[k].eop;
        pat[k][1][0] = inp[k].eop + inp[k].e0*cCos(n,off ) + inp[k].e1*cSin(n,off );
        pat[k][0][1] = inp[k].eop + inp[k].e0*cCos(n,offp) + inp[k].e1*cSin(n,offp);
        pat[k][1][1] = inp[k].a11[off];
    }

    ////////////////////////////////////////////////////////////////////////////////
    //   Stream Out Bi-3
    ////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////
    // Hybrid method, cubic construction same with regular patch (start)
    ///////////////////////////////////////////////////////////////////

    float4 bc[4][4];

    bc[0][3]=pat[3][0][0]; bc[1][3]=pat[3][0][1]; bc[2][3]=pat[2][1][0]; bc[3][3]=pat[2][0][0];
    bc[0][2]=pat[3][1][0]; bc[1][2]=pat[3][1][1]; bc[2][2]=pat[2][1][1]; bc[3][2]=pat[2][0][1];   
    bc[0][1]=pat[0][0][1]; bc[1][1]=pat[0][1][1]; bc[2][1]=pat[1][1][1]; bc[3][1]=pat[1][1][0];
    bc[0][0]=pat[0][0][0]; bc[1][0]=pat[0][1][0]; bc[2][0]=pat[1][0][1]; bc[3][0]=pat[1][0][0];

    ///////////////////////////////////////////////////////////
    // degree lowering from 5 to 3 for b11
    ///////////////////////////////////////////////////////////

    // Stream out
    [unroll]
    for (i=0; i<4; i++) 
        for (j=0; j<4; j++)
            STREAM_OUT(bc[i][j]);
}

////////////////////////////////////////////////////////////////////////////////

    [maxvertexcount(19)]
void GSConstruction3( triangleadj VS_OUTPUT_irreg input[6], inout PointStream<GS_OUTPUT> Stream, uint p_ID : SV_PrimitiveID )
{   //can handle 1024 dwords for outstream ==> less than 128 points (pos, color) to evaluate
    uint i, k;

    GS_OUTPUT output; 

    const uint num_sides = 3;

    VS_OUTPUT_irreg inp[num_sides];
    for (i = 0; i < num_sides; ++i)
        inp[i] = input[i];

    uint pID = p_ID + gStart[num_sides];

    uint rot_packed = gOffsetData.Load(int3(pID, 0, 0));
    uint rot_off[num_sides];
    [unroll]
    for (i = 0; i < num_sides; ++i)
        rot_off[i] = (rot_packed >> (NBITS_OFFSET*i)) & NBITS_OFFSET_MASK;

    Gon4Patch pat;

    // Compute the corner tangential control points
    [unroll] // unroll needed no matter what
    for (k = 0; k < num_sides; ++k) {
        const uint km = (k+num_sides-1) % num_sides;

        const uint n = inp[k].n;
        const uint off  = rot_off[k];
        const uint offp = (off+1)%n;

        pat[k ].b300 = inp[k].eop;
        pat[k ].b210 = inp[k].eop + inp[k].e0*cCos(n,off ) + inp[k].e1*cSin(n,off );
        pat[km].b120 = inp[k].eop + inp[k].e0*cCos(n,offp) + inp[k].e1*cSin(n,offp);
    }

    // Compute the edges.
    Point b004 = float4(0,0,0,0);
    const float mu = (1 - cCos(num_sides,1));
    const float mu_inv = 1 / mu;
    const float gamma = 2*mu;
    const float weight = WEIGHT[num_sides]; //gamma-1;

    [unroll]
    for (k = 0; k < num_sides; ++k) {
        const uint kp = (k+1) % num_sides;
        const uint n0 = inp[k ].n;
        const uint n1 = inp[kp].n;

        const float c0 = cCos(n0,1);
        const float c1 = cCos(n1,1);
        const float ss_inv = 1 / (cSin(n0,1) + cSin(n1,1));

        const uint off0 = rot_off[k ];
        const uint off1 = rot_off[kp];

        const float4 b310 = (3.0/4.0) * pat[k].b210 + (1.0/4.0) * pat[k ].b300;
        const float4 b130 = (3.0/4.0) * pat[k].b120 + (1.0/4.0) * pat[kp].b300;
        const float4 a11  = inp[k ].a11[ off0         ];
        const float4 a21  = inp[kp].a11[ off1         ];
        const float4 a11m = inp[k ].a11[(off0-1+n0)%n0];
        const float4 a21m = inp[kp].a11[(off1+1   )%n1];

        pat[k].b211 = b310 + ( (6.0/24.0) * mu_inv * (      1+c0) ) * (pat[k].b120 - pat[k].b210)
                           + ( (3.0/24.0) * mu_inv * (gamma-1-c1) ) * (pat[k].b210 - pat[k].b300)
                           + ( (3.0/ 8.0) * mu_inv * ss_inv ) * (a11 - a11m);

        pat[k].b121 = b130 + ( (6.0/24.0) * mu_inv * (      1+c1) ) * (pat[k].b210 - pat[k ].b120)
                           + ( (3.0/24.0) * mu_inv * (gamma-1-c0) ) * (pat[k].b120 - pat[kp].b300)
                           + ( (3.0/ 8.0) * mu_inv * ss_inv ) * (a21 - a21m);

        // b004 += pat[k].b300 + (3.0/2.0) * (a11 + a11m + a21 + a21m) + 9 * a11;
        // Equivalent to equation below due to the contribution of a11 and
        // the previous iteration's a21.
        b004 += weight*pat[k].b300 + (3.0/2.0) * (a11m + a21m) + 12 * a11;
    }
    b004 /= (15+weight)*num_sides;

    // Compute b112 for 4-sided extraordinary patches
    [unroll]
    for (k = 0; k < num_sides; ++k) {
        const uint km  = (k+(num_sides-1)) % num_sides;
        const uint kp  = (k+1) % num_sides;
        pat[k].b112 =
              (3.0/2.0 ) * b004
            - (1.0/12.0) * pat[km].b300
            - (1.0/24.0) * ( pat[km].b210 + pat[kp].b120 )
            - (1.0/6.0 ) * ( pat[km].b211 + pat[kp].b121 );
    }

    ////////////////////////////////////////////////////////////////////////////////
    //   Stream Out
    ////////////////////////////////////////////////////////////////////////////////

    STREAM_OUT(b004);
    [unroll]
    for (k = 0; k < num_sides; ++k) {
        STREAM_OUT(pat[k].b300);
        STREAM_OUT(pat[k].b210);
        STREAM_OUT(pat[k].b120);
        STREAM_OUT(pat[k].b211);
        STREAM_OUT(pat[k].b121);
        STREAM_OUT(pat[k].b112);
    }
}

////////////////////////////////////////////////////////////////////////////////

    [maxvertexcount(25)]
void GSConstruction4( triangleadj VS_OUTPUT_irreg input[6], inout PointStream<GS_OUTPUT> Stream, uint p_ID : SV_PrimitiveID )
{   //can handle 1024 dwords for outstream ==> less than 128 points (pos, color) to evaluate
    uint i, k;

    GS_OUTPUT output; 

    const uint num_sides = 4;

    VS_OUTPUT_irreg inp[num_sides];
    for (i = 0; i < num_sides; ++i)
        inp[i] = input[i];

    uint pID = p_ID + gStart[num_sides];

    uint rot_packed = gOffsetData.Load(int3(pID, 0, 0));
    uint rot_off[num_sides];
    [unroll]
    for (i = 0; i < num_sides; ++i)
        rot_off[i] = (rot_packed >> (NBITS_OFFSET*i)) & NBITS_OFFSET_MASK;

    Gon4Patch pat;

    // Compute the corner tangential control points
    [unroll] // unroll needed no matter what
    for (k = 0; k < num_sides; ++k) {
        const uint km = (k+num_sides-1) % num_sides;

        const uint n = inp[k].n;
        const uint off  = rot_off[k];
        const uint offp = (off+1)%n;

        pat[k ].b300 = inp[k].eop;
        pat[k ].b210 = inp[k].eop + inp[k].e0*cCos(n,off ) + inp[k].e1*cSin(n,off );
        pat[km].b120 = inp[k].eop + inp[k].e0*cCos(n,offp) + inp[k].e1*cSin(n,offp);
    }

    // Compute the edges.
    Point b004 = float4(0,0,0,0);
    const float mu = (1 - cCos(num_sides,1));
    const float mu_inv = 1 / mu;
    const float gamma = 2*mu;
    const float weight = WEIGHT[num_sides]; //gamma-1;

    [unroll]
    for (k = 0; k < num_sides; ++k) {
        const uint kp = (k+1) % num_sides;
        const uint n0 = inp[k ].n;
        const uint n1 = inp[kp].n;

        const float c0 = cCos(n0,1);
        const float c1 = cCos(n1,1);
        const float ss_inv = 1 / (cSin(n0,1) + cSin(n1,1));

        const uint off0 = rot_off[k ];
        const uint off1 = rot_off[kp];

        const float4 b310 = (3.0/4.0) * pat[k].b210 + (1.0/4.0) * pat[k ].b300;
        const float4 b130 = (3.0/4.0) * pat[k].b120 + (1.0/4.0) * pat[kp].b300;
        const float4 a11  = inp[k ].a11[ off0         ];
        const float4 a21  = inp[kp].a11[ off1         ];
        const float4 a11m = inp[k ].a11[(off0-1+n0)%n0];
        const float4 a21m = inp[kp].a11[(off1+1   )%n1];

        pat[k].b211 = b310 + ( (6.0/24.0) * mu_inv * (      1+c0) ) * (pat[k].b120 - pat[k].b210)
                           + ( (3.0/24.0) * mu_inv * (gamma-1-c1) ) * (pat[k].b210 - pat[k].b300)
                           + ( (3.0/ 8.0) * mu_inv * ss_inv ) * (a11 - a11m);

        pat[k].b121 = b130 + ( (6.0/24.0) * mu_inv * (      1+c1) ) * (pat[k].b210 - pat[k ].b120)
                           + ( (3.0/24.0) * mu_inv * (gamma-1-c0) ) * (pat[k].b120 - pat[kp].b300)
                           + ( (3.0/ 8.0) * mu_inv * ss_inv ) * (a21 - a21m);

        // b004 += pat[k].b300 + (3.0/2.0) * (a11 + a11m + a21 + a21m) + 9 * a11;
        // Equivalent to equation below due to the contribution of a11 and
        // the previous iteration's a21.
        b004 += weight*pat[k].b300 + (3.0/2.0) * (a11m + a21m) + 12 * a11;
    }
    b004 /= (15+weight)*num_sides;

    // Compute b112 for 4-sided extraordinary patches
    [unroll]
    for (k = 0; k < num_sides; ++k) {
        const uint km  = (k+num_sides-1) % num_sides;
        const uint km2 = (k+num_sides-2) % num_sides;
        const uint kp  = (k+1) % num_sides;
        const uint kp2 = (k+2) % num_sides;
        pat[k].b112 = b004
            + (3.0/16.0) * (pat[k ].b211 + pat[k ].b121 - pat[kp ].b121 - pat[km ].b211)
            + (1.0/16.0) * (pat[kp].b211 + pat[km].b121 - pat[kp2].b211 - pat[km2].b121);
    }

    ////////////////////////////////////////////////////////////////////////////////
    //   Stream Out
    ////////////////////////////////////////////////////////////////////////////////

    STREAM_OUT(b004);
    [unroll]
    for (k = 0; k < num_sides; ++k) {
        STREAM_OUT(pat[k].b300);
        STREAM_OUT(pat[k].b210);
        STREAM_OUT(pat[k].b120);
        STREAM_OUT(pat[k].b211);
        STREAM_OUT(pat[k].b121);
        STREAM_OUT(pat[k].b112);
    }
}

////////////////////////////////////////////////////////////////////////////////

    [maxvertexcount(31)]
void GSConstruction5( triangleadj VS_OUTPUT_irreg input[6], inout PointStream<GS_OUTPUT> Stream, uint p_ID : SV_PrimitiveID )
{   //can handle 1024 dwords for outstream ==> less than 128 points (pos, color) to evaluate
    uint i, k;

    GS_OUTPUT output; 

    const uint num_sides = 5;

    VS_OUTPUT_irreg inp[num_sides];
    for (i = 0; i < num_sides; ++i)
        inp[i] = input[i];

    uint pID = p_ID + gStart[num_sides];

    uint rot_packed = gOffsetData.Load(int3(pID, 0, 0));
    uint rot_off[num_sides];
    [unroll]
    for (i = 0; i < num_sides; ++i)
        rot_off[i] = (rot_packed >> (NBITS_OFFSET*i)) & NBITS_OFFSET_MASK;

    Gon5Patch pat;

    // Compute the corner tangential control points
    [unroll] // unroll needed no matter what
    for (k = 0; k < num_sides; ++k) {
        const uint km = (k+num_sides-1) % num_sides;

        const uint n = inp[k].n;
        const uint off  = rot_off[k];
        const uint offp = (off+1)%n;

        pat[k ].b300 = inp[k].eop;
        pat[k ].b210 = inp[k].eop + inp[k].e0*cCos(n,off ) + inp[k].e1*cSin(n,off );
        pat[km].b120 = inp[k].eop + inp[k].e0*cCos(n,offp) + inp[k].e1*cSin(n,offp);
    }

    // Compute the edges.
    Point b004 = float4(0,0,0,0);
    const float mu = (1 - cCos(num_sides,1));
    const float mu_inv = 1 / mu;
    const float gamma = 2*mu;
    const float weight = WEIGHT[num_sides]; //gamma-1;

    [unroll]
    for (k = 0; k < num_sides; ++k) {
        const uint kp = (k+1) % num_sides;
        const uint n0 = inp[k ].n;
        const uint n1 = inp[kp].n;

        const float c0 = cCos(n0,1);
        const float c1 = cCos(n1,1);
        const float ss_inv = 1 / (cSin(n0,1) + cSin(n1,1));

        const uint off0 = rot_off[k ];
        const uint off1 = rot_off[kp];

        const float4 b310 = (3.0/4.0) * pat[k].b210 + (1.0/4.0) * pat[k ].b300;
        const float4 b130 = (3.0/4.0) * pat[k].b120 + (1.0/4.0) * pat[kp].b300;
        const float4 a11  = inp[k ].a11[ off0         ];
        const float4 a21  = inp[kp].a11[ off1         ];
        const float4 a11m = inp[k ].a11[(off0-1+n0)%n0];
        const float4 a21m = inp[kp].a11[(off1+1   )%n1];

        pat[k].b211 = b310 + ( (6.0/24.0) * mu_inv * (      1+c0) ) * (pat[k].b120 - pat[k].b210)
                           + ( (3.0/24.0) * mu_inv * (gamma-1-c1) ) * (pat[k].b210 - pat[k].b300)
                           + ( (3.0/ 8.0) * mu_inv * ss_inv ) * (a11 - a11m);

        pat[k].b121 = b130 + ( (6.0/24.0) * mu_inv * (      1+c1) ) * (pat[k].b210 - pat[k ].b120)
                           + ( (3.0/24.0) * mu_inv * (gamma-1-c0) ) * (pat[k].b120 - pat[kp].b300)
                           + ( (3.0/ 8.0) * mu_inv * ss_inv ) * (a21 - a21m);

        // b004 += pat[k].b300 + (3.0/2.0) * (a11 + a11m + a21 + a21m) + 9 * a11;
        // Equivalent to equation below due to the contribution of a11 and
        // the previous iteration's a21.
        b004 += weight*pat[k].b300 + (3.0/2.0) * (a11m + a21m) + 12 * a11;
    }
    b004 /= (15+weight)*num_sides;

    float4 b202[num_sides];
    const float cN = cCos(num_sides,1);
    const float k2 = float(0.5) / (1 - cN);
    const float k1 = 1 - 2*k2;
    for (k = 0; k < num_sides; ++k) {
        const uint km = (k+num_sides-1) % num_sides;
        b202[k] = (k1*k1 + float(2.0/4.0)*k1*k2) * pat[k].b300
                +          float(3.0/4.0)*k1*k2  * (pat[k].b210 + pat[km].b120)
                +                           k2  * (pat[k].b211 + pat[km].b121);
    }
    const float c2N = cCos(num_sides,2);
    for (k = 0; k < num_sides; ++k) {
        const uint kp1 = (k+1) % num_sides;
        const uint kp2 = (k+2) % num_sides;
        const uint kp3 = (k+3) % num_sides;
        const uint kp4 = (k+4) % num_sides;
        pat[k].b112 = (1-cN) * (b004 +
                float(1.0/5.0) * (
                    b202[kp3]
                    -4*c2N     * ( b202[k  ] + b202[kp1] )
                    -4*c2N*c2N * ( b202[kp2] + b202[kp4] )
                    )
                );
    }

/*
    // Compute b112 for odd-sided extraordinary patches
    float4 b202[num_sides];
    float4 r11[num_sides];
    float4 e0 = float4(0,0,0,0);
    float4 e1 = float4(0,0,0,0);
    const float cN = cCos(num_sides,1);
    const float k2 = float(0.5) / (1 - cN);
    const float k1 = 1 - 2*k2;

    [unroll]
    for (k = 0; k < num_sides; ++k) {
        const uint km = (k+num_sides-1) % num_sides;
        //b202[k] = k1*k1 *  pat[k].b400
        //        + k1*k2 * (pat[k].b310 + pat[km].b130)
        //        +    k2 * (pat[k].b211 + pat[km].b121);
        // Equivalent to the formula above, but using variables that
        // actually exist.
        b202[k] = (k1*k1 + float(2.0/4.0)*k1*k2) * pat[k].b300
                +          float(3.0/4.0)*k1*k2  * (pat[k].b210 + pat[km].b120)
                +                           k2  * (pat[k].b211 + pat[km].b121);

        e0 += b202[k] * cCos(num_sides,k);
        e1 += b202[k] * cSin(num_sides,k);
    }
    e0 /= num_sides;
    e1 /= num_sides;

    [unroll]
    for (k = 0; k < num_sides; ++k) {
        const float4 b103 = b004 + e0*cCos(num_sides,k) + e1*cSin(num_sides,k);
        r11[k] = (1-cN) * b103 + cN * b202[k];
    }

    [unroll]
    for (k = 0; k < num_sides; ++k) {
        pat[k].b112 = float4(0,0,0,0);
        for (j = 0; j < num_sides; ++j) {
            const float sgn = ( ((num_sides+k-j)%num_sides) & 1 ) ? -1 : 1;
            pat[k].b112 += sgn * r11[j];
        }
    }
*/
    ////////////////////////////////////////////////////////////////////////////////
    //   Stream Out
    ////////////////////////////////////////////////////////////////////////////////

    STREAM_OUT(b004);
    [unroll]
    for (k = 0; k < num_sides; ++k) {
        STREAM_OUT(pat[k].b300);
        STREAM_OUT(pat[k].b210);
        STREAM_OUT(pat[k].b120);
        STREAM_OUT(pat[k].b211);
        STREAM_OUT(pat[k].b121);
        STREAM_OUT(pat[k].b112);
    }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// EVALUATION PASS
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

float4
compute_lighting(float3 pos, float3 normal, float4 color)
{
/*
    // Silver 
    float4 gAmbientMtrl = float4(0.19225f, 0.19225f, 0.19225f, 1.0f);
    float4 gAmbientLight = float4(1.0f, 1.0f, 1.0f, 0.0f);
    float4 gDiffuseMtrl = float4(0.50754f, 0.50754f, 0.50754f, 1.0f); 
    float4 gDiffuseLight = float4(1.0f, 1.0f, 1.0f, 1.0f);
    float4 gSpecularMtrl = float4(0.50823f, 0.50823f, 0.50823f, 1.0f);
    float4 gSpecularLight = float4(1.0f, 1.0f, 1.0f, 1.0f);
    float  gSpecularPower = 20;
*/

    float4 gAmbientMtrl = float4(0.2f, 0.2f, 0.2f, 1.0f);
    float4 gDiffuseMtrl = color; // float4(0.55f, 0.47f, 0.14f, 1.0f);
    float4 gSpecularMtrl = float4(1.0f, 1.0f, 1.0f, 1.0f);
    float  gSpecularPower = 30;

    float4 gAmbientLight = float4(1.0f, 1.0f, 1.0f, 1.0f);
    float4 gDiffuseLight = float4(1.0f, 1.0f, 1.0f, 1.0f);
    float4 gSpecularLight = float4(1.0f, 1.0f, 1.0f, 1.0f);
 
    // Compute the vector from the vertex to the eye position.
    float3 toEye = normalize(gEyePosW - pos);

    // Compute the reflection vector.
    float3 r1 = reflect(-gLightVecW1, normal);

    // Determine how much (if any) specular light makes it into the eye.
    float t1  = pow(max(dot(r1, toEye), 0.0f), gSpecularPower);

    // Determine the diffuse light intensity that strikes the vertex.
    float s1 = max(dot(gLightVecW1, normal), 0.0f);

    // Compute the ambient, diffuse and specular terms separatly. 
    float3 spec = t1*(gSpecularMtrl*gSpecularLight).rgb;
    float3 diffuse = s1*(gDiffuseMtrl*gDiffuseLight).rgb;

    float3 ambient = gAmbientMtrl*gAmbientLight;

    float4 finalColor = float4 (ambient+diffuse+spec, 1); 
    //float4 finalColor = ambient+diffuse+spec; 

    return finalColor;
}

////////////////////////////////////////////////////////////

float4
compute_color(uint CASE, float4 pos, float3 normal) {
    float4 color = (gFlags & (1 << BIT_GROUP_COLOR)) ? gDiffuseMtrl[CASE] : gDiffuseMtrl[DEFAULT];
#if (PHONG)
    return color;
#else
    return compute_lighting(pos, normal, color);
#endif
}

////////////////////////////////////////////////////////////////////////////////

VS_OUTPUT_eval VSEvaluation_cubic( uint TYPE, Buffer<float4> ControlPoints, VS_INPUT_eval input,  uint vID : SV_InstanceID )
{
    uint i, j;
    VS_OUTPUT_eval output;
    float4 pos;
    float3 normal;

    pos = ControlPoints.Load(MAX_CONTROL_POINTS[TYPE]*vID);
//  uint type = (uint)pos.w;

    float uc[6], vc[6], ud[6], vd[6];

    float4 su=float4(0,0,0,0);
    float4 sv=float4(0,0,0,0);
    float3 p =float3(0,0,0);
    //uint beginIndex=MAX_CONTROL_POINTS[TYPE]*input.patchID;         
    uint beginIndex=MAX_CONTROL_POINTS[TYPE]*vID;    

    float um1 = (1-input.pos.x);
    float um2 = um1*um1;
    float um3 = um2*um1;

    float u1 = input.pos.x;
    float u2 = u1*u1;
    float u3 = u2*u1;

    float vm1 = (1-input.pos.y);
    float vm2 = vm1*vm1;
    float vm3 = vm2*vm1;

    float v1 = input.pos.y;
    float v2 = v1*v1;
    float v3 = v2*v1;

    uc[0] = um3           ;
    uc[1] = um2 * u1 * 3.0;
    uc[2] = um1 * u2 * 3.0;
    uc[3] =       u3      ;

    vc[0] = vm3           ;
    vc[1] = vm2 * v1 * 3.0;
    vc[2] = vm1 * v2 * 3.0;
    vc[3] =       v3      ;

    ud[0] = um2                   *-3.0;
    ud[1] = um1 * (1 - 3*u1)      * 3.0;
    ud[2] =       (2 - 3*u1) * u1 * 3.0;
    ud[3] =                    u2 * 3.0;

    vd[0] = vm2                   *-3.0;
    vd[1] = vm1 * (1 - 3*v1)      * 3.0;
    vd[2] =       (2 - 3*v1) * v1 * 3.0;
    vd[3] =                    v2 * 3.0;

    float4 cpts[4][4];

    [unroll]
    for (i=0; i<4; i++) 
        for (j=0; j<4; j++)  {
            cpts[i][j] = ControlPoints.Load(beginIndex + i*4 + j); 
            p  +=   uc[i]*vc[j]*cpts[i][j];
        }

    pos = float4(p[0], p[1], p[2], 1.0);


    [unroll]
    for (i=0; i<4; i++) 
        for (j=0; j<4; j++) {
            su +=  ud[i]*vc[j]*cpts[i][j]; 
            sv +=  uc[i]*vd[j]*cpts[i][j]; 
        }

    [branch]
    if (TYPE == POLAR && v1 == 1.0) {
        su = cpts[0][2] - cpts[3][3];
        sv = cpts[3][2] - cpts[3][3];
    }

    normal=cross(su, sv);  

    pos = mul( pos, matWorldViewProj ); 
    normal = normalize(mul(normal, (float3x3)matWorldView)); 

    output.pos=pos;
    output.normal = normal;
    output.color = compute_color(TYPE, pos, normal);

    return output;
}

////////////////////////////////////////

VS_OUTPUT_eval VSEvaluation_reg(VS_INPUT_eval input,  uint vID : SV_InstanceID )
{
    return VSEvaluation_cubic( REG, gControlPoints_reg, input, vID );
}

VS_OUTPUT_eval VSEvaluation_polar(VS_INPUT_eval input,  uint vID : SV_InstanceID )
{
    VS_INPUT_eval inp = input;
    // Map the unit triangle domain to the unit quad domain.
    inp.pos.x = (inp.pos.y == 1.0) ? 0 : (input.pos.x / (1-input.pos.y));
    return VSEvaluation_cubic( POLAR, gControlPoints_polar, inp, vID );
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// Helper function for evaluation
struct TriEval {
    float4 P, Pu, Pv;
    float4 Puu, Puv, Pvv; // Typically used for bump mapping
};
inline TriEval
evaluate_triangle( Buffer<float4> ControlPoints, uint vID, uint num_sides, uint pnum,
                   float u, float v, float w )
{
    TriEval output;

    uint beginIndex = MAX_CONTROL_POINTS[num_sides]*vID;    

    //           (004)
    //           14   <-- (u,v,w) = (0,0,1)
    //          12 13
    //          9 10 11
    //         5  6  7  8
    // (300)  0  1  2  3  4 (030)
    //
    // Points stored in the order:
    // 004 ... ( 300 210 120 211 121 112 ) x num_sides

    // Load the control points needed for computation.
    const uint next = (pnum + 1            ) % num_sides;
    const uint prev = (pnum - 1 + num_sides) % num_sides;
    const float4 b004  = ControlPoints.Load(beginIndex             );
    const float4 b300  = ControlPoints.Load(beginIndex + 6*pnum + 1);
    const float4 b210  = ControlPoints.Load(beginIndex + 6*pnum + 2);
    const float4 b120  = ControlPoints.Load(beginIndex + 6*pnum + 3);
    const float4 b211  = ControlPoints.Load(beginIndex + 6*pnum + 4);
    const float4 b121  = ControlPoints.Load(beginIndex + 6*pnum + 5);
    const float4 b112  = ControlPoints.Load(beginIndex + 6*pnum + 6);

    const float4 bm120 = ControlPoints.Load(beginIndex + 6*prev + 3);
    const float4 bm121 = ControlPoints.Load(beginIndex + 6*prev + 5);
    const float4 bm112 = ControlPoints.Load(beginIndex + 6*prev + 6);

    const float4 bp300 = ControlPoints.Load(beginIndex + 6*next + 1);
    const float4 bp210 = ControlPoints.Load(beginIndex + 6*next + 2);
    const float4 bp211 = ControlPoints.Load(beginIndex + 6*next + 4);
    const float4 bp112 = ControlPoints.Load(beginIndex + 6*next + 6);

    const float4 bm130 = 0.25 *  b300 + 0.75 * bm120;
    const float4 bp310 = 0.25 * bp300 + 0.75 * bp210;

    // Compute the control points of the patch to be evaluated.
    const float cN = cCos(num_sides, 1);
    const float k2 = (0.5 / (1 - cN));
    const float k1 = (1.0 - 2 * k2);
    float4 rc[15];

    rc[ 0] = b300;
    rc[ 1] = (0.25) * b300 + (0.75) * b210;
    rc[ 2] = (0.5 ) * b210 + (0.5 ) * b120;
    rc[ 3] = (0.75) * b120 + (0.25) * bp300;
    rc[ 4] = bp300;

    rc[ 5] = k1 * rc[ 0] + k2 * (rc[ 1] + bm130);
    rc[ 6] = b211;
    rc[ 7] = b121;
    rc[ 8] = k1 * rc[ 4] + k2 * (rc[ 3] + bp310);
    rc[ 9] = k1 * rc[ 5] + k2 * (rc[ 6] + bm121);
    rc[10] = b112;
    rc[11] = k1 * rc[ 8] + k2 * (rc[ 7] + bp211);
    rc[12] = k1 * rc[ 9] + k2 * (rc[10] + bm112);
    rc[13] = k1 * rc[11] + k2 * (rc[10] + bp112);
    rc[14] = b004;

            // Tianyun's deCasteljau
            float4 tmp[3][10], Ps, Pss, Pt, Ptt, Pst;
            tmp[0][0] = u*rc[ 0] + v*rc[ 1]+ w*rc[ 5];
            tmp[0][1] = u*rc[ 1] + v*rc[ 2]+ w*rc[ 6];
            tmp[0][2] = u*rc[ 2] + v*rc[ 3]+ w*rc[ 7];
            tmp[0][3] = u*rc[ 3] + v*rc[ 4]+ w*rc[ 8];
            tmp[0][4] = u*rc[ 5] + v*rc[ 6]+ w*rc[ 9];
            tmp[0][5] = u*rc[ 6] + v*rc[ 7]+ w*rc[10];
            tmp[0][6] = u*rc[ 7] + v*rc[ 8]+ w*rc[11];
            tmp[0][7] = u*rc[ 9] + v*rc[10]+ w*rc[12];
            tmp[0][8] = u*rc[10] + v*rc[11]+ w*rc[13];
            tmp[0][9] = u*rc[12] + v*rc[13]+ w*rc[14];

            tmp[1][0] = u*tmp[0][0] + v*tmp[0][1]+ w*tmp[0][4];
            tmp[1][1] = u*tmp[0][1] + v*tmp[0][2]+ w*tmp[0][5];
            tmp[1][2] = u*tmp[0][2] + v*tmp[0][3]+ w*tmp[0][6];
            tmp[1][3] = u*tmp[0][4] + v*tmp[0][5]+ w*tmp[0][7];
            tmp[1][4] = u*tmp[0][5] + v*tmp[0][6]+ w*tmp[0][8];
            tmp[1][5] = u*tmp[0][7] + v*tmp[0][8]+ w*tmp[0][9];


            tmp[2][0] = u*tmp[1][0] + v*tmp[1][1]+ w*tmp[1][3];
            tmp[2][1] = u*tmp[1][1] + v*tmp[1][2]+ w*tmp[1][4];
            tmp[2][2] = u*tmp[1][3] + v*tmp[1][4]+ w*tmp[1][5];

            output.P  = u*tmp[2][0]+v*tmp[2][1]+w*tmp[2][2];

            output.Pu = 4*(tmp[2][0] - tmp[2][2]);
            output.Pv = 4*(tmp[2][1] - tmp[2][2]);
            output.Puu = 12*(tmp[1][0] - 2*tmp[1][3]);
            output.Puu = output.Puu + 12*tmp[1][5];
            output.Pvv = 12*(tmp[1][2] - 2*tmp[1][4]);
            output.Pvv = output.Pvv + 12*tmp[1][5];
            output.Puv = 12*(tmp[1][1] - tmp[1][4] - tmp[1][3] + tmp[1][5]);

    // In the boundary case for the sake of watertightness.
    [branch]
    if (w == 0) {
        // The v-coordinate is what needs to be used to evaluate at the
        // boundary. Since u=1-v, it is also used here.
        output.P = ( (  u*u*u)*b300 + (  v*v*v)*bp300 )
                 + ( (3*u*u*v)*b210 + (3*u*v*v)*b120 );
    }

    return output;
}

////////////////////////////////////////////////////////////////////////////////

VS_OUTPUT_eval VSEvaluation3( VS_INPUT_eval input,  uint vID : SV_InstanceID )
{
    VS_OUTPUT_eval output;

    const uint num_sides = 3;
    float u, v, w;

    // Figure out which sector needs to be evaluated.
    // For the triangle furthest away from (u,v,w) = (0,0,1), the transformation is
    //  [ u' ]   [ 1 0 -1 ] [ u ]
    //  [ v' ] = [ 0 1 -1 ] [ v ]
    //  [ w' ]   [ 0 0  3 ] [ w ]
    uint pnum;
    [branch]
    if (input.pos.z <= input.pos.y && input.pos.z < input.pos.x) {
        u = input.pos.x - input.pos.z;
        v = input.pos.y - input.pos.z;
        w = 3*input.pos.z;
        pnum = 0;
    }
    else if (input.pos.x <= input.pos.y) {
        u = input.pos.y - input.pos.x;
        v = input.pos.z - input.pos.x;
        w = 3*input.pos.x;
        pnum = 1;
    }
    else {
        u = input.pos.z - input.pos.y;
        v = input.pos.x - input.pos.y;
        w = 3*input.pos.y;
        pnum = 2;
    }

    float4 pos;
    float3 normal;

    // Evaluate the patch
    TriEval res = evaluate_triangle( gControlPoints3, vID, num_sides, pnum, u, v, w );
    pos = float4(res.P.xyz,1);
    normal = cross(res.Pu.xyz, res.Pv.xyz);  

    // Compute the final projected positions and normals
    pos = mul( pos, matWorldViewProj ); 
    normal = normalize(mul(normal, (float3x3)matWorldView)); 

    output.pos = pos;
    output.normal = normal;
    output.color = compute_color(num_sides, pos, normal);

    return output;
}

////////////////////////////////////////////////////////////

VS_OUTPUT_eval VSEvaluation4( VS_INPUT_eval input,  uint vID : SV_InstanceID )
{
    VS_OUTPUT_eval output;

    const uint num_sides = 4;
    float u, v, w;

    // Figure out which sector needs to be evaluated.
    uint pnum;
    [branch]
    if (input.pos.y <= input.pos.x && input.pos.x + input.pos.y <= 1) { //T1
        u = 1.0 - input.pos.x - input.pos.y;
        v = input.pos.x - input.pos.y;
        w = 2 * input.pos.y;
        pnum = 0;
    }
    else if (input.pos.y <= input.pos.x && input.pos.x + input.pos.y > 1) { //T2
        u = input.pos.x - input.pos.y;
        v = input.pos.x + input.pos.y - 1;
        w = 2 * (1 - input.pos.x);
        pnum = 1;
    }
    else if (input.pos.y > input.pos.x && input.pos.x + input.pos.y <= 1) { //T4
        u = input.pos.y - input.pos.x;
        v = 1 - input.pos.y - input.pos.x;
        w = 2 * input.pos.x;
        pnum = 3;
    }
    else { //T3
        u = input.pos.x + input.pos.y - 1;
        v = input.pos.y - input.pos.x;
        w = 2 * (1 - input.pos.y);
        pnum = 2;
    }

    float4 pos;
    float3 normal;

    // Evaluate the patch
    TriEval res = evaluate_triangle( gControlPoints4, vID, num_sides, pnum, u, v, w );
    pos = float4(res.P.xyz,1);
    normal = cross(res.Pu.xyz, res.Pv.xyz);  

    // Compute the final projected positions and normals
    pos = mul( pos, matWorldViewProj ); 
    normal = normalize(mul(normal, (float3x3)matWorldView)); 

    output.pos = pos;
    output.normal = normal;
    output.color = compute_color(num_sides, pos, normal);

    return output;
}

////////////////////////////////////////////////////////////

VS_OUTPUT_eval VSEvaluation5( VS_INPUT_eval input,  uint iID : SV_InstanceID )
{
    VS_OUTPUT_eval output;

    const uint num_sides = 5;
    float u = input.pos.x;
    float v = input.pos.y;
    float w = input.pos.z;

    const uint vID = iID / num_sides;
    const uint pnum = iID % num_sides;

    float4 pos;
    float3 normal;

    // Evaluate the patch
    TriEval res = evaluate_triangle( gControlPoints5, vID, num_sides, pnum, u, v, w );
    pos = float4(res.P.xyz,1);
    normal = cross(res.Pu.xyz, res.Pv.xyz);  

    // Compute the final projected positions and normals
    pos = mul( pos, matWorldViewProj ); 
    normal = normalize(mul(normal, (float3x3)matWorldView)); 

    output.pos = pos;
    output.normal = normal;
    output.color = compute_color(num_sides, pos, normal);

    return output;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#if (PHONG)

// For Phong shading
float4 PSScenemain_phong( VS_OUTPUT_eval input )
{             
    return compute_lighting(input.pos, input.normal, input.color);
}

////////////////////////////////////////
#else // => (!PHONG)
////////////////////////////////////////

// For Gouraud shading
float4 PSScenemain_gouraud( VS_OUTPUT_eval input )
{
    return input.color;
}

#endif

////////////////////////////////////////

float4 PSScenemain( VS_OUTPUT_eval input ) : SV_TARGET
{
#if (PHONG)
        return PSScenemain_phong(input);
#else
        return PSScenemain_gouraud(input);
#endif
}

//////////////////////////////////////////////////////////////////////////////////
// For rendering the control mesh
//////////////////////////////////////////////////////////////////////////////////

// Vertex shader
VS_OUTPUT_mesh VSMesh( VS_INPUT_mesh input, uint vID : SV_VertexID )
{
    VS_OUTPUT_mesh output;
    output.pos    = float4( gVertexLocation.Load(int3(vID,gRowTex,0)).xyz, 1);
    output.pos    = mul( output.pos, matWorldViewProj ); 
    output.normal = float3(0, 0, 0);
//  output.normal = normalize(mul(normal, (float3x3)matWorldView)); 

    [branch]
    if (gFlags & (1 << BIT_DARK_BG))
        output.color = float4(1.0, 1.0, 1.0, 1.0);
    else
        output.color = float4(0.0, 0.0, 0.0, 1.0);

    return output;
}

// Pixel shader
float4 PSMesh( VS_OUTPUT_mesh input ) : SV_TARGET
{
    return input.color;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

technique10 NgonPass1
{
    pass p0
    {
        SetVertexShader( CompileShader( vs_4_0, VSConstruction_irreg() ) );
        SetGeometryShader( ConstructGSWithSO( CompileShader( gs_4_0,  GSConstruction_polar()), "POSITION.xyzw" ) );
        SetPixelShader( NULL);
    }  
    pass p1
    {
        SetVertexShader( CompileShader( vs_4_0, VSEvaluation_polar() ) );
        SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_4_0, PSScenemain() ));

        SetDepthStencilState( RenderWithStencilState, 0 );
    }
}

technique10 NgonPass2
{
    pass p0
    {
        SetVertexShader( CompileShader( vs_4_0, VSConstruction_reg() ) );
        SetGeometryShader( ConstructGSWithSO( CompileShader( gs_4_0,  GSConstruction_reg()), "POSITION.xyzw" ) );
        SetPixelShader( NULL);
    }  
    pass p1
    {
        SetVertexShader( CompileShader( vs_4_0, VSEvaluation_reg() ) );
        SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_4_0, PSScenemain() ));

        SetDepthStencilState( RenderWithStencilState, 0 );
    }
}

technique10 NgonPass3
{
    pass p0
    {
        SetVertexShader( CompileShader( vs_4_0, VSConstruction_irreg() ) );
        SetGeometryShader( ConstructGSWithSO( CompileShader( gs_4_0,  GSConstruction3()), "POSITION.xyzw" ) );
        SetPixelShader( NULL);
    }  
    pass p1
    {
        SetVertexShader( CompileShader( vs_4_0, VSEvaluation3() ) );
        SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_4_0, PSScenemain() ));

        SetDepthStencilState( RenderWithStencilState, 0 );
    }
}

technique10 NgonPass4
{
    pass p0
    {
        SetVertexShader( CompileShader( vs_4_0, VSConstruction_irreg() ) );
        SetGeometryShader( ConstructGSWithSO( CompileShader( gs_4_0,  GSConstruction4()), "POSITION.xyzw" ) );
        SetPixelShader( NULL);
    }  
    pass p1
    {
        SetVertexShader( CompileShader( vs_4_0, VSEvaluation4() ) );
        SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_4_0, PSScenemain() ));

        SetDepthStencilState( RenderWithStencilState, 0 );
    }
}

technique10 NgonPass5
{
    pass p0
    {
        SetVertexShader( CompileShader( vs_4_0, VSConstruction_irreg() ) );
        SetGeometryShader( ConstructGSWithSO( CompileShader( gs_4_0,  GSConstruction5()), "POSITION.xyzw" ) );
        SetPixelShader( NULL);
    }  
    pass p1
    {
        SetVertexShader( CompileShader( vs_4_0, VSEvaluation5() ) );
        SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_4_0, PSScenemain() ));

        SetDepthStencilState( RenderWithStencilState, 0 );
    }
}

// For rendering the input mesh
technique10 MeshPass
{
    pass p0
    {
        SetVertexShader( CompileShader( vs_4_0, VSMesh() ) );
        SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_4_0, PSMesh() ));

        SetDepthStencilState( RenderWithStencilState, 0 );
    }  
}

