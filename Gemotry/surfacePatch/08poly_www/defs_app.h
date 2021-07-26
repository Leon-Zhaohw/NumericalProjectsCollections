
#ifndef DEFS_APP_H
#define DEFS_APP_H

////////////////////////////////////////////////////////////

ID3D10EffectMatrixVariable*   g_pMatWorldView;
ID3D10EffectMatrixVariable*   g_pMatWorldViewProj;
ID3D10EffectVectorVariable*   g_pLightVecW1;
ID3D10EffectVectorVariable*   g_pEyePosW;
ID3D10EffectScalarVariable*   g_pFlags;
ID3D10EffectScalarVariable*   g_pRowTex;
ID3D10EffectScalarVariable*   g_pStart;

CModelViewerCamera      g_Camera;               // A model viewing camera

ID3DX10Font*            g_pFont10 = NULL;       
ID3DX10Sprite*          g_pSprite10 = NULL;
CDXUTTextHelper*        g_pTxtHelper = NULL;

int cVertIndex=0;
int cIndexIndex=0;

// number of regular patches
int g_numFaces;
int g_start[MAX_FACE_SIZE+1] = {0};
int g_num[MAX_FACE_SIZE+1] = {0};

ostream& operator << ( ostream& out, const D3DXVECTOR3& v )
{
    out<<"{ "
        << v.x << ", "
        << v.y << ", "
        << v.z << " }";

    return out;
}

REAL xmax=-10000.0;
REAL xmin=10000.0;
REAL ymax=-10000.0;
REAL ymin=10000.0;
REAL zmax=-10000.0;
REAL zmin=10000.0;


/////////////////////////////////
// 1 for x, 2 for y 3 for z axis
/////////////////////////////////
unsigned int longAxis=0;
REAL maxAmount=0;


//
// Vertex format to use in this example
//
#pragma pack( push, 1 )
struct ConstructionVertexFormat
{
    UINT32 valence;
};
struct EvaluationVertexFormat
{
    UINT32 patchID;
    D3DXVECTOR3 position;
};
//struct ControlMeshFormat
//{
//    unsigned int patchID;
//};
#pragma pack (pop)


bool bDarkBackground = false; // Background color choice
bool bTessChange = false; // Change Tessellation Factor
bool bMorphing   = false; // Morphing -> rotation around axis
bool bFrog       = false; // Load frog animation
bool bAnimation  = false; // Animation
bool bText       = true ; // Text up top including the frame rate
bool bFullScreen = false; // Full screen
bool bShowMesh   = false; // Display control mesh

int g_numFrames = 1; // The number of frames in the pre-scripted animation (set in initialize())

//short: 16 bits, int: 32 bits

vector<unsigned short int> vIndices;
vector<ConstructionVertexFormat> vVertices;
vector<unsigned short int> vMeshIndices;

Polyhedron P; 
Polyhedron pFrog[100];


//
// D3D10 input element descriptions corresponding to our VertexFormat
//
D3D10_INPUT_ELEMENT_DESC ConstructionVertexFormatElements[] = 
{
    { "BLENDINDICES", 0, DXGI_FORMAT_R32_UINT, 0, 0, D3D10_INPUT_PER_VERTEX_DATA, 0 },
};

D3D10_INPUT_ELEMENT_DESC EvaluationVertexFormatElements[] = 
{
    { "BLENDINDICES", 0, DXGI_FORMAT_R32_UINT, 0, 0, D3D10_INPUT_PER_VERTEX_DATA, 0 },
    { "SV_POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 4, D3D10_INPUT_PER_VERTEX_DATA, 0 },
};

//
// Structure to keep track of application state
//
struct AppState
{
    ID3D10Effect*           pEffect10;
    
    // Shader passes
    PassInfo NgonPass[MAX_FACE_SIZE+1];

    // Mesh rendering pass
    ID3D10EffectTechnique *pMeshTechnique;

    // Streamout related stuff
    BufferResource ControlPoints[MAX_FACE_SIZE+1];

    // Layouts
    ID3D10InputLayout*  pConstructionLayout;
    ID3D10InputLayout*  pEvaluationLayout;

    // Vertex and index buffers
    ID3D10Buffer*  pVertexBuffer;
    ID3D10Buffer*  pIndexBuffer;
    ID3D10Buffer*  pEvaluationVertexBuffer;
    ID3D10Buffer*  pEvaluationIndexBuffer;
    ID3D10Buffer*  pMeshIndexBuffer;

    // vertex location Texture
    Texture2DResource VertexLocation;

    // ring index Texture
    Texture2DResource RingIndex;

    // offset data
    Texture2DResource OffsetData; 

    unsigned int  nValidEvaluationIndices;
    unsigned int  nValidFaces;

    D3DXVECTOR3   cameraEye;
    D3DXVECTOR3   cameraAt;
    D3DXVECTOR3   cameraUp;

    unsigned int  TF;
    unsigned int  RowTex;
    UINT32        flags;
};

//--------------------------------------------------------------------------------------
// Forward declarations 
//--------------------------------------------------------------------------------------
bool    CALLBACK ModifyDeviceSettings( DXUTDeviceSettings* pDeviceSettings, void* pUserContext );
void    CALLBACK OnFrameMove( double fTime, float fElapsedTime, void* pUserContext );
LRESULT CALLBACK MsgProc( HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam, bool* pbNoFurtherProcessing, void* pUserContext );
HRESULT CALLBACK OnD3D10SwapChainResized( ID3D10Device* pd3dDevice, IDXGISwapChain *pSwapChain, const DXGI_SURFACE_DESC* pBackBufferSurfaceDesc, void* pUserContext );
void    CALLBACK OnD3D10FrameRender( ID3D10Device* pd3dDevice, double fTime, float fElapsedTime, void* pUserContext );
void    CALLBACK OnD3D10ReleasingSwapChain( void* pUserContext );
void    CALLBACK OnD3D10DestroyDevice( void* pUserContext );

HRESULT tessellation(AppState &state);
void CALLBACK RenderFrame(
        ID3D10Device * pd3dDevice,    // Direct3D device pointer
        DOUBLE fTime,                 // Time since the application started
        FLOAT fElapsedTime,           // Time since the last call to RenderFrame()
        void* pUserContext );

#endif // DEFS_APP_H

