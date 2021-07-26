
// This file is for making wrapper functions around the creation of DirectX
// resources.
//
// NOTE:
// - If any resource needs to be re-created, use the Recreate* versions of
//   the functions.

#ifndef UTILS_H
#define UTILS_H

#include <vector>
using namespace std;

template <typename ResourceType>
inline void RELEASE(ResourceType *&res) { if (res) res->Release(); }

////////////////////////////////////////////////////////////////////////////////
// Viewable resources
////////////////////////////////////////////////////////////////////////////////

// The following class encapsulates "named" resources (directly viewable in
// the shader).

template <typename ResourceType>
struct NamedResource
{
    ResourceType*                          res;
    ID3D10ShaderResourceView*             view;
    ID3D10EffectShaderResourceVariable*    var;
};

// Texture-specific
typedef NamedResource<ID3D10Texture2D> Texture2DResource;

// Buffer-specific
typedef NamedResource<ID3D10Buffer   > BufferResource;

////////////////////////////////////////////////////////////////////////////////
// Keep track of resources
////////////////////////////////////////////////////////////////////////////////

// Store all the created resources to be released later.
vector<IUnknown*> gResources;

// Add a resource for automatic cleanup.
inline void TrackResource(IUnknown* res) { gResources.push_back(res); }
template <typename T>
inline void TrackResource(T* res) { TrackResource(dynamic_cast<IUnknown*>(res)); }

// Release nr resources starting at the location of res.
void ReleaseResources(IUnknown* res, int numResources) {
    assert(numResources > 0);
    vector<IUnknown*>::iterator iter = find(gResources.begin(), gResources.end(), res);
    if (iter == gResources.end()) return; // nothing to release
    vector<IUnknown*>::iterator i = iter;
    for(; i < iter+numResources; ++i)
        (*i)->Release();
    gResources.erase(iter, iter+numResources);
}
template <typename T>
inline void ReleaseResources(T* res, int nr) { ReleaseResources(dynamic_cast<IUnknown*>(res), nr); }

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// CREATE A MUTABLE TEXTURE
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
// Call before filling the 2D texture. Pointer to texture returned in data.
HRESULT
BEGIN_CreateMutableTexture2D(ID3D10Effect* effect, const string ShaderResourceName,
        UINT width, UINT height, DXGI_FORMAT format, // D3D10_USAGE usage, // UINT CPUAccessFlags,
        Texture2DResource &tex, void **data)
{
    ID3D10Texture2D*                        &pTex = tex.res;
    ID3D10ShaderResourceView*          &pTexRView = tex.view;
    ID3D10EffectShaderResourceVariable* &pTexRVar = tex.var;

    assert(tex.res == NULL && tex.view == NULL && tex.var == NULL);

    HRESULT hr = S_OK;

    // Associate the texture with a shader variable
    pTexRVar = effect->GetVariableByName(ShaderResourceName.c_str())->AsShaderResource();

    // Set up the texture parameters
    D3D10_TEXTURE2D_DESC desc;
    desc.Width = width;
    desc.Height = height;
    desc.MipLevels = 1;
    desc.ArraySize = 1;
    desc.Format = format;
    desc.SampleDesc.Count = 1;
    desc.Usage = D3D10_USAGE_DYNAMIC;
    desc.BindFlags = D3D10_BIND_SHADER_RESOURCE;
    desc.CPUAccessFlags = D3D10_CPU_ACCESS_WRITE;
    desc.MiscFlags = 0;
    hr = DXUTGetD3D10Device()->CreateTexture2D( &desc, NULL, &pTex);

    RETURN_IF_FAILED(hr);
    TrackResource(pTex); // track for automatic release

    // Create the shader-resource view
    D3D10_SHADER_RESOURCE_VIEW_DESC srv;
    ZeroMemory( &srv, sizeof(srv) );
    srv.Format = desc.Format;
    srv.ViewDimension = D3D10_SRV_DIMENSION_TEXTURE2D;
    srv.Texture2D.MostDetailedMip = 0;
    srv.Texture2D.MipLevels = 1;
    hr = DXUTGetD3D10Device()->CreateShaderResourceView( pTex, &srv, &pTexRView );

    RETURN_IF_FAILED(hr);
    TrackResource(pTexRView); // track for automatic release
    pTexRVar->SetResource(pTexRView);

    D3D10_MAPPED_TEXTURE2D vlTex;
    hr = pTex->Map( D3D10CalcSubresource(0, 0, 1), D3D10_MAP_WRITE_DISCARD, 0, &vlTex );

    RETURN_IF_FAILED(hr);
    *data = vlTex.pData;

    return hr;
}

// Call after filling the texture.
void
END_CreateMutableTexture2D(Texture2DResource &tex)
{
    tex.res->Unmap( D3D10CalcSubresource(0, 0, 1) );
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// CREATE AN IMMUTABLE TEXTURE
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
// Call before filling the 2D texture. Pointer to texture returned in data.
HRESULT
CreateTexture2D(ID3D10Effect* effect, const string ShaderResourceName,
        UINT width, UINT height, DXGI_FORMAT format, D3D10_USAGE usage, //UINT CPUAccessFlags,
        Texture2DResource &tex, void *data)
{
    ID3D10Texture2D*                        &pTex = tex.res;
    ID3D10ShaderResourceView*          &pTexRView = tex.view;
    ID3D10EffectShaderResourceVariable* &pTexRVar = tex.var;

    assert(tex.res == NULL && tex.view == NULL && tex.var == NULL);

    HRESULT hr = S_OK;

    // Associate the texture with a shader variable
    pTexRVar = effect->GetVariableByName(ShaderResourceName.c_str())->AsShaderResource();

    // Set up the texture parameters
    D3D10_TEXTURE2D_DESC desc;
    desc.Width = width;
    desc.Height = height;
    desc.MipLevels = 1;
    desc.ArraySize = 1;
    desc.Format = format;
    desc.SampleDesc.Count = 1;
    desc.Usage = usage; // D3D10_USAGE_DEFAULT, D3D10_USAGE_IMMUTABLE, D3D10_USAGE_DYNAMIC, D3D10_USAGE_STAGING;
    desc.BindFlags = D3D10_BIND_SHADER_RESOURCE;
    desc.CPUAccessFlags = 0; //CPUAccessFlags; // 0, D3D10_CPU_ACCESS_READ, D3D10_CPU_ACCESS_WRITE;
    desc.MiscFlags = 0;
    D3D10_SUBRESOURCE_DATA texData = {0};
    texData.pSysMem = data;
    texData.SysMemPitch = width * 4;  /// <====== NOT WORKING. CHANGE DATA STRUCTURE TO VERTEX/INDEX/CONSTANT BUFFER?
    hr = DXUTGetD3D10Device()->CreateTexture2D( &desc, &texData, &pTex);

    RETURN_IF_FAILED(hr);
    TrackResource(pTex); // track for automatic release

    // Create the shader-resource view
    D3D10_SHADER_RESOURCE_VIEW_DESC srv;
    ZeroMemory( &srv, sizeof(srv) );
    srv.Format = desc.Format;
    srv.ViewDimension = D3D10_SRV_DIMENSION_TEXTURE2D;
    srv.Texture2D.MostDetailedMip = 0;
    srv.Texture2D.MipLevels = 1;
    hr = DXUTGetD3D10Device()->CreateShaderResourceView( pTex, &srv, &pTexRView );

    RETURN_IF_FAILED(hr);
    TrackResource(pTexRView); // track for automatic release
    pTexRVar->SetResource(pTexRView);

    return hr;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// CREATE A BUFFER
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

HRESULT
CreateUnnamedImmutableBuffer(UINT bindFlags, UINT size, ID3D10Buffer* &pBuf, void *data)
{
    assert(pBuf == NULL);

    // Set up the buffer parameters
    D3D10_BUFFER_DESC bufDesc = {0};
    bufDesc.BindFlags            = bindFlags; // D3D10_BIND_VERTEX_BUFFER or D3D10_BIND_INDEX_BUFFER
    bufDesc.ByteWidth            = size;      // size in bytes
    bufDesc.CPUAccessFlags       = 0;         // No CPU access required
    bufDesc.Usage                = D3D10_USAGE_IMMUTABLE; // Can't be changed after creation

    D3D10_SUBRESOURCE_DATA bufData = {0};
    bufData.pSysMem = data;

    HRESULT hr = DXUTGetD3D10Device()->CreateBuffer( &bufDesc, &bufData, &pBuf);
    RETURN_IF_FAILED(hr);
    TrackResource(pBuf); // track for automatic release
    return hr;
}

////////////////////////////////////////////////////////////////////////////////

inline HRESULT
RecreateUnnamedImmutableBuffer(UINT bindFlags, UINT size, ID3D10Buffer* &pBuf, void *data)
{
    ReleaseResources(pBuf, 1);
    pBuf = NULL;
    return CreateUnnamedImmutableBuffer(bindFlags, size, pBuf, data);
}

////////////////////////////////////////////////////////////////////////////////

// Create a named buffer with a resource view.
HRESULT
CreateNamedBufferWithView(ID3D10Effect* effect, const string ShaderResourceName,
        UINT bindFlags, UINT size, UINT cpu_access_flags, D3D10_USAGE usage,
        DXGI_FORMAT format, UINT element_width,
        BufferResource &buf,
        void *data = NULL)
{
    ID3D10Buffer*                           &pBuf = buf.res;
    ID3D10ShaderResourceView*          &pBufRView = buf.view;
    ID3D10EffectShaderResourceVariable* &pBufRVar = buf.var;

    assert(buf.res == NULL && buf.view == NULL && buf.var == NULL);

    // Associate the texture with a shader variable
    pBufRVar = effect->GetVariableByName(ShaderResourceName.c_str())->AsShaderResource();

    // Set up the buffer parameters
    D3D10_BUFFER_DESC bufDesc = {0};
    bufDesc.BindFlags            = bindFlags; // D3D10_BIND_VERTEX_BUFFER or D3D10_BIND_INDEX_BUFFER
    bufDesc.ByteWidth            = size;      // size in bytes
    bufDesc.CPUAccessFlags       = cpu_access_flags; // D3D10_CPU_ACCESS_WRITE | D3D10_CPU_ACCESS_READ
    bufDesc.Usage                = usage; // D3D10_USAGE_DEFAULT, D3D10_USAGE_IMMUTABLE, D3D10_USAGE_DYNAMIC, D3D10_USAGE_STAGING.

    D3D10_SUBRESOURCE_DATA bufData = {0};
    bufData.pSysMem = data;

    HRESULT hr = DXUTGetD3D10Device()->CreateBuffer( &bufDesc, data ? &bufData : NULL, &pBuf);

    RETURN_IF_FAILED(hr);
    TrackResource(pBuf); // track for automatic release

    // Create the shader-resource view
    D3D10_SHADER_RESOURCE_VIEW_DESC SRVDesc;
    ZeroMemory( &SRVDesc, sizeof(SRVDesc) );
    SRVDesc.Format = format;
    SRVDesc.ViewDimension = D3D10_SRV_DIMENSION_BUFFER;
    SRVDesc.Buffer.ElementOffset = 0;
    SRVDesc.Buffer.ElementWidth = element_width;
    hr = DXUTGetD3D10Device()->CreateShaderResourceView( pBuf, &SRVDesc, &pBufRView);

    RETURN_IF_FAILED(hr);
    TrackResource(pBufRView); // track for automatic release
    pBufRVar->SetResource(pBufRView);

    return hr;
}

////////////////////////////////////////////////////////////////////////////////

void
ReleaseTrackedResources()
{
    for (unsigned i = 0; i < gResources.size(); ++i)
        RELEASE( gResources[i] );
}

////////////////////////////////////////////////////////////////////////////////
// This class is for handling the construction/evaluation passes.
////////////////////////////////////////////////////////////////////////////////

struct AppState;

// This stores information for two-pass construction and evaluation.
struct PassInfo {
    ID3D10EffectTechnique *technique;
    BufferResource *streamout;
    int num_sides;

    static vector<ID3D10ShaderResourceView*> resource_views;

    ////////////////////////////////////////////////////////////
    // Constructor/destructor
    PassInfo() { technique = NULL; streamout = NULL; num_sides = 0; }

    ////////////////////////////////////////////////////////////
    // Initialize the pass
    HRESULT
    initialize(struct AppState &state, const string TechniqueName, BufferResource *Streamout, int NumSides);

    ////////////////////////////////////////////////////////////
    // Set up the passes
    // pass: 0 => construction, 1 => evaluation
    static void setup(UINT pass, ID3D10Device *pd3dDevice, struct AppState &state);
    static void clean_up(UINT pass, ID3D10Device *pd3dDevice, struct AppState &state);

    ////////////////////////////////////////////////////////////
    // Execute the pass
    // pass: 0 => construction, 1 => evaluation
    HRESULT
    execute(UINT pass, ID3D10Device *pd3dDevice, struct AppState &state,
            UINT num_patches, UINT num_instances, UINT start_facet = 0);
};

////////////////////////////////////////////////////////////////////////////////

#endif // UTILS_H

