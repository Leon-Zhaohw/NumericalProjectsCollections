
#ifndef RENDER_H
#define RENDER_H

/////////////////////////////////////////////////////////////////////////////
// Morphing
/////////////////////////////////////////////////////////////////////////////

int iMorph=0;                   // Rotation amount
bool bMorphIncrease=true;       // iMorph increase or decrease

void morphing(AppState &state)
{
    
    D3D10_MAPPED_TEXTURE2D vlTex;
    state.VertexLocation.res->Map( D3D10CalcSubresource(0, 0, 1), D3D10_MAP_WRITE_DISCARD, 0, &vlTex );

    float* vLocation = (float*) vlTex.pData;

    // Modify 100 models for Frog animation
#if 0
    if(bFrog)
        for(int k=0; k<100; k++){

            Vertex_iterator vv=pFrog[k].vertices_begin();

            for(int i=0; i<pFrog[k].size_of_vertices(); ++i){

                float angle=PI*iMorph*0.002/maxAmount;

                // rotation around X-axis
                if(longAxis==1)
                {
                    vLocation[k*MAX_ALL + i*4+1]= cos(angle*vv->point().x())*(vv->point().y())-sin(angle*vv->point().x())*(vv->point().z());
                    vLocation[k*MAX_ALL + i*4+2]= sin(angle*vv->point().x())*(vv->point().y())+cos(angle*vv->point().x())*(vv->point().z());
                }        
                
                // rotation around Y-axis
                else if(longAxis==2)
                {
                    vLocation[k*MAX_ALL + i*4  ]= cos(angle*vv->point().y())*(vv->point().x())+sin(angle*vv->point().y())*(vv->point().z());
                    vLocation[k*MAX_ALL + i*4+2]=-sin(angle*vv->point().y())*(vv->point().x())+cos(angle*vv->point().y())*(vv->point().z());
                }
                
                // rotation around Z-axis
                else
                {
                    vLocation[k*MAX_ALL + i*4  ]= cos(angle*vv->point().z())*(vv->point().x())-sin(angle*vv->point().z())*(vv->point().y());
                    vLocation[k*MAX_ALL + i*4+1]= sin(angle*vv->point().z())*(vv->point().x())+cos(angle*vv->point().z())*(vv->point().y());
                }
                
                vv++;
            }
        }
    else
    {
#endif
        Vertex_iterator vv=P.vertices_begin();

        for(int i = 0; i < (int)P.size_of_vertices(); ++i){

            float angle = float(PI*iMorph*0.002/maxAmount);

            // rotation around X-axis
            if(longAxis==1)
            {
                vLocation[i*4+0]= vv->point().x();
                vLocation[i*4+1]= cos(angle*vv->point().x())*(vv->point().y())-sin(angle*vv->point().x())*(vv->point().z());
                vLocation[i*4+2]= sin(angle*vv->point().x())*(vv->point().y())+cos(angle*vv->point().x())*(vv->point().z());
            }        
            
            // rotation around Y-axis
            else if(longAxis==2)
            {
                vLocation[i*4  ]= cos(angle*vv->point().y())*(vv->point().x())+sin(angle*vv->point().y())*(vv->point().z());
                vLocation[i*4+1]= vv->point().y();
                vLocation[i*4+2]=-sin(angle*vv->point().y())*(vv->point().x())+cos(angle*vv->point().y())*(vv->point().z());
            }
            
            // rotation around Z-axis
            else
            {
                vLocation[i*4  ]= cos(angle*vv->point().z())*(vv->point().x())-sin(angle*vv->point().z())*(vv->point().y());
                vLocation[i*4+1]= sin(angle*vv->point().z())*(vv->point().x())+cos(angle*vv->point().z())*(vv->point().y());
                vLocation[i*4+2]= vv->point().z();
            }
            
            vv++;
        }
//  }

    state.VertexLocation.res->Unmap( D3D10CalcSubresource(0, 0, 1) );
    
    if(bMorphIncrease)iMorph+=2;
    else iMorph-=2;

    if(iMorph==100)
        bMorphIncrease=false;
    else if(iMorph==-100)
        bMorphIncrease=true;
}

///////////////////////////////////////////////////////////////////////////////
// RenderFrame
//
// This function is called by DXUT each time the frame needs rendering.
// Generally DXUT tries to maintain a constant frame rate by calling this
// function at regular intervals.
//
void CALLBACK RenderFrame(
        ID3D10Device * pd3dDevice,    // Direct3D device pointer
        DOUBLE fTime,                 // Time since the application started
        FLOAT fElapsedTime,           // Time since the last call to RenderFrame()
        void* pUserContext )          // Unused in this application
{
    HRESULT hr = S_OK;
    AppState &state = *(AppState*)pUserContext;

    // Clear the render target to a predefined background color.
    //
    const FLOAT WHITE[4] = { 1.0f, 1.0f, 1.0f, 0.0f };
    const FLOAT BLACK[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
    const FLOAT *rgbaBackground =  bDarkBackground ? BLACK : WHITE;
    pd3dDevice->ClearRenderTargetView( DXUTGetD3D10RenderTargetView(), rgbaBackground );
    pd3dDevice->ClearDepthStencilView( DXUTGetD3D10DepthStencilView(), D3D10_CLEAR_DEPTH|D3D10_CLEAR_STENCIL, 1.0f, 0 );
    // clear the depth stencil
    //ID3D10DepthStencilView* pDSV = DXUTGetD3D10DepthStencilView();
    //pd3dDevice->ClearDepthStencilView( pDSV, D3D10_CLEAR_DEPTH, 1.0, 0 ); 
    //ID3D10RenderTargetView* pRTV = NULL;


    // Morphing
    if(bMorphing) morphing(state);

    // Change the tessellation buffer if the tessellation factor has changed
    if(bTessChange) hr = tessellation(state);

    // Set shader constants for lighting computation
    D3DXMATRIX matWorld         = *g_Camera.GetWorldMatrix();
    D3DXMATRIX matProj          = *g_Camera.GetProjMatrix();
    D3DXMATRIX matView          = *g_Camera.GetViewMatrix();
    D3DXMATRIX matWorldView     = matWorld*matView;
    D3DXMATRIX matWorldViewProj = matWorldView*matProj;
    g_pMatWorldView->SetMatrix( (float*)&matWorldView );
    g_pMatWorldViewProj->SetMatrix( (float*)&matWorldViewProj );
    D3DXVECTOR3 viewLightDir1=D3DXVECTOR3(0,0,-1);
    g_pLightVecW1->SetFloatVector( (float*)&viewLightDir1 );
    D3DXVECTOR3 vecAt(0.0f,0.0f,10.0f);
    D3DXVec3TransformNormal( &vecAt, &vecAt, &matView ); g_pEyePosW->SetFloatVector((float*)vecAt);

    // Set shader constants for Frog Animation 
    if (bAnimation)
		state.RowTex = (state.RowTex+1) % g_numFrames;
    g_pRowTex->SetInt( state.RowTex );

    // Construct and evaluate the patches
    for (UINT p = 0; p < NUM_PASSES; ++p) {
        // Set up the current pass
        PassInfo::setup(p, pd3dDevice, state);

        // Execute the passes for each case
        for (int i = 0; i <= MAX_TYPE; ++i) {
//if (i != REG) continue;
            if (g_num[i] == 0) continue;
            state.NgonPass[i].execute(p, pd3dDevice, state, g_num[i], NUM_OF_INSTANCES, g_start[i]);
        }

        // Clean up
        PassInfo::clean_up(p, pd3dDevice, state);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Show the input mesh
    ////////////////////////////////////////////////////////////////////////////////

    if (bShowMesh) {
        ID3D10Buffer* targets[1] = { state.pVertexBuffer };
        UINT stride[1] = { sizeof(ConstructionVertexFormat) };
        UINT offset[1] = { 0 };

        // Set up input buffers
        pd3dDevice->IASetInputLayout( state.pConstructionLayout );
        pd3dDevice->IASetVertexBuffers( 0, 1, targets, stride, offset );
        pd3dDevice->IASetIndexBuffer( state.pMeshIndexBuffer, DXGI_FORMAT_R16_UINT, 0 );
        pd3dDevice->IASetPrimitiveTopology( D3D10_PRIMITIVE_TOPOLOGY_LINELIST ); 

        D3D10_TECHNIQUE_DESC techDesc;
        state.pMeshTechnique->GetDesc( &techDesc );
        assert( techDesc.Passes == 1 ); // there should be only 1 pass
        state.pMeshTechnique->GetPassByIndex( 0 )->Apply(0);

        pd3dDevice->DrawIndexed((int)vMeshIndices.size(), 0, 0);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MISC
    ////////////////////////////////////////////////////////////////////////////////

    if (bText) RenderText();
}

////////////////////////////////////////////////////////////////////////////////
// PassInfo class to do the two-pass approach (declared in util.h)
////////////////////////////////////////////////////////////////////////////////

// Initialize static vector
vector<ID3D10ShaderResourceView*> PassInfo::resource_views;

////////////////////////////////////////////////////////////
// Initialize the pass
HRESULT
PassInfo::initialize(struct AppState &state, const string TechniqueName, BufferResource *Streamout, int NumSides)
{
    num_sides = NumSides;
    streamout = Streamout;
    resource_views.push_back(streamout->view);

    // Initialize the technique
    technique = state.pEffect10->GetTechniqueByName( TechniqueName.c_str() );

    return S_OK;
}

////////////////////////////////////////////////////////////
// Set up the passes
// pass: 0 => construction, 1 => evaluation
void
PassInfo::setup(UINT pass, ID3D10Device *pd3dDevice, struct AppState &state)
{
    assert(0 <= pass && pass < NUM_PASSES);

    static ID3D10DepthStencilState* pDSOld = NULL;
    static ID3D10DepthStencilState* pDSState = NULL;
    static UINT OldStencil = 0;

    ID3D10Buffer* targets[1];
    UINT stride[1];
    UINT offset[1] = { 0 };

    // For the construction pass
    if (pass == 0) {
        // Set up input buffers
        targets[0] = state.pVertexBuffer;
        stride[0] = sizeof(ConstructionVertexFormat);
        pd3dDevice->IASetInputLayout( state.pConstructionLayout );
        pd3dDevice->IASetVertexBuffers( 0, 1, targets, stride, offset );
        pd3dDevice->IASetIndexBuffer( state.pIndexBuffer, DXGI_FORMAT_R16_UINT, 0 );
        pd3dDevice->IASetPrimitiveTopology( D3D10_PRIMITIVE_TOPOLOGY_TRIANGLELIST_ADJ );  

        // Disable rasterization for this pass
        pd3dDevice->OMGetDepthStencilState( &pDSOld, &OldStencil );
        D3D10_DEPTH_STENCIL_DESC dsdesc = {0};
        pd3dDevice->CreateDepthStencilState( &dsdesc, &pDSState );
        pd3dDevice->OMSetDepthStencilState( pDSState, 0 );
    }

    // For the evaluation pass
    else if (pass == 1) {
        // Get back to normal
        targets[0] = NULL;
        pd3dDevice->SOSetTargets( 1, targets, offset ); 
        // re-enable rasterization
        pd3dDevice->OMSetDepthStencilState( pDSOld, OldStencil );
        // Cleanup stencil state resources
        RELEASE(pDSOld);
        RELEASE(pDSState);
        OldStencil = 0;

        stride[0] = sizeof(EvaluationVertexFormat) ;
        pd3dDevice->IASetInputLayout( state.pEvaluationLayout );
        pd3dDevice->IASetVertexBuffers( 0, 1, &state.pEvaluationVertexBuffer, stride, offset );
        pd3dDevice->IASetIndexBuffer( state.pEvaluationIndexBuffer, DXGI_FORMAT_R16_UINT, 0 );
        pd3dDevice->IASetPrimitiveTopology( D3D10_PRIMITIVE_TOPOLOGY_TRIANGLELIST ); 

        pd3dDevice->VSSetShaderResources( 0, (UINT)PassInfo::resource_views.size(),
                                                  &PassInfo::resource_views[0] );
//      pd3dDevice->VSSetConstantBuffers( 0, (UINT)PassInfo::resource_views.size(),
//                                                &PassInfo::resource_views[0] );
    }
}

////////////////////////////////////////////////////////////
// Clean up the passes
// pass: 0 => construction, 1 => evaluation
void
PassInfo::clean_up(UINT pass, ID3D10Device *pd3dDevice, struct AppState &state)
{
    // nothing to be done
}

////////////////////////////////////////////////////////////
// Execute the passes
// pass: 0 => construction, 1 => evaluation
HRESULT
PassInfo::execute(UINT pass, ID3D10Device *pd3dDevice, struct AppState &state,
        UINT num_patches, UINT num_instances, UINT start_facet)
{
    assert(0 <= pass && pass < NUM_PASSES);

    // Nothing to render if no indices.
    if (num_patches == 0) return S_OK;

    HRESULT hr = S_OK;

    D3D10_TECHNIQUE_DESC techDesc;
    technique->GetDesc( &techDesc );
    assert( techDesc.Passes == NUM_PASSES ); // there should be two passes

    // Construction
    if (pass == 0) {
        ID3D10Buffer* targets[1] = { streamout->res };
        UINT offset[1] = { 0 };
        pd3dDevice->SOSetTargets( 1, targets, offset );
        technique->GetPassByIndex( 0 )->Apply(0);
        pd3dDevice->DrawIndexedInstanced( num_patches*6, num_instances, start_facet*6, 0, 0 );
    }

    // Evaluation
    else {
        technique->GetPassByIndex( 1 )->Apply(0);
        if (num_sides == 3) {
            pd3dDevice->DrawIndexedInstanced( state.nValidEvaluationIndices/2, num_patches*num_instances, 0, 0, 0 );
        }
        else if (num_sides == 4) {
            pd3dDevice->DrawIndexedInstanced( state.nValidEvaluationIndices  , num_patches*num_instances, 0, 0, 0 );
        }
        else {
            pd3dDevice->DrawIndexedInstanced( state.nValidEvaluationIndices/2, num_patches*num_instances*num_sides, 0, 0, 0 );
        }
    }

    return hr;
}

////////////////////////////////////////////////////////////////////////////////

#endif // RENDER_H

