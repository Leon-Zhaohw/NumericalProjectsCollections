
#ifndef EVENT_HANDLERS_H
#define EVENT_HANDLERS_H

////////////////////////////////////////////////////////////////////////////////

#include "DXUT.h"

//--------------------------------------------------------------------------------------
// Handle updates to the scene.  This is called regardless of which D3D API is used
//--------------------------------------------------------------------------------------
void CALLBACK OnFrameMove( double fTime, float fElapsedTime, void* pUserContext )
{
    // Update the camera's position based on user input \    
    //cout<<"moving frame..."<<fElapsedTime<<endl;
    g_Camera.FrameMove( fElapsedTime );
}

void CALLBACK OnKeyboardEvent(
        UINT nChar,
        bool bKeyDown,
        bool bAltDown,
        void* pUserContext
        )
{
    AppState* pState = (AppState*)pUserContext;

    if( bKeyDown )
    {
        switch( nChar ) {
            case '9':
                pState->TF = 5;
                bTessChange=true;
                break;
                
            case '0':
                pState->TF = 1;
                bTessChange=true;
                break;
            case '1':
                pState->TF = 2;
                bTessChange=true;
                break;
            case '2':
                pState->TF = 4;
                bTessChange=true;
                break;
            case '3':
                pState->TF = 8;
                bTessChange=true;
                break;
            case '4':
                pState->TF = 16;
                bTessChange=true;
                break;
            case '5':
                pState->TF = 32;
                bTessChange=true;
                break;
            case '6':
                pState->TF = 64;
                bTessChange=true;
                break;
            case '7':
                pState->TF = 128;
                bTessChange=true;
                break;
            case 'A':
                if (bAnimation) bAnimation=false ;
                else bAnimation=true;
                break;
            case 'G':
                toggle_bit(pState->flags, BIT_GROUP_COLOR);
                g_pFlags->SetInt( pState->flags );
                cout << "Group color is "
                    << ( get_bit(pState->flags, BIT_GROUP_COLOR) ? "ON" : "OFF" )
                    << endl;
                break;
            case 'B':
                bDarkBackground = !bDarkBackground;
                toggle_bit(pState->flags, BIT_DARK_BG);
                g_pFlags->SetInt( pState->flags );
                cout << "Background is "
                    << ( get_bit(pState->flags, BIT_DARK_BG) ? "DARK" : "LIGHT" )
                    << endl;
                break;
                break;
            case 'T':
                bText = !bText;
                break;
            case 'F':
                DXUTToggleFullScreen();
                bFullScreen = !bFullScreen;
                break;
            case 'C':
                bShowMesh = !bShowMesh;
                break;
            case 'R':
				pState->RowTex = 0;
                break;
        }
    }
}

//--------------------------------------------------------------------------------------
// Render the help and statistics text
//--------------------------------------------------------------------------------------
void RenderText()
{
    g_pTxtHelper->Begin();
    g_pTxtHelper->SetInsertionPos( 2, 0 );
    g_pTxtHelper->SetForegroundColor( D3DXCOLOR( 0.0f, 0.5f, 0.0f, 1.0f ) );
    g_pTxtHelper->DrawTextLine( DXUTGetFrameStats( DXUTIsVsyncEnabled() ) );  
    g_pTxtHelper->DrawTextLine( DXUTGetDeviceStats() );
    //g_pTxtHelper->SetForegroundColor( D3DXCOLOR( 0.0f, 0.0f, 0.0f, 1.0f ) );
    //g_pTxtHelper->DrawFormattedTextLine( L"Frame per second: %f", fps );
    g_pTxtHelper->End();
}
//--------------------------------------------------------------------------------------
// Handle messages to the application
//--------------------------------------------------------------------------------------
LRESULT CALLBACK MsgProc( HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam, bool* pbNoFurtherProcessing, void* pUserContext )
{
    //cout<<"handling msgs..."<<endl;
    // Pass all remaining windows messages to camera so it can respond to user input
    g_Camera.HandleMessages( hWnd, uMsg, wParam, lParam );

    return 0;
}
//--------------------------------------------------------------------------------------
// Create any D3D10 resources that depend on the back buffer
//--------------------------------------------------------------------------------------
HRESULT CALLBACK OnD3D10SwapChainResized( ID3D10Device* pd3dDevice, IDXGISwapChain *pSwapChain, const DXGI_SURFACE_DESC* pBackBufferSurfaceDesc, void* pUserContext )
{
    HRESULT hr = S_OK;

    // Setup the camera's projection parameters
    float fAspectRatio = pBackBufferSurfaceDesc->Width / (FLOAT)pBackBufferSurfaceDesc->Height;
    //cout<<"ratio: "<<fAspectRatio<<endl;
    g_Camera.SetProjParams( D3DX_PI/4, fAspectRatio, 0.05f, 800.0f );

    g_Camera.SetWindow( pBackBufferSurfaceDesc->Width, pBackBufferSurfaceDesc->Height );
    g_Camera.SetButtonMasks( MOUSE_LEFT_BUTTON, MOUSE_WHEEL, MOUSE_MIDDLE_BUTTON );

    return hr;
}
//--------------------------------------------------------------------------------------
// Release D3D10 resources created in OnD3D10ResizedSwapChain 
//--------------------------------------------------------------------------------------
void CALLBACK OnD3D10ReleasingSwapChain( void* pUserContext )
{
}
//--------------------------------------------------------------------------------------
// Called right before creating a D3D9 or D3D10 device, allowing the app to modify the device settings as needed
//--------------------------------------------------------------------------------------
bool CALLBACK ModifyDeviceSettings( DXUTDeviceSettings* pDeviceSettings, void* pUserContext )
{
    // For the first device created if its a REF device, optionally display a warning dialog box
    static bool s_bFirstTime = true;
    if( s_bFirstTime )
    {
        s_bFirstTime = false;
        if( (DXUT_D3D9_DEVICE == pDeviceSettings->ver && pDeviceSettings->d3d9.DeviceType == D3DDEVTYPE_REF) ||
                (DXUT_D3D10_DEVICE == pDeviceSettings->ver && pDeviceSettings->d3d10.DriverType == D3D10_DRIVER_TYPE_REFERENCE) )
            DXUTDisplaySwitchingToREFWarning( pDeviceSettings->ver );
    }

    return true;
}

#endif // EVENT_HANDLERS_H

