
#ifndef INITIALIZE_H
#define INITIALIZE_H

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// Animation
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void LoadAnimationData(float *vLocation, const string &pref)
{
    // Load files for pre-scripted frog animation
    if(bFrog)
    {
        OUT1_BEGIN("Loading frog files for animation...");
        for(int i=0; i<FROG_FRAMES; ++i)
        {
            string prefix = pref + "/monsterfrog_mapped_";
            string iFileName = str_index(prefix.c_str(), i, ".off", 2);

            ifstream iImage(iFileName.c_str(), ios::in);
            if (!iImage) FATAL_ERROR("Unable to open '" << iFileName << "'\n");
            iImage >> pFrog[i];
            if (!iImage) FATAL_ERROR("Error reading model from '" << iFileName << "'\n");
            
            Vertex_iterator viter=pFrog[i].vertices_begin();
            for(int j=0; j<(int)pFrog[i].size_of_vertices(); ++j)
            {
                vLocation[i*MAX_ALL*4 + j*4  ]=viter->point().x();
                vLocation[i*MAX_ALL*4 + j*4+1]=viter->point().y();
                vLocation[i*MAX_ALL*4 + j*4+2]=viter->point().z();
                vLocation[i*MAX_ALL*4 + j*4+3]=(float)viter->vertex_degree();
                viter++;
            }
        }
        OUT1_DONE();
    }

    // Set up twisting animation for other loaded models
    else {
        for (int iMorph = 0; iMorph < g_numFrames; ++iMorph) {
            Vertex_iterator vv = P.vertices_begin();

            // Do one back and forth swing for a full cycle
            float angFact = (iMorph <   MORPH_FRAMES/4) ? (float)iMorph
                          : (iMorph < 3*MORPH_FRAMES/4) ? (MORPH_FRAMES/2) - (float)iMorph 
                          :                               (float)(iMorph - MORPH_FRAMES);

            for(int i = 0; i < (int)P.size_of_vertices(); ++i){

                float angle = float(PI*angFact*0.005/maxAmount);

                // rotation around X-axis
                if(longAxis==1)
                {
                    vLocation[iMorph*MAX_ALL*4 + i*4 + 0]= vv->point().x();
                    vLocation[iMorph*MAX_ALL*4 + i*4 + 1]= cos(angle*vv->point().x())*(vv->point().y())-sin(angle*vv->point().x())*(vv->point().z());
                    vLocation[iMorph*MAX_ALL*4 + i*4 + 2]= sin(angle*vv->point().x())*(vv->point().y())+cos(angle*vv->point().x())*(vv->point().z());
                }        
                
                // rotation around Y-axis
                else if(longAxis==2)
                {
                    vLocation[iMorph*MAX_ALL*4 + i*4 + 0]= cos(angle*vv->point().y())*(vv->point().x())+sin(angle*vv->point().y())*(vv->point().z());
                    vLocation[iMorph*MAX_ALL*4 + i*4 + 1]= vv->point().y();
                    vLocation[iMorph*MAX_ALL*4 + i*4 + 2]=-sin(angle*vv->point().y())*(vv->point().x())+cos(angle*vv->point().y())*(vv->point().z());
                }
                
                // rotation around Z-axis
                else
                {
                    vLocation[iMorph*MAX_ALL*4 + i*4 + 0]= cos(angle*vv->point().z())*(vv->point().x())-sin(angle*vv->point().z())*(vv->point().y());
                    vLocation[iMorph*MAX_ALL*4 + i*4 + 1]= sin(angle*vv->point().z())*(vv->point().x())+cos(angle*vv->point().z())*(vv->point().y());
                    vLocation[iMorph*MAX_ALL*4 + i*4 + 2]= vv->point().z();
                }
                
                vv++;
            }

#if 0
            // Output the intermediate vertices so that they may be used to create a separate .bv file and viewed offline.
            if (iMorph == MORPH_FRAMES/4) {
                const char filename[] = "morphed_verts.out";
                OUT1_BEGIN("Saving morphed vertices for frame " << iMorph << " to '" << filename << "'...");
                ofstream fout(filename);
                fout << P.size_of_vertices() << endl;
                for (int i = 0; i < (int)P.size_of_vertices(); ++i) {
                    fout << vLocation[iMorph*MAX_ALL*4 + i*4 + 0] << " "
                        << vLocation[iMorph*MAX_ALL*4 + i*4 + 1] << " "
                        << vLocation[iMorph*MAX_ALL*4 + i*4 + 2] << endl;
                }
                fout.close();
                OUT1_DONE();
            }
#endif
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Tessellation 
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
HRESULT
tessellation(AppState &state)
{
    OUT1_BEGIN("Tessellation table...");

    ////////////////////////////////////////////////////////////
    // Create the new tessellation
    ////////////////////////////////////////////////////////////

    const int N = state.TF;
 
    vector<EvaluationVertexFormat> quadVertices((N+1)*(N+1));
    vector<unsigned short int> quadIndices;
    quadIndices.reserve(N*N*6);

    const float scale = 1.0f / N;
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            quadVertices[(N+1)*i+j].position = D3DXVECTOR3(i*scale, j*scale, (N-i-j)*scale);
            quadVertices[(N+1)*i+j].patchID = 0;
        }
    }

    // Half of the quad, which is a triangle
    //     --> i
    //   | . . . .
    // j | . . .
    //   V . .
    //     .
    // Go along diagonal triangle strips (actually enumerating triangles, though).
    // First two strips enumerated as i,j:
    //     0,1  0,0  1,0
    //     0,2  0,1  1,1  1,0  2,0
    // k goes top-left to bottom-right, and l goes bottom-left to top-right.
    //     i,j = l,k-l alternating with l,k-l-1 and ending with k,0
    const int factor = N + 1;
    for (int k = 1; k <= N; ++k) {
        for (int l = 0; l < k; ++l) {
            quadIndices.push_back((unsigned short)( l   *factor + (k-l  )));
            quadIndices.push_back((unsigned short)( l   *factor + (k-l-1)));
            quadIndices.push_back((unsigned short)((l+1)*factor + (k-l-1)));

            if (l < k-1) {
                quadIndices.push_back((unsigned short)((l+1)*factor + (k-l-1)));
                quadIndices.push_back((unsigned short)( l   *factor + (k-l-1)));
                quadIndices.push_back((unsigned short)((l+1)*factor + (k-l-2)));
            }
        }
    }

    // The other half of the quad
    //     --> i
    //   |       .
    // j |     . .
    //   V   . . .
    //     . . . .
    // Go along diagonal triangle strips (actually enumerating triangles, though).
    // First two and last strips enumerated as i,j:
    //     N,N-1  N,N  N-1,N
    //     N,N-2  N,N-1  N-1,N-1  N-1,N  N-2,N
    //     ...
    //     N,0  N,1  N-1,1  N-1,2 ... 1,N  0,N
    // k goes bottom-right to top-left, and l goes top-right to bottom-left.
    //     i,j = l,N+k-l alternating with l,N+k-l+1 and ending with k,N
    for (int k = N-1; k >= 0; --k) {
        for (int l = N; l > k; --l) {
            quadIndices.push_back((unsigned short)( l   *factor + (N+k-l  )));
            quadIndices.push_back((unsigned short)( l   *factor + (N+k-l+1)));
            quadIndices.push_back((unsigned short)((l-1)*factor + (N+k-l+1)));

            if (l > k+1) {
                quadIndices.push_back((unsigned short)((l-1)*factor + (N+k-l+1)));
                quadIndices.push_back((unsigned short)( l   *factor + (N+k-l+1)));
                quadIndices.push_back((unsigned short)((l-1)*factor + (N+k-l+2)));
            }
        }
    }

    ////////////////////////////////////////////////////////////
    // Store the quad tessellation back into the vertex buffer.
    ////////////////////////////////////////////////////////////
 
    HRESULT hr = S_OK;

    // Create the vertex buffer (and release the previous one)
    hr = RecreateUnnamedImmutableBuffer(D3D10_BIND_VERTEX_BUFFER,
            (UINT)quadVertices.size()*sizeof(EvaluationVertexFormat),
            state.pEvaluationVertexBuffer, (void *)&quadVertices[0]);
    RETURN_IF_FAILED(hr);

    // Create the index buffer (and release the previous one)
    hr = RecreateUnnamedImmutableBuffer(D3D10_BIND_INDEX_BUFFER,
            (UINT)quadIndices.size()*2,
            state.pEvaluationIndexBuffer, (void *)&quadIndices[0]);
    RETURN_IF_FAILED(hr);

    state.nValidEvaluationIndices = 6 * N * N;

    ////////////////////////////////////////////////////////////

    bTessChange = false;
    OUT1("Tessellation Factor = " << N);

    OUT1_DONE();

    return hr;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Pre-process mesh -- convert polygons with even number of sides >4 to
// quads.
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
preprocess_mesh(PolygonMesh &m)
{
    OUT1_BEGIN("CHECKING and modifying the control mesh as needed...");

    ////////////////////////////////////////////////////////////
    // TEMPORARY HACK
    // Replace polar structures with n-gons (this is a hack around
    // Blender's current inability to model n-gons).

    bool flag = true;
    while (flag) {
        flag = false;
        for (PolygonMesh::Vertex_iterator vi = m.vertices_begin(); vi != m.vertices_end(); ++vi) {
            if (vi->flat_polar() && (int)vi->vertex_degree() <= MAX_FACE_SIZE) {
                OUT1("Detected flat polar structure. Converting to " << vi->vertex_degree() << "-gon.");
                //  NOT SAFE FOR A TETRAHEDON OR ANY GENERAL 2-LAYER DIPOLAR OBJECT.
                //  (That's because it would end up in two facets sharing
                //  the same vertices and having opposite orientations.)
                m.erase_center_vertex(vi->halfedge());
                flag = true;
                break;
            }
        }
    }
    // END TEMPORARY HACK
    ////////////////////////////////////////////////////////////

    // Convert even-sided extraordinary n-gons to quads.
    for (PolygonMesh::Facet_iterator fi = m.facets_begin(); fi != m.facets_end(); ++fi) {
        const int num_sides = (int)fi->facet_degree();
        if (num_sides != 4 && num_sides % 2 == 0) {
            OUT1("Detected " << num_sides << "-gon. Converting to quads.");

            // Figure out what the new point should be (averaging could
            // result in too much flattenning).
            RealPoint p(0,0,0);
            PolygonMesh::Halfedge_handle h = fi->halfedge();
            for(int i = 0; i < num_sides; ++i) {
                p += KToRealPoint(h->vertex()->point());
                h = h->next();
            }
            p /= num_sides;

            // Create the central vertex
            h = m.create_center_vertex(fi->halfedge());
            h->vertex()->point() = RealToKPoint(p);

            // Merge the resulting triangles into quads
            for(int i = 0; i < num_sides/2; ++i) {
                h = m.join_facet(h);
                h = h->next()->next()->opposite()->prev();
            }
        }
    }

    // For diagnostics, return the valence stats on the mesh.
    vector<int> vals(500, 0);
    int max_val = 0;
    for (PolygonMesh::Vertex_iterator vi = m.vertices_begin(); vi != m.vertices_end(); ++vi) {
        int v = (int)vi->halfedge()->vertex_degree();
        assert(v >= 0);
        if (v > max_val) max_val = v;
        if (v >= (int)vals.size()) vals.resize(v+1);
        ++vals[v];
    }

    OUT1_BEGIN("Number of vertices by valence:");
    for (int i = 0; i <= max_val; ++i) {
        if (vals[i] > 0) OUT1("Valence " << i << ": " << vals[i]);
    }
    OUT1_DONE();
    if (max_val > MAX_VALENCE)
        FATAL_ERROR("MAX_VALENCE (" << MAX_VALENCE << ") exceeded.");

    // The mesh is not allowed to have borders.
    if (!m.is_closed())
        FATAL_ERROR("Mesh has a border.");

    OUT1_DONE();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// General initialization
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

HRESULT
initialize( AppState &state, const string &filename )
{
    HRESULT hr = S_OK;

    state.TF=2;
    state.RowTex=0;

    /////////////////////////////////////////////////////////////////////
    // Do basic DirectX 10 setup
    /////////////////////////////////////////////////////////////////////

    OUT1_BEGIN("Doing basic DirectX setup (create window, setup callbacks, and load shader)...");

    if( !DXUTIsD3D10Available() )
        return DXUTERR_NODIRECT3D; // re-use an existing error since defining our own is a pain

    //hr = DXUTInit();
    hr = DXUTInit( true, true, NULL );
    hr = FAILED(hr) ? hr : DXUTCreateWindow(L"N-gon patches" );
    // Create the DirectX device
    if (bFullScreen)
        hr = FAILED(hr) ? hr : DXUTCreateDevice( false, 0, 0 );
    else
        hr = FAILED(hr) ? hr : DXUTCreateDevice( true, 800, 800 );
    RETURN_IF_FAILED(hr);

    //
    // Make sure we support Direct3D 10, and create a D3D10 device that will render
    // on to the application window.
    //
    D3DX10CreateFont( DXUTGetD3D10Device(), 15, 0, FW_BOLD, 1, FALSE, DEFAULT_CHARSET, 
            OUT_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH | FF_DONTCARE, 
            L"Arial", &g_pFont10 ) ;
    TrackResource(g_pFont10); // track for automatic release
    D3DX10CreateSprite( DXUTGetD3D10Device(), 512, &g_pSprite10 );
    TrackResource(g_pSprite10); // track for automatic release
    g_pTxtHelper = new CDXUTTextHelper( NULL, NULL, g_pFont10, g_pSprite10, 15 );

    RETURN_IF_FAILED(hr);

    // Set a callback that DXUT can call when the frame needs rendering
    DXUTSetCallbackDeviceChanging( ModifyDeviceSettings );
    DXUTSetCallbackMsgProc( MsgProc );
    DXUTSetCallbackFrameMove( OnFrameMove );
    DXUTSetCallbackD3D10SwapChainResized( OnD3D10SwapChainResized );
    //DXUTSetCallbackD3D10SwapChainReleasing( OnD3D10ReleasingSwapChain );
    //DXUTSetCallbackD3D10DeviceDestroyed( OnD3D10DestroyDevice );
    DXUTSetCallbackD3D10FrameRender( &RenderFrame, (void*)&state );
    DXUTSetCallbackKeyboard( OnKeyboardEvent, (void*)&state );

    // Setup the camera's view parameters
    D3DXVECTOR3 vecEye(0.0f, 0.0f, -8.0f);
    //D3DXVECTOR3 vecAt(0.5*(xmin+xmax),0.5*(ymin+ymax),0.5*(zmin+zmax));
    D3DXVECTOR3 vecAt(0, 0, 0);
    g_Camera.SetViewParams( &vecEye, &vecAt );
    //the near plane cannot start with 0.0f
    g_Camera.SetProjParams( D3DX_PI/4, 1.0, 0.05f, 800.0f );

    /////////////////////////////////////////////////////////////////////
    // Load the shader and make some basic variable correspondences
    /////////////////////////////////////////////////////////////////////

    //
    // Create the resources that we'll
    // use later for drawing
    //
    DWORD dwShaderFlags = D3D10_SHADER_ENABLE_STRICTNESS;
#if defined( DEBUG ) || defined( _DEBUG )
    // Set the D3D10_SHADER_DEBUG flag to embed debug information in the shaders.
    // Setting this flag improves the shader debugging experience, but still allows 
    // the shaders to be optimized and to run exactly the way they will run in 
    // the release configuration of this program.
    dwShaderFlags |= D3D10_SHADER_DEBUG;
#endif
    ID3D10Blob* pBlob = NULL;
    hr=D3DX10CreateEffectFromFile( L"surface.fx", NULL, NULL, "fx_4_0", dwShaderFlags, 0, DXUTGetD3D10Device(), NULL, NULL, &state.pEffect10, &pBlob, NULL  );
    TrackResource(state.pEffect10); // track for automatic release

    LPVOID l_pError = NULL;
    if( pBlob )
    {
        l_pError = pBlob->GetBufferPointer();
        // then cast to    a char* to see it in the locals window
        OutputDebugString( (LPCWSTR)pBlob->GetBufferPointer() );
		cout << "SHADER MESSAGES:\n";
		const int SIZE = 1000 ;
		char buf[SIZE+1];
		int l = (int)strlen((char*)l_pError);
		int s = SIZE < l ? SIZE : l;
		memcpy(buf, l_pError, s);
		buf[s] = (char)0;
		cout << (char*)buf << endl;
        pBlob->Release();
    }
    if( FAILED( hr ) )
    {
        MessageBox( NULL, L"The FX file cannot be located or there was a shader compilation error.", L"Error", MB_OK );
        return hr;
    }

    // Set up variable correspondences
    g_pRowTex = state.pEffect10->GetVariableByName( "gRowTex" )->AsScalar();
    g_pStart  = state.pEffect10->GetVariableByName( "gStart"  )->AsScalar();

    // Obtain the parameter handles
    g_pLightVecW1       = state.pEffect10->GetVariableByName( "gLightVecW1"      )->AsVector();
    g_pMatWorldView     = state.pEffect10->GetVariableByName( "matWorldView"     )->AsMatrix();
    g_pMatWorldViewProj = state.pEffect10->GetVariableByName( "matWorldViewProj" )->AsMatrix();
    g_pEyePosW          = state.pEffect10->GetVariableByName( "gEyePosW"         )->AsVector();
    g_pFlags            = state.pEffect10->GetVariableByName( "gFlags"           )->AsScalar();

    // Set the initial settings
    state.flags = 0;
    set_bit(state.flags, BIT_GROUP_COLOR);
    //set_bit(state.flags, BIT_DARK_BG);
    g_pFlags->SetInt( state.flags );

    OUT1_DONE();

    /////////////////////////////////////////////////////////////////////
    // Set up the buffers
    /////////////////////////////////////////////////////////////////////

    OUT1_BEGIN("SETTING UP buffers...");

    // Domain tessellation
    hr = tessellation(state);
    RETURN_IF_FAILED(hr);

    /////////////////////////////////////////////////////////////////////

    g_numFrames = bFrog ? FROG_FRAMES : MORPH_FRAMES;

    // Vertex Location + Valence Definition Begin : Texture
    float* vLocation;
    hr = BEGIN_CreateMutableTexture2D(state.pEffect10, "gVertexLocation",
            MAX_ALL, g_numFrames, DXGI_FORMAT_R32G32B32A32_FLOAT,
            state.VertexLocation,
            (void **)&vLocation);
    RETURN_IF_FAILED(hr);

    // Offset Definition Begin : Texture
    UINT32* pOffset;
    hr = BEGIN_CreateMutableTexture2D(state.pEffect10, "gOffsetData",
            MAX_ALL, 1, DXGI_FORMAT_R32_UINT,
            state.OffsetData,
            (void **)&pOffset);
    RETURN_IF_FAILED(hr);
    
    // 1-ring Index Definition Begin : Texture
    float* pRIndex;
    hr = BEGIN_CreateMutableTexture2D(state.pEffect10, "gRingIndex",
            64, MAX_ALL, DXGI_FORMAT_R32_FLOAT,
            state.RingIndex,
            (void **)&pRIndex);
    RETURN_IF_FAILED(hr);

    ////////////////////////////////////////////////////////////
    // Load files
    ////////////////////////////////////////////////////////////
    string inFileName = (bFrog) ? filename + "/monsterfrog_mapped_00.off" : filename;
    ifstream infile(inFileName.c_str(), ios::in);
    if (!infile) { FATAL_ERROR("Input OFF File couldn't be opened"); } 
    infile >> P; 
    if (!infile) { FATAL_ERROR("Error reading mesh from the input OFF file"); } 
    //if (!P.is_pure_quad()) { FATAL_ERROR("Input Mesh must be quadrilateral"); }
    preprocess_mesh(P);
    P.normalize_border();

    ////////////////////////////////////////////////////////////
    // Reorder all the facets
    ////////////////////////////////////////////////////////////

    vector<Halfedge_handle> facet_list;
    Facet_iterator fIter = P.facets_begin();

    // Mark all the polar vertices (needed to check whether facets are
    // polar).
    for (Vertex_iterator vv = P.vertices_begin(); vv != P.vertices_end(); ++vv) {
        vv->is_polar() = vv->is_polar_vertex_with_regular_links(1);   
    }

    ////////////////////
    // Count each type
    g_numFaces = 0;
    for(int i = 0; i < (int)P.size_of_facets(); ++i, ++fIter)
    {
        Halfedge_handle h = fIter->halfedge();
        const int num_sides = (int)h->facet_degree();

        // It is assumed that there are no n-gons for n even and n!=4.
        assert(num_sides == 4 || num_sides % 2 != 0);

        // Ignore facets that are too large
        if (num_sides > MAX_FACE_SIZE) continue;

        ++g_numFaces;

        // Count the number of irregular facets
        if (fIter->facet_degree() == 4
                && h                        ->vertex_degree() == 4
                && h->next()                ->vertex_degree() == 4
                && h->next()->next()        ->vertex_degree() == 4
                && h->next()->next()->next()->vertex_degree() == 4
//              && h                        ->opposite()->facet_degree() == 4
//              && h->next()                ->opposite()->facet_degree() == 4
//              && h->next()->next()        ->opposite()->facet_degree() == 4
//              && h->next()->next()->next()->opposite()->facet_degree() == 4
                )
        {
            fIter->type() = REG;
        }

        // Polar case
        else if ( fIter->is_polar_face(&h) )
        {
            fIter->type() = POLAR;
        }

        // n-gon case
        else
        {
            fIter->type() = num_sides;
        }

        ++g_num[fIter->type()];
    }

    // Compute the start indices for all the irregular facets and set them
    // as the shader variables.
    for (int i = 0; i <= MAX_TYPE; ++i) {
        g_start[i] = g_start[i-1] + g_num[i-1];
    }
    g_pStart->SetIntArray( g_start, 0, MAX_TYPE+1 );

    ////////////////////
    // Put Regular+Irregular patches in order of type
    fIter = P.facets_begin();
    facet_list.resize(g_numFaces);
    int typeCounter[MAX_TYPE+1] = {0};
    for(int i = 0; i < (int)P.size_of_facets(); ++i, ++fIter)
    {
        Halfedge_handle h = fIter->facet_begin();
        const int num_sides = (int)h->facet_degree();

        // Ignore facets that are too large
        if (num_sides > MAX_FACE_SIZE) continue;

        int type = fIter->type();

        // Set h to the polar edge if the facet is polar
        if (fIter->is_polar_face(&h)) {
            h = h->next();
            assert(h->prev()->vertex()->is_polar());
        }

        facet_list[g_start[type] + typeCounter[type]] = h;
        ++typeCounter[type];
    }

    // sanity check: counters should match #elements
    for (int i = 0; i <= MAX_TYPE; ++i)
        assert(g_num[i] == typeCounter[i]);

    OUT1("Total number of facets to be processed: " << g_numFaces);
    OUT1("    Number of regular facets: " << g_num[REG]);
    OUT1("    Number of polar facets:   " << g_num[POLAR]);
    for (int i = 3; i <= MAX_TYPE; ++i) {
        OUT1("    Number of irregular " << i << "-gons: " << g_num[i]);
    }
    OUT1("Total number of facets (including those ignored): " << P.size_of_facets());
    
    ////////////////////////////////////////////////////////////
    // Vertex buffer
    ////////////////////////////////////////////////////////////

    vVertices.resize(P.size_of_vertices());
    vIndices.resize(P.size_of_facets()*6);

    //setup vertex buffer and index buffer
    Halfedge_handle h1[MAX_VALENCE], h2[MAX_VALENCE];
    int vCounter, numFacet, numVert;
    vCounter = 0;
    Vertex_iterator vv = P.vertices_begin();
    
    // vertex location
    do {
        vv->t_index = vCounter;

        if (vv->vertex_degree() > MAX_VALENCE) WARN("MAX_VALENCE is greater than "<<MAX_VALENCE<<": "<<vv->vertex_degree());
        Point p = vv->point();
        if (p.x() > xmax) xmax = p.x();
        if (p.x() < xmin) xmin = p.x();
        if (p.y() > ymax) ymax = p.y();
        if (p.y() < ymin) ymin = p.y();
        if (p.z() > zmax) zmax = p.z();
        if (p.z() < zmin) zmin = p.z();

        vCounter++;

    } while ( ++vv != P.vertices_end() );

    vv=P.vertices_begin();
    vCounter=0;
    float avgx = (xmax + xmin) / 2;
    float avgy = (ymax + ymin) / 2;
    float avgz = (zmax + zmin) / 2;
    do {
        // Shifting points to the origin
        Point &p=vv->point();
        p = Point(p.x() - avgx, p.y() - avgy, p.z() - avgz);

        // Save Vertex Location to vertex texture
        vLocation[vCounter*4  ] = p.x();
        vLocation[vCounter*4+1] = p.y();
        vLocation[vCounter*4+2] = p.z();
        vLocation[vCounter*4+3] = (float)vv->vertex_degree();
        vCounter++;

    } while ( ++vv != P.vertices_end() );

    // Find longest axis for Morphing
    if((xmax-xmin)>(ymax-ymin))
    {
        if((xmax-xmin)>(zmax-zmin)){
            longAxis=1; maxAmount=0.5f*(xmax-xmin);
        }
        else {
            longAxis=3; maxAmount=0.5f*(zmax-zmin);
        }
    }
    else
    {
        if((ymax-ymin)>(zmax-zmin)){
            longAxis=2; maxAmount=0.5f*(ymax-ymin);
        }
        else {
            longAxis=3; maxAmount=0.5f*(zmax-zmin);
        }
    }

    // Load all the frames of the animation
    LoadAnimationData(vLocation, filename);

    END_CreateMutableTexture2D(state.VertexLocation);

    ////////////////////////////////////////////////////////////
    // Rotational offset buffer
    ////////////////////////////////////////////////////////////

    for(int i = 0; i < (int)P.size_of_facets(); ++i)
    {
        Halfedge_handle h = facet_list[i];
        int num_sides = (int)h->facet_degree();
        pOffset[i] = 0;
        for(int j = 0; j < num_sides; ++j)
        {
            // The half-edge itself is along the v-coordinate, not the
            // u-coordinate. Need the rotation index of the u-coordinate.
            int r = (int)(h->next()->opposite()->rotation_index_of_edge());
            assert(0 <= r && r < (int)h->vertex_degree() && h->vertex_degree() <= 1<<NBITS_OFFSET);
            pOffset[i] |= (r << (NBITS_OFFSET*j));
            h = h->next();
        }
    }

    END_CreateMutableTexture2D(state.OffsetData);

    ////////////////////////////////////////////////////////////
    // Ring Index buffer
    ////////////////////////////////////////////////////////////

    // save index of 1-ring vertices in the index texture per vertex
    Vertex_iterator vIter = P.vertices_begin();

    int rIndexStart=0;
    for(int i=0; i<(int)P.size_of_vertices(); ++i)
    {
        
        HV_circulator hvC = vIter->vertex_begin();
        //cout<<vIter->vertex_begin()->vertex()->t_index<<" -> "<<rIndexStart<<" : ";
        
        // save ring index
        
        for(int j=0; j<(int)vIter->vertex_begin()->vertex_degree(); ++j)
        {
            pRIndex[i*64+j*3  ]=(float)(hvC->opposite()->vertex()->t_index);
            pRIndex[i*64+j*3+1]=(float)(hvC->opposite()->next()->vertex()->t_index);
            pRIndex[i*64+j*3+2]=(float)(hvC->opposite()->prev()->prev()->opposite()->vertex()->t_index);
            --hvC;
        }
         
        vVertices[i].valence = (int)vIter->vertex_begin()->vertex_degree();
        
        // assertion for max valence is smaller or equal to MAX VALENCE
        assert(vVertices[i].valence <= MAX_VALENCE);
        
        ++vIter;
    }

    END_CreateMutableTexture2D(state.RingIndex);
    // end of ring index

    // Create the vertex buffer
    hr = CreateUnnamedImmutableBuffer(D3D10_BIND_VERTEX_BUFFER, (UINT)vVertices.size()*sizeof(vVertices[0]),
            state.pVertexBuffer, (void *)&vVertices[0]);
    RETURN_IF_FAILED(hr);

    //////////////////////////////////////////////////////////////
    // Facet indices (input to the geometry shader)
    //////////////////////////////////////////////////////////////

    numFacet=0; numVert=0; vCounter=0;
    for (numFacet = 0; numFacet < (int)facet_list.size(); ++numFacet) {
        Halfedge_handle hf = facet_list[numFacet];

        // Pick the starting halfedge carefully in the polar case.
//      PolygonMesh::Halfedge_const_handle hf;
//      if (f->type() == POLAR) { // Polar case
//          bool test = f->is_polar_face(&hf);
//          assert(test && "fIter->is_polar_face();");
//          assert(hf->vertex()->is_polar());
//          hf = hf->next(); // the polar vertex will be enumerated last
//      }
//      else { // N-gon case
//          hf = f->halfedge();
//      }

        // Copy the vertex indices of the n-gon.
        for (int i = 0; i < (int)hf->facet_degree(); ++i, hf = hf->next()) {
            vIndices[6*numFacet+i] = hf->vertex()->t_index;
        }
        // Set the rest of the vertex indices to the first one.
        for (int i = (int)hf->facet_degree(); i < 6; ++i) {
            vIndices[6*numFacet+i] = vIndices[6*numFacet];
        }

        //if (numFacet >= g_num[REG]) {
        //    vIndices[6*numFacet  ]=0; vIndices[6*numFacet+1]=0; vIndices[6*numFacet+2]=0;
        //    vIndices[6*numFacet+3]=0; vIndices[6*numFacet+4]=0; vIndices[6*numFacet+5]=0; 
        //}
    }

    // Create the index buffer
    hr = CreateUnnamedImmutableBuffer(D3D10_BIND_INDEX_BUFFER, (UINT)vIndices.size()*sizeof(vIndices[0]),
            state.pIndexBuffer, (void *)&vIndices[0]);
    RETURN_IF_FAILED(hr);

    ////////////////////////////////////////////////////////////
    // Mesh Indices (for mesh display)
    ////////////////////////////////////////////////////////////

    numFacet=0; numVert=0; vCounter=0;

    vMeshIndices.reserve(2*P.size_of_halfedges()); // line list - 2 indices per halfedge
    for (numFacet = 0; numFacet < (int)facet_list.size(); ++numFacet) {
        Halfedge_handle hf = facet_list[numFacet];

        // Copy the vertex indices of the n-gon.
        for (int i = 0; i < (int)hf->facet_degree(); ++i, hf = hf->next()) {
            vMeshIndices.push_back( hf->vertex()->t_index );
            vMeshIndices.push_back( hf->next()->vertex()->t_index );
        }
    }

    // Create the index buffer
    hr = CreateUnnamedImmutableBuffer(D3D10_BIND_INDEX_BUFFER, (UINT)vMeshIndices.size()*sizeof(vMeshIndices[0]),
            state.pMeshIndexBuffer, (void *)&vMeshIndices[0]);
    RETURN_IF_FAILED(hr);

    ////////////////////////////////////////////////////////////
    // Control Mesh (Stream out)
    ////////////////////////////////////////////////////////////

    UINT element_width;

    // Create a different buffer for each patch type.
    for (int i = 0; i <= MAX_TYPE; ++i) {
        if (g_num[i] == 0) continue;

        element_width = ( MAX_CONTROL_POINTS[i]*(g_num[i]+1)*NUM_OF_INSTANCES );
        assert(element_width > 0);
        hr = CreateNamedBufferWithView(state.pEffect10, str_index("gControlPoints", i),
                D3D10_BIND_STREAM_OUTPUT | D3D10_BIND_SHADER_RESOURCE,
                element_width * sizeof(D3DXVECTOR4),    // Size
                0,                                      // No CPU access required
                D3D10_USAGE_DEFAULT,                    // Can't be changed after creation
                DXGI_FORMAT_R32G32B32A32_FLOAT,         // Format
                element_width,                          // Data element width
                state.ControlPoints[i],
                NULL); // no data

        RETURN_IF_FAILED(hr);
    }

    ////////////////////////////////////////////////////////////
    // Set up the shader techniques
    ////////////////////////////////////////////////////////////

    OUT1("Set up Shader Techniques...");

    int VALID_INDEX = -1;

    // Get techniques from the shader
    for (int i = 0; i <= MAX_TYPE; ++i) {
		// Do not load the shaders that will not get used. For some reason
        // NOT using loaded shaders results in a HUGE performance hit.
        if (g_num[i] == 0) continue;
        VALID_INDEX = i;

        int num_sides = NUM_SIDES[i];
        state.NgonPass[i].initialize(state, str_index("NgonPass", i), &state.ControlPoints[i], num_sides);
    }
    assert(VALID_INDEX != -1);

    // Create the layouts
    D3D10_PASS_DESC PassDesc;
    // Construction pass layout
    hr = state.NgonPass[VALID_INDEX].technique->GetPassByIndex( 0 )->GetDesc( &PassDesc );
    RETURN_IF_FAILED(hr);
    hr = DXUTGetD3D10Device()->CreateInputLayout( ConstructionVertexFormatElements, sizeof(ConstructionVertexFormatElements)/sizeof(ConstructionVertexFormatElements[0]), 
            PassDesc.pIAInputSignature, PassDesc.IAInputSignatureSize, &state.pConstructionLayout ) ;
    RETURN_IF_FAILED(hr);
    TrackResource(state.pConstructionLayout); // track for automatic release
    // Evaluation pass layout
    hr = state.NgonPass[VALID_INDEX].technique->GetPassByIndex( 1 )->GetDesc( &PassDesc );

    RETURN_IF_FAILED(hr);
    hr = DXUTGetD3D10Device()->CreateInputLayout( EvaluationVertexFormatElements, sizeof(EvaluationVertexFormatElements)/sizeof(EvaluationVertexFormatElements[0]), 
            PassDesc.pIAInputSignature, PassDesc.IAInputSignatureSize, &state.pEvaluationLayout ) ;
    RETURN_IF_FAILED(hr);
    TrackResource(state.pEvaluationLayout); // track for automatic release

    // Mesh rendering pass (same input layout as construction pass)
    state.pMeshTechnique = state.pEffect10->GetTechniqueByName( "MeshPass" );

    OUT1_DONE();

    return hr;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Cleanup
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
clean_up( AppState &state, HRESULT hr )
{
    OUT1_BEGIN("Cleaning up...");

    //
    // Cleanup
    //
    if( DXUTGetD3D10Device() ) DXUTGetD3D10Device()->ClearState();

    // Release all the tracked resources (including those
    // created/initialized by functions in util.h).
    ReleaseTrackedResources();

    SAFE_DELETE( g_pTxtHelper );

    DXUTShutdown( hr );

    OUT1_DONE();
}

#endif // INITIALIZE_H

