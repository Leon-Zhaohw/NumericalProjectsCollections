/*
This file is part of Cubica.
 
Cubica is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Cubica is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Cubica.  If not, see <http://www.gnu.org/licenses/>.
*/
// PIXIE.h: interface for the RENDER_MAN class.
//
//////////////////////////////////////////////////////////////////////

#include "RENDER_MAN.h"

#define USING_TEAPOT 0

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////////////
RENDER_MAN::RENDER_MAN()
{
#if USING_RENDERMAN
  _qualityRender = true;
  _commandLine = false;
  _boxOn = false;
  _blurOn = false;
  _white[0] = 1.0;
  _white[1] = 1.0;
  _white[2] = 1.0;
  _red[0] = 1.0;
  _red[1] = 0.0;
  _red[2] = 0.0;

  //_liverRed[0] = 152.0 / 255.0;
  //_liverRed[0] = 200.0 / 255.0 / 2;
  //_liverRed[1] = 46.0 / 255.0 / 2;
  //_liverRed[2] = 56.0 / 255.0 / 2;
  _liverRed[0] = 200.0 / 255.0 / 4;
  _liverRed[1] = 46.0 / 255.0 / 4;
  _liverRed[2] = 56.0 / 255.0 / 4;

  //_spinY = -22.5;
  _spinX = 0;
  _spinY = -90.0;
  _spinZ = 180;

  //_xRes = 640/4;
  //_yRes = 480/4;
  //_xRes = 640/2;
  //_yRes = 480/2;
  _xRes = 640;
  _yRes = 480;
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// turn on ray tracing
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::enableRaytracing()
{
#if USING_RENDERMAN
  // turn on ray tracing and set the bias
  RtInt trace = 1;
  RiAttribute("visibility", "trace", (RtPointer)&trace, RI_NULL);
  //RtFloat bias = 0.01;
  RtFloat bias = 0.01;
  RiAttribute("trace", "bias", (RtPointer)&bias, RI_NULL);

  RtInt transmission = 1; 
  //RtString transmissionhitmode = "primitive";
  // this kills subsurface
  RiAttribute("visibility", "int transmission", (RtPointer)&transmission, RI_NULL);
  //RiAttribute("shade", "transmissionhitmode", (RtPointer)&transmissionhitmode, RI_NULL);

  RtInt diffuse = 1;
  RiAttribute("visibility", "int diffuse", (RtPointer)&diffuse, RI_NULL);
  // turns on specular reflections
  RtInt specular = 1;
  RiAttribute("visibility", "int specular",(RtPointer)&specular, RI_NULL);
#endif 
}

///////////////////////////////////////////////////////////////////////////////////////
// turn on ambient occlustion
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::enableAmbientOcclusion()
{
#if USING_RENDERMAN
  RiShadingInterpolation("smooth");
  RiShadingRate(0.01f);
  RtInt trace = 1;
  RtFloat bias = 0.005;
  //RtInt samples = 100;
  RtInt samples = 10;
  RiAttribute("visibility", "trace", (RtPointer)&trace, RI_NULL);
  RiAttribute("trace", "bias", (RtPointer)&bias, RI_NULL);
  //RiSurface("occsurf2", "samples", (RtPointer)&samples, RI_NULL);
  RiSurface("ambientocclusion", "samples", (RtPointer)&samples, RI_NULL);
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// turn on a normal light
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::enableExplicitLight()
{
#if USING_RENDERMAN
  RtFloat intensity0 = 1.0f;
  //RtFloat from[] = {10.0f, 10.0f, 0.0f};
  //RtFloat to[] = {-10.0f, -10.0f, 0.0f};
  RtFloat from0[] = {1.0f, 1.0f, 0.0f};
  RtFloat to0[] = {-1.0f, -1.0f, 0.0f};
  RiLightSource("distantlight", 
                "intensity", (RtPointer)&intensity0,
                "from", (RtPointer)from0,
                "to", (RtPointer)to0,
                "lightcolor", (RtPointer)_white,
                RI_NULL);

  /*
  // DEBUG LIGHT
  RtFloat from2[] = {-1.0f, -1.0f, 0.0f};
  RtFloat to2[] = {1.0f, 1.0f, 0.0f};
  RiLightSource("distantlight", 
                "intensity", (RtPointer)&intensity0,
                "from", (RtPointer)from2,
                "to", (RtPointer)to2,
                "lightcolor", (RtPointer)_white,
                RI_NULL);
                */

  RtFloat intensity1 = 10.0f;
  RtFloat from1[] = {3.0f, 3.0f, -3.0f};
  RiLightSource("pointlight", 
  //RiLightSource("pointlight_rts", 
                "intensity", (RtPointer)&intensity1,
                "from", (RtPointer)from1,
                "lightcolor", (RtPointer)_white,
                RI_NULL);
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// turn on raytraced environment map
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::enableLiverSkydome()
{
#if USING_RENDERMAN
  RtString envmap = "./shaders/paris-cropped-blurred-pixar.exr";
  //RtString envmap = "./shaders/paris-cropped-blurred-new-pixar.exr";
  RtColor filter = {0.5, 0.5, 0.5};
  RiDeclare("envmap", "uniform string");
  RiDeclare("filter", "uniform color");
  RiLightSource("./shaders/environmentlight",
                "envmap", &envmap, 
                "filter", &filter,
                RI_NULL);
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// turn on raytraced light
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::enableRaytracedLiverLight()
{
#if USING_RENDERMAN
  RtFloat intensity1 = 4.0f;
  RtFloat from1[] = {2.5f, 3.0f, 3.0f};
  RtColor skinColor = {0.1, 0.1, 0.1};
  /*
  // this is the low quality light
  RiLightSource("pointlight_rts", 
                "intensity", (RtPointer)&intensity1,
                "from", (RtPointer)from1,
                "lightcolor", (RtPointer)skinColor,
                RI_NULL);
                */

  // this is the high quality light
  RtFloat intensity2 = 2.0f;
  RtFloat from2[] = {3.0f, 3.0f, -3.0f};
  RiLightSource("pointlight_rts", 
                "intensity", (RtPointer)&intensity2,
                "from", (RtPointer)from2,
                "lightcolor", (RtPointer)skinColor,
                RI_NULL);

  enableLiverSkydome();

  /*
  //RtString envmap = "./shaders/paris-pixar.exr";
  //RtString envmap = "./shaders/paris-latlong-pixar.exr";
  //RtString envmap = "./shaders/paris-cross-pixar.exr";
  //RtString envmap = "./shaders/paris-cropped-pixar.exr";
  RtString envmap = "./shaders/paris-cropped-blurred-pixar.exr";
  //RtColor filter = {1, 1, 1};
  RtColor filter = {0.5, 0.5, 0.5};
  RiDeclare("envmap", "uniform string");
  RiDeclare("filter", "uniform color");
  RiLightSource("./shaders/environmentlight",
                "envmap", &envmap, 
                "filter", &filter,
                RI_NULL);
                */
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// turn on raytraced light
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::enableRaytracedDragonLight()
{
#if USING_RENDERMAN
  RtFloat intensity0 = 1.5f;
  RtFloat from0[] = {-2.0f, 2.0f, 2.0f};
  RiLightSource("pointlight_rts", 
                "intensity", (RtPointer)&intensity0,
                "from", (RtPointer)from0,
                "lightcolor", (RtPointer)_white,
                RI_NULL);

  /*
  RtFloat intensity1 = 1.0f;
  RtFloat from1[] = {-2.5f, -3.0f, 3.0f};
  RiLightSource("pointlight_rts", 
                "intensity", (RtPointer)&intensity1,
                "from", (RtPointer)from1,
                "lightcolor", (RtPointer)_white,
                RI_NULL);
                */

  RtFloat intensity2 = 1.0f;
  RtFloat from2[] = {3.0f, 3.0f, 3.0f};
  RiLightSource("pointlight_rts", 
                "intensity", (RtPointer)&intensity2,
                "from", (RtPointer)from2,
                "lightcolor", (RtPointer)_white,
                RI_NULL);
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// turn on raytraced light
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::enableRaytracedLight()
{
#if USING_RENDERMAN
  RtFloat intensity0 = 1.0f;
  //RtFloat from[] = {10.0f, 10.0f, 0.0f};
  //RtFloat to[] = {-10.0f, -10.0f, 0.0f};
  RtFloat from0[] = {1.0f, 1.0f, 0.0f};
  RtFloat to0[] = {-1.0f, -1.0f, 0.0f};
  /*
  RiLightSource("distantlight_rts", 
                "intensity", (RtPointer)&intensity0,
                "from", (RtPointer)from0,
                "to", (RtPointer)to0,
                "lightcolor", (RtPointer)_white,
                RI_NULL);
                */

  //RtFloat intensity1 = 10.0f;
  //RtFloat intensity1 = 5.0f;
  RtFloat intensity1 = 7.5f;
  RtFloat from1[] = {-2.5f, -3.0f, -3.0f};
  RtColor skinColor = {1.0, 0.941, 0.847};
  RtFloat samples = 10;
  RiLightSource("pointlight_rts", 
                "samples", (RtPointer)&samples,
                "intensity", (RtPointer)&intensity1,
                "from", (RtPointer)from1,
                "lightcolor", (RtPointer)skinColor,
                RI_NULL);

  RtFloat intensity2 = 2.0f;
  RtFloat from2[] = {3.0f, 3.0f, -3.0f};
  RiLightSource("pointlight_rts", 
                "intensity", (RtPointer)&intensity2,
                "samples", (RtPointer)&samples,
                "from", (RtPointer)from2,
                "lightcolor", (RtPointer)skinColor,
                RI_NULL);

  // Doesn't work for subsurface
#if 0      
  RiTransformBegin();
  {
    // key light
    RiTransformBegin();
    {
      RtMatrix matrix = {
        { 0.813101, 0, 0.582123, 0  },
        {-0.415911, 0.699663, 0.580938, 0  },
        {0.40729, 0.714473, -0.568897, 0  },
        {368.851, 630.315, -556.897, 1}};
      RiTransform(matrix);

      RtFloat coneAngle = 0.698132;
      RtFloat penumbraAngle = 0.20944;
      RiDeclare("coneAngle", "uniform float");
      RiDeclare("penumbraAngle", "uniform float");
      RiLightSource("key_spotLightShape1_",
                    "coneAngle", (RtPointer)&coneAngle,
                    "penumbraAngle", (RtPointer)&penumbraAngle,
                    RI_NULL);
    }
    RiTransformEnd();

    // back light
    RiTransformBegin();
    {
      RtMatrix matrix = {
        {-0.910684, 0, -0.413104, 0 }, 
        {  0.0517757, 0.992115, -0.114139, 0  },
        {-0.409847, 0.125333, 0.903503, 0  },
        {-236.715, -2.54367, 506.986, 1}};

      RiTransform(matrix);
      RtFloat coneAngle = 0.698132;
      RtFloat penumbraAngle = 0.20944;
      RiDeclare("coneAngle", "uniform float");
      RiDeclare("penumbraAngle", "uniform float");
      RiLightSource("back_spotLight2Shape_",
                    "coneAngle", (RtPointer)&coneAngle,
                    "penumbraAngle", (RtPointer)&penumbraAngle,
                    RI_NULL);
    }
    RiTransformEnd();

    // fill light
    RiTransformBegin();
    {
      RtMatrix matrix = {
        {0.855364, 0, -0.518027, 0},
        {-0.134075, 0.965926, -0.221385, 0  },
        {-0.500376, -0.258819, -0.826218, 0},  
        {-391.132, -260.053, -648.299, 1}};
      RiTransform(matrix);
      RtFloat coneAngle = 0.698132;
      RtFloat penumbraAngle = 0.20944;
      RiDeclare("coneAngle", "uniform float");
      RiDeclare("penumbraAngle", "uniform float");
      RiLightSource("fill_spotLightShape3_",
                    "coneAngle", (RtPointer)&coneAngle,
                    "penumbraAngle", (RtPointer)&penumbraAngle,
                    RI_NULL);
    }
    RiTransformEnd();

    /*
    // eye light 1
    RiTransformBegin();
    {
      RtMatrix matrix = {
      {0.766044, -4.33681e-19, -0.642788, 0},
      {0.00673113, 0.999945, 0.00802185, 0},  
      {-0.642752, 0.0104718, -0.766002, 0},
      {-183.693, -49.1112, -278.43, 1}};
      RiTransform(matrix);
      RtFloat coneAngle = 0.698132;
      RtFloat penumbraAngle = 0.0;
      RiDeclare("coneAngle", "uniform float");
      RiDeclare("penumbraAngle", "uniform float");
      RiLightSource("Eye_HL_spotLightShape5_",
                    "coneAngle", (RtPointer)&coneAngle,
                    "penumbraAngle", (RtPointer)&penumbraAngle,
                    RI_NULL);
    }
    RiTransformEnd();

    // eye light 2
    RiTransformBegin();
    {
      RtMatrix matrix = {
    	  {0.924546, 0, 0.38107, 0},
        {0.0199437, 0.99863, -0.048387, 0 },
        {0.380548, -0.052336, -0.923279, 0 }, 
        {121.625, -42.8316, -406.604, 1}};
      RiTransform(matrix);
      RtFloat coneAngle = 0.698132;
      RtFloat penumbraAngle = 0.0;
      RiDeclare("coneAngle", "uniform float");
      RiDeclare("penumbraAngle", "uniform float");
      RiLightSource("Eye_HL_spotLightShape4_",
                    "coneAngle", (RtPointer)&coneAngle,
                    "penumbraAngle", (RtPointer)&penumbraAngle,
                    RI_NULL);
    }
    RiTransformEnd();
    */
  }
  RiTransformEnd();
#endif
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// turn on an environment light
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::enableEnvironmentLight()
{
#if USING_RENDERMAN
  // squash settings
  float samples = 128;
  float Kenv = 1;
  //float Kenv = 5; // for squashing?
  float Kocc = 1;
  RiLightSource("envlight2", 
  //RiLightSource("envlightocclusion", 
                "float samples", (RtPointer)&samples,
                "string envmap", (RtPointer)&_envmap,
                "float Kenv", (RtPointer)&Kenv,
                "float Kocc", (RtPointer)&Kocc,
                RI_NULL);
  /*
  float samples = 10;
  float Kenv = 1;
  float Kocc = 1;
  RiLightSource("envlight2", 
  //RiLightSource("envlightocclusion", 
                "float samples", (RtPointer)&samples,
                "string envmap", (RtPointer)&envmap,
                "float Kenv", (RtPointer)&Kenv,
                "float Kocc", (RtPointer)&Kocc,
                RI_NULL);
                */
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// turn on a sky sphere
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::enableSkySphere()
{
#if USING_RENDERMAN
  // sky sphere
  RiAttributeBegin();
    RtInt notrace = 0;
    RiAttribute("visibility", "trace", (RtPointer)&notrace, RI_NULL);
    RtString transparent = "transparent";
    RiAttribute("visibility", "transmission", (RtPointer)&transparent, RI_NULL);
    //RiScale(100,100,100);
    RiScale(30,30,30);
    RiSurface("skysphere", "string mapname", (RtPointer)&_envmap, RI_NULL);
    RiSphere(1.0, -1, 1, 360, RI_NULL);
  RiAttributeEnd();
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// turn on a solid color sky sphere
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::enableSolidColorSky()
{
#if USING_RENDERMAN
  // sky sphere
  RiAttributeBegin();
    RtInt notrace = 0;
    RiAttribute("visibility", "trace", (RtPointer)&notrace, RI_NULL);
    RtString transparent = "transparent";
    RiAttribute("visibility", "transmission", (RtPointer)&transparent, RI_NULL);
    RiScale(30,30,30);
    RiSurface("solidsky", RI_NULL);
    RiSphere(1.0, -1, 1, 360, RI_NULL);
  RiAttributeEnd();
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// enable a displacement shader
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::enableDisplacementShader(float frequency, float amplitude)
{
#if USING_RENDERMAN
  RtFloat amp = amplitude;
  RtFloat freq = frequency;
  //RiDisplacement("displaceNoise", RI_MAXAMP, (RtPointer)&maxAmp, RI_NULL);
  RiDisplacement("displaceNoise", 
                 "float freq", (RtPointer)&freq,
                 "float amp", (RtPointer)&amp, RI_NULL);

  RtFloat sphere = amp * 0.5;
  RtString space("shader");
  RiAttribute("displacementbound", "sphere", (RtPointer)&sphere, "space", (RtPointer)&space, RI_NULL);
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// disable the displacement shader
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::disableDisplacementShader()
{
#if USING_RENDERMAN
  RiDisplacement("", RI_NULL);
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// enable a floor at y = 0
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::enableFloor()
{
#if USING_RENDERMAN
  RiTransformBegin();
    RiRotate(90.0f, 0.0f, 0.0f, 1.0f);
    RiRotate(90.0f, 0.0f, 1.0f, 0.0f);
    RiDisk(0.0f, 100.0f, 360.0f, RI_NULL);
  RiTransformEnd();
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// draw a glass cube 
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::cube()
{
#if USING_RENDERMAN
  // total faces in mesh
  RtInt nfaces = 1;
  
  // number of vertices in each face
  RtInt nvertices[] = {4,4,4,4,4};

  // vertex indices - array length should be sum of all elements in nvertices
  RtInt ntags = 1;
  RtPoint P[] = {
    {-0.5, -0.5, 0}, {0.5, -0.5, 0}, {0.5, 0.5, 0}, {-0.5, 0.5, 0}
  };
  RtPoint N[] = {
    {0,0,-1},{0,0,-1},{0,0,-1},{0,0,-1}
  };
  RtInt vertices[] = {
    3,2,1,0
  };

  RiTransformBegin();
  RiTranslate(0.5, 0.5,0);
  
  // front
  RiTransformBegin();
    RiPointsPolygons(nfaces, nvertices, vertices, RI_P, P, RI_N, N, RI_NULL);
  RiTransformEnd();

  // back 
  RiTransformBegin();
    RiTranslate(0, 0, 1);
    RiRotate(180,0,1,0);
    RiPointsPolygons(nfaces, nvertices, vertices, RI_P, P, RI_NULL);
  RiTransformEnd();

  RiTransformBegin();
    RiTranslate(-0.5,0,0.5);
    RiRotate(90,0,1,0);
    RiPointsPolygons(nfaces, nvertices, vertices, RI_P, P, RI_NULL);
  RiTransformEnd();

  RiTransformBegin();
    RiTranslate(0.5,0,0.5);
    RiRotate(-90,0,1,0);
    RiPointsPolygons(nfaces, nvertices, vertices, RI_P, P, RI_NULL);
  RiTransformEnd();

  RiTransformBegin();
    RiTranslate(0,0.5,0.5);
    RiRotate(90,1,0,0);
    RiPointsPolygons(nfaces, nvertices, vertices, RI_P, P, RI_NULL);
  RiTransformEnd();

  RiTransformBegin();
    RiTranslate(0,-0.5,0.5);
    RiRotate(-90,1,0,0);
    RiPointsPolygons(nfaces, nvertices, vertices, RI_P, P, RI_NULL);
  RiTransformEnd();
  RiTransformEnd();
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// Draw a set of points
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawPoints(vector<VEC3F*>& points, bool sliced)
{
#if USING_RENDERMAN
  _red[0] = 2;
  RiColor(_red);
  for (int x = 0; x < points.size(); x++)
  {
    //if ((*points[x])[0] < 0.5)
    //if ((*points[x])[0] >= 0.45)
    {
      RiTransformBegin();
        VEC3F& point = *points[x];
        RiTranslate(point[0], point[1], point[2]);
        int res = 32;
        float scaling = 0.15 / res;
        RiScale(scaling, scaling, scaling);
        RiSphere(1.0f, -1.0f, 1.0f, 360.0f, RI_NULL);
      RiTransformEnd();
    }
  }
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// bake out the point cloud for the liver
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::bakeLiverSubsurface(string filename, TET_MESH* tetMesh)
{
#if USING_RENDERMAN
  cout << " Baking out subsurface for " << filename.c_str() << endl;

  // view settings
  float tiltCameraDown = 0;
  float zoomOut = 1.0;
  float cameraY = 0.5;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);

  string ribName = filename + string(".bake.rib");
  
  // zoomed in version
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    string bakedFilename = filename + string(".bake.tif");
    sprintf(buffer, "%s", bakedFilename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGBA, RI_NULL);
    RtFloat fov = fieldOfView;

    //RiPixelSamples(4,4);
    //RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiShadingInterpolation("smooth");
    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    RiFormat(_xRes, _yRes,1.0);
    //RiFormat(_xRes / 2, _yRes / 2,1.0);
    //RiFormat(640, 480,1.0);

    RiDisplayChannel("float _area", RI_NULL);
    RiDisplayChannel("color _radiance_t", RI_NULL);

    RiWorldBegin();
      RtInt zero = 0;
      enableRaytracing();
      RiAttribute("cull", "hidden", (RtPointer)&zero, RI_NULL);
      RiAttribute("cull", "backfacing", (RtPointer)&zero, RI_NULL);
      RiAttribute("dice", "rasterorient", (RtPointer)&zero, RI_NULL);

      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiTranslate(-0.5, -0.5, -0.5);

        enableRaytracedLiverLight();

        string ptcFilename = filename + string(".ptc");
        sprintf(buffer, "%s", ptcFilename.c_str());
        RtString pointCloud = buffer;
        RtString channels = "_area,_radiance_t";
        RiAttributeBegin();
        RiDeclare("filename", "uniform string");
        RiDeclare("displaychannels", "uniform string");
        RiSurface("./shaders/bake_radiance_t", 
                  (RtToken)"filename", (RtPointer)&pointCloud,
                  (RtToken)"displaychannels", (RtPointer)&channels,
                  RI_NULL);
        RiSides(1);

        //RiColor(_red);
        RiColor(_liverRed);
        tetMesh->drawEmbeddingToRenderMan();
        RiAttributeEnd();
      RiTransformEnd();
    RiWorldEnd();
  RiEnd();
  string prman = string("prman -progress ") + ribName;
  system(prman.c_str());
  string rmRib = string("rm ") + ribName;
  system(rmRib.c_str());
  string rmTiff = string("rm ") + bakedFilename;
  system(rmTiff.c_str());
  
  // Run diffusion on the point cloud
  string blurredPtc = filename + string(".sss.ptc");
  string ptfilter = string("ptfilter -filter ssdiffusion -partial 1 ");
  ptfilter += ptcFilename + string(" ");
  ptfilter += blurredPtc;
  cout << " Running " << ptfilter << endl;
  system(ptfilter.c_str());
  cout << " Point cloud filename: " << blurredPtc << endl;

  // remove the unblurred point cloud
  string rmPtc = string("rm ") + ptcFilename;
  system(rmPtc.c_str());
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// bake out the point cloud for the head
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::bakeHeadSubsurface(string filename, TET_MESH* tetMesh)
{
#if USING_RENDERMAN
  cout << " Baking out subsurface for " << filename.c_str() << endl;

  // view settings
  float tiltCameraDown = 11;
  float zoomOut = 0.6;
  float cameraY = 0.55;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);

  string ribName = filename + string(".bake.rib");
  
  // zoomed in version
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    string bakedFilename = filename + string(".bake.tif");
    sprintf(buffer, "%s", bakedFilename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGBA, RI_NULL);
    RtFloat fov = fieldOfView;

    //RiPixelSamples(4,4);
    //RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiShadingInterpolation("smooth");
    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    //RiFormat(640, 480,1.0);
    //RiFormat(320, 240,1.0);
    RiFormat(_xRes, _yRes,1.0);
    //RiFormat(4 * _xRes, 4 * _yRes,1.0);

    RiDisplayChannel("float _area", RI_NULL);
    RiDisplayChannel("color _radiance_t", RI_NULL);

    //RiShadingRate(RI_INFINITY);
    //RiShadingRate(100);

    RiWorldBegin();
      RtInt zero = 0;
      enableRaytracing();
      RiAttribute("cull", "hidden", (RtPointer)&zero, RI_NULL);
      RiAttribute("cull", "backfacing", (RtPointer)&zero, RI_NULL);
      RiAttribute("dice", "rasterorient", (RtPointer)&zero, RI_NULL);

      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiTranslate(-0.5, -0.5, -0.5);

        enableRaytracedLight();

        string ptcFilename = filename + string(".ptc");
        sprintf(buffer, "%s", ptcFilename.c_str());
        RtString pointCloud = buffer;
        RtString channels = "_area,_radiance_t";
        RiAttributeBegin();
        RiDeclare("filename", "uniform string");
        RiDeclare("displaychannels", "uniform string");
        RiSurface("./shaders/bake_radiance_t", 
                  (RtToken)"filename", (RtPointer)&pointCloud,
                  (RtToken)"displaychannels", (RtPointer)&channels,
                  RI_NULL);
        RiSides(1);

#if USING_TEAPOT        
        cout << " TEAPOT TRANSFORM ENABLED" << endl;
        RiTranslate(0.25, 0.65, 0.0);
        RiRotate(90,1,0,0);
        RiScale(0.1, 0.1, 0.1);
        RtColor superWhite = {10,10,10};
        RiColor(superWhite);
        RiGeometry("teapot", RI_NULL);
#else
        tetMesh->drawHeadEmbeddingToRenderMan(true);
#endif
        RiAttributeEnd();
      RiTransformEnd();
    RiWorldEnd();
  RiEnd();
  string prman = string("prman -progress ") + ribName;
  system(prman.c_str());
  string rmRib = string("rm ") + ribName;
  system(rmRib.c_str());
  //string rmTiff = string("rm ") + bakedFilename;
  //system(rmTiff.c_str());
  //string gzip = string("gzip -f ") + string(filename);
  //system(gzip.c_str());
  
  // Run diffusion on the point cloud
  string blurredPtc = filename + string(".sss.ptc");
  string ptfilter = string("ptfilter -filter ssdiffusion -partial 1 ");
  ptfilter += ptcFilename + string(" ");
  ptfilter += blurredPtc;
  cout << " Running " << ptfilter << endl;
  system(ptfilter.c_str());
  cout << " Point cloud filename: " << blurredPtc << endl;

  // remove the unblurred point cloud
  string rmPtc = string("rm ") + ptcFilename;
  system(rmPtc.c_str());
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// render a frame of the head
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawHeadWithBee(string filename, TET_MESH* tetMesh, OBJ* bee)
{
#if USING_RENDERMAN
  //RtString spath[] = {"./::/hydra_nfs/tedkim/svnrepos/shaders/"};
  //RiOption("searchpath", "shader", (RtPointer)spath, RI_NULL);

  // do a fast preview?
  bool fastRender = true;

  // bake a point cloud for the head
  TIMER headTimer;
  if (!fastRender)
    bakeHeadSubsurface(filename, tetMesh);
  cout << "Bake took " << headTimer.timing() << " s (" << headTimer.timing() / 60 << " min) " << endl;
  cout << " BAKING IS DISABLED" << endl;

  cout << " Rendering image " << filename.c_str() << endl;

  // view settings
  float tiltCameraDown = 11;
  float zoomOut = 0.6;
  float cameraY = 0.55;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);

  string ribName = filename + string(".rib");
  
  // zoomed in version
  TIMER renderTimer;
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);
    RtFloat fov = fieldOfView;

    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    RiFormat(_xRes, _yRes, 1.0);
    //RiFormat(320, 240,1.0);
    //RiFormat(160, 120,1.0);

    RiWorldBegin();
      RtColor white = {1.0, 1.0, 1.0};
      RiImager("background", "background", (RtPointer)white, RI_NULL);

      //enableRaytracedLight();
      //enableRaytracing();
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiTranslate(-0.5, -0.5, -0.5);
        enableRaytracedLight();
        enableRaytracing();

        drawObjWithMaterials(bee);

        RiColor(white);
        if (fastRender)
          tetMesh->drawHeadEmbeddingToRenderMan(false, false, filename);
        else
          tetMesh->drawHeadEmbeddingToRenderMan(false, true, filename);

      RiTransformEnd();
    RiWorldEnd();
  RiEnd();

  string prman = string("prman -progress ") + ribName;
  system(prman.c_str());
  cout << "Render took " << renderTimer.timing() << " s (" << renderTimer.timing() / 60 << " min) " << endl;

  TIMER cleanupTimer;
  string rmRib = string("rm ") + ribName;
  system(rmRib.c_str());
  string rmPtc = string("rm ") + filename + string(".sss.ptc");
  system(rmPtc.c_str());
  string rmBaketif = string("rm ") + filename + string(".bake.tif");
  system(rmBaketif.c_str());
  cout << "Cleanup took " << cleanupTimer.timing() << " s (" << cleanupTimer.timing() / 60 << " min) " << endl;

  //string gzip = string("gzip -f ") + string(filename);
  //system(gzip.c_str());
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// render a frame of the head
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawHeadEmbedding(string filename, TET_MESH* tetMesh)
{
#if USING_RENDERMAN
  //RtString spath[] = {"./::/hydra_nfs/tedkim/svnrepos/shaders/"};
  //RiOption("searchpath", "shader", (RtPointer)spath, RI_NULL);

  // do a fast preview?
  bool fastRender = true;

  // bake a point cloud for the head
  TIMER headTimer;
  if (!fastRender)
    bakeHeadSubsurface(filename, tetMesh);
  cout << "Bake took " << headTimer.timing() << " s (" << headTimer.timing() / 60 << " min) " << endl;
  cout << " BAKING IS DISABLED" << endl;

  cout << " Rendering image " << filename.c_str() << endl;

  // view settings
  float tiltCameraDown = 11;
  float zoomOut = 0.6;
  float cameraY = 0.55;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);

  string ribName = filename + string(".rib");
  
  // zoomed in version
  TIMER renderTimer;
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);
    RtFloat fov = fieldOfView;

    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
    
    cout << " Rendering at res: " << _xRes << " " << _yRes << endl; 
    RiFormat(_xRes, _yRes, 1.0);
    //RiFormat(320, 240,1.0);
    //RiFormat(160, 120,1.0);

    RiWorldBegin();
      RtColor white = {1.0, 1.0, 1.0};
      //RiImager("background", "background", (RtPointer)white, RI_NULL);

      //enableRaytracedLight();
      //enableRaytracing();
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiTranslate(-0.5, -0.5, -0.5);
        enableRaytracedLight();
        enableRaytracing();

#if USING_TEAPOT        
        sprintf(buffer, "%s.sss.ptc", filename.c_str());
        RtString sssFile = buffer;
        RiDeclare("unitlength", "uniform float");
        //RtFloat unitlength = 0.001;
        RtFloat unitlength = 0.01;
        RiSurface("./shaders/render_ssdiffusion", 
                  (RtToken)"filename", (RtPointer)&sssFile, 
                  (RtToken)"unitlength", (RtPointer)&unitlength,
                  RI_NULL);
        cout << " TEAPOT TRANSFORM ENABLED" << endl;
        RiTranslate(0.25, 0.65, 0.0);
        RiRotate(90,1,0,0);
        RiScale(0.1, 0.1, 0.1);
        RtColor superWhite = {10,10,10};
        RiColor(superWhite);
        RiGeometry("teapot", RI_NULL);
#else
        if (fastRender)
          tetMesh->drawHeadEmbeddingToRenderMan(false, false, filename);
        else
          tetMesh->drawHeadEmbeddingToRenderMan(false, true, filename);
#endif
      RiTransformEnd();
    RiWorldEnd();
  RiEnd();

  string prman = string("prman -progress ") + ribName;
  system(prman.c_str());
  cout << "Render took " << renderTimer.timing() << " s (" << renderTimer.timing() / 60 << " min) " << endl;

  TIMER cleanupTimer;
  string rmRib = string("rm ") + ribName;
  system(rmRib.c_str());
  string rmPtc = string("rm ") + filename + string(".sss.ptc");
  system(rmPtc.c_str());
  string rmBaketif = string("rm ") + filename + string(".bake.tif");
  system(rmBaketif.c_str());
  cout << "Cleanup took " << cleanupTimer.timing() << " s (" << cleanupTimer.timing() / 60 << " min) " << endl;

  //string gzip = string("gzip -f ") + string(filename);
  //system(gzip.c_str());
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// render a frame of the Menger membrane
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawPendulum(string filename, TET_MESH* tetMesh, vector<VEC3F*>& points, bool sliced, bool embedded)
{
#if USING_RENDERMAN
  cout << " Rendering image " << filename.c_str() << endl;

  // view settings
  //float tiltCameraDown = 11;
  float tiltCameraDown = 8;
  float zoomOut = 1.1f;
  float cameraY = 0.55;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);
  string ribName = filename + string(".rib");
  
  // zoomed in version
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);

    RtFloat fov = fieldOfView;
    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    RiFormat(640, 480,1.0);
    enableExplicitLight();
    RiWorldBegin();
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiColor(_white);
        RiSurface("matte", RI_NULL);

        // TODO: DRAW TRIANGLES HERE
        RiTranslate(-0.5, -0.5, -0.5);
        if (sliced)
          tetMesh->drawZSliceToRenderMan();
        else if (embedded)
          tetMesh->drawHeadEmbeddingToRenderMan();
        else
          tetMesh->drawSurfaceToRenderMan();

        RiColor(_red);
        drawPoints(points, sliced);
      RiTransformEnd();
    RiWorldEnd();
  RiEnd();

  string prman = string("prman -progress ") + ribName;
  system(prman.c_str());
  string rmRib = string("rm ") + ribName;
  system(rmRib.c_str());
  //string gzip = string("gzip -f ") + string(filename);
  //system(gzip.c_str());
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// render a frame of the Menger membrane
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawSkinningBegin(string filename)
{
#if USING_RENDERMAN
  // view settings
  //float tiltCameraDown = 11;
  float tiltCameraDown = 8;
  //float zoomOut = 0.25f;
  float zoomOut = 2.0;
  //float zoomOut = 4.0;
  float cameraY = 0.5;
  float cameraX = 0.75;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);
  string ribName = filename + string(".rib");
  
  // zoomed in version
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);

    RtFloat fov = fieldOfView;
    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    RiFormat(_xRes, _yRes, 1.0);
    enableExplicitLight();
    RiWorldBegin();
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiTranslate(-0.5, -0.5, -0.5);
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// render a frame of the Menger membrane
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawBeamFinalBegin(string filename, TET_MESH* tetMesh)
{
#if USING_RENDERMAN
  // view settings
  //float tiltCameraDown = 11;
  float tiltCameraDown = 8;
  //float zoomOut = 0.25f;
  float zoomOut = 1.1f;
  float cameraY = 0.55;
  float cameraX = 0.4;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);
  string ribName = filename + string(".rib");
  
  // zoomed in version
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);

    RtFloat fov = fieldOfView;
    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    RiFormat(_xRes, _yRes, 1.0);
    //enableExplicitLight();
    /*
    RtFloat intensity1 = 10.0f;
    RtFloat from1[] = {3.0f, 3.0f, -3.0f};
    RiLightSource("pointlight_rts", 
                "intensity", (RtPointer)&intensity1,
                "from", (RtPointer)from1,
                "lightcolor", (RtPointer)_white,
                RI_NULL);
                */
    RiWorldBegin();
      RtColor white = {1.0, 1.0, 1.0};
      //RiImager("background", "background", (RtPointer)white, RI_NULL);
      RiTransformBegin();

        RiRotate(-45, 0, 1,0);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiTranslate(-0.5, -0.5, -0.5);
        /*
        {
          string ptcFilename = filename + string(".sss.ptc");
          sprintf(buffer, "%s", ptcFilename.c_str());
          RtString sssFile = buffer;
          RtFloat unitlength = 0.001;
          RtFloat Ks = 1.0;

          RiDeclare("unitlength", "uniform float");
          RiDeclare("Ks", "uniform float");
          RiDeclare("ptcname", "uniform string");
          RiSurface("./shaders/livershader", 
              (RtToken)"ptcname", (RtPointer)&sssFile,
              (RtToken)"unitlength", (RtPointer)&unitlength,
              (RtToken)"Ks", (RtPointer)&Ks,
              RI_NULL);
        }
        */
        enableRaytracedLiverLight();
        enableRaytracing();
  // turn on ray tracing and set the bias
  RtInt trace = 1;
  RiAttribute("visibility", "trace", (RtPointer)&trace, RI_NULL);
  RtFloat bias = 0.001;
  RiAttribute("trace", "bias", (RtPointer)&bias, RI_NULL);

  RtInt transmission = 1; 
  //RtString transmissionhitmode = "primitive";
  // this kills subsurface
  RiAttribute("visibility", "int transmission", (RtPointer)&transmission, RI_NULL);
  //RiAttribute("shade", "transmissionhitmode", (RtPointer)&transmissionhitmode, RI_NULL);

  RtInt diffuse = 1;
  RiAttribute("visibility", "int diffuse", (RtPointer)&diffuse, RI_NULL);
  // turns on specular reflections
  RtInt specular = 1;
  RiAttribute("visibility", "int specular",(RtPointer)&specular, RI_NULL); 
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// render a frame of the Menger membrane
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawPendulumFinalBegin(string filename, TET_MESH* tetMesh)
{
#if USING_RENDERMAN
  // view settings
  //float tiltCameraDown = 11;
  float tiltCameraDown = 8;
  //float zoomOut = 0.25f;
  float zoomOut = 1.1f;
  float cameraY = 0.55;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);
  string ribName = filename + string(".rib");
  
  // zoomed in version
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);

    RtFloat fov = fieldOfView;
    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    RiFormat(_xRes, _yRes, 1.0);
    //enableExplicitLight();
    /*
    RtFloat intensity1 = 10.0f;
    RtFloat from1[] = {3.0f, 3.0f, -3.0f};
    RiLightSource("pointlight_rts", 
                "intensity", (RtPointer)&intensity1,
                "from", (RtPointer)from1,
                "lightcolor", (RtPointer)_white,
                RI_NULL);
                */
    RiWorldBegin();
      RtColor white = {1.0, 1.0, 1.0};
      //RiImager("background", "background", (RtPointer)white, RI_NULL);
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiTranslate(-0.5, -0.5, -0.5);
        /*
        {
          string ptcFilename = filename + string(".sss.ptc");
          sprintf(buffer, "%s", ptcFilename.c_str());
          RtString sssFile = buffer;
          RtFloat unitlength = 0.001;
          RtFloat Ks = 1.0;

          RiDeclare("unitlength", "uniform float");
          RiDeclare("Ks", "uniform float");
          RiDeclare("ptcname", "uniform string");
          RiSurface("./shaders/livershader", 
              (RtToken)"ptcname", (RtPointer)&sssFile,
              (RtToken)"unitlength", (RtPointer)&unitlength,
              (RtToken)"Ks", (RtPointer)&Ks,
              RI_NULL);
        }
        */
        enableRaytracedLiverLight();
        enableRaytracing();
  // turn on ray tracing and set the bias
  RtInt trace = 1;
  RiAttribute("visibility", "trace", (RtPointer)&trace, RI_NULL);
  RtFloat bias = 0.001;
  RiAttribute("trace", "bias", (RtPointer)&bias, RI_NULL);

  RtInt transmission = 1; 
  //RtString transmissionhitmode = "primitive";
  // this kills subsurface
  RiAttribute("visibility", "int transmission", (RtPointer)&transmission, RI_NULL);
  //RiAttribute("shade", "transmissionhitmode", (RtPointer)&transmissionhitmode, RI_NULL);

  RtInt diffuse = 1;
  RiAttribute("visibility", "int diffuse", (RtPointer)&diffuse, RI_NULL);
  // turns on specular reflections
  RtInt specular = 1;
  RiAttribute("visibility", "int specular",(RtPointer)&specular, RI_NULL); 

#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// render a frame of the Menger membrane
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawPendulumBegin(string filename)
{
#if USING_RENDERMAN
  // view settings
  //float tiltCameraDown = 11;
  float tiltCameraDown = 8;
  //float zoomOut = 0.25f;
  float zoomOut = 1.1f;
  float cameraY = 0.55;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);
  string ribName = filename + string(".rib");
  
  // zoomed in version
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);

    RtFloat fov = fieldOfView;
    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    RiFormat(_xRes, _yRes, 1.0);
    enableExplicitLight();
    RiWorldBegin();
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiTranslate(-0.5, -0.5, -0.5);
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// begin drawing the dragon
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawDragonBegin(string filename)
{
#if USING_RENDERMAN
  // view settings
  float tiltCameraDown = 8;
  float zoomOut = 2.1f;
  float cameraY = 0.4;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);
  string ribName = filename + string(".rib");
  
  // zoomed in version
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);

    RtFloat fov = fieldOfView;
    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    RiFormat(640, 480,1.0);
    enableExplicitLight();
    RiWorldBegin();
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiColor(_white);
        RiSurface("matte", RI_NULL);

        RiTranslate(-0.5, -0.5, -0.5);
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// finish drawing
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawEnd(string filename)
{
#if USING_RENDERMAN
  string ribName = filename + string(".rib");

      RiTransformEnd();
    RiWorldEnd();
  RiEnd();

  cout << " Drawing frame " << filename.c_str() << endl;
  string prman = string("prman -progress ") + ribName;
  system(prman.c_str());
  string rmRib = string("rm ") + ribName;
  system(rmRib.c_str());
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// render a frame of the Menger membrane
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawDragon(string filename, TET_MESH* tetMesh, vector<VEC3F*>& points, bool sliced, bool embedded)
{
#if USING_RENDERMAN
  cout << " Rendering image " << filename.c_str() << endl;

  // view settings
  //float tiltCameraDown = 11;
  float tiltCameraDown = 8;
  float zoomOut = 2.1f;
  //float cameraY = 0.55;
  float cameraY = 0.4;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);
  string ribName = filename + string(".rib");
  
  // zoomed in version
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);

    RtFloat fov = fieldOfView;
    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    RiFormat(640, 480,1.0);
    enableExplicitLight();
    RiWorldBegin();
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiColor(_white);
        RiSurface("matte", RI_NULL);

        // TODO: DRAW TRIANGLES HERE
        RiTranslate(-0.5, -0.5, -0.5);
        tetMesh->drawSurfaceToRenderMan();

        RiColor(_red);
        drawPoints(points, sliced);
      RiTransformEnd();
    RiWorldEnd();
  RiEnd();

  string prman = string("prman -progress ") + ribName;
  system(prman.c_str());
  string rmRib = string("rm ") + ribName;
  system(rmRib.c_str());
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// render a frame of the Menger membrane
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawHeadBegin(string filename)
{
#if USING_RENDERMAN
  // view settings
  float tiltCameraDown = 11;
  float zoomOut = 0.6;
  //float zoomOut = 2;
  float cameraY = 0.55;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);

  string ribName = filename + string(".rib");
  
  // zoomed in version
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);
    RtFloat fov = fieldOfView;

    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    RiFormat(640, 480,1.0);
    //RiFormat(320, 240,1.0);
    //RiFormat(160, 120,1.0);
    enableExplicitLight();
    RiWorldBegin();
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiColor(_white);
        RiSurface("matte", RI_NULL);
        
        RiTranslate(-0.5, -0.5, -0.5);
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// render a frame of the Menger membrane
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawTetMesh(string filename, TET_MESH* tetMesh, vector<VEC3F*>& points, bool sliced, bool embedded)
{
#if USING_RENDERMAN
  cout << " Rendering image " << filename.c_str() << endl;

  // view settings
  //float spinModelY = 22.5;
  //float spinModelY = 90.0;
  //float spinModelZ = 180;
  float tiltCameraDown = 11;
  float zoomOut = 0.6;
  float cameraY = 0.55;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);

  string ribName = filename + string(".rib");
  
  // zoomed in version
  //RiBegin("launch:prman? -t -ctrl");
  //RiBegin("launch:prman?");
  //RiBegin(RI_NULL);
  //RiBegin("output.rib");
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);

    //RiClipping(100.0 * RI_EPSILON, RI_INFINITY);
    RtFloat fov = fieldOfView;

    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    //RtInt progress = 1;
    //RtInt progress = 0;
    //RiOption("statistics", "progress", (RtPointer)&progress, RI_NULL);
    //RiFormat(640, 480,1.0);
    RiFormat(1280 , 960 ,1.0);

    enableExplicitLight();

    RiWorldBegin();
      RtColor white = {1.0, 1.0, 1.0};
      //RiImager("background", "background", (RtPointer)white, RI_NULL);
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        //RiRotate(spinModelY, 0, 1,0);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        //RiRotate(spinModelZ, 0, 0,1);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiColor(_white);
        RiSurface("matte", RI_NULL);

        // TODO: DRAW TRIANGLES HERE
        RiTranslate(-0.5, -0.5, -0.5);
        if (sliced)
          tetMesh->drawZSliceToRenderMan();
        else if (embedded)
          tetMesh->drawHeadEmbeddingToRenderMan();
        else
        {
          tetMesh->drawSurfaceToRenderMan();
          RtColor black = {0,0,0};
          RiColor(black);
          tetMesh->drawOutlinesToRenderMan();
        }
        //tetMesh->drawExhaustiveToRenderMan();

        //drawAxes();

        RiColor(_red);
        drawPoints(points, sliced);

        //RiSphere(1.0f, -1.0f, 1.0f, 360.0f, RI_NULL);
      RiTransformEnd();
    RiWorldEnd();
  RiEnd();

  //string prman = string("prman -progress ") + ribName;
  string prman = string("prman ") + ribName;
  system(prman.c_str());
  string rmRib = string("rm ") + ribName;
  system(rmRib.c_str());
  //string gzip = string("gzip -f ") + string(filename);
  //system(gzip.c_str());
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// render a frame of the Menger membrane
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawObjs(string filename, vector<OBJ*>& objs)
{
#if USING_RENDERMAN
  cout << " Rendering image " << filename.c_str() << endl;

  // view settings
  //float spinModelY = 22.5;
  //float spinModelY = 90.0;
  //float spinModelZ = 180;
  float tiltCameraDown = 11;
  float zoomOut = 0.6;
  float cameraY = 0.55;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);

  string ribName = filename + string(".rib");
  
  // zoomed in version
  //RiBegin("launch:prman? -t -ctrl");
  //RiBegin("launch:prman?");
  //RiBegin(RI_NULL);
  //RiBegin("output.rib");
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);

    //RiClipping(100.0 * RI_EPSILON, RI_INFINITY);
    RtFloat fov = fieldOfView;

    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    //RtInt progress = 1;
    //RtInt progress = 0;
    //RiOption("statistics", "progress", (RtPointer)&progress, RI_NULL);
    RiFormat(640, 480,1.0);
    //RiFormat(1280 , 960 ,1.0);

    enableExplicitLight();

    RiWorldBegin();
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        //RiRotate(spinModelY, 0, 1,0);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        //RiRotate(spinModelZ, 0, 0,1);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiColor(_white);
        RiSurface("matte", RI_NULL);

        // TODO: DRAW TRIANGLES HERE
        RiTranslate(-0.5, -0.5, -0.5);

        for (int x = 0; x < objs.size(); x++)
          drawObj(objs[x]);

        //RiSphere(1.0f, -1.0f, 1.0f, 360.0f, RI_NULL);
      RiTransformEnd();
    RiWorldEnd();
  RiEnd();

  string prman = string("prman -progress ") + ribName;
  system(prman.c_str());
  string rmRib = string("rm ") + ribName;
  system(rmRib.c_str());
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// draw tree mesh surface
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawTreeSurface(string filename, TET_MESH* tetMesh)
{
#if USING_RENDERMAN
  cout << " Rendering image " << filename.c_str() << endl;

  // view settings
  float tiltCameraDown = 11;
  float zoomOut = 1.0;
  float cameraY = 0.55;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  // zoomed in version
  RiBegin("output.rib");
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);

    RtFloat fov = fieldOfView;

    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    RtInt progress = 0;
    RiFormat(640, 480,1.0);

    enableExplicitLight();

    RiWorldBegin();
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiColor(_white);
        RiSurface("matte", RI_NULL);

        // TODO: DRAW TRIANGLES HERE
        RiTranslate(-0.5, -0.5, -0.5);
        tetMesh->drawSurfaceToRenderMan();
      RiTransformEnd();
    RiWorldEnd();
  RiEnd();

  string prman = string("prman output.rib");
  system(prman.c_str());
  string rmRib = string("rm output.rib");
  system(rmRib.c_str());
#endif
}

#if 0
///////////////////////////////////////////////////////////////////////////////////////
// draw embedded tree mesh
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawEmbeddedTree(string filename, PARTITIONED_TET_MESH* tetMesh)
{
#if USING_RENDERMAN
  cout << " Rendering image " << filename.c_str() << endl;

  // view settings
  float tiltCameraDown = 11;
  float zoomOut = 1.0;
  float cameraY = 0.55;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  /*
  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);
  */

  // zoomed in version
  RiBegin("output.rib");
  //RiBegin(RI_NULL);
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);

    RtFloat fov = fieldOfView;

    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    //RtInt progress = 1;
    RtInt progress = 0;
    RiFormat(640, 480,1.0);

    enableExplicitLight();

    RiWorldBegin();
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiColor(_white);
        RiSurface("matte", RI_NULL);

        // TODO: DRAW TRIANGLES HERE
        RiTranslate(-0.5, -0.5, -0.5);
        tetMesh->drawEmbeddingToRenderMan();
      RiTransformEnd();
    RiWorldEnd();
  RiEnd();

  //string prman = string("prman -progress output.rib");
  string prman = string("prman output.rib");
  system(prman.c_str());
  string rmRib = string("rm output.rib");
  system(rmRib.c_str());
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// render a frame of the Menger membrane
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::debugPartitionedEmbedding(string filename, PARTITIONED_TET_MESH* tetMesh)
{
#if USING_RENDERMAN
  cout << " Rendering image " << filename.c_str() << endl;

  // view settings
  float tiltCameraDown = 11;
  float zoomOut = 1.0;
  float cameraY = 0.55;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);

  // zoomed in version
  RiBegin("output.rib");
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);

    RtFloat fov = fieldOfView;

    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    //RtInt progress = 1;
    RtInt progress = 0;
    RiFormat(640, 480,1.0);

    enableExplicitLight();

    RiWorldBegin();
      RiTransformBegin();
        RiTranslate(-0.125, 0.0, 0.0);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiColor(_white);
        RiSurface("matte", RI_NULL);

        // TODO: DRAW TRIANGLES HERE
        RiTranslate(-0.5, -0.5, -0.5);
        tetMesh->drawSurfaceToRenderMan();
      RiTransformEnd();

      RiTransformBegin();
        RiTranslate(0.125, 0.0, 0.0);
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiColor(_white);
        RiSurface("matte", RI_NULL);

        // TODO: DRAW TRIANGLES HERE
        RiTranslate(-0.5, -0.5, -0.5);
        tetMesh->drawEmbeddingToRenderMan();
      RiTransformEnd();
    RiWorldEnd();
  RiEnd();

  //string prman = string("prman -progress output.rib");
  string prman = string("prman output.rib");
  system(prman.c_str());
  string rmRib = string("rm output.rib");
  system(rmRib.c_str());
#endif
}
#endif

///////////////////////////////////////////////////////////////////////////////////////
// render a frame of the Menger membrane
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::debugEmbedding(string filename, TET_MESH* tetMesh)
{
#if USING_RENDERMAN
  cout << " Rendering image " << filename.c_str() << endl;

  // view settings
  //float spinModelY = 22.5;
  //float spinModelY = 90.0;
  //float spinModelZ = 180;
  float tiltCameraDown = 11;
  float zoomOut = 1.0;
  //float zoomOut = 20.0;
  float cameraY = 0.55;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);

  // zoomed in version
  //RiBegin(RI_NULL);
  RiBegin("output.rib");
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);

    //RiClipping(RI_EPSILON, RI_INFINITY);
    RtFloat fov = fieldOfView;

    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    //RtInt progress = 1;
    RtInt progress = 0;
    //RiOption("statistics", "progress", (RtPointer)&progress, RI_NULL);
    RiFormat(640, 480,1.0);
    //RiFormat(1280 , 960 ,1.0);

    enableExplicitLight();

    RiWorldBegin();
      RiTransformBegin();
        RiTranslate(-0.125, 0.0, 0.0);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        //RiRotate(spinModelY, 0, 1,0);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        //RiRotate(spinModelZ, 0, 0,1);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiColor(_white);
        RiSurface("matte", RI_NULL);

        // TODO: DRAW TRIANGLES HERE
        RiTranslate(-0.5, -0.5, -0.5);
        tetMesh->drawSurfaceToRenderMan();
      RiTransformEnd();

      RiTransformBegin();
        RiTranslate(0.125, 0.0, 0.0);
        RiTranslate(0.5, 0.5, 0.5);
        //RiRotate(spinModelY, 0, 1,0);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        //RiRotate(spinModelZ, 0, 0,1);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiColor(_white);
        RiSurface("matte", RI_NULL);

        // TODO: DRAW TRIANGLES HERE
        RiTranslate(-0.5, -0.5, -0.5);
        //tetMesh->drawHeadEmbeddingToRenderMan();
        tetMesh->drawEmbeddingToRenderMan();
      RiTransformEnd();
    RiWorldEnd();
  RiEnd();

  string prman = string("prman -progress output.rib");
  system(prman.c_str());
  string rmRib = string("rm output.rib");
  system(rmRib.c_str());
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawAxes(Real size)
{
#if USING_RENDERMAN
  RtString linear = "linear";
  RtString nonperiodic = "nonperiodic";
  RtInt points = 2;

  RtInt nvertices[1];
  nvertices[0] = 2;

  RtString constantwidth = "constantwidth";
  RiDeclare("constantwidth", "uniform float");
  RtFloat width = 0.1;

  RtPoint P[2];
  P[0][0] = 0.0;
  P[0][1] = 0.0;
  P[0][2] = 0.0;

  RtColor superRed = {10,0,0};
  RiColor(superRed);
  P[1][0] = size;
  P[1][1] = 0.0;
  P[1][2] = 0.0;
  RiCurves("linear", 1, nvertices, "nonperiodic", RI_P, &P, constantwidth, &width, RI_NULL);

  RtColor superGreen = {0,10,0};
  RiColor(superGreen);
  P[1][0] = 0.0;
  P[1][1] = size;
  P[1][2] = 0.0;
  RiCurves("linear", 1, nvertices, "nonperiodic", RI_P, &P, constantwidth, &width, RI_NULL);

  RtColor superBlue = {0,0,10};
  RiColor(superBlue);
  P[1][0] = 0.0;
  P[1][1] = 0.0;
  P[1][2] = size;
  RiCurves("linear", 1, nvertices, "nonperiodic", RI_P, &P, constantwidth, &width, RI_NULL);
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
// Note that these materials are only intended for the bee model
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawObjWithMaterials(OBJ* obj)
{
#if USING_RENDERMAN
  int x;
  vector<VEC3F>& vertices = obj->vertices;
  vector<VEC3F>& normals = obj->normals;
  vector<OBJ::Face>& faces = obj->faces;

  if (normals.size() != vertices.size())
  {
    cout << " Normals have not been computed! Call ComputeVertexNormals! " << endl;
    cout << endl;
  }

#ifndef _WIN32
  vector<int> materialStarts = obj->materialStarts();
  materialStarts.push_back(faces.size());

  RtColor colors[7];
  RtColor yellow = {1,1,0};
  RtColor black = {0,0,0};
  RtColor gray = {0.5, 0.5, 0.5};

  for (int i = 0; i < 7; i++)
    colors[i][0] = colors[i][1] = colors[i][2] = 0;

   // 1 2 3 2 1 4
  colors[0][0] = 0;
  colors[0][1] = 0;
  colors[0][2] = 0;

  colors[1][0] = 1;
  colors[1][1] = 1;
  colors[1][2] = 0;

  colors[2][0] = 0;
  colors[2][1] = 0;
  colors[2][2] = 0;

  colors[3][0] = 0.5;
  colors[3][1] = 0.5;
  colors[3][2] = 0.5;

  colors[4][0] = 0.5;
  colors[4][1] = 0.5;
  colors[4][2] = 0.5;

  colors[5][0] = 1;
  colors[5][1] = 1;
  colors[5][2] = 0;

  for (int i = 0; i < materialStarts.size() - 1; i++)
  {
    int nfaces = materialStarts[i + 1] - materialStarts[i];
    int start = materialStarts[i];
    cout << " Material start: " << start << endl;

    // init point locations and normals
    RtPoint* P = new RtPoint[vertices.size()];
    RtPoint* N = new RtPoint[normals.size()];
    for (x = 0; x < vertices.size(); x++)
    {
      P[x][0] = vertices[x][0];
      P[x][1] = vertices[x][1];
      P[x][2] = vertices[x][2];

      N[x][0] = normals[x][0];
      N[x][1] = normals[x][1];
      N[x][2] = normals[x][2];

      if (isnan(N[x][0]) || isnan(N[x][1]) || isnan(N[x][2]))
        N[x][0] = N[x][1] = N[x][2] = 0.0;
    }

    // init all to triangles
    RtInt* nvertices = new RtInt[nfaces];
    for (x = 0; x < nfaces; x++)
      nvertices[x] = 3;

    // init faces
    RtInt* rifaces = new RtInt[3 * nfaces];
    for (x = 0; x < nfaces; x++)
    {
      OBJ::Face face = faces[start + x];
      rifaces[x * 3]     = face.vertices[0];
      rifaces[x * 3 + 1] = face.vertices[1];
      rifaces[x * 3 + 2] = face.vertices[2];
    }
  
    RiColor(colors[i]);

    RiPointsPolygons(nfaces, nvertices, rifaces, RI_P, P, RI_N, (RtPointer)N, RI_NULL);
    delete[] rifaces;

    delete[] nvertices;
    delete[] P;
    delete[] N;
  }
#endif
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawObj(OBJ* obj)
{
#if USING_RENDERMAN
  int x;
  vector<VEC3F>& vertices = obj->vertices;
  vector<VEC3F>& normals = obj->normals;
  vector<OBJ::Face>& faces = obj->faces;
  int nfaces = faces.size();

  if (normals.size() != vertices.size())
  {
    cout << " Normals have not been computed! Call ComputeVertexNormals! " << endl;
    cout << endl;
  }

#ifndef _WIN32
  // init point locations and normals
  RtPoint* P = new RtPoint[vertices.size()];
  RtPoint* N = new RtPoint[normals.size()];
  for (x = 0; x < vertices.size(); x++)
  {
    P[x][0] = vertices[x][0];
    P[x][1] = vertices[x][1];
    P[x][2] = vertices[x][2];

    N[x][0] = normals[x][0];
    N[x][1] = normals[x][1];
    N[x][2] = normals[x][2];

    if (isnan(N[x][0]) || isnan(N[x][1]) || isnan(N[x][2]))
      N[x][0] = N[x][1] = N[x][2] = 0.0;
  }

  // init all to triangles
  RtInt* nvertices = new RtInt[nfaces];
  for (x = 0; x < nfaces; x++)
    nvertices[x] = 3;

  // init faces
  RtInt* rifaces = new RtInt[3 * nfaces];
  for (x = 0; x < nfaces; x++)
  {
    OBJ::Face face = faces[x];
    rifaces[x * 3]     = face.vertices[0];
    rifaces[x * 3 + 1] = face.vertices[1];
    rifaces[x * 3 + 2] = face.vertices[2];
  }
  
  RiPointsPolygons(nfaces, nvertices, rifaces, RI_P, P, RI_N, (RtPointer)N, RI_NULL);
  delete[] rifaces;

  delete[] nvertices;
  delete[] P;
  delete[] N;
#endif
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawCompression(string filename, TET_MESH* tetMesh, vector<VEC3F*>& points)
{
#if USING_RENDERMAN
  //cout << " BAKING HAS BEEN DISABLED" << endl;
  if (_qualityRender)
    bakeLiverSubsurface(filename, tetMesh);

  cout << " Rendering image " << filename.c_str() << endl;

  // view settings
  float tiltCameraDown = 0;
  float zoomOut = 1.5;
  float cameraY = 0.5;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);

  string ribName = filename + string(".rib");
  
  // zoomed in version
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);
    RtFloat fov = fieldOfView;

    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    RiFormat(_xRes, _yRes ,1.0);
    //RiFormat(640, 480,1.0);
    //
    if (!_qualityRender)
      enableExplicitLight();

    RiWorldBegin();
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiColor(_white);
        string ptcFilename = filename + string(".sss.ptc");
        sprintf(buffer, "%s", ptcFilename.c_str());
        RtString sssFile = buffer;
        RtFloat unitlength = 0.001;
        RtFloat Ks = 1.0;

        if (_qualityRender)
        {
          RiDeclare("unitlength", "uniform float");
          RiDeclare("Ks", "uniform float");
          RiDeclare("ptcname", "uniform string");
          RiSurface("./shaders/livershader", 
              (RtToken)"ptcname", (RtPointer)&sssFile,
              (RtToken)"unitlength", (RtPointer)&unitlength,
              (RtToken)"Ks", (RtPointer)&Ks,
              RI_NULL);
        }
        else
          RiSurface("matte", RI_NULL);

        if (_qualityRender)
          RiColor(_liverRed);
        else
          RiColor(_red);

        RiTranslate(-0.5, -0.5, -0.5);
        if (_qualityRender)
        {
          enableRaytracedLiverLight();
          enableRaytracing();
        }
        // TODO: DRAW TRIANGLES HERE
        //tetMesh->drawEmbeddingToRenderMan();
        tetMesh->drawSurfaceToRenderMan();
        drawPoints(points, false);

        if (_qualityRender)
        {
          RtString envmap = "./shaders/paris-cropped-blurred-pixar.exr";
          RiDeclare("envmap", "uniform string");
          RiSurface("./shaders/raymetal", 
                    "envmap", &envmap, RI_NULL);
        }
      RiTransformEnd();
    RiWorldEnd();
  RiEnd();

  string prman = string("prman -progress ") + ribName;
  system(prman.c_str());
  string rmRib = string("rm ") + ribName;
  system(rmRib.c_str());
  if (_qualityRender)
  {
    string rmPtc = string("rm ") + filename + string(".sss.ptc");
    system(rmPtc.c_str());
  }
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawLiver(string filename, TET_MESH* tetMesh, VEC3F translation, MATRIX3 rotation, bool grabbed, OBJ& tray, OBJ* openClamp, OBJ* closedClamp, VEC3F* secondClamp)
{
#if USING_RENDERMAN
  //cout << " BAKING HAS BEEN DISABLED" << endl;
  if (_qualityRender)
    bakeLiverSubsurface(filename, tetMesh);

  cout << " Rendering image " << filename.c_str() << endl;

  // view settings
  float tiltCameraDown = 0;
  float zoomOut = 1.0;
  float cameraY = 0.5;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);

  string ribName = filename + string(".rib");
  
  // zoomed in version
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);
    RtFloat fov = fieldOfView;

    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    RiFormat(_xRes, _yRes ,1.0);
    //RiFormat(640, 480,1.0);
    //
    if (!_qualityRender)
      enableExplicitLight();

    RiWorldBegin();
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiColor(_white);
        string ptcFilename = filename + string(".sss.ptc");
        sprintf(buffer, "%s", ptcFilename.c_str());
        RtString sssFile = buffer;
        RtFloat unitlength = 0.001;
        RtFloat Ks = 1.0;

        if (_qualityRender)
        {
          RiDeclare("unitlength", "uniform float");
          RiDeclare("Ks", "uniform float");
          RiDeclare("ptcname", "uniform string");
          RiSurface("./shaders/livershader", 
              (RtToken)"ptcname", (RtPointer)&sssFile,
              (RtToken)"unitlength", (RtPointer)&unitlength,
              (RtToken)"Ks", (RtPointer)&Ks,
              RI_NULL);
        }
        else
          RiSurface("matte", RI_NULL);

        if (_qualityRender)
          RiColor(_liverRed);
        else
          RiColor(_red);

        RiTranslate(-0.5, -0.5, -0.5);
        if (_qualityRender)
        {
          enableRaytracedLiverLight();
          enableRaytracing();
        }
        // TODO: DRAW TRIANGLES HERE
        tetMesh->drawEmbeddingToRenderMan();
        //tetMesh->drawSurfaceToRenderMan();

        /*
        RtColor blue;
        blue[0] = 0.0;
        blue[1] = 0.0;
        blue[2] = 1.0;
        RiColor(blue);
        tetMesh->drawSurfaceToRenderMan();
        */

        if (_qualityRender)
        {
          RtString envmap = "./shaders/paris-cropped-blurred-pixar.exr";
          //RtString envmap = "./shaders/paris-cropped-blurred-new-pixar.exr";
          RiDeclare("envmap", "uniform string");
          RiSurface("./shaders/raymetal", 
                    "envmap", &envmap, RI_NULL);
        }

        // draw the tray
        RiColor(_white);
        drawObj(&tray);

        // draw the other clamp
        RiTransformBegin();
          if (secondClamp != NULL)
          {
            RiTranslate((*secondClamp)[0], (*secondClamp)[1], (*secondClamp)[2]);
            RiRotate(45, 0,1,0);
            RiRotate(45, 1,0,0);
            RiRotate(180, 0,1,0);
            //RiTranslate(0.0, 0.025, 0.0);
            drawObj(closedClamp);
          }
        RiTransformEnd();

        if (grabbed)
          RiColor(_red);
        else
        {
          RtColor blue;
          blue[0] = 0.0;
          blue[1] = 0.0;
          blue[2] = 1.0;
          RiColor(blue);
        }

        // draw the forcep
        RiTransformBegin();
          RtMatrix transform;
          for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
              transform[i][j] = 0.0;

          transform[3][0] = translation[0];
          transform[3][1] = translation[1];
          transform[3][2] = translation[2];
          transform[3][3] = 1.0;

          // OpenGL and RenderMan use different handed matrices
          transform[0][0] = rotation(0,0);
          transform[0][1] = rotation(0,1);
          transform[0][2] = -rotation(0,2);
          transform[0][3] = 0.0;
          transform[1][0] = rotation(1,0);
          transform[1][1] = rotation(1,1);
          transform[1][2] = -rotation(1,2);
          transform[1][3] = 0.0;
          transform[2][0] = -rotation(2,0);
          transform[2][1] = -rotation(2,1);
          transform[2][2] = rotation(2,2);
          transform[2][3] = 0.0;

          // Must concat or it will blow the camera matrix away!
          RiConcatTransform(transform);
          RiRotate(90, 0,0,1);

          //RiTranslate(0.02, 0.0, 0.0);
          //RiTranslate(0.0, 0.025, 0.0);
          // been using this one
          //RiTranslate(-0.025, 0.0, 0.0);
          if (grabbed)
            drawObj(closedClamp);
          else
            drawObj(openClamp);

        RiTransformEnd();
      RiTransformEnd();
    RiWorldEnd();
  RiEnd();

  string prman = string("prman -progress ") + ribName;
  system(prman.c_str());
  string rmRib = string("rm ") + ribName;
  system(rmRib.c_str());
  if (_qualityRender)
  {
    string rmPtc = string("rm ") + filename + string(".sss.ptc");
    system(rmPtc.c_str());
  }
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void RENDER_MAN::drawRaytraced(string filename, TET_MESH* tetMesh)
{
#if USING_RENDERMAN
  // view settings
  //float tiltCameraDown = 11;
  float tiltCameraDown = 8;
  float zoomOut = 0.75f;
  float cameraY = 0.55;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);
  string ribName = filename + string(".rib");
  
  // zoomed in version
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);

    RtFloat fov = fieldOfView;
    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    RiFormat(_xRes, _yRes, 1.0);
    RiWorldBegin();
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiTranslate(-0.5, -0.5, -0.5);

        enableRaytracedDragonLight();
        enableRaytracing();

        RtColor blue;
        blue[0] = 0.0;
        blue[1] = 0.0;
        blue[2] = 0.0;
        RiColor(blue);
        RiSurface("matte", RI_NULL);
        //tetMesh->drawOutlinesToRenderMan();

        RiColor(_white);
        RiSurface("matte", RI_NULL);
        tetMesh->drawExhaustiveToRenderMan();
      RiTransformEnd();
    RiWorldEnd();
  RiEnd();

  cout << " Drawing frame " << filename.c_str() << endl;
  string prman = string("prman -progress ") + ribName;
  system(prman.c_str());
  string rmRib = string("rm ") + ribName;
  system(rmRib.c_str());
#endif
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//void RENDER_MAN::drawColored(string filename, TET_MESH* tetMesh)
void RENDER_MAN::drawColored(string filename, SUBSPACE_TET_MESH* tetMesh)
{
#if USING_RENDERMAN
  // view settings
  //float tiltCameraDown = 11;
  float tiltCameraDown = 8;
  float zoomOut = 0.75f;
  float cameraY = 0.55;
  float cameraX = 0.5;
  float fieldOfView = 20.0f;

  RtString str = "gzip";
  RiOption("rib", "compression", &str, RI_NULL);
  string ribName = filename + string(".rib");
  
  // zoomed in version
  RiBegin((char*)ribName.c_str());
    char buffer[256];
    sprintf(buffer, "%s", filename.c_str());
    RtString outputFile = buffer;
    RtString imageType = "file";
    RiDisplay((char*)outputFile, (char*)imageType, RI_RGB, RI_NULL);

    RtFloat fov = fieldOfView;
    RiPixelSamples(4,4);
    RiPixelFilter(RiCatmullRomFilter, 3,3);

    RiProjection(RI_PERSPECTIVE, "fov", &fov, RI_NULL);
    RtInt eyesplits = 10;
    RiOption("limits", "eyesplits", (RtPointer)&eyesplits, RI_NULL);
    RiTranslate(-cameraX, -cameraY, zoomOut);
    RiRotate(-tiltCameraDown, 1,0,0);
     
    RiFormat(_xRes, _yRes, 1.0);
    RiWorldBegin();
      RiTransformBegin();
        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinX, 1, 0,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinY, 0, 1,0);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiRotate(_spinZ, 0, 0,1);
        RiTranslate(-0.5, -0.5, -0.5);

        RiTranslate(0.5, 0.5, 0.5);
        RiScale(0.5, 0.5, 0.5);

        RiTranslate(-0.5, -0.5, -0.5);

        enableRaytracedDragonLight();
        enableRaytracing();

        RtColor lines;
        //lines[0] = 0.25;
        //lines[1] = 0.25;
        //lines[2] = 0.25;
        lines[0] = 1; lines[1] = 1; lines[2] = 1;
        lines[0] = 0; lines[1] = 0; lines[2] = 0;
        RiColor(lines);
        RiSurface("matte", RI_NULL);
        tetMesh->drawOutlinesToRenderMan();

        RiColor(_white);
        RiSurface("matte", RI_NULL);
        tetMesh->drawExhaustiveToRenderMan();
        //tetMesh->drawColoredToRenderMan();
        //tetMesh->drawKeyTetsToRenderMan();
        //tetMesh->drawSurfaceForcesToRenderMan();
        drawPoints(tetMesh->surfaceVertices(), false);
      RiTransformEnd();
    RiWorldEnd();
  RiEnd();

  cout << " Drawing frame " << filename.c_str() << endl;
  string prman = string("prman -progress ") + ribName;
  system(prman.c_str());
  string rmRib = string("rm ") + ribName;
  system(rmRib.c_str());
#endif
}

