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
// PIXIE.h: interface for the PIXIE class.
//
//////////////////////////////////////////////////////////////////////

#ifndef RENDER_MAN_H
#define RENDER_MAN_H

#include <SETTINGS.h>
#include <string>
#include <TET_MESH.h>
#include <PARTITIONED_TET_MESH.h>
#include <SUBSPACE_TET_MESH.h>
#include <OBJ.h>
#include <TIMER.h>

#ifdef USING_RENDERMAN
#include <ri.h>
#endif

using namespace std;

//////////////////////////////////////////////////////////////////////
// Interface to Pixie RenderMan implementation
//////////////////////////////////////////////////////////////////////
class RENDER_MAN {

public:
  RENDER_MAN();
  ~RENDER_MAN() {};

  void drawTetMesh(string filename, TET_MESH* tetMesh, vector<VEC3F*>& points, bool sliced = false, bool embedded = false);
  void drawHeadEmbedding(string filename, TET_MESH* tetMesh);
  void drawHeadWithBee(string filename, TET_MESH* tetMesh, OBJ* bee);
  void drawObjs(string filename, vector<OBJ*>& objs);
  void debugEmbedding(string filename, TET_MESH* tetMesh);
  //void debugPartitionedEmbedding(string filename, PARTITIONED_TET_MESH* tetMesh);
  void drawPendulum(string filename, TET_MESH* tetMesh, vector<VEC3F*>& points, bool sliced = false, bool embedded = false);
  void drawDragon(string filename, TET_MESH* tetMesh, vector<VEC3F*>& points, bool sliced = false, bool embedded = false);

  void drawLiver(string filename, TET_MESH* tetMesh, VEC3F translation, MATRIX3 rotation, bool grabbed, OBJ& tray, OBJ* openClamp, OBJ* closedClamp, VEC3F* secondClamp);
  void drawCompression(string filename, TET_MESH* tetMesh, vector<VEC3F*>& points);

  void drawRaytraced(string filename, TET_MESH* tetMesh);
  //void drawColored(string filename, TET_MESH* tetMesh);
  void drawColored(string filename, SUBSPACE_TET_MESH* tetMesh);

  void drawHeadBegin(string filename);
  void drawDragonBegin(string filename);
  void drawPendulumBegin(string filename);
  void drawPendulumFinalBegin(string filename, TET_MESH* tetMesh);
  void drawBeamFinalBegin(string filename, TET_MESH* tetMesh);
  void drawSkinningBegin(string filename);
  void drawEnd(string filename);

  void drawTreeSurface(string filename, TET_MESH* tetMesh);
  //void drawEmbeddedTree(string filename, PARTITIONED_TET_MESH* tetMesh);

  float& xSpinAngle() { return _spinX; };
  float& ySpinAngle() { return _spinY; };
  float& zSpinAngle() { return _spinZ; };

  bool& qualityRender() { return _qualityRender; };
  int& xRes() { return _xRes; };
  int& yRes() { return _yRes; };

private:
  void bakeHeadSubsurface(string filename, TET_MESH* tetMesh);
  void bakeLiverSubsurface(string filename, TET_MESH* tetMesh);

  bool  _qualityRender;
  bool  _commandLine;
  bool  _boxOn;
  bool  _blurOn;
#if USING_RENDERMAN
  RtColor _cyan;
  RtColor _red;
  RtColor _liverRed;
  RtColor _green;
  RtColor _blue;
  RtColor _white;
  RtString _envmap;
  RtString _raytrace;
#endif

  void enableRaytracing();
  void enableAmbientOcclusion();
  void enableExplicitLight();
  void enableRaytracedLight();
  void enableRaytracedLiverLight();
  void enableRaytracedDragonLight();
  void enableLiverSkydome();
  void enableEnvironmentLight();
  void enableSkySphere();
  void enableSolidColorSky();
  void enableDisplacementShader(float frequency, float amplitude);
  void disableDisplacementShader();
  void enableFloor();
public:
  void cube();

private:  
  void drawPoints(vector<VEC3F*>& points, bool sliced);
  void drawObj(OBJ* obj);
  void drawObjWithMaterials(OBJ* obj);

  void drawAxes(Real size = 1.0);

  float _spinX;
  float _spinY;
  float _spinZ;

  int _xRes;
  int _yRes;
};

#endif
