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
#include "ISO_STUFFER.h"

//////////////////////////////////////////////////////////////////////
// Draw all the tets as triangles
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawTets0()
{
  for (int x = 0; x < _tets0.size(); x++)
    _tets0[x].drawTriangles();
}

//////////////////////////////////////////////////////////////////////
// Draw all the tets as triangles
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawTets1()
{
  for (int x = 0; x < _tets1.size(); x++)
    _tets1[x].drawTriangles();
}

//////////////////////////////////////////////////////////////////////
// Draw all the tets as triangles
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawTets2()
{
  for (int x = 0; x < _tets2.size(); x++)
    _tets2[x].drawTriangles();
}

//////////////////////////////////////////////////////////////////////
// Draw all the tets as triangles
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawTets3()
{
  for (int x = 0; x < _tets3.size(); x++)
    _tets3[x].drawTriangles();
}

//////////////////////////////////////////////////////////////////////
// Draw all the tets as triangles
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawTets4()
{
  for (int x = 0; x < _tets4.size(); x++)
    _tets4[x].drawTriangles();
}

//////////////////////////////////////////////////////////////////////
// Draw all the tets as triangles
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawTets5()
{
  for (int x = 0; x < _tets5.size(); x++)
    _tets5[x].drawTriangles();
}

//////////////////////////////////////////////////////////////////////
// Generate all possible tets in the grid (for debugging)
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateType0Tets()
{
  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
      {
        int index = x + y * _xRes + z * _slabSize;

        VEC3F* vert000 = corner000(x,y,z);
        VEC3F* vert100 = corner100(x,y,z);
        VEC3F* vert010 = corner010(x,y,z);
        VEC3F* vert110 = corner110(x,y,z);
        VEC3F* vert001 = corner001(x,y,z);
        VEC3F* vert101 = corner101(x,y,z);
        VEC3F* vert011 = corner011(x,y,z);
        VEC3F* vert111 = corner111(x,y,z);

        VEC3F* center = centerVec(x,y,z);
        VEC3F* left = centerVec(x-1,y,z);
        VEC3F* near = centerVec(x,y,z-1);

        // type 0 tet
        // the first from the left in Fig. 2 of Isosurface Stuffing
        TET type0(center, left, vert011, vert010);
        _tets0.push_back(type0);

        // type 0 reflected about y
        TET type0Reflected(center, left, vert000, vert001);
        _tets0.push_back(type0Reflected);
      }
}

//////////////////////////////////////////////////////////////////////
// Generate all possible tets in the grid (for debugging)
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateType1Tets()
{
  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
      {
        int index = x + y * _xRes + z * _slabSize;

        VEC3F* vert000 = corner000(x,y,z);
        VEC3F* vert100 = corner100(x,y,z);
        VEC3F* vert010 = corner010(x,y,z);
        VEC3F* vert110 = corner110(x,y,z);
        VEC3F* vert001 = corner001(x,y,z);
        VEC3F* vert101 = corner101(x,y,z);
        VEC3F* vert011 = corner011(x,y,z);
        VEC3F* vert111 = corner111(x,y,z);

        VEC3F* center = centerVec(x,y,z);
        VEC3F* left = centerVec(x-1,y,z);
        VEC3F* near = centerVec(x,y,z-1);

        // type 1 tet
        // the second from the left in Fig. 2 of Isosurface Stuffing
        TET type1(center, left, vert001, vert011);
        _tets1.push_back(type1);

        // type 1 reflected about y
        TET type1Reflected(center, left, vert010, vert000);
        _tets1.push_back(type1Reflected);
      }
}

//////////////////////////////////////////////////////////////////////
// Generate all possible tets in the grid (for debugging)
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateType2Tets()
{
  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
      {
        int index = x + y * _xRes + z * _slabSize;

        VEC3F* vert000 = corner000(x,y,z);
        VEC3F* vert100 = corner100(x,y,z);
        VEC3F* vert010 = corner010(x,y,z);
        VEC3F* vert110 = corner110(x,y,z);
        VEC3F* vert001 = corner001(x,y,z);
        VEC3F* vert101 = corner101(x,y,z);
        VEC3F* vert011 = corner011(x,y,z);
        VEC3F* vert111 = corner111(x,y,z);

        VEC3F* center = centerVec(x,y,z);
        VEC3F* left = centerVec(x-1,y,z);
        VEC3F* near = centerVec(x,y,z-1);

        // type 2 tet
        // rightmost in Fig. 2 of Isosurface Stuffing
        TET type2(center, near, vert000, vert010);
        _tets2.push_back(type2);

        // type 2 tet reflecte about x
        // rightmost in Fig. 2 of Isosurface Stuffing
        TET type2Reflected(center, near, vert110, vert100);
        _tets2.push_back(type2Reflected);
      }
}

//////////////////////////////////////////////////////////////////////
// Generate all possible tets in the grid (for debugging)
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateType3Tets()
{
  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
      {
        int index = x + y * _xRes + z * _slabSize;

        VEC3F* vert000 = corner000(x,y,z);
        VEC3F* vert100 = corner100(x,y,z);
        VEC3F* vert010 = corner010(x,y,z);
        VEC3F* vert110 = corner110(x,y,z);
        VEC3F* vert001 = corner001(x,y,z);
        VEC3F* vert101 = corner101(x,y,z);
        VEC3F* vert011 = corner011(x,y,z);
        VEC3F* vert111 = corner111(x,y,z);

        VEC3F* center = centerVec(x,y,z);
        VEC3F* left = centerVec(x-1,y,z);
        VEC3F* near = centerVec(x,y,z-1);
        VEC3F* below = centerVec(x,y-1,z);

        // type 3 tet
        TET type3(center, below, vert001, vert000);
        _tets3.push_back(type3);

        // type 3 tet reflected
        TET type3Reflected(center, below, vert100, vert101);
        _tets3.push_back(type3Reflected);
      }
}

//////////////////////////////////////////////////////////////////////
// Generate all possible tets in the grid (for debugging)
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateType4Tets()
{
  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
      {
        int index = x + y * _xRes + z * _slabSize;

        VEC3F* vert000 = corner000(x,y,z);
        VEC3F* vert100 = corner100(x,y,z);
        VEC3F* vert010 = corner010(x,y,z);
        VEC3F* vert110 = corner110(x,y,z);
        VEC3F* vert001 = corner001(x,y,z);
        VEC3F* vert101 = corner101(x,y,z);
        VEC3F* vert011 = corner011(x,y,z);
        VEC3F* vert111 = corner111(x,y,z);

        VEC3F* center = centerVec(x,y,z);
        VEC3F* left = centerVec(x-1,y,z);
        VEC3F* near = centerVec(x,y,z-1);
        VEC3F* below = centerVec(x,y-1,z);

        // type 4 tet
        TET type4(center, below, vert000, vert100);
        _tets4.push_back(type4);

        // type 4 tet reflected
        TET type4Reflected(center, below, vert101, vert001);
        _tets4.push_back(type4Reflected);
      }
}

//////////////////////////////////////////////////////////////////////
// Generate all possible tets in the grid (for debugging)
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateType5Tets()
{
  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
      {
        int index = x + y * _xRes + z * _slabSize;

        VEC3F* vert000 = corner000(x,y,z);
        VEC3F* vert100 = corner100(x,y,z);
        VEC3F* vert010 = corner010(x,y,z);
        VEC3F* vert110 = corner110(x,y,z);
        VEC3F* vert001 = corner001(x,y,z);
        VEC3F* vert101 = corner101(x,y,z);
        VEC3F* vert011 = corner011(x,y,z);
        VEC3F* vert111 = corner111(x,y,z);

        VEC3F* center = centerVec(x,y,z);
        VEC3F* left = centerVec(x-1,y,z);
        VEC3F* near = centerVec(x,y,z-1);
        VEC3F* below = centerVec(x,y-1,z);

        // type 4 tet
        TET type5(center, near, vert010, vert110);
        _tets5.push_back(type5);

        // type 4 tet reflected
        TET type5Reflected(center, near, vert100, vert000);
        _tets5.push_back(type5Reflected);
      }
}

//////////////////////////////////////////////////////////////////////
// generate the final tets
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCase4Only()
{
  // clear previous final tets
  _finalTets.clear();

  // for each inside tet
  for (int x = 0; x < _insideTets.size(); x++)
  {
    // count number of inside/outside/snapped points
    int insidePoints[] = {-1,-1,-1,-1};
    int outsidePoints[] = {-1,-1,-1,-1};
    int snappedPoints[] = {-1,-1,-1,-1};
    int insideTotal = 0;
    int outsideTotal = 0;
    int snappedTotal = 0;
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexIDs[_insideTets[x]->vertices[y]];
      map<int,bool>::iterator finder = _snapped.find(vertexID);
      if (finder != _snapped.end())
      {
        snappedPoints[y] = vertexID;
        snappedTotal++;
      }
      else if (_inside[vertexID])
      {
        insidePoints[y] = vertexID;
        insideTotal++;
      }
      else
      {
        outsidePoints[y] = vertexID;
        outsideTotal++;
      }
    }

    // count number of edge cuts
    VEC3F* cuts[] = {NULL, NULL, NULL, NULL, NULL, NULL};
    pair<VEC3F*, VEC3F*> edges[6];
    edges[0] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[1]);
    edges[1] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[2]);
    edges[2] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[3]);
    edges[3] = createEdge(_insideTets[x]->vertices[1],
                          _insideTets[x]->vertices[2]);
    edges[4] = createEdge(_insideTets[x]->vertices[2],
                          _insideTets[x]->vertices[3]);
    edges[5] = createEdge(_insideTets[x]->vertices[3],
                          _insideTets[x]->vertices[1]);

    // look for a cut along each edge
    int cutsTotal = 0;
    for (int y = 0; y < 6; y++)
    {
      map<pair<VEC3F*,VEC3F*>, VEC3F>::iterator cutFinder;
      cutFinder = _cutPoints.find(edges[y]);
      if (cutFinder != _cutPoints.end())
      {
        cuts[y] = &(cutFinder->second);
        cutsTotal++;
      }
    }
    
    // case 0-3 in Figure 3 of paper
    if (outsideTotal == 0)
    {
      _finalTets.push_back(_insideTets[x]);
      continue;
    }

    // case 4
    if (outsideTotal == 1 && cutsTotal == 1)
    {
      generateCase4(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// generate the final tets
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCase5Only()
{
  // clear previous final tets
  _finalTets.clear();

  // for each inside tet
  for (int x = 0; x < _insideTets.size(); x++)
  {
    // count number of inside/outside/snapped points
    int insidePoints[] = {-1,-1,-1,-1};
    int outsidePoints[] = {-1,-1,-1,-1};
    int snappedPoints[] = {-1,-1,-1,-1};
    int insideTotal = 0;
    int outsideTotal = 0;
    int snappedTotal = 0;
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexIDs[_insideTets[x]->vertices[y]];
      map<int,bool>::iterator finder = _snapped.find(vertexID);
      if (finder != _snapped.end())
      {
        snappedPoints[y] = vertexID;
        snappedTotal++;
      }
      else if (_inside[vertexID])
      {
        insidePoints[y] = vertexID;
        insideTotal++;
      }
      else
      {
        outsidePoints[y] = vertexID;
        outsideTotal++;
      }
    }

    // count number of edge cuts
    VEC3F* cuts[] = {NULL, NULL, NULL, NULL, NULL, NULL};
    pair<VEC3F*, VEC3F*> edges[6];
    edges[0] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[1]);
    edges[1] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[2]);
    edges[2] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[3]);
    edges[3] = createEdge(_insideTets[x]->vertices[1],
                          _insideTets[x]->vertices[2]);
    edges[4] = createEdge(_insideTets[x]->vertices[2],
                          _insideTets[x]->vertices[3]);
    edges[5] = createEdge(_insideTets[x]->vertices[3],
                          _insideTets[x]->vertices[1]);

    // look for a cut along each edge
    int cutsTotal = 0;
    for (int y = 0; y < 6; y++)
    {
      map<pair<VEC3F*,VEC3F*>, VEC3F>::iterator cutFinder;
      cutFinder = _cutPoints.find(edges[y]);
      if (cutFinder != _cutPoints.end())
      {
        cuts[y] = &(cutFinder->second);
        cutsTotal++;
      }
    }
    
    // case 0-3 in Figure 3 of paper
    if (outsideTotal == 0)
    {
      _finalTets.push_back(_insideTets[x]);
      continue;
    }
    // case 4
    if (outsideTotal == 1 && cutsTotal == 1)
    {
      generateCase4(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
    // case 5
    if (outsideTotal == 2 && cutsTotal == 2)
    {
      generateCase5(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// generate the final tets
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCase6Only()
{
  // clear previous final tets
  _finalTets.clear();

  // for each inside tet
  for (int x = 0; x < _insideTets.size(); x++)
  {
    // count number of inside/outside/snapped points
    int insidePoints[] = {-1,-1,-1,-1};
    int outsidePoints[] = {-1,-1,-1,-1};
    int snappedPoints[] = {-1,-1,-1,-1};
    int insideTotal = 0;
    int outsideTotal = 0;
    int snappedTotal = 0;
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexIDs[_insideTets[x]->vertices[y]];
      map<int,bool>::iterator finder = _snapped.find(vertexID);
      if (finder != _snapped.end())
      {
        snappedPoints[y] = vertexID;
        snappedTotal++;
      }
      else if (_inside[vertexID])
      {
        insidePoints[y] = vertexID;
        insideTotal++;
      }
      else
      {
        outsidePoints[y] = vertexID;
        outsideTotal++;
      }
    }

    // count number of edge cuts
    VEC3F* cuts[] = {NULL, NULL, NULL, NULL, NULL, NULL};
    pair<VEC3F*, VEC3F*> edges[6];
    edges[0] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[1]);
    edges[1] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[2]);
    edges[2] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[3]);
    edges[3] = createEdge(_insideTets[x]->vertices[1],
                          _insideTets[x]->vertices[2]);
    edges[4] = createEdge(_insideTets[x]->vertices[2],
                          _insideTets[x]->vertices[3]);
    edges[5] = createEdge(_insideTets[x]->vertices[3],
                          _insideTets[x]->vertices[1]);

    // look for a cut along each edge
    int cutsTotal = 0;
    for (int y = 0; y < 6; y++)
    {
      map<pair<VEC3F*,VEC3F*>, VEC3F>::iterator cutFinder;
      cutFinder = _cutPoints.find(edges[y]);
      if (cutFinder != _cutPoints.end())
      {
        cuts[y] = &(cutFinder->second);
        cutsTotal++;
      }
    }
    
    // case 0-3 in Figure 3 of paper
    if (outsideTotal == 0)
    {
      _finalTets.push_back(_insideTets[x]);
      continue;
    }

    // case 4
    if (outsideTotal == 1 && cutsTotal == 1)
    {
      generateCase4(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
    // case 5
    if (outsideTotal == 2 && cutsTotal == 2)
    {
      generateCase5(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
    // case 6
    if (insideTotal == 1 && cutsTotal == 3 && outsideTotal == 3)
    {
      generateCase6(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }

  }
}

//////////////////////////////////////////////////////////////////////
// generate the final tets
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCase1Only()
{
  // clear previous final tets
  _finalTets.clear();

  // for each inside tet
  for (int x = 0; x < _insideTets.size(); x++)
  {
    // count number of inside/outside/snapped points
    int insidePoints[] = {-1,-1,-1,-1};
    int outsidePoints[] = {-1,-1,-1,-1};
    int snappedPoints[] = {-1,-1,-1,-1};
    int insideTotal = 0;
    int outsideTotal = 0;
    int snappedTotal = 0;
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexIDs[_insideTets[x]->vertices[y]];
      map<int,bool>::iterator finder = _snapped.find(vertexID);
      if (finder != _snapped.end())
      {
        snappedPoints[y] = vertexID;
        snappedTotal++;
      }
      else if (_inside[vertexID])
      {
        insidePoints[y] = vertexID;
        insideTotal++;
      }
      else
      {
        outsidePoints[y] = vertexID;
        outsideTotal++;
      }
    }

    // count number of edge cuts
    VEC3F* cuts[] = {NULL, NULL, NULL, NULL, NULL, NULL};
    pair<VEC3F*, VEC3F*> edges[6];
    edges[0] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[1]);
    edges[1] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[2]);
    edges[2] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[3]);
    edges[3] = createEdge(_insideTets[x]->vertices[1],
                          _insideTets[x]->vertices[2]);
    edges[4] = createEdge(_insideTets[x]->vertices[2],
                          _insideTets[x]->vertices[3]);
    edges[5] = createEdge(_insideTets[x]->vertices[3],
                          _insideTets[x]->vertices[1]);

    // look for a cut along each edge
    int cutsTotal = 0;
    for (int y = 0; y < 6; y++)
    {
      map<pair<VEC3F*,VEC3F*>, VEC3F>::iterator cutFinder;
      cutFinder = _cutPoints.find(edges[y]);
      if (cutFinder != _cutPoints.end())
      {
        cuts[y] = &(cutFinder->second);
        cutsTotal++;
      }
    }
    
    // case 0-3 in Figure 3 of paper
    if (outsideTotal == 0)
    {
      _finalTets.push_back(_insideTets[x]);
      continue;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// generate the final tets
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCase7Only()
{
  // clear previous final tets
  _finalTets.clear();

  // for each inside tet
  for (int x = 0; x < _insideTets.size(); x++)
  {
    // count number of inside/outside/snapped points
    int insidePoints[] = {-1,-1,-1,-1};
    int outsidePoints[] = {-1,-1,-1,-1};
    int snappedPoints[] = {-1,-1,-1,-1};
    int insideTotal = 0;
    int outsideTotal = 0;
    int snappedTotal = 0;
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexIDs[_insideTets[x]->vertices[y]];
      map<int,bool>::iterator finder = _snapped.find(vertexID);
      if (finder != _snapped.end())
      {
        snappedPoints[y] = vertexID;
        snappedTotal++;
      }
      else if (_inside[vertexID])
      {
        insidePoints[y] = vertexID;
        insideTotal++;
      }
      else
      {
        outsidePoints[y] = vertexID;
        outsideTotal++;
      }
    }

    // count number of edge cuts
    VEC3F* cuts[] = {NULL, NULL, NULL, NULL, NULL, NULL};
    pair<VEC3F*, VEC3F*> edges[6];
    edges[0] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[1]);
    edges[1] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[2]);
    edges[2] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[3]);
    edges[3] = createEdge(_insideTets[x]->vertices[1],
                          _insideTets[x]->vertices[2]);
    edges[4] = createEdge(_insideTets[x]->vertices[2],
                          _insideTets[x]->vertices[3]);
    edges[5] = createEdge(_insideTets[x]->vertices[3],
                          _insideTets[x]->vertices[1]);

    // look for a cut along each edge
    int cutsTotal = 0;
    for (int y = 0; y < 6; y++)
    {
      map<pair<VEC3F*,VEC3F*>, VEC3F>::iterator cutFinder;
      cutFinder = _cutPoints.find(edges[y]);
      if (cutFinder != _cutPoints.end())
      {
        cuts[y] = &(cutFinder->second);
        cutsTotal++;
      }
    }
    
    // case 0-3 in Figure 3 of paper
    if (outsideTotal == 0)
    {
      _finalTets.push_back(_insideTets[x]);
      continue;
    }

    // case 4
    if (outsideTotal == 1 && cutsTotal == 1)
    {
      generateCase4(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
    // case 5
    if (outsideTotal == 2 && cutsTotal == 2)
    {
      generateCase5(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
    // case 6
    if (insideTotal == 1 && cutsTotal == 3 && outsideTotal == 3)
    {
      generateCase6(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
    // cases 7 and 9
    if (outsideTotal == 1 && snappedTotal == 1)
    {
      // get the snapped vertex
      VEC3F* snappedVertex = NULL;
      for (int y = 0; y < 4; y++)
        if (snappedPoints[y] > 0)
          snappedVertex = &_vertices[snappedPoints[y]];

      // get the outside vertex
      VEC3F* outsideVertex = NULL;
      for (int y = 0; y < 4; y++)
        if (outsidePoints[y] > 0)
          outsideVertex = &_vertices[outsidePoints[y]];

      // make an edge
      EDGE edge = createEdge(snappedVertex, outsideVertex);

      // get the edge color
      EDGE_TYPE edgeType = _edgeTypes[edge];

      if (edgeType == RED)
        generateCase7(_insideTets[x], snappedPoints, insidePoints, outsidePoints, cuts, edges);

      continue;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// generate the final tets
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCase9Only()
{
  // clear previous final tets
  _finalTets.clear();

  // for each inside tet
  for (int x = 0; x < _insideTets.size(); x++)
  {
    // count number of inside/outside/snapped points
    int insidePoints[] = {-1,-1,-1,-1};
    int outsidePoints[] = {-1,-1,-1,-1};
    int snappedPoints[] = {-1,-1,-1,-1};
    int insideTotal = 0;
    int outsideTotal = 0;
    int snappedTotal = 0;
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexIDs[_insideTets[x]->vertices[y]];
      map<int,bool>::iterator finder = _snapped.find(vertexID);
      if (finder != _snapped.end())
      {
        snappedPoints[y] = vertexID;
        snappedTotal++;
      }
      else if (_inside[vertexID])
      {
        insidePoints[y] = vertexID;
        insideTotal++;
      }
      else
      {
        outsidePoints[y] = vertexID;
        outsideTotal++;
      }
    }

    // count number of edge cuts
    VEC3F* cuts[] = {NULL, NULL, NULL, NULL, NULL, NULL};
    pair<VEC3F*, VEC3F*> edges[6];
    edges[0] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[1]);
    edges[1] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[2]);
    edges[2] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[3]);
    edges[3] = createEdge(_insideTets[x]->vertices[1],
                          _insideTets[x]->vertices[2]);
    edges[4] = createEdge(_insideTets[x]->vertices[2],
                          _insideTets[x]->vertices[3]);
    edges[5] = createEdge(_insideTets[x]->vertices[3],
                          _insideTets[x]->vertices[1]);

    // look for a cut along each edge
    int cutsTotal = 0;
    for (int y = 0; y < 6; y++)
    {
      map<pair<VEC3F*,VEC3F*>, VEC3F>::iterator cutFinder;
      cutFinder = _cutPoints.find(edges[y]);
      if (cutFinder != _cutPoints.end())
      {
        cuts[y] = &(cutFinder->second);
        cutsTotal++;
      }
    }
    
    // case 0-3 in Figure 3 of paper
    if (outsideTotal == 0)
    {
      _finalTets.push_back(_insideTets[x]);
      continue;
    }

    // case 4
    if (outsideTotal == 1 && cutsTotal == 1)
    {
      generateCase4(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
    // case 5
    if (outsideTotal == 2 && cutsTotal == 2)
    {
      generateCase5(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
    // case 6
    if (insideTotal == 1 && cutsTotal == 3 && outsideTotal == 3)
    {
      generateCase6(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
    // cases 7 and 9
    if (outsideTotal == 1 && snappedTotal == 1)
    {
      // get the snapped vertex
      VEC3F* snappedVertex = NULL;
      for (int y = 0; y < 4; y++)
        if (snappedPoints[y] > 0)
          snappedVertex = &_vertices[snappedPoints[y]];

      // get the outside vertex
      VEC3F* outsideVertex = NULL;
      for (int y = 0; y < 4; y++)
        if (outsidePoints[y] > 0)
          outsideVertex = &_vertices[outsidePoints[y]];

      // make an edge
      EDGE edge = createEdge(snappedVertex, outsideVertex);

      // get the edge color
      EDGE_TYPE edgeType = _edgeTypes[edge];

      if (edgeType == RED)
        generateCase7(_insideTets[x], snappedPoints, insidePoints, outsidePoints, cuts, edges);
      else
        generateCase9(_insideTets[x], insidePoints, outsidePoints, cuts, edges);

      continue;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// generate the final tets
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCase10Only()
{
  // clear previous final tets
  _finalTets.clear();

  // for each inside tet
  for (int x = 0; x < _insideTets.size(); x++)
  {
    // count number of inside/outside/snapped points
    int insidePoints[] = {-1,-1,-1,-1};
    int outsidePoints[] = {-1,-1,-1,-1};
    int snappedPoints[] = {-1,-1,-1,-1};
    int insideTotal = 0;
    int outsideTotal = 0;
    int snappedTotal = 0;
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexIDs[_insideTets[x]->vertices[y]];
      map<int,bool>::iterator finder = _snapped.find(vertexID);
      if (finder != _snapped.end())
      {
        snappedPoints[y] = vertexID;
        snappedTotal++;
      }
      else if (_inside[vertexID])
      {
        insidePoints[y] = vertexID;
        insideTotal++;
      }
      else
      {
        outsidePoints[y] = vertexID;
        outsideTotal++;
      }
    }

    // count number of edge cuts
    VEC3F* cuts[] = {NULL, NULL, NULL, NULL, NULL, NULL};
    pair<VEC3F*, VEC3F*> edges[6];
    edges[0] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[1]);
    edges[1] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[2]);
    edges[2] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[3]);
    edges[3] = createEdge(_insideTets[x]->vertices[1],
                          _insideTets[x]->vertices[2]);
    edges[4] = createEdge(_insideTets[x]->vertices[2],
                          _insideTets[x]->vertices[3]);
    edges[5] = createEdge(_insideTets[x]->vertices[3],
                          _insideTets[x]->vertices[1]);

    // look for a cut along each edge
    int cutsTotal = 0;
    for (int y = 0; y < 6; y++)
    {
      map<pair<VEC3F*,VEC3F*>, VEC3F>::iterator cutFinder;
      cutFinder = _cutPoints.find(edges[y]);
      if (cutFinder != _cutPoints.end())
      {
        cuts[y] = &(cutFinder->second);
        cutsTotal++;
      }
    }
    
    // case 0-3 in Figure 3 of paper
    if (outsideTotal == 0)
    {
      _finalTets.push_back(_insideTets[x]);
      continue;
    }

    // case 4
    if (outsideTotal == 1 && cutsTotal == 1)
    {
      generateCase4(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
    // case 5
    if (outsideTotal == 2 && cutsTotal == 2)
    {
      generateCase5(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
    // case 6
    if (insideTotal == 1 && cutsTotal == 3 && outsideTotal == 3)
    {
      generateCase6(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
    // cases 7 and 9
    if (outsideTotal == 1 && snappedTotal == 1)
    {
      // get the snapped vertex
      VEC3F* snappedVertex = NULL;
      for (int y = 0; y < 4; y++)
        if (snappedPoints[y] > 0)
          snappedVertex = &_vertices[snappedPoints[y]];

      // get the outside vertex
      VEC3F* outsideVertex = NULL;
      for (int y = 0; y < 4; y++)
        if (outsidePoints[y] > 0)
          outsideVertex = &_vertices[outsidePoints[y]];

      // make an edge
      EDGE edge = createEdge(snappedVertex, outsideVertex);

      // get the edge color
      EDGE_TYPE edgeType = _edgeTypes[edge];

      if (edgeType == RED)
        generateCase7(_insideTets[x], snappedPoints, insidePoints, outsidePoints, cuts, edges);
      else
        generateCase9(_insideTets[x], insidePoints, outsidePoints, cuts, edges);

      continue;
    }
    // case 10
    if (outsideTotal == 1 && insideTotal == 3)
      {
        generateCase10(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
        continue;
      }
  }
}

//////////////////////////////////////////////////////////////////////
// load in the inside/outside info for faster debugging
//////////////////////////////////////////////////////////////////////
bool ISO_STUFFER::loadInsideOutside(const char* filename)
{
  _inside.clear();

  FILE* file = fopen(filename, "rb");
  if (file == NULL)
  {
    cout << " INSIDE CACHE LOAD FAILED " << endl;
    return false;
  }

  int size;
  fread((void*)&size, sizeof(int), 1, file);
  for (int x = 0; x < size; x++)
  {
    int index;
    bool insideOutside;
    fread((void*)&index, sizeof(int), 1, file);
    fread((void*)&insideOutside, sizeof(bool), 1, file);

    _inside[index] = insideOutside;
  }

  fclose(file);
  cout << " Loaded cached inside/outside info from " << filename << endl;
  return true;
}

//////////////////////////////////////////////////////////////////////
// load in the distances for faster debugging
//////////////////////////////////////////////////////////////////////
bool ISO_STUFFER::loadDistances(const char* filename)
{
  _distances.clear();

  FILE* file = fopen(filename, "rb");
  if (file == NULL)
  {
    cout << " CACHE LOAD FAILED " << endl;
    return false;
  }

  int size;
  fread((void*)&size, sizeof(int), 1, file);
  for (int x = 0; x < size; x++)
  {
    int index;
    float distance;
    fread((void*)&index, sizeof(int), 1, file);
    fread((void*)&distance, sizeof(Real), 1, file);

    _distances[index] = distance;
  }

  fclose(file);
  cout << " Loaded cached distances from " << filename << endl;
  return true;
}

//////////////////////////////////////////////////////////////////////
// save out inside/outside info for faster debugging
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::saveInsideOutside(const char* filename)
{
  FILE* file = fopen(filename, "wb");

  if (file == NULL)
    cout << " CACHE SAVE FAILED " << endl;

  int size = _inside.size();
  fwrite((void*)&size, sizeof(int), 1, file);

  map<int,bool>::iterator i;
  for (i = _inside.begin(); i != _inside.end(); ++i)
  {
    int index = i->first;
    bool inside = i->second;
    fwrite((void*)&index, sizeof(int), 1, file);
    fwrite((void*)&inside, sizeof(bool), 1, file);
  }
  fclose(file);
  cout << " Saved cached inside/outside to " << filename << endl;
}

//////////////////////////////////////////////////////////////////////
// save out distances for faster debugging
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::saveDistances(const char* filename)
{
  FILE* file = fopen(filename, "wb");

  if (file == NULL)
    cout << " CACHE SAVE FAILED " << endl;

  int size = _distances.size();
  fwrite((void*)&size, sizeof(int), 1, file);

  map<int,Real>::iterator i;
  for (i = _distances.begin(); i != _distances.end(); ++i)
  {
    int index = i->first;
    float distance = i->second;
    fwrite((void*)&index, sizeof(int), 1, file);
    fwrite((void*)&distance, sizeof(Real), 1, file);
  }
  fclose(file);
  cout << " Saved cached distances to " << filename << endl;
}
