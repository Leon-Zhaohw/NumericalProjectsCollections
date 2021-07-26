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
// GLSL_TET_MESH.cpp
//
//////////////////////////////////////////////////////////////////////

#include <GLSL_TET_MESH.h>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <algorithm>

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
GLSL_TET_MESH::GLSL_TET_MESH(SUBSPACE_TET_MESH *tetMesh, OBJ *embeddedMesh) :
  _tetMesh(tetMesh),
  _embeddedMesh(embeddedMesh)
{
  
  _surfaceVertVboID = 0;
  _surfaceVertNeighbor0VboID = 0;
  _surfaceVertNeighbor1VboID = 0;
  _surfaceVertSize = 0;
 
  _n1RestTexID = 0;
  _n2RestTexID = 0; // */
 
  constructTetMeshVbo();
  constructUBasisTex();
  constructEmbeddedMeshBuffers(embeddedMesh);
  constructTetFaceInfo();

  _tetVertexShader = new GLSL_SHADER();
  if (_tetVertexShader != NULL)
  {
    _tetVertexShader->attachVert("src/rendering/glsl/tetvertex.vert");
    _tetVertexShader->attachFrag("src/rendering/glsl/tetvertex.frag");
    _tetVertexShader->setProgramConst("Q_SIZE", _UmatrixPaddedWidth/4, 1);
    _tetVertexShader->setProgramConst("Q_STEP", (float)(_UbasisTexWidth/4), 1);
    _tetVertexShader->compile();
  }
  _genRotationShader = new GLSL_SHADER();
  if (_genRotationShader != NULL)
  {
    _genRotationShader->attachVert("src/rendering/glsl/tetgenrot.vert");
    _genRotationShader->attachFrag("src/rendering/glsl/tetgenrot.frag");
    _genRotationShader->compile();
  }
  _drawBasesShader = new GLSL_SHADER();
  if (_drawBasesShader != NULL)
  {
    _drawBasesShader->attachVert("src/rendering/glsl/drawbases1.vert");
    _drawBasesShader->attachGeom("src/rendering/glsl/drawbases1.geom");
    _drawBasesShader->attachFrag("src/rendering/glsl/drawbases.frag");
    _drawBasesShader->setGeomInput(GL_POINTS);
    _drawBasesShader->setGeomOutput(GL_TRIANGLE_STRIP);
    _drawBasesShader->setGeomVerticesOut(12);
    _drawBasesShader->compile();
  }  
  _drawNeighborsShader = new GLSL_SHADER();
  if (_drawNeighborsShader != NULL)
  {
    _drawNeighborsShader->attachVert("src/rendering/glsl/drawbases2.vert");
    _drawNeighborsShader->attachGeom("src/rendering/glsl/drawbases2.geom");
    _drawNeighborsShader->attachFrag("src/rendering/glsl/drawbases.frag");
    _drawNeighborsShader->setGeomInput(GL_POINTS);
    _drawNeighborsShader->setGeomOutput(GL_LINE_STRIP);
    _drawNeighborsShader->setGeomVerticesOut(12);
    _drawNeighborsShader->compile();
  }  
  _drawFacesShader = new GLSL_SHADER();
  if (_drawFacesShader != NULL)
  {
    _drawFacesShader->attachVert("src/rendering/glsl/drawfaces.vert");
    _drawFacesShader->attachGeom("src/rendering/glsl/drawfaces.geom");
    _drawFacesShader->attachFrag("src/rendering/glsl/drawfaces.frag");
    _drawFacesShader->setGeomInput(GL_POINTS);
    _drawFacesShader->setGeomOutput(GL_TRIANGLE_STRIP);
    _drawFacesShader->setGeomVerticesOut(6);
    _drawFacesShader->compile();
  }  
  _drawAOVShader = new GLSL_SHADER();
  if (_drawAOVShader != NULL)
  {
    _drawAOVShader->attachVert("src/rendering/glsl/aov1.vert");
    _drawAOVShader->attachGeom("src/rendering/glsl/aov1.geom");
    _drawAOVShader->attachFrag("src/rendering/glsl/aov1.frag");
    _drawAOVShader->setGeomInput(GL_POINTS);
    _drawAOVShader->setGeomOutput(GL_TRIANGLE_STRIP);
    _drawAOVShader->setGeomVerticesOut(14);
    _drawAOVShader->compile();
  }  
  _drawAOV2Shader = new GLSL_SHADER();
  if (_drawAOV2Shader != NULL)
  {
    _drawAOV2Shader->attachVert("src/rendering/glsl/aov4.vert");
    _drawAOV2Shader->attachGeom("src/rendering/glsl/aov2.geom");
    _drawAOV2Shader->attachFrag("src/rendering/glsl/aov1.frag");
    _drawAOV2Shader->setGeomInput(GL_TRIANGLES);
    _drawAOV2Shader->setGeomOutput(GL_TRIANGLE_STRIP);
    _drawAOV2Shader->setGeomVerticesOut(16);
    _drawAOV2Shader->compile();
  }
  _drawAOV3Shader = new GLSL_SHADER();
  if (_drawAOV3Shader != NULL)
  {
    _drawAOV3Shader->attachVert("src/rendering/glsl/aov2.vert");
    _drawAOV3Shader->attachGeom("src/rendering/glsl/aov3.geom");
    _drawAOV3Shader->attachFrag("src/rendering/glsl/aov3.frag");
    _drawAOV3Shader->setGeomInput(GL_TRIANGLES);
    _drawAOV3Shader->setGeomOutput(GL_LINE_STRIP);
    _drawAOV3Shader->setGeomVerticesOut(26);
    _drawAOV3Shader->compile();
  } 

  _renderFromTextureTest = new GLSL_SHADER();
  if (_renderFromTextureTest != NULL)
  {
    _renderFromTextureTest->attachVert("src/rendering/glsl/embedded17.vert");
    _renderFromTextureTest->attachFrag("src/rendering/glsl/embedded17.frag");
    _renderFromTextureTest->compile();
  } 


}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
GLSL_TET_MESH::~GLSL_TET_MESH()
{
  if (_tetMeshVboID != 0)  {
    glDeleteBuffers(1, &_tetMeshVboID);
    _tetMeshVboID = 0;
  }
  if (_tetMeshInfoVboID != 0)  {
    glDeleteBuffers(1, &_tetMeshInfoVboID);
    _tetMeshInfoVboID = 0;
  }
  if (_tetFaceInfo0VboID != 0)  {
    glDeleteBuffers(1, &_tetFaceInfo0VboID);
    _tetFaceInfo0VboID = 0;
  }
  if (_tetFaceInfo1VboID != 0)  {
    glDeleteBuffers(1, &_tetFaceInfo1VboID);
    _tetFaceInfo1VboID = 0;
  }
  if (_embeddedMeshVboID != 0)  {
    glDeleteBuffers(1, &_embeddedMeshVboID);
    _embeddedMeshVboID = 0;
  }
  if (_embeddedNormalVboID != 0)  {
    glDeleteBuffers(1, &_embeddedNormalVboID);
    _embeddedNormalVboID = 0;
  }
  if (_embeddedVertLoc0VboID != 0)  {
    glDeleteBuffers(1, &_embeddedVertLoc0VboID);
    _embeddedVertLoc0VboID = 0;
  }
  if (_embeddedVertLoc1VboID != 0)  {
    glDeleteBuffers(1, &_embeddedVertLoc1VboID);
    _embeddedVertLoc1VboID = 0;
  }
  if (_embeddedVertLoc2VboID != 0)  {
    glDeleteBuffers(1, &_embeddedVertLoc2VboID);
    _embeddedVertLoc2VboID = 0;
  }
  if (_embeddedMeshIboID != 0)  {
    glDeleteBuffers(1, &_embeddedMeshIboID);
    _embeddedMeshIboID = 0;
  }

  if (_UcoordVboID != 0)  {
    glDeleteBuffers(1, &_UcoordVboID);
    _UcoordVboID = 0;
  }
  if (_UbasisTexID != 0)  {
    glDeleteTextures(1, &_UbasisTexID);
    _UbasisTexID = 0;
  }

  if (_surfaceVertVboID != 0)  {
    glDeleteBuffers(1, &_surfaceVertVboID);
    _surfaceVertVboID = 0;
  }
  if (_surfaceVertNeighbor0VboID != 0)  {
    glDeleteBuffers(1, &_surfaceVertNeighbor0VboID);
    _surfaceVertNeighbor0VboID = 0;
  }
  if (_surfaceVertNeighbor1VboID != 0)  {
    glDeleteBuffers(1, &_surfaceVertNeighbor1VboID);
    _surfaceVertNeighbor1VboID = 0;
  }
  if (_n1RestTexID != 0)  {
    glDeleteTextures(1, &_n1RestTexID);
    _n1RestTexID = 0;
  }
  if (_n2RestTexID != 0)  {
    glDeleteTextures(1, &_n2RestTexID);
    _n2RestTexID = 0;
  } // */

  if (_tetVertexShader != NULL)  {
    delete _tetVertexShader;
    _tetVertexShader = NULL;
  }
  if (_genRotationShader != NULL)  {
    delete _genRotationShader;
    _genRotationShader = NULL;
  }
  if (_drawFacesShader != NULL)  {
    delete _drawFacesShader;
    _drawFacesShader = NULL;
  }
  if (_drawBasesShader != NULL)  {
    delete _drawBasesShader;
    _drawBasesShader = NULL;
  } 
  if (_drawNeighborsShader != NULL)  {
    delete _drawNeighborsShader;
    _drawNeighborsShader = NULL;
  }  
  if (_drawAOVShader != NULL)  {
    delete _drawAOVShader;
    _drawAOVShader = NULL;
  }   
  if (_drawAOV2Shader != NULL)  {
    delete _drawAOV2Shader;
    _drawAOV2Shader = NULL;
  }  
  if (_drawAOV3Shader != NULL)  {
    delete _drawAOV3Shader;
    _drawAOV3Shader = NULL;
  }  
  if (_renderFromTextureTest != NULL)  {
    delete _renderFromTextureTest;
    _renderFromTextureTest = NULL;
  }   

  delete _Dv;
  delete _Df;
  delete _Q;
  _embeddedMesh = NULL;

}

//////////////////////////////////////////////////////////////////////
// Functions to generate random numbers 

// Normal distribution with mean mu, standard deviation sigma
float NormalRandom(float mu, float sigma)
{
  double var = sigma*sigma;
  double x, p, q;

  while(1)
  {
    x = (double)(rand() % 123987456) / 123987456.0;

    // Shift x from (0.0, 1.0) to (-10s, 10s)
    x = sigma * 20.0f * (x - 0.5f);
    // Find probability of x given normal distribution N(mu,var)
    p = exp(-((x-mu)*(x-mu))/(2.0f*var)) / sqrt(2.0f*M_PI*var);
    // Get uniform random to determine if x is our chosen value  
    q = (double)(rand() % 123987456) / 123987456.0;
    if (q < p)
      return (float)x;
  }
}

// Chi Squared with k degrees of freedom
float ChiSquareRandom(int k)
{
  float Q = 0.0f;
  float x = 0.0f;
  for (int i = 0; i < k; i++)
  {
    x = NormalRandom(0.0f, 1.0f);
    Q += x*x;
  }
  return Q;
}

// Student's t distribution with v degrees of freedom
float tRandom(int v)
{
  float df = (float)v;
  float Z = NormalRandom(0.0f,1.0f);
  float V = 0.0f;
  while (V <= 0.0f)
    V = ChiSquareRandom(v);

  Z = Z * sqrt(df / V);
  return Z;
}

//////////////////////////////////////////////////////////////////////
// Helper functions
Real clamp(Real value, Real min, Real max)
{
  if (value < min)
    return min;
  if (value > max)
    return max;
  return value;
}

Real dot(VEC3F& a, VEC3F& b)
{
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}


//////////////////////////////////////////////////////////////////////
// Construct a vector of
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH::constructDeltaVector(OBJ *embeddedMesh, Real maxDelta,
                                         Real threshold, int ra_max )
{
  MATRIX Q;
  MATRIX B;
  MATRIX Df;
  MATRIX Dv;
  bool isQ, isB, isDf, isDv;

  string filename = _tetMesh->filename();
  filename += string(".dface");
  isDf = Df.read(filename.c_str());
    
  filename = _tetMesh->filename();
  filename += string(".bmat");
  isB = B.read(filename.c_str());

  filename = _tetMesh->filename();
  filename += string(".dvert");
  isDv = Dv.read(filename.c_str());
 
  filename = _tetMesh->filename();
  filename += string(".qmat");
  isQ = Q.read(filename.c_str());
  
  filename = _tetMesh->filename();
  filename += string(".trainingset");
  FILE* file = fopen(filename.c_str(), "rb");

  if (file == NULL)
  {
    cout << " Tried to open training set file " << filename.c_str() << endl;
    return;
  }

  cout << " Found an training set file: " << filename.c_str() << endl;

  // read dimensions
  int size;
  fread((void*)&size, sizeof(int), 1, file);

  vector<VECTOR> trainingForces;
  vector<VECTOR> trainingQs;

  // read in the forces
  for (int x = 0; x < size; x++)
  {
    VECTOR sample(file);
    trainingForces.push_back(sample);
  }

  // read in the qs
  for (int x = 0; x < size; x++)
  {
    VECTOR sample(file);
    trainingQs.push_back(sample);
  }

  fclose(file);

  int numPoses = trainingQs.size();
  int embeddedVertSize = embeddedMesh->vertices.size();
  vector<VEC3> vertices;
  for (int i = 0; i < embeddedVertSize; i++)
    vertices.push_back(embeddedMesh->vertices[i]);

  vector<VEC3F>& baryEmbeddings = _tetMesh->barycentricEmbeddings();
  vector<int>& tetEmbeddings = _tetMesh->tetEmbeddings();
  
  int numFaces = embeddedMesh->faces.size();
  if (!isQ)   Q.resizeAndWipe(_tetMesh->q().size(), numPoses);
  if (!isDf)  Df.resizeAndWipe(numFaces, numPoses);
  if (!isDv)  Dv.resizeAndWipe(embeddedVertSize, numPoses);

  cout << " Number of poses in training set is: " << numPoses << endl;
  cout << " ---------------------------------------------------------------------------------------------------| " << endl;
  
  for (int i = 0; i < numPoses; i++)
  {
    cout << " Pose " << (i+1) << " / " << numPoses << ": ";
    flush(cout);
    VECTOR &q = _tetMesh->q();
    cout << ".";
    if (!isQ)
    {
      for (int j = 0; j < q.size(); j++)
      {
        q[j] = 0.025 * trainingQs[i][j];
        Q(j,i) = q[j];
        //cout << Q(i,j) << " ";
      }    
      //cout << endl;
    }

    cout << ".";
    flush(cout);
    _tetMesh->updateSurfaceMesh();

    vector<TET>& tets = _tetMesh->tets();

    cout << ".";
    flush(cout);
    for (int x = 0; x < embeddedVertSize; x++)
    {
      TET& tet = tets[tetEmbeddings[x]];
      VEC3F p = (baryEmbeddings[x][0] * (*tet.vertices[0])) + (baryEmbeddings[x][1] * (*tet.vertices[1])) + 
                (baryEmbeddings[x][2] * (*tet.vertices[2])) + ((1.0 - baryEmbeddings[x][0] - baryEmbeddings[x][1] - baryEmbeddings[x][2]) * (*tet.vertices[3]));
      embeddedMesh->vertices[x] = p;
    }

    embeddedMesh->ComputeVertexNormals();

    // Draw it here to see if it's correct.
    
    cout << ".";
    flush(cout);
    vector<OBJ::Face>& faces = embeddedMesh->faces;
    int step = (int)((Real)numFaces / 80.0);
    if (!isDf)
    {
      for (int f = 0; f < numFaces; f++)
      {
        if (f%step == 0)
        {
          cout << ".";
          flush(cout);
        }
        TRIANGLE tri(&embeddedMesh->vertices[faces[f].vertices[0]], &embeddedMesh->vertices[faces[f].vertices[1]], &embeddedMesh->vertices[faces[f].vertices[2]]);
        Df(f,i) = calculateDelta(tri, maxDelta, threshold, *embeddedMesh);
      }
    }
    cout << "!" << endl;
    flush(cout);
  }

  filename = _tetMesh->filename();
  filename += string(".qmat");
  Q.write(filename.c_str());

  cout << " Finished calculating Df, saving.." << endl;
  flush(cout);
 
  filename = _tetMesh->filename();
  filename += string(".dface");
  Df.write(filename.c_str());

  cout << " constructing Dv from Df; ";
  flush(cout);

  if (!isDv)
  {
    Dv.clear();
    for (int p = 0; p < numPoses; p++)
    {
      for (int f = 0; f < numFaces; f++)
      {
        int index = embeddedMesh->faces[f].vertices[0];
        if (Dv(index,p) < Df(f,p))
          Dv(index,p) = Df(f,p);
        index = embeddedMesh->faces[f].vertices[1];
        if (Dv(index,p) < Df(f,p))
          Dv(index,p) = Df(f,p);
        index = embeddedMesh->faces[f].vertices[2];
        if (Dv(index,p) < Df(f,p))
          Dv(index,p) = Df(f,p);
      }
    }
  }

  filename = _tetMesh->filename();
  filename += string(".dvert");
  Dv.write(filename.c_str());

  cout << "performing reduction on Dv; " << endl;
  flush(cout);

  MATRIX A(Dv.rows(), Dv.cols());
  VECTOR impt;
  Dv.PCA(A,impt);
  Real hi = impt(0);
  Real factor = 0.01;
  int ra = 1;
  int max = (ra_max < A.rows()) ? ra_max : A.rows();

  for (int i = 1; i < max; i++)
  {
    if (impt(i) > hi*factor)
    {
      ra = i+1;
    }
  }
  cout << "Impt: ";
  for (int i = 0; i < impt.size(); i++)
  {
    cout << impt(i) << "   ";
  }
  cout << endl;  flush(cout);


  cout << " A -> " << A.rows() << " x " << A.cols() << endl; flush(cout);
  MATRIX A_ra = A.getSubmatrix(0, A.rows(), 0, ra);
  cout << " A_ra = A reduced-> " << A_ra.rows() << " x " << A_ra.cols() << endl; flush(cout);
  MATRIX tempMatrix = A_ra.transpose() * Dv;
  cout << " M = A_ra^T * Dv -> " << tempMatrix.rows() << " x " << tempMatrix.cols() << endl; flush(cout);
  cout << " Q -> " << Q.rows() << " x " << Q.cols() << endl; flush(cout);
  MATRIX Qinv;
  Q.pseudoInverse(Qinv);
  cout << " Qinv -> " << Qinv.rows() << " x " << Qinv.cols() << endl; flush(cout);
  cout << " B = M * Qinv "; flush(cout);
  B = tempMatrix * Qinv;
  cout << " B -> " << B.rows() << " x " << B.cols() << endl; flush(cout);

  //Udelta = A_ra * B;     // This Udelta is the equivalent of the U basis matrix.
                         // Udelta * q = delta, vector for the pose defined by q
 
  filename = _tetMesh->filename();
  filename += string(".bmat");
  B.write(filename.c_str());

  _Dv = new MATRIX(Dv);
  _Df = new MATRIX(Df);
  _Q = new MATRIX(Q);
  _embeddedMesh = embeddedMesh;
  _prevPose = -10;
 
  cout << " complete!" << endl;

  //return tempMatrix;
}

//////////////////////////////////////////////////////////////////////
// Calculate a nice delta, given the triangle, a max delta, and threshold percent
//////////////////////////////////////////////////////////////////////
Real GLSL_TET_MESH::calculateDelta(TRIANGLE &T, Real delta_max, Real threshold, OBJ &obj)
{
  // Construct m[4] of prism P, given T and delta_max
  int vertsize = obj.vertices.size();
  vector<VEC3F> m;
  VEC3F a = *T.vertex(1) - *T.vertex(0);
  VEC3F b = *T.vertex(2) - *T.vertex(1);
  VEC3F c = *T.vertex(0) - *T.vertex(2);
  VEC3F tempCross = cross(b,a);   unitize(tempCross);
  VEC3F m3 = tempCross;

  tempCross = cross(a,m3);    unitize(tempCross);
  m.push_back(delta_max*tempCross);

  tempCross = cross(b,m3);    unitize(tempCross);
  m.push_back(delta_max*tempCross);

  tempCross = cross(c,m3);    unitize(tempCross);
  m.push_back(delta_max*tempCross);

  m.push_back(delta_max*m3);

  // Make a histogram that stores the AO contributions of points in certain distance ranges from T
  int bins = 25;
  int myBin;
  Real *H = new Real[bins];
  for (int j = 0; j < bins; j++)
    H[j] = 0.0;
 
  Real *AO = new Real[vertsize];
  Real *dist = new Real[vertsize];
  Real totAO = 0.0;
  Real curAO = 0.0;

  // Calculate the total AO and distances of the points in range of T 
  for (int i = 0; i < vertsize; i++)
  {
    VEC3F &x = obj.vertices[i];
    AO[i] = 0.0;
    dist[i] = 0.0;

    //  if (x is NOT inside prism(m[4],T)) go to next i
    if (isPointInVolume(x, T, m))
    {
      VEC3F &n = obj.normals[i];
      AO[i] = clamp(calculateAO(x, n, T), 0.0, 1.0);
      dist[i] = calculatePointTriDistance(x, T, m, delta_max);
      AO[i] *= (1.0 - (dist[i] / delta_max));

      myBin = (int)((Real)bins * (dist[i] / delta_max));
      myBin = (myBin < bins) ? myBin : (bins-1);
      H[myBin] += AO[i];
      totAO += AO[i];
    }
  }

  // Step from furthest distance bin to nearer bins, subtracting AO contributions along the way
  Real myDelta = 0.0;
  curAO = totAO;
  
  for (int x = bins-1; x >= 0; x--)
  {
    curAO -= H[x];

    // Once we drop below the requested threshold, pick the distance of 
    // the bin one step up (i.e. conservative)
    if (curAO < (threshold*totAO))
    {
      myDelta = ((Real)(x+1) / (Real)bins) * delta_max;
      break;
    }
  }

  // Clean up arrays we put on the heap
  delete H;
  delete AO;
  delete dist;

  return myDelta;
}

//////////////////////////////////////////////////////////////////////
// Determine whether point x is inside the triangular prism volume defined by t and m
//////////////////////////////////////////////////////////////////////
bool GLSL_TET_MESH::isPointInVolume(VEC3F& x, TRIANGLE& t, vector<VEC3F>& m)
{
  if (m.size() != 4)
    return false;

  VEC3F q;
    
  // Check the triangle plane (bottom)
  q = x - (*t.vertex(0));
  if (dot(q,m[3]) > 0.0)
    return false;

  // Check the top
  q = q + m[3];
  if (dot(q,m[3]) < 0.0)
    return false;

  // Check the three sides of the prism volume
  for (int i = 0; i < 3; i++)
  {
    q = x - (*t.vertex(i)) + m[i];
    if (dot(q,m[i]) < 0.0)
      return false;
  }

  // Got here so point is inside the volume
  return true;
}


//////////////////////////////////////////////////////////////////////
// Calculate the form factor accessibility that the triangle t projects onto the 
// point x, given a surface normal n (t clipped to the positive half-space)
//////////////////////////////////////////////////////////////////////
Real GLSL_TET_MESH::calculateAO(VEC3F& x, VEC3F& n, TRIANGLE& tri, bool debugmode)
{
  // Something is wrong in here.
  // Could be clipping, or form-factor calculation. I'm thinking clipping.
  //
  // !! Create a version of the box scene and test. !!
  
  Real epsilon = 0.001;
  // Clip to positive half space of tangent plane at this point
  Real dots[3];
  Real pLen;
  int pdots = 3;
  VEC3F P[4];
  P[0] = (*tri.vertex(0)) - x;
  P[1] = (*tri.vertex(1)) - x;
  P[2] = (*tri.vertex(2)) - x;
  P[3] = x - x;
  VEC3F N = n;
  unitize(N);
  if (debugmode)
  {
    cout << "calculateAO debug print:" << endl;
    cout << " P0: " << P[0];
    cout << " P1: " << P[1];
    cout << " P2: " << P[2];
    cout << " P3: " << P[3] << endl;
    cout << " N:  " << N << endl;
  }

  // Determine how many, and which points are in the positive half-space
  dots[0] = dot(P[0], N);
  pLen = dot(P[0],P[0]);
  if (pLen > 0.0)
    dots[0] = dots[0] / sqrt(pLen);
  if (dots[0] < -epsilon)
    pdots -= 1;
  
  dots[1] = dot(P[1], N);
  pLen = dot(P[1],P[1]);
  if (pLen > 0.0)
    dots[1] = dots[1] / sqrt(pLen);
  if (dots[1] < -epsilon)
    pdots -= 1;
  
  dots[2] = dot(P[2], N);
  pLen = dot(P[2],P[2]);
  if (pLen > 0.0)
    dots[2] = dots[2] / sqrt(pLen);
  if (dots[2] < -epsilon)
    pdots -= 1;

  if (debugmode)
    cout << "pdots: " << pdots << endl;
  
  // pdots is the number of vertices in the positive half-space
  if (pdots == 0) 
    return 0.0;

  // Clip the negative points to the positive half-space
  bool isTriangle = true;
  if (pdots == 1) 
  {    
    // Clip into triangle
    VEC3F edge;
    Real t;
    if (dots[0] >= -epsilon)
    {
     // color = vec3(1.0,0.0,0.0);
      edge = P[0] - P[1];
      t = (-dot(P[1], N)) / (epsilon + dot(edge, N));
      P[1] += (t * edge);
      edge = P[0] - P[2];
      t = (-dot(P[2],N)) / (epsilon + dot(edge, N));
      P[2] += (t * edge);
    } 
    else if (dots[1] >= -epsilon)
    {
      //color = vec3(1.0,1.0,0.0);

      edge = P[1] - P[0];
      t = (-dot(P[0],N)) / (epsilon + dot(edge, N));
      P[0] += (t * edge);
      edge = P[1] - P[2];
      t = (-dot(P[2],N)) / (epsilon + dot(edge, N));
      P[2] += (t * edge);
    } 
    else // if (dots[2] >= 0.0)
    {
      //color = vec3(1.0,0.0,1.0);
      edge = P[2] - P[0];
      t = (-dot(P[0],N)) / (epsilon + dot(edge, N));
      P[0] += (t * edge);
      edge = P[2] - P[1];
      t = (-dot(P[1],N)) / (epsilon + dot(edge, N));
      P[1] += (t * edge);
    }
  }

  // Clipped region will be a quad if two of three verts is negative
  if (pdots == 2) 
  {
    isTriangle = false;    
    // Clip into quad
    VEC3F edge;
    Real t;
    if (dots[0] < -epsilon)
    {
      //color = vec3(0.0,1.0,1.0);
      edge = P[2] - P[0];
      t = (-dot(P[0],N)) / (epsilon + dot(edge, N));
      P[3] = P[0] + (t * edge);
      edge = P[1] - P[0];
      t = (-dot(P[0],N)) / (epsilon + dot(edge, N));
      P[0] += (t * edge);
    }
    else if (dots[1] < -epsilon)
    {
      //color = vec3(0.0,0.0,1.0);
      P[3] = P[0];
      edge = P[0] - P[1];
      t = (-dot(P[1],N)) / (epsilon + dot(edge, N));
      P[0] = P[1] + (t * edge);
      edge = P[2] - P[1];
      t = (-dot(P[1],N)) / (epsilon + dot(edge, N));
      P[1] += (t * edge);
    } 
    else 
    {
      //color = vec3(0.5,0.5,0.5);
      edge = P[0] - P[2];
      t = (-dot(P[2],N)) / (epsilon + dot(edge, N));
      P[3] = P[2] + (t * edge);
      edge = P[1] - P[2];
      t = (-dot(P[2],N)) / (epsilon + dot(edge, N));
      P[2] += (t * edge);
    }
  }

  if (debugmode)
  {
    cout << "After clipping:" << endl;
    cout << " P0: " << P[0];
    cout << " P1: " << P[1];
    cout << " P2: " << P[2];
    cout << " P3: " << P[3] << endl;
  }

  if (dot(P[0],P[0]) > 0.0)
      unitize(P[0]);
  if (dot(P[1],P[1]) > 0.0)
      unitize(P[1]);
  if (dot(P[2],P[2]) > 0.0)
      unitize(P[2]);
  if (dot(P[3],P[3]) > 0.0)
      unitize(P[3]);

  // Calculate ambient occlusion accessibility using clipped polygon
  Real AOp = 0.0;
  Real contrib;
  Real PidotPj;
  VEC3F PixPj;
  if (debugmode) cout << " AO = ";
  if (isTriangle)
  {
    PidotPj = clamp(dot(P[0],P[1]),-1.0, 1.0);
    PixPj = cross(P[1],P[0]);
    if (dot(PixPj,PixPj) > 0.0)
      unitize(PixPj);
    contrib = acos(PidotPj) * dot(N, PixPj);
    if (debugmode) cout << contrib << " + ";
    AOp += contrib;
  
    PidotPj = clamp(dot(P[1],P[2]),-1.0, 1.0);
    PixPj = cross(P[2],P[1]);
    if (dot(PixPj,PixPj) > 0.0)
      unitize(PixPj);
    contrib = acos(PidotPj) * dot(N, PixPj);
    if (debugmode) cout << contrib << " + ";
    AOp += contrib;

    PidotPj = clamp(dot(P[2],P[0]),-1.0, 1.0);
    PixPj = cross(P[0],P[2]);
    if (dot(PixPj,PixPj) > 0.0)
      unitize(PixPj);
    contrib = acos(PidotPj) * dot(N, PixPj);
    if (debugmode) cout << contrib << " = ";
    AOp += contrib;
      
  } else {
    PidotPj = clamp(dot(P[0],P[1]),-1.0, 1.0);
    PixPj = cross(P[1],P[0]);
    if (dot(PixPj,PixPj) > 0.0)
      unitize(PixPj);
    contrib = acos(PidotPj) * dot(N, PixPj);
    if (debugmode) cout << contrib << " + ";
    AOp += contrib;

    PidotPj = clamp(dot(P[1],P[2]),-1.0, 1.0);
    PixPj = cross(P[2],P[1]);
    if (dot(PixPj,PixPj) > 0.0)
      unitize(PixPj);
    contrib = acos(PidotPj) * dot(N, PixPj);
    if (debugmode) cout << contrib << " + ";
    AOp += contrib;

    PidotPj = clamp(dot(P[2],P[3]),-1.0, 1.0);
    PixPj = cross(P[3],P[2]);
    if (dot(PixPj,PixPj) > 0.0)
      unitize(PixPj);
    contrib = acos(PidotPj) * dot(N, PixPj);
    if (debugmode) cout << contrib << " + ";
    AOp += contrib;

    PidotPj = clamp(dot(P[3],P[0]),-1.0, 1.0);
    PixPj = cross(P[0],P[3]);
    if (dot(PixPj,PixPj) > 0.0)
      unitize(PixPj);
    contrib = acos(PidotPj) * dot(N, PixPj);
    if (debugmode) cout << contrib << " = ";
    AOp += contrib;
    
  }

  if (debugmode) cout << AOp << " / 2pi = ";

  // Divide by 2pi to normalize it.
  AOp *= 0.159154943091895;  //   = 1.0 / (2.0 * M_PI);

  if (debugmode) cout << AOp << endl;
  
  return AOp;
}

//////////////////////////////////////////////////////////////////////
// Pseudo distance between t and x
//////////////////////////////////////////////////////////////////////
Real GLSL_TET_MESH::calculatePointTriDistance(VEC3F& x, TRIANGLE& t, vector<VEC3F>& m, Real delta)
{
  if (m.size() != 4)
    return delta;

  Real dist = 1.0;
  Real g;
  for (int i = 0; i < 4; i++)
  {
    VEC3F p = x - (*t.vertex(i%3));
    g = 1.0 + dot(p,m[i])/delta;
    g = (g < 0.0) ? 0.0 : ((g > 1.0) ? 1.0 : g);
    dist = dist * g;
  }
  
  // dist holds 1.0 near t, linearily falling to 0.0 outward to a distance of delta.
  // Want to remap this to be 0.0 at t, delta at delta.
  return delta*(1.0 - dist);
}

//////////////////////////////////////////////////////////////////////
// Draw method to display an embedded mesh
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH::debugDraw(int pose, bool drawdeltas)
{
  vector<OBJ::Face>& faces = _embeddedMesh->faces;
  int numFaces = faces.size();
  int vertSize = _embeddedMesh->vertices.size();

  Real delta = 0.25;
  
  // Draw it here to see if it's correct.
  if (_prevPose != pose)
  {
    vector<VEC3F>& baryEmbeddings = _tetMesh->barycentricEmbeddings();
    vector<int>& tetEmbeddings = _tetMesh->tetEmbeddings();
   
    VECTOR &q = _tetMesh->q();
    MATRIX &Q = *_Q;
    for (int j = 0; j < q.size(); j++)
    {
      q[j] = (pose >= 0) ? Q(j,pose) : 0.0;
      //q[j] *= 0.1;
    }    
 
    _tetMesh->updateSurfaceMesh();
    //_tetMesh->drawSurfaceFaces();

    vector<TET>& tets = _tetMesh->tets();

    if ((int)_embeddedMesh->texcoords.size() != vertSize)
      _embeddedMesh->texcoords.resize(vertSize);
    
    for (int x = 0; x < vertSize; x++)
    {
      TET& tet = tets[tetEmbeddings[x]];
      VEC3F p = (baryEmbeddings[x][0] * (*tet.vertices[0])) + (baryEmbeddings[x][1] * (*tet.vertices[1])) + 
                (baryEmbeddings[x][2] * (*tet.vertices[2])) + ((1.0 - baryEmbeddings[x][0] - baryEmbeddings[x][1] - baryEmbeddings[x][2]) * (*tet.vertices[3]));
      _embeddedMesh->vertices[x] = p;
      _embeddedMesh->texcoords[x][0] = _embeddedMesh->texcoords[x][1] = 0.0;
    }

    _embeddedMesh->ComputeVertexNormals();

    for (int f = 0; f < numFaces; f++)
    {
      //TRIANGLE T(&_embeddedMesh->vertices[faces[f].vertices[0]], &_embeddedMesh->vertices[faces[f].vertices[1]], &_embeddedMesh->vertices[faces[f].vertices[2]]);
      VEC3F &v0 = _embeddedMesh->vertices[faces[f].vertices[0]];
      VEC3F &v1 = _embeddedMesh->vertices[faces[f].vertices[1]];
      VEC3F &v2 = _embeddedMesh->vertices[faces[f].vertices[2]];
      TRIANGLE T(&v0,&v1,&v2);
      
      vector<VEC3F> m;                    // m appears to be calculated correctly
      VEC3F a = v1 - v0;
      VEC3F b = v2 - v1;
      VEC3F c = v0 - v2;
      VEC3F tempCross = cross(b,a);   unitize(tempCross);
      VEC3F m3 = tempCross;

      tempCross = cross(a,m3);    unitize(tempCross);
      m.push_back(delta*tempCross);

      tempCross = cross(b,m3);    unitize(tempCross);
      m.push_back(delta*tempCross);

      tempCross = cross(c,m3);    unitize(tempCross);
      m.push_back(delta*tempCross);

      m.push_back(delta*m3);

      // Calculate the total AO and distances of the points in range of T 
      /*for (int i = 0; i < vertSize; i++)
      {
        VEC3F &x = _embeddedMesh->vertices[i];

        //  if (x is NOT inside prism(m[4],T)) go to next i
        if (isPointInVolume(x, T, m))                           // isPointInVolume appears to be correct
        {
          VEC3F &n = _embeddedMesh->normals[i];                 
          Real AO = clamp(calculateAO(x, n, T), 0.0, 1.0); 
          Real dist = calculatePointTriDistance(x, T, m, delta);    
          AO *= (1.0 - (dist / delta));
          
          _embeddedMesh->texcoords[i][0] += AO;
        }
      }*/
    }

    _prevPose = pose;
  }

  /////////
  
  MATRIX &Df = *_Df;
    
  glBegin(GL_TRIANGLES);
  for (int f = 0; f < numFaces; f++)
  {
    float z;
    if (pose < 0 || pose >= 100)
    {
      z = 1.0f;
    } else {
      Real faceDelta = Df(f, pose);
      z = 4.0f * (float)faceDelta;
    }

    VEC3F &v0 = _embeddedMesh->vertices[faces[f].vertices[0]];
    VEC3F &v1 = _embeddedMesh->vertices[faces[f].vertices[1]];
    VEC3F &v2 = _embeddedMesh->vertices[faces[f].vertices[2]];

    //float z = _embeddedMesh->texcoords[faces[f].vertices[0]][0];
    glColor3f(z,z,z);
    glVertex3f((GLfloat)v0[0], (GLfloat)v0[1], (GLfloat)v0[2]);

    //z = _embeddedMesh->texcoords[faces[f].vertices[1]][0];
    glColor3f(z,z,z);
    glVertex3f((GLfloat)v1[0], (GLfloat)v1[1], (GLfloat)v1[2]);

    //z = _embeddedMesh->texcoords[faces[f].vertices[2]][0];
    glColor3f(z,z,z);
    glVertex3f((GLfloat)v2[0], (GLfloat)v2[1], (GLfloat)v2[2]);
  }
  glEnd();

  if (drawdeltas == false)
    return;

  
  glColor3f(1.0f,1.0f,1.0f);
  glBegin(GL_LINES);
  for (int f = 0; f < numFaces; f++)
  {
    if (f%50 != 0)
      continue;

    //TRIANGLE tri(&_embeddedMesh->vertices[faces[f].vertices[0]], &_embeddedMesh->vertices[faces[f].vertices[1]], &_embeddedMesh->vertices[faces[f].vertices[2]]);
    VEC3F &v0 = _embeddedMesh->vertices[faces[f].vertices[0]];
    VEC3F &v1 = _embeddedMesh->vertices[faces[f].vertices[1]];
    VEC3F &v2 = _embeddedMesh->vertices[faces[f].vertices[2]];

    vector<VEC3F> m;
    VEC3F a = v1 - v0;
    VEC3F b = v2 - v1;
    VEC3F c = v0 - v2;
    VEC3F tempCross = cross(b,a);   unitize(tempCross);
    VEC3F m3 = tempCross;
    tempCross = cross(a,m3);    unitize(tempCross);
    m.push_back(tempCross);
    tempCross = cross(b,m3);    unitize(tempCross);
    m.push_back(tempCross);
    tempCross = cross(c,m3);    unitize(tempCross);
    m.push_back(tempCross);
    m.push_back(m3);

    float z;
    if (pose < 0 || pose >= 100)
    {
      z = 0.5f;
    } else {
      Real faceDelta = Df(f, pose);
      z = (float)faceDelta;
    }
         
    VEC3F start = (v0 + v1 + v2) / 3.0;
    start = start - 0.01 * m[3];
        
    //VEC3F n = -m3;    unitize(n);
    //VEC3F end = start + (((Real)z) * n);

    //glColor3f((z)*2.0f,(z)*2.0f,(z)*2.0f);
    
    VEC3F end;
    end = start - 0.02*m[3];    
    glColor3f(1.0f,0.0f,0.0f);
    glVertex3f((GLfloat)start[0], (GLfloat)start[1], (GLfloat)start[2]);
    glVertex3f((GLfloat)end[0], (GLfloat)end[1], (GLfloat)end[2]);

    end = start + 0.02*m[0];   
    glColor3f(0.0f,1.0f,0.0f); 
    glVertex3f((GLfloat)start[0], (GLfloat)start[1], (GLfloat)start[2]);
    glVertex3f((GLfloat)end[0], (GLfloat)end[1], (GLfloat)end[2]);

    end = start + 0.02*m[1];   
    glColor3f(0.0f,0.0f,1.0f); 
    glVertex3f((GLfloat)start[0], (GLfloat)start[1], (GLfloat)start[2]);
    glVertex3f((GLfloat)end[0], (GLfloat)end[1], (GLfloat)end[2]);

    end = start + 0.02*m[2];
    glColor3f(0.8f, 0.8f, 0.0f); 
    glVertex3f((GLfloat)start[0], (GLfloat)start[1], (GLfloat)start[2]);
    glVertex3f((GLfloat)end[0], (GLfloat)end[1], (GLfloat)end[2]);

    glColor3f(0.0f,0.0f,0.0f);
    start = v0 - 0.01*m[3];
    end = v1 - 0.01*m[3];
    glVertex3f((GLfloat)start[0], (GLfloat)start[1], (GLfloat)start[2]);
    glVertex3f((GLfloat)end[0], (GLfloat)end[1], (GLfloat)end[2]);

    start = v1 - 0.01*m[3];
    end = v2 - 0.01*m[3];
    glVertex3f((GLfloat)start[0], (GLfloat)start[1], (GLfloat)start[2]);
    glVertex3f((GLfloat)end[0], (GLfloat)end[1], (GLfloat)end[2]);

    start = v2 - 0.01*m[3];
    end = v0 - 0.01*m[3];
    glVertex3f((GLfloat)start[0], (GLfloat)start[1], (GLfloat)start[2]);
    glVertex3f((GLfloat)end[0], (GLfloat)end[1], (GLfloat)end[2]);

  }
  glEnd(); // */


}

//////////////////////////////////////////////////////////////////////
// Draw method to display an embedded mesh
//////////////////////////////////////////////////////////////////////
bool GLSL_TET_MESH::drawEmbedded()
{
  vector<OBJ::Face>& faces = _embeddedMesh->faces;
  if ( _embeddedMesh == NULL )
  {
    cout << "NULL mesh" << endl;
    return false;
  }
  int numFaces = faces.size();
  int vertSize = _embeddedMesh->vertices.size();

  Real delta = 0.25;
  
  vector<VEC3F>& baryEmbeddings = _tetMesh->barycentricEmbeddings();
  vector<int>& tetEmbeddings = _tetMesh->tetEmbeddings();
 
  _tetMesh->updateSurfaceMesh();
  //_tetMesh->drawSurfaceFaces();

  vector<TET>& tets = _tetMesh->tets();

  if ((int)_embeddedMesh->texcoords.size() != vertSize)
    _embeddedMesh->texcoords.resize(vertSize);
  
  for (int x = 0; x < vertSize; x++)
  {
    TET& tet = tets[tetEmbeddings[x]];
    VEC3F p = (baryEmbeddings[x][0] * (*tet.vertices[0]))
            + (baryEmbeddings[x][1] * (*tet.vertices[1]))
            + (baryEmbeddings[x][2] * (*tet.vertices[2]))
            + ((1.0 - baryEmbeddings[x][0] - baryEmbeddings[x][1]
                    - baryEmbeddings[x][2]) * (*tet.vertices[3]));
    _embeddedMesh->vertices[x] = p;
    _embeddedMesh->texcoords[x][0] = _embeddedMesh->texcoords[x][1] = 0.0;
  }

  _embeddedMesh->ComputeVertexNormals();

  for (int f = 0; f < numFaces; f++)
  {
    VEC3F &v0 = _embeddedMesh->vertices[faces[f].vertices[0]];
    VEC3F &v1 = _embeddedMesh->vertices[faces[f].vertices[1]];
    VEC3F &v2 = _embeddedMesh->vertices[faces[f].vertices[2]];
    TRIANGLE T(&v0,&v1,&v2);
    
    vector<VEC3F> m;                    // m appears to be calculated correctly
    VEC3F a = v1 - v0;
    VEC3F b = v2 - v1;
    VEC3F c = v0 - v2;
    VEC3F tempCross = cross(b,a);   unitize(tempCross);
    VEC3F m3 = tempCross;

    tempCross = cross(a,m3);    unitize(tempCross);
    m.push_back(delta*tempCross);

    tempCross = cross(b,m3);    unitize(tempCross);
    m.push_back(delta*tempCross);

    tempCross = cross(c,m3);    unitize(tempCross);
    m.push_back(delta*tempCross);

    m.push_back(delta*m3);

    // Calculate the total AO and distances of the points in range of T 
    /*for (int i = 0; i < vertSize; i++)
    {
      VEC3F &x = _embeddedMesh->vertices[i];

      //  if (x is NOT inside prism(m[4],T)) go to next i
      // isPointInVolume appears to be correct
      if (isPointInVolume(x, T, m))
      {
        VEC3F &n = _embeddedMesh->normals[i];                 
        Real AO = clamp(calculateAO(x, n, T), 0.0, 1.0); 
        Real dist = calculatePointTriDistance(x, T, m, delta);    
        AO *= (1.0 - (dist / delta));
        
        _embeddedMesh->texcoords[i][0] += AO;
      }
    }*/
  }

  /////////
  
  glBegin(GL_TRIANGLES);
  for (int f = 0; f < numFaces; f++)
  {
    float z;

    VEC3F &v0 = _embeddedMesh->vertices[faces[f].vertices[0]];
    VEC3F &v1 = _embeddedMesh->vertices[faces[f].vertices[1]];
    VEC3F &v2 = _embeddedMesh->vertices[faces[f].vertices[2]];

    //float z = _embeddedMesh->texcoords[faces[f].vertices[0]][0];
    glVertex3f((GLfloat)v0[0], (GLfloat)v0[1], (GLfloat)v0[2]);

    //z = _embeddedMesh->texcoords[faces[f].vertices[1]][0];
    glVertex3f((GLfloat)v1[0], (GLfloat)v1[1], (GLfloat)v1[2]);

    //z = _embeddedMesh->texcoords[faces[f].vertices[2]][0];
    glVertex3f((GLfloat)v2[0], (GLfloat)v2[1], (GLfloat)v2[2]);
  }
  glEnd();

  return true;
}

//////////////////////////////////////////////////////////////////////
// Numerically testing several 'known' values to make sure the functions are correct
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH::debugTestDeltaFunctions()
{
  VEC3F v0(0.0,0.0,0.0);
  VEC3F v1(2.0,0.0,0.0);
  VEC3F v2(0.0,2.0,0.0);

  TRIANGLE t(&v0,&v1,&v2);
  vector<VEC3F> m;
  VEC3F a = v1 - v0;
  VEC3F b = v2 - v1;
  VEC3F c = v0 - v2;
  VEC3F tx = cross(b,a);    unitize(tx);
  VEC3F m3 = tx;
  tx = cross(a,m3);   unitize(tx);
  m.push_back(tx);
  tx = cross(b,m3);   unitize(tx);
  m.push_back(tx);
  tx = cross(c,m3);   unitize(tx);
  m.push_back(tx);
  m.push_back(m3);
  
  VEC3F x;

  cout << " == TEST : isPointInVolume == " << endl;

  x[0] = x[1] = 0.667;
  x[2] = 0.0;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = 0.5;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = 0.99;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = 1.0;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = 1.01;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = -0.5;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;

  x[0] = x[1] = -1.0;
  x[2] = 0.0;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = 0.5;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = 0.99;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = 1.0;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = 1.01;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = -0.5;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;

  x[0] = x[1] = -1.01;
  x[2] = 0.0;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = 0.5;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = 0.99;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = 1.0;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = 1.01;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = -0.5;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;

  x[0] = x[1] = 2.0;
  x[2] = 0.0;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = 0.5;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = 0.99;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = 1.0;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = 1.01;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;
  x[2] = -0.5;  cout << " x = " << x << " -> ";
  if (isPointInVolume(x,t,m))    cout << "Yes" << endl;  else     cout << "No" << endl;

  cout << " == TEST : calculatePointTriDistance == " << endl;

  Real dist;
  x[0] = x[1] = 0.0;
  x[2] = 0.0; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;
  x[2] = 0.5; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;
  x[2] = 1.0; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;

  x[0] = x[1] = 1.0;
  x[2] = 0.0; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;
  x[2] = 0.5; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;
  x[2] = 1.0; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;

  x[0] = x[1] = -1.0;
  x[2] = 0.0; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;
  x[2] = 0.5; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;
  x[2] = 1.0; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;

  x[0] = 0.0; x[1] = 2.0;
  x[2] = 0.0; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;
  x[2] = 0.5; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;
  x[2] = 1.0; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;

  x[0] = 2.0; x[1] = 0.0;
  x[2] = 0.0; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;
  x[2] = 0.5; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;
  x[2] = 1.0; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;

  x[0] = -0.5; x[1] = 4.0;
  x[2] = 0.0; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;
  x[2] = 0.5; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;
  x[2] = 1.0; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;

  x[0] = 2.5; x[1] = -0.5;
  x[2] = 0.0; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;
  x[2] = 0.5; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;
  x[2] = 1.0; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;

  x[0] = x[1] = 1.0 + sqrt(0.125);
  x[2] = 0.0; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;
  x[2] = 0.5; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;
  x[2] = 1.0; cout << " x = " << x << " -> ";
  dist = calculatePointTriDistance(x, t, m, 1.0);
  cout << " dist: " << dist << endl;

  cout << " == TEST : calculateAO == " << endl;

  Real AO;
  VEC3F n(0.0,0.0,1.0);
  x[0] = x[1] = x[2] = 0.0;

  v0[0] = 1.0; v0[1] = 0.0; v0[2] = 0.0;
  v1[0] = 0.0; v1[1] = 1.0; v1[2] = 0.0;
  v2[0] = 0.0; v2[1] = 0.0; v2[2] = 1.0;
  AO = calculateAO(x, n, t, true);
  cout << "AO : " << AO << endl;

  v0[0] = 1.0; v0[1] = 0.0; v0[2] = 0.0;
  v1[0] = 0.0; v1[1] = 0.0; v1[2] = 1.0;
  v2[0] = 0.0; v2[1] = 1.0; v2[2] = 0.0;
  AO = calculateAO(x, n, t, true);
  cout << "AO : " << AO << endl;

  v0[0] = 2.0; v0[1] = 0.0; v0[2] = -1.0;
  v1[0] = 0.0; v1[1] = 0.0; v1[2] = 1.0;
  v2[0] = 0.0; v2[1] = 2.0; v2[2] = -1.0;
  AO = calculateAO(x, n, t, true);
  cout << "AO : " << AO << endl;
  cout << t;

  v0[0] = 1.0; v0[1] = 0.0; v0[2] = 0.0;
  v1[0] = 0.0; v1[1] = -0.5; v1[2] = 1.0;
  v2[0] = 0.0; v2[1] = 0.5; v2[2] = 1.0;
  AO = calculateAO(x, n, t, true);
  cout << "AO : " << AO << endl;
 
  v0[0] = 1.0; v0[1] = 0.0; v0[2] = -0.1;
  v1[0] = 0.0; v1[1] = -0.5; v1[2] = 1.0;
  v2[0] = 0.0; v2[1] = 0.5; v2[2] = 1.0;
  AO = calculateAO(x, n, t, true);
  cout << "AO : " << AO << endl;

  v0[0] = 1.0; v0[1] = 0.0; v0[2] = -0.2;
  v1[0] = 0.0; v1[1] = -0.5; v1[2] = 1.0;
  v2[0] = 0.0; v2[1] = 0.5; v2[2] = 1.0;
  AO = calculateAO(x, n, t, true);
  cout << "AO : " << AO << endl;

  v0[0] = 1.0; v0[1] = 0.0; v0[2] = -0.5;
  v1[0] = 0.0; v1[1] = -0.5; v1[2] = 1.0;
  v2[0] = 0.0; v2[1] = 0.5; v2[2] = 1.0;
  AO = calculateAO(x, n, t, true);
  cout << "AO : " << AO << endl;

  v0[0] = 1.0; v0[1] = 0.0; v0[2] = -1.0;
  v1[0] = 0.0; v1[1] = -0.5; v1[2] = 1.0;
  v2[0] = 0.0; v2[1] = 0.5; v2[2] = 1.0;
  AO = calculateAO(x, n, t, true);
  cout << "AO : " << AO << endl;
 
  v0[0] = 1.0; v0[1] = 0.0; v0[2] = -0.01;
  v1[0] = 0.0; v1[1] = 0.0; v1[2] = 1.0;
  v2[0] = 0.0; v2[1] = 1.0; v2[2] = 0.0;
  AO = calculateAO(x, n, t, true);
  cout << "AO : " << AO << endl;


}

//////////////////////////////////////////////////////////////////////
// Construct and initialize tetMesh vertex buffer object and related variables
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH::constructTetMeshVbo()
{
  // Copy rest pose vertex data to float array
  vector<VEC3F>& tetVerts = _tetMesh->vertices();
  _tetMeshSize = tetVerts.size();
  GLfloat *vertices = new GLfloat[_tetMeshSize * 3];

  if (vertices == NULL)
  {
    cout << "Could not allocate memory for vertex array" << endl;
    return;
  } 
  
  for (GLsizei x = 0; x < _tetMeshSize; x++)
  {
    /*deformX = 0.0f;
    deformY = 0.0f;
    deformZ = 0.0f;

    for (int y = 0; y < tempU.cols(); y++)
    {
      deformX += (GLfloat)(tempU.data()[(3*x)*tempU.cols() + y]);
      deformY += (GLfloat)(tempU.data()[(3*x+1)*tempU.cols() + y]);
      deformZ += (GLfloat)(tempU.data()[(3*x+2)*tempU.cols() + y]);
    }*/

    vertices[x*3]   = (GLfloat)tetVerts[x][0];// + mult*deformX;
    vertices[x*3+1] = (GLfloat)tetVerts[x][1];// + mult*deformY;
    vertices[x*3+2] = (GLfloat)tetVerts[x][2];// + mult*deformZ;
  }

  // Generate a buffer object for the mesh VBO
  glGenBuffers(1, &_tetMeshVboID);
  if (_tetMeshVboID == 0)
  {
    cout << "Could not generate a GL buffer object" << endl;
    delete vertices;
    return;
  }

  // Attempt to copy mesh data to GPU.
  glBindBuffer(GL_ARRAY_BUFFER, _tetMeshVboID);
  glBufferData(GL_ARRAY_BUFFER,
               _tetMeshSize * 3 * sizeof(GLfloat),
               (GLvoid*)vertices, GL_STATIC_DRAW);

  // Unbind mesh VBO
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  delete vertices;
  //vertices = NULL;
}

//////////////////////////////////////////////////////////////////////
// Copy the U basis matrix to a texture, and set the textures dimensions
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH::constructUBasisTex()
{
// Copy U basis matrix to float array
  MATRIX& U = _tetMesh->U();
  int Uwidth = U.cols();
  if ((Uwidth % 4) > 0)
    Uwidth += 4 - (Uwidth % 4);
  int Uheight = U.rows();
  int numEntries = Uwidth * Uheight;

  _UmatrixPaddedWidth = Uwidth;

  _UbasisTexWidth = (int)sqrt((float)numEntries);
  if ((_UbasisTexWidth % (3*Uwidth)) > 0)
    _UbasisTexWidth += ((3*Uwidth) - (_UbasisTexWidth % (3*Uwidth)));

  //cout << "Texwidth: " << texwidth << endl;
 
  _UbasisTexHeight = numEntries / _UbasisTexWidth;
  if ((numEntries % _UbasisTexWidth) > 0)
    _UbasisTexHeight += 1;

  int texelCount = _UbasisTexWidth * _UbasisTexHeight;

  GLfloat *texcoord = new GLfloat[(U.rows()/3)*2];
  if (texcoord == NULL)
  {
    cout << "Could not allocate memory for texcoord array" << endl;
  }

  GLfloat *Utex = new GLfloat[texelCount];
  if (Utex == NULL)
  {
    cout << "Could not allocate memory for matrix array" << endl;
    return;
  } 
 
  for (int x = 0; x < texelCount; x++)
  {
    Utex[x] = 0.0f;
  }

  int texrow = 0;
  int texcol = 0;
  for (int y = 0; y < U.rows(); y++)
  {
    if (y%3 == 0)
    {
      // Every third row in U corresponds to the start of the x-coord row of the
      // vertex. Construct texcoords to point to this.
      texcoord[(y/3)*2]   = ((GLfloat)(texcol) + 0.5f) / (GLfloat)_UbasisTexWidth;
      texcoord[(y/3)*2+1] = ((GLfloat)(texrow) + 0.5f) / (GLfloat)_UbasisTexHeight;

      /*if (y<60)
      {
        cout << "  TexCoord for " << (y/3) << ": ";
        cout <<  texcoord[(y/3)*2] << " " <<  texcoord[(y/3)*2+1] << endl;
      }*/
    }
    for (int x = 0; x < Uwidth; x++)
    {
      if (x < U.cols())
      {
        Utex[texrow*_UbasisTexWidth + texcol] = (GLfloat)(U.data()[y*U.cols() + x]);
      } else {
        Utex[texrow*_UbasisTexWidth + texcol] = 0.0f;
      }
      
      texcol++;
      if (texcol >= _UbasisTexWidth)
      {
        texrow++;
        texcol = 0;
      }
    }
  }

  /* // For debugging; clear U, and set texel 0, 5, 10 to red, green, blue (respectively)
  for (int x = 0; x < texelCount; x++)
  {
    Utex[x] = 0.0f;
  }

  for (int x = 0; x < (texelCount/60); x++)
  {
    Utex[x*60+0] = 1.0f;
    Utex[x*60+1] = 0.0f;
    Utex[x*60+2] = 0.0f;
    Utex[x*60+3] = 0.5f;

    Utex[x*60+20] = 0.0f;
    Utex[x*60+21] = 1.0f;
    Utex[x*60+22] = 0.0f;
    Utex[x*60+23] = 0.25f;

    Utex[x*60+40] = 0.0f;
    Utex[x*60+41] = 0.0f;
    Utex[x*60+42] = 1.0f;
    Utex[x*60+43] = 0.75f;
  }// */

  // Create GL texture object to hold matrix (assumes float and NPOT texture support)
  glGenTextures(1, &_UbasisTexID);
  if (_UbasisTexID == 0)
  {
    cout << "Could not allocate memory for texture object" << endl;
    delete Utex;
    delete texcoord;
    return;
  }

  glBindTexture(GL_TEXTURE_2D, _UbasisTexID);
   
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB,
               _UbasisTexWidth/4, _UbasisTexHeight,
               0, GL_RGBA, GL_FLOAT, Utex);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

  cout << " U matrix packed into " << _UbasisTexWidth << "x" << _UbasisTexHeight << endl;
		
  glBindTexture(GL_TEXTURE_2D, 0);
		
  delete Utex;
  //Utex = NULL;
    
  // Generate a buffer object for the mesh VBO
  glGenBuffers(1, &_UcoordVboID);
  if (_UcoordVboID == 0)
  {
    cout << "Could not generate a GL buffer object" << endl;
    delete texcoord;
    return;
  }

   // Attempt to copy texcoords to GPU.
  glBindBuffer(GL_ARRAY_BUFFER, _UcoordVboID);
  glBufferData(GL_ARRAY_BUFFER, (U.rows()/3) * 2 * sizeof(GLfloat),
               (GLvoid*)texcoord, GL_STATIC_DRAW);
		
  // Unbind VBO
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  delete texcoord;
  //texcoord = NULL;

}


//////////////////////////////////////////////////////////////////////
// Construct the VBOs and textures needed to render the given embedded obj
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH::constructEmbeddedMeshBuffers(OBJ *embeddedMesh)
{
  vector<OBJ::Face>& faces = embeddedMesh->faces;
  int numFaces = faces.size();
  _embeddedMeshIboSize = numFaces * 3;
  _embeddedVertSize = embeddedMesh->vertices.size();

  cout << " Embedded mesh has " << _embeddedVertSize;
  cout << " vertices and " << numFaces;
  cout << " faces." << endl;

  GLuint* faceIndices = new GLuint[_embeddedMeshIboSize];
  for (int y = 0; y < numFaces; y++)
  {
    OBJ::Face face = faces[y];
    faceIndices[y*3]   = face.vertices[0];
    faceIndices[y*3+1] = face.vertices[1];
    faceIndices[y*3+2] = face.vertices[2];
  }

  glGenBuffers(1, &_embeddedMeshIboID);
  if (_embeddedMeshIboID == 0)
  {
    cout << "Could not generate a GL buffer object" << endl;
    delete faceIndices;
    return;
  }

  // Attempt to copy indices to GPU.
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _embeddedMeshIboID);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
               _embeddedMeshIboSize * sizeof(GLuint),
               (GLvoid*)faceIndices, GL_STATIC_DRAW);
		
  // Unbind surface face IBO
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  delete faceIndices;

  vector<VEC3F>& baryEmbeddings = _tetMesh->barycentricEmbeddings();
  vector<int>& tetEmbeddings = _tetMesh->tetEmbeddings();
  vector<TET>& tets = _tetMesh->tets();
  
  GLfloat *barycoords = new GLfloat[_embeddedVertSize * 3];
  if (barycoords == NULL)
  {
    cout << "Error occurred allocating barycentric coords array" << endl;
    return;
  }
  GLfloat *normals = new GLfloat[_embeddedVertSize * 3];
  if (normals == NULL)
  {
    cout << "Error occurred allocating embedded normal array" << endl;
    delete barycoords;
    return;
  }
  GLfloat *vertLoc0 = new GLfloat[_embeddedVertSize * 4];
  if (vertLoc0 == NULL)
  {
    cout << "Error occurred allocating embedded vertloc array" << endl;
    delete barycoords;
    delete normals;
    return;
  }
  GLfloat *vertLoc1 = new GLfloat[_embeddedVertSize * 4];
  if (vertLoc1 == NULL)
  {
    cout << "Error occurred allocating embedded vertloc array" << endl;
    delete barycoords;
    delete normals;
    delete vertLoc0;
    return;
  }
  GLfloat *vertLoc2 = new GLfloat[_embeddedVertSize * 2];
  if (vertLoc2 == NULL)
  {
    cout << "Error occurred allocating embedded vertloc array" << endl;
    delete barycoords;
    delete normals;
    delete vertLoc0;
    delete vertLoc1;
    return;
  }

  embeddedMesh->ComputeVertexNormals();
  //embeddedMesh->SmoothVertexNormals();
  /*for (int x = 0; x < 10; x++)
    embeddedMesh->SmoothVertexNormals();// */

  int order[4];

  _tetVertexTexWidth = (GLsizei)(sqrt(_tetMeshSize)) + 1;
  _tetVertexTexHeight = _tetMeshSize / _tetVertexTexWidth;
  if ((_tetMeshSize % _tetVertexTexWidth) > 0)
    _tetVertexTexHeight += 1;

  _embeddedVertexTexWidth = (GLsizei)(sqrt(_embeddedVertSize)) + 1;
  _embeddedVertexTexHeight = _embeddedVertSize / _embeddedVertexTexWidth;
  if ((_embeddedVertSize % _embeddedVertexTexWidth) > 0)
    _embeddedVertexTexWidth += 1;

  for (int x = 0; x < _embeddedVertSize; x++)
  {
    TET& tet = tets[tetEmbeddings[x]];
    VEC3F vertnorm = embeddedMesh->normals[x];

    orderTet(tet, vertnorm, order);

    /*VEC3F& v0 = *(tet.vertices[order[0]]);
    VEC3F& v1 = *(tet.vertices[order[1]]);
    VEC3F& v2 = *(tet.vertices[order[2]]);
    //VEC3F& v3 = *(tet.vertices[order[3]]); */
    
    

    vertLoc0[x*4]   = (GLfloat)(((double)(_tetMesh->vertexID(tet.vertices[order[0]]) % _tetVertexTexWidth) 
                          + (0.5)) / (double)_tetVertexTexWidth); 
    vertLoc0[x*4+1] = (GLfloat)(((double)(_tetMesh->vertexID(tet.vertices[order[0]]) / _tetVertexTexWidth) 
                          + (0.5)) / (double)_tetVertexTexHeight); 
    
    vertLoc0[x*4+2] = (GLfloat)(((double)(_tetMesh->vertexID(tet.vertices[order[1]]) % _tetVertexTexWidth) 
                          + (0.5)) / (double)_tetVertexTexWidth); 
    vertLoc0[x*4+3] = (GLfloat)(((double)(_tetMesh->vertexID(tet.vertices[order[1]]) / _tetVertexTexWidth) 
                          + (0.5)) / (double)_tetVertexTexHeight); 
    
    vertLoc1[x*4]   = (GLfloat)(((double)(_tetMesh->vertexID(tet.vertices[order[2]]) % _tetVertexTexWidth) 
                          + (0.5)) / (double)_tetVertexTexWidth); 
    vertLoc1[x*4+1] = (GLfloat)(((double)(_tetMesh->vertexID(tet.vertices[order[2]]) / _tetVertexTexWidth) 
                          + (0.5)) / (double)_tetVertexTexHeight); 

    vertLoc1[x*4+2] = (GLfloat)(((double)(_tetMesh->vertexID(tet.vertices[order[3]]) % _tetVertexTexWidth) 
                          + (0.5)) / (double)_tetVertexTexWidth); 
    vertLoc1[x*4+3] = (GLfloat)(((double)(_tetMesh->vertexID(tet.vertices[order[3]]) / _tetVertexTexWidth) 
                          + (0.5)) / (double)_tetVertexTexHeight);

    vertLoc2[x*2]   = (GLfloat)(((double)(x % _embeddedVertexTexWidth) + (0.5)) / (double)_embeddedVertexTexWidth); 
    vertLoc2[x*2+1] = (GLfloat)(((double)(x / _embeddedVertexTexWidth) + (0.5)) / (double)_embeddedVertexTexHeight); 

    
    barycoords[x*3]   = (GLfloat)((order[0] == 3)
                  ? (1.0 - baryEmbeddings[x][0] - baryEmbeddings[x][1] - baryEmbeddings[x][2])
                  : baryEmbeddings[x][order[0]]);
    barycoords[x*3+1] = (GLfloat)((order[1] == 3)
                  ? (1.0 - baryEmbeddings[x][0] - baryEmbeddings[x][1] - baryEmbeddings[x][2])
                  : baryEmbeddings[x][order[1]]);
    barycoords[x*3+2] = (GLfloat)((order[2] == 3)
                  ? (1.0 - baryEmbeddings[x][0] - baryEmbeddings[x][1] - baryEmbeddings[x][2])
                  : baryEmbeddings[x][order[2]]);

    /*VEC3F norm(0.0,0.0,0.0);
    norm += (Real)(barycoords[x*3])   * tetNormals[_tetMesh->vertexID(tet.vertices[order[0]])];
    norm += (Real)(barycoords[x*3+1]) * tetNormals[_tetMesh->vertexID(tet.vertices[order[1]])];
    norm += (Real)(barycoords[x*3+2]) * tetNormals[_tetMesh->vertexID(tet.vertices[order[2]])];
    unitize(norm);
 
    VEC3F tan(0.0,0.0,0.0);
    tan += (Real)(barycoords[x*3])   * tetTangents[_tetMesh->vertexID(tet.vertices[order[0]])];
    tan += (Real)(barycoords[x*3+1]) * tetTangents[_tetMesh->vertexID(tet.vertices[order[1]])];
    tan += (Real)(barycoords[x*3+2]) * tetTangents[_tetMesh->vertexID(tet.vertices[order[2]])];
    unitize(tan);

    VEC3F bitan = cross(norm, tan);  unitize(bitan);

    MATRIX3 basis;
    basis(0,0) = tan[0];   basis(0,1) = tan[1];   basis(0,2) = tan[2];
    basis(1,0) = bitan[0]; basis(1,1) = bitan[1]; basis(1,2) = bitan[2];
    basis(2,0) = norm[0];  basis(2,1) = norm[1];  basis(2,2) = norm[2];
    MATRIX3 invertedBasis = basis.inverse();
    
    normals[x*3]   = invertedBasis(0,0) * vertnorm[0] + invertedBasis(0,1) * vertnorm[1] + invertedBasis(0,2) * vertnorm[2];
    normals[x*3+1] = invertedBasis(1,0) * vertnorm[0] + invertedBasis(1,1) * vertnorm[1] + invertedBasis(1,2) * vertnorm[2];
    normals[x*3+2] = invertedBasis(2,0) * vertnorm[0] + invertedBasis(2,1) * vertnorm[1] + invertedBasis(2,2) * vertnorm[2];
    // */

    normals[x*3]   = vertnorm[0];
    normals[x*3+1] = vertnorm[1];
    normals[x*3+2] = vertnorm[2];
  }

  cout << " got here 3 " << endl;
  /*delete tetNormals;
  delete tetTangents;*/

  // Generate a buffer object for the mesh VBO
  glGenBuffers(1, &_embeddedMeshVboID);
  if (_embeddedMeshVboID == 0)
  {
    cout << "Could not generate a GL buffer object" << endl;
    delete barycoords;
    delete normals;
    delete vertLoc0;
    delete vertLoc1;
    delete vertLoc2;
    return;
  }

  // Attempt to copy mesh data to GPU.
  glBindBuffer(GL_ARRAY_BUFFER, _embeddedMeshVboID);
  glBufferData(GL_ARRAY_BUFFER, _embeddedVertSize * 3 * sizeof(GLfloat),
               (GLvoid*)barycoords, GL_STATIC_DRAW);

  delete barycoords;
		
  // Generate a buffer object for the mesh VBO
  glGenBuffers(1, &_embeddedNormalVboID);
  if (_embeddedNormalVboID == 0)
  {
    cout << "Could not generate a GL buffer object" << endl;
    delete normals;
    delete vertLoc0;
    delete vertLoc1;
    delete vertLoc2;
    return;
  }

  // Attempt to copy mesh data to GPU.
  glBindBuffer(GL_ARRAY_BUFFER, _embeddedNormalVboID);
  glBufferData(GL_ARRAY_BUFFER, _embeddedVertSize * 3 * sizeof(GLfloat),
               (GLvoid*)normals, GL_STATIC_DRAW);
	
  delete normals;
  
  glGenBuffers(1, &_embeddedVertLoc0VboID);
  if (_embeddedVertLoc0VboID == 0)
  {
    cout << "Could not generate a GL buffer object" << endl;
    delete vertLoc0;
    delete vertLoc1;
    delete vertLoc2;
    return;
  }

  // Attempt to copy mesh data to GPU.
  glBindBuffer(GL_ARRAY_BUFFER, _embeddedVertLoc0VboID);
  glBufferData(GL_ARRAY_BUFFER, _embeddedVertSize * 4 * sizeof(GLfloat),
               (GLvoid*)vertLoc0, GL_STATIC_DRAW);

  delete vertLoc0;
	
  glGenBuffers(1, &_embeddedVertLoc1VboID);
  if (_embeddedVertLoc1VboID == 0)
  {
    cout << "Could not generate a GL buffer object" << endl;
    delete vertLoc1;
    delete vertLoc2;
    return;
  }

  // Attempt to copy mesh data to GPU.
  glBindBuffer(GL_ARRAY_BUFFER, _embeddedVertLoc1VboID);
  glBufferData(GL_ARRAY_BUFFER, _embeddedVertSize * 4 * sizeof(GLfloat),
               (GLvoid*)vertLoc1, GL_STATIC_DRAW);
	
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  delete vertLoc1;

  glGenBuffers(1, &_embeddedVertLoc2VboID);
  if (_embeddedVertLoc2VboID == 0)
  {
    cout << "Could not generate a GL buffer object" << endl;
    delete vertLoc2;
    return;
  }

  // Attempt to copy mesh data to GPU.
  glBindBuffer(GL_ARRAY_BUFFER, _embeddedVertLoc2VboID);
  glBufferData(GL_ARRAY_BUFFER, _embeddedVertSize * 2 * sizeof(GLfloat),
               (GLvoid*)vertLoc2, GL_STATIC_DRAW);
	
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  delete vertLoc2;

}

//////////////////////////////////////////////////////////////////////
// Construct the face info used to generate normal and tangent textures (int the instances)
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH::constructTetFaceInfo()
{
  GLfloat *vertexinfo = new GLfloat[_tetMeshSize * 3];

  if (vertexinfo == NULL)
  {
    cout << "Could not allocate memory for tex coord array" << endl;
    return;
  } 
 
  for (GLsizei x = 0; x < _tetMeshSize; x++)
  {
    // tex location in _tetVertexTexID of this vertex
    vertexinfo[x*3]   = (GLfloat)(x % _tetVertexTexWidth) / (GLfloat)_tetVertexTexWidth;
    vertexinfo[x*3+1] = (GLfloat)(x / _tetVertexTexWidth) / (GLfloat)_tetVertexTexHeight;
    // 1.0 if not constrained, 0.0 otherwise
    if ((int)x < _tetMesh->unconstrainedNodes())
      vertexinfo[x*3+2] = 1.0f;
    else 
      vertexinfo[x*3+2] = 0.0f;
  }

  // Generate a buffer object for the mesh VBO
  glGenBuffers(1, &_tetMeshInfoVboID);
  if (_tetMeshInfoVboID == 0)
  {
    cout << "Could not generate a GL buffer object" << endl;
    delete vertexinfo;
    return;
  }

  // Attempt to copy mesh data to GPU.
  glBindBuffer(GL_ARRAY_BUFFER, _tetMeshInfoVboID);
  glBufferData(GL_ARRAY_BUFFER, _tetMeshSize * 3 * sizeof(GLfloat),
               (GLvoid*)vertexinfo, GL_STATIC_DRAW);
		
  // Unbind mesh VBO
  glBindBuffer(GL_ARRAY_BUFFER, 0);
		
  delete vertexinfo;

  vector<pair<int,int> >& surfaceFaces = _tetMesh->surfaceFaces();
  vector<TET>& tets = _tetMesh->tets();

  map< VEC3F*, vector<int> > neighbors;

  vector<VEC3F*>& surfaceVertices = _tetMesh->surfaceVertices();
  
  for (int x = 0; x < (int)surfaceVertices.size(); x++)
  {
    VEC3F* vertex = surfaceVertices[x];
    vector<VEC3F*> oneRing;
    _tetMesh->oneRing(vertex, oneRing);
    
    for (int y = 0; y < (int)oneRing.size(); y++)
    {
      if (_tetMesh->isOnSurface(oneRing[y]))
      {
        neighbors[vertex].push_back(_tetMesh->vertexID(oneRing[y]));
      }
    }
  }

  /*vector<VEC3F>& allvertices = _tetMesh->vertices();
  
  for (int x = 0; x < allvertices.size(); x++)
  {
    VEC3F* vertex = &allvertices[x];
    vector<VEC3F*> oneRing;
    _tetMesh->oneRing(vertex, oneRing);
    
    for (int y = 0; y < oneRing.size(); y++)
    {
      //if (_tetMesh->isOnSurface(oneRing[y]))
      //{
        neighbors[vertex].push_back(_tetMesh->vertexID(oneRing[y]));
      //}
    }
  }*/

  _tetFaceSize = surfaceFaces.size();

  GLfloat *faceinfo0 = new GLfloat[_tetFaceSize * 4];
  GLfloat *faceinfo1 = new GLfloat[_tetFaceSize * 2];

  for (int x = 0; x < _tetFaceSize; x++)
  {
    TRIANGLE triangle = tets[surfaceFaces[x].first].face(surfaceFaces[x].second);
    
    VEC3F up;   up[0] = 0.0;  up[1] = 0.0;  up[2] = 1.0;
    VEC3F n = triangle.normal();
    VEC3F t = cross(up, n);   unitize(t);

    VEC3F& v0 = *(triangle.vertex(0));
    VEC3F& v1 = *(triangle.vertex(1));
    VEC3F& v2 = *(triangle.vertex(2));
    VEC3F a = v1 - v0;    unitize(a);
    VEC3F b = v2 - v1;    unitize(b);
    VEC3F c = v0 - v2;    unitize(c);

    Real dots[3];
    dots[0] = a[0]*t[0] + a[1]*t[1] + a[2]*t[2];
    dots[1] = b[0]*t[0] + b[1]*t[1] + b[2]*t[2];
    dots[2] = c[0]*t[0] + c[1]*t[1] + c[2]*t[2];

    int closest = 0;
    for (int y = 1; y < 3; y++)
    {
      if (dots[y] > dots[closest])
        closest = y;
    }
  
    faceinfo0[x*4]   = (GLfloat)((double)(_tetMesh->vertexID(triangle.vertex((0+closest)%3)) % _tetVertexTexWidth) / (double)_tetVertexTexWidth);
    faceinfo0[x*4+1] = (GLfloat)((double)(_tetMesh->vertexID(triangle.vertex((0+closest)%3)) / _tetVertexTexWidth) / (double)_tetVertexTexHeight);
    
    faceinfo0[x*4+2] = (GLfloat)((double)(_tetMesh->vertexID(triangle.vertex((1+closest)%3)) % _tetVertexTexWidth) / (double)_tetVertexTexWidth);
    faceinfo0[x*4+3] = (GLfloat)((double)(_tetMesh->vertexID(triangle.vertex((1+closest)%3)) / _tetVertexTexWidth) / (double)_tetVertexTexHeight);
    
    faceinfo1[x*2]   = (GLfloat)((double)(_tetMesh->vertexID(triangle.vertex((2+closest)%3)) % _tetVertexTexWidth) / (double)_tetVertexTexWidth);
    faceinfo1[x*2+1] = (GLfloat)((double)(_tetMesh->vertexID(triangle.vertex((2+closest)%3)) / _tetVertexTexWidth) / (double)_tetVertexTexHeight);
 
  }

  _surfaceVertSize = neighbors.size();
  GLfloat *surfvert = new GLfloat[_surfaceVertSize * 2];
  GLfloat *surfneighbors0 = new GLfloat[_surfaceVertSize * 4];
  GLfloat *surfneighbors1 = new GLfloat[_surfaceVertSize * 4];

  if (surfvert == NULL)  {
    cout << "Could not allocate memory for surface vertex tex coord array" << endl;
    return;
  }
  if (surfneighbors0 == NULL)  {
    cout << "Could not allocate memory for surface vertex neighbor tex coord array" << endl;
    delete surfvert;
    return;
  }
  if (surfneighbors1 == NULL)  {
    cout << "Could not allocate memory for surface vertex neighbor tex coord array" << endl;
    delete surfvert;
    delete surfneighbors0;
    return;
  } // */

  map< VEC3F*, vector<int> >::iterator it;
  int i = 0;
  for (it = neighbors.begin(); it != neighbors.end(); it++, i++)
  {
    VEC3F *vertex = it->first;
    int order[4];
    orderNeighbors(vertex, it->second, order);

    /*if (it->second.size() > 3)
    {
      order[0] = it->second[0];
      order[1] = it->second[1];
      order[2] = it->second[2];
      order[3] = it->second[3];
    } else {
      order[0] = it->second[0];
      order[1] = it->second[0];
      order[2] = it->second[1];
      order[3] = it->second[2];
    }*/

    surfvert[i*2]   = (GLfloat)((double)(_tetMesh->vertexID(vertex) % _tetVertexTexWidth) / (double)_tetVertexTexHeight);
    surfvert[i*2+1] = (GLfloat)((double)(_tetMesh->vertexID(vertex) / _tetVertexTexWidth) / (double)_tetVertexTexHeight);
    
    surfneighbors0[i*4]   = (GLfloat)((double)(order[0] % _tetVertexTexWidth) / (double)_tetVertexTexHeight);
    surfneighbors0[i*4+1] = (GLfloat)((double)(order[0] / _tetVertexTexWidth) / (double)_tetVertexTexHeight);

    surfneighbors0[i*4+2] = (GLfloat)((double)(order[1] % _tetVertexTexWidth) / (double)_tetVertexTexHeight);
    surfneighbors0[i*4+3] = (GLfloat)((double)(order[1] / _tetVertexTexWidth) / (double)_tetVertexTexHeight);

    surfneighbors1[i*4]   = (GLfloat)((double)(order[2] % _tetVertexTexWidth) / (double)_tetVertexTexHeight);
    surfneighbors1[i*4+1] = (GLfloat)((double)(order[2] / _tetVertexTexWidth) / (double)_tetVertexTexHeight);

    surfneighbors1[i*4+2] = (GLfloat)((double)(order[3] % _tetVertexTexWidth) / (double)_tetVertexTexHeight);
    surfneighbors1[i*4+3] = (GLfloat)((double)(order[3] / _tetVertexTexWidth) / (double)_tetVertexTexHeight);
    // */
  }

  // Generate a buffer object for the VBO
  glGenBuffers(1, &_surfaceVertVboID);
  glGenBuffers(1, &_surfaceVertNeighbor0VboID);
  glGenBuffers(1, &_surfaceVertNeighbor1VboID);

  if (_surfaceVertVboID == 0 || _surfaceVertNeighbor0VboID == 0
   || _surfaceVertNeighbor1VboID == 0)
  {
    cout << "Could not generate a GL buffer object" << endl;
    delete surfvert;
    delete surfneighbors0;
    delete surfneighbors1;
    return;
  }

  cout << " Surface Vert size is " << _surfaceVertSize << endl;
 
  glBindBuffer(GL_ARRAY_BUFFER, _surfaceVertVboID);
  glBufferData(GL_ARRAY_BUFFER, _surfaceVertSize * 2 * sizeof(GLfloat),
               (GLvoid*)surfvert, GL_STATIC_DRAW);
  
  glBindBuffer(GL_ARRAY_BUFFER, _surfaceVertNeighbor0VboID);
  glBufferData(GL_ARRAY_BUFFER, _surfaceVertSize * 4 * sizeof(GLfloat),
               (GLvoid*)surfneighbors0, GL_STATIC_DRAW);
  
  glBindBuffer(GL_ARRAY_BUFFER, _surfaceVertNeighbor1VboID);
  glBufferData(GL_ARRAY_BUFFER, _surfaceVertSize * 4 * sizeof(GLfloat),
               (GLvoid*)surfneighbors1, GL_STATIC_DRAW);

  delete surfvert;
  delete surfneighbors0;
  delete surfneighbors1;

  genRestFrame();
  
  // Generate a buffer object for the VBO
  glGenBuffers(1, &_tetFaceInfo0VboID);
  if (_tetFaceInfo0VboID == 0)
  {
    cout << "Could not generate a GL buffer object" << endl;
    delete faceinfo0;
    delete faceinfo1;
    return;
  }

  // Attempt to copy data to GPU.
  glBindBuffer(GL_ARRAY_BUFFER, _tetFaceInfo0VboID);
  glBufferData(GL_ARRAY_BUFFER, _tetFaceSize * 4 * sizeof(GLfloat),
               (GLvoid*)faceinfo0, GL_STATIC_DRAW);

  delete faceinfo0;

  glGenBuffers(1, &_tetFaceInfo1VboID);
  if (_tetFaceInfo1VboID == 0)
  {
    cout << "Could not generate a GL buffer object" << endl;
    delete faceinfo1;
    return;
  }

  // Attempt to copy data to GPU.
  glBindBuffer(GL_ARRAY_BUFFER, _tetFaceInfo1VboID);
  glBufferData(GL_ARRAY_BUFFER, _tetFaceSize * 2 * sizeof(GLfloat),
               (GLvoid*)faceinfo1, GL_STATIC_DRAW);
		
  // Unbind VBO
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  delete faceinfo1;
}

//////////////////////////////////////////////////////////////////////
// Given a TET, stores in order[4] the ordered indices corresponding to
// the vertices in tet, such that cross((v1-v0),(v2-v0)) is the outward
// normal of the tet
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH::orderTet(TET& tet, VEC3F& normal, int *order)
{
  //TET& tet = tets[tetEmbeddings[x]];
  VEC3F& v0 = *(tet.vertices[0]);
  VEC3F& v1 = *(tet.vertices[1]);
  VEC3F& v2 = *(tet.vertices[2]);
  VEC3F& v3 = *(tet.vertices[3]);

  double dots[4];
  unitize(normal);

  VEC3F& n = tet.face(0).normal();
  dots[0] = (double)(n[0]*normal[0] + n[1]*normal[1] + n[2]*normal[2]);

  n = tet.face(1).normal();
  dots[1] = (double)(n[0]*normal[0] + n[1]*normal[1] + n[2]*normal[2]);

  n = tet.face(2).normal();
  dots[2] = (double)(n[0]*normal[0] + n[1]*normal[1] + n[2]*normal[2]);

  n = tet.face(3).normal();
  dots[3] = (double)(n[0]*normal[0] + n[1]*normal[1] + n[2]*normal[2]);

  int largest = 0;
  for (int x = 1; x < 4; x++)
  {
    if (dots[x] > dots[largest])
      largest = x;
  }

  VEC3F a, b, c;
  switch (largest) {
    case 0:
      a = v1 - v0;    unitize(a);
      b = v3 - v1;    unitize(b);
      c = v0 - v3;    unitize(c);
      break;
    case 1:
      a = v2 - v0;    unitize(a);
      b = v1 - v2;    unitize(b);
      c = v0 - v1;    unitize(c);
      break;
    case 2:
      a = v2 - v3;    unitize(a);
      b = v0 - v2;    unitize(b);
      c = v3 - v0;    unitize(c);
      break;
    case 3:
      a = v2 - v1;    unitize(a);
      b = v3 - v2;    unitize(b);
      c = v1 - v3;    unitize(c);
      break;
    default:
      cout << " The index variable 'largest' does not have a valid value (in orderTet)" << endl;
      return;
      break;
  }

  VEC3F up;   up[0] = 0.0;  up[1] = 0.0;  up[2] = 1.0;
  n = tet.face(largest).normal();
  VEC3F t = cross(up, n);   unitize(t);

  dots[0] = a[0]*t[0] + a[1]*t[1] + a[2]*t[2];
  dots[1] = b[0]*t[0] + b[1]*t[1] + b[2]*t[2];
  dots[2] = c[0]*t[0] + c[1]*t[1] + c[2]*t[2];

  int closest = 0;
  for (int x = 1; x < 3; x++)
  {
    if (dots[x] > dots[closest])
      closest = x;
  }

  if (closest == 0) // if (a[2] <= b[2] && a[2] <= c[2])
  {
    switch (largest) {
      case 0:
        order[0] = 0;  order[1] = 1;  order[2] = 3;  order[3] = 2;
        break;
      case 1:
        order[0] = 0;  order[1] = 2;  order[2] = 1;  order[3] = 3;
        break;
      case 2:
        order[0] = 3;  order[1] = 2;  order[2] = 0;  order[3] = 1;
        break;
      case 3:
        order[0] = 1;  order[1] = 2;  order[2] = 3;  order[3] = 0;
        break;
      default:
        // Should never get here
        break;
    }
  }
  else if (closest == 1) // if (b[2] <= a[2] && b[2] <= c[2])
  {
    switch (largest) {
      case 0:
        order[0] = 1;  order[1] = 3;  order[2] = 0;  order[3] = 2;
        break;
      case 1:
        order[0] = 2;  order[1] = 1;  order[2] = 0;  order[3] = 3;
        break;
      case 2:
        order[0] = 2;  order[1] = 0;  order[2] = 3;  order[3] = 1;
        break;
      case 3:
        order[0] = 2;  order[1] = 3;  order[2] = 1;  order[3] = 0;
        break;
      default:
        // Should never get here
        break;
    }
  } 
  else  // //if (c[2] <= a[2] && c[2] <= b[2])
  {
    switch (largest) {
      case 0:
        order[0] = 3;  order[1] = 0;  order[2] = 1;  order[3] = 2;
        break;
      case 1:
        order[0] = 1;  order[1] = 0;  order[2] = 2;  order[3] = 3;
        break;
      case 2:
        order[0] = 0;  order[1] = 3;  order[2] = 2;  order[3] = 1;
        break;
      case 3:
        order[0] = 3;  order[1] = 1;  order[2] = 2;  order[3] = 0;
        break;
      default:
        // Should never get here
        break;
    }
  }

  //cout << " Order: (" << order[0] << "," << order[1] << "," << order[2] << "," << order[3] << ")" << endl;
}


///////
bool compdist( pair< pair<int,int>, float > a,  pair< pair<int,int>, float > b )
{
  return (a.second < b.second);
}
bool compscore( pair< vector<int>, float> a,  pair< vector<int>, float> b )
{
  // want largest
  return (a.second > b.second);
}

//////////////////////////////////////////////////////////////////////
// Given array of neighbor indices, and a center vertex, generate an order of
// 4 s.t. 0,1,2 are close together, 3 is ~90 deg away from them.
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH::orderNeighbors(VEC3F* center, vector<int>& neighbors, int *order)
{
  int n = neighbors.size();
  //cout << " " << n;
  //int numPairs = (n * (n - 1)) / 2;
  vector< pair<VEC3F,int> > vecgrpList;
  vector< pair< pair<int,int>, float > > distances;
  
  // Only two neighbors; grab one and duplicate it, stick other one as 'ortho' vector
  if (n == 2)
  {
    cout << endl << " Found vertex with only 2 neighbors!" << endl;
    order[0] = order[1] = order[2] = neighbors[0];
    order[3] = neighbors[1];
    return;
  }

 // float xz_dist;//   = sqrt( x*x + z*z )
 // float latitude;//  = atan2( xz_dist, y )
 // float longitude;// = atan2( x, z )

  for (int x = 0; x < n; x++)
  {
    VEC3F v = (*_tetMesh->vertices(neighbors[x])) - (*center);
    unitize(v);
    pair<VEC3F,int> newpair(v,x);
    vecgrpList.push_back(newpair);
   // xz_dist = sqrt((float)(v[0]*v[0] + v[2]*v[2]));
   // latitude = atan2((float)xz_dist, (float)v[1]);
   // longitude = (atan2((float)v[0], (float)v[2])*sin(latitude)) / 3.14159f;
   // latitude = latitude / 3.14159f;
    //cout << " > " << x << ": (" << v[0] << "," << v[1] << "," << v[2] << ")" << endl;
  }

  for (int x = 0; x < n; x++)
  {
    for (int y = x+1; y < n; y++)
    {
      pair<int,int> p(x,y);
      VEC3F& v0 = (vecgrpList[x].first);
      VEC3F& v1 = (vecgrpList[y].first);
      float d = 0.5f * (1.0f - (v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2]));
      pair< pair<int,int>, float > newdist(p,d);
      distances.push_back(newdist);
    }
  }

  sort (distances.begin(), distances.end(), compdist);

  // If only three neighbors, duplicate one, and order
  if (n == 3)
  {
    // First two distance comparisons will have a common vector between them. Duplicate the other one.
    pair<int,int>& p1 = distances[0].first;
    pair<int,int>& p2 = distances[1].first;
    if (p1.first == p2.first || p1.first == p2.second)
    {
      order[0] = order[1] = neighbors[p1.second];
      order[2] = neighbors[p1.first];
    } else {
      order[0] = order[1] = neighbors[p1.first];
      order[2] = neighbors[p1.second];
    }
    // The uncommon vector from p2 is the 'ortho' vector
    if (p2.first == p1.first || p2.first == p1.second)
    {
      order[3] = neighbors[p2.second];
    } else {
      order[3] = neighbors[p2.first];
    }

    return;
  }

  if (n > 4)
  {
    order[0] = neighbors[1];
    order[1] = neighbors[2];
    order[2] = neighbors[2];
    order[3] = neighbors[4]; // */
  } else {
    order[0] = neighbors[1];
    order[1] = neighbors[2];
    order[2] = neighbors[2];
    order[3] = neighbors[3]; // */
  }

  /*vector< pair< vector<int>, float> > scoreList;
  //cout << " > n == " << n << endl;

  for (int a = 0; a < n; a++)
  {
    for (int b = a+1; b < n; b++)
    {
      for (int c = b+1; c < n; c++)
      {
        for (int d = c+1; d < n; d++)
        {
          //unsigned int id = (1 << a) | (1 << b) | (1 << c) | (1 << d);
          VEC3F& va = (vecgrpList[a].first);
          VEC3F& vb = (vecgrpList[b].first);
          VEC3F& vc = (vecgrpList[c].first);
          VEC3F& vd = (vecgrpList[d].first);

          float score = 4.0f - (float)(abs(va*vb) + abs(va*vc) + abs(vb*vd) + abs(vc*vd) + va*vd + vb*vc);
          
          //cout << " id: " << id << "  score: " << score << endl;

          vector<int> indices;
          indices.push_back(vecgrpList[a].second);
          indices.push_back(vecgrpList[b].second);
          indices.push_back(vecgrpList[c].second);
          indices.push_back(vecgrpList[d].second);

          pair< vector<int>, float> scorePair(indices,score);
          
          scoreList.push_back(scorePair);
        } 
      }
    }
  }
  //cout <<  " -- < ";
 
  sort (scoreList.begin(), scoreList.end(), compscore);

  vector<int>& best = scoreList[0].first;

  //cout << " Best has score of " << scoreList[0].second << ", " << neighbors[best[0]] << ", " << neighbors[best[1]] << ", " << neighbors[best[2]] << ", " << neighbors[best[3]] << endl;

  order[0] = neighbors[best[0]];
  order[1] = neighbors[best[1]];
  order[2] = neighbors[best[2]];
  order[3] = neighbors[best[3]];// */
    
  //cout << " > ";
  /*bool stop = false;
  vector< pair< pair<int,int>, float > >::iterator it;
  vector<int> groupCount;
  groupCount.resize(n);
  for (int x = 0; x < n; x++)
    groupCount[x] = 1;
  vector<int> cluster;

  //for (it = distances.begin(); it != distances.end(); it++)
  //{
  //  cout << " > " << it->first.first << "," << it->first.second << ": " << it->second << endl;
  //}

  for (it = distances.begin(); it != distances.end() && !stop; it++)
  {
    //cout << " > VecGroup List: ";
    //for (int x = 0; x < n; x++)
    //{
    //  cout << vecgrpList[x].second << " ";
    //}
    //cout << endl;

    int grp1 = vecgrpList[it->first.first].second;
    int grp2 = vecgrpList[it->first.second].second;

    //cout << " Looking at " << grp1 << " and " << grp2 << ":";

    if (grp1 == grp2)
    {
      cout << " Found a cycle; this should not happen." << endl;
      stop = true;      
    } 
    else
    {
      //cout << " Union; ";

      // grp1 and grp2 are connected, so take the smaller grp number and use that for all of them
      int newgrp = (grp1 < grp2) ? grp1 : grp2;
      int change = (grp1 < grp2) ? grp2 : grp1;
      //cout << "make all " << change << " into " << newgrp << endl;

      //cout << " G" << newgrp << " > " << groupCount[newgrp] << " & ";
      //cout << "G" << change << " > " << groupCount[change] << " => ";
      
      groupCount[newgrp] += groupCount[change];
      groupCount[change] = 0;

      //cout << "G" << newgrp << " > " << groupCount[newgrp] << endl;

      for (int x = 0; x < n; x++)
      {
        if (vecgrpList[x].second == change)
          vecgrpList[x].second = newgrp;
      }

      // Three in a group means we found three vectors that are near each other.
      if (groupCount[newgrp] == 3)
      {
        //cout << " Found three vectors near each other: ";
        stop = true;

        for (int x = 0; x < n; x++)
        {
          if (vecgrpList[x].second == newgrp)
          {
            //cout << x << " ";
            cluster.push_back(x);
          }
        }
        //cout << endl;
      } 
      else if (groupCount[newgrp] == 4)
      {
        //cout << " Found four vectors near each other: ";
        stop = true;

        // Four near each other means that two pairs were joined. We want the three closest, so ideally
        // we'd like pair 1, which is closer, and then the nearest vector from pair 2. As a compromise,
        // we take the middle two (shared edge between clusters) and then a vector from one of the original
        // pairs. They're in "numerical" order, which doesn't mean anything, so this is effectively random.
        for (int x = 0; x < n; x++)
        {
          if (vecgrpList[x].second == newgrp)
          {
            if (x == it->first.first || x == it->first.second)
            {
              //cout << x << " ";
              cluster.insert(cluster.begin(), x);
            }
            else
            {
              //cout << x << " ";
              cluster.push_back(x);
            }
          }
        }
        // print out cluster (for debugging)
        //for (int x = 0; x < 4; x++) cout << cluster[x] << " ";
        //cout << endl;
      }
    }
  }

  // This shouldn't happen. Cluster size of < 3 should have been caught before now.
  if (cluster.size() < 3)
  {
    cout << " huh? " << endl << endl;
    return;
  }

  VEC3F avg(0,0,0);
  //cout << "using " << vecgrpList[cluster[0]].first << ", ";
  //cout << vecgrpList[cluster[1]].first << ", ";
  //cout << vecgrpList[cluster[2]].first << endl;

  avg = vecgrpList[cluster[0]].first + vecgrpList[cluster[1]].first + vecgrpList[cluster[2]].first;
  //cout << " avg: " << avg << endl;
  unitize(avg);

  distances.clear();

  for (int x = 3; x < cluster.size(); x++)
  {
    pair<int,int> p(-1,cluster[x]);
    VEC3F& v0 = (vecgrpList[cluster[x]].first);
    float d = abs(v0[0]*avg[0] + v0[1]*avg[1] + v0[2]*avg[2]);
    pair< pair<int,int>, float > newdist(p,d);
    distances.push_back(newdist);
  }

  sort (distances.begin(), distances.end(), comp);

  //for (it = distances.begin(); it != distances.end(); it++)
  //{
  //  cout << " > " << it->first.first << "," << it->first.second << ": " << it->second << endl;
  //}

  //cout << " Avg: (" << avg[0] << "," << avg[1] << "," << avg[2] << ")" << endl;

  // First three are the first three clustered vectors
  order[0] = neighbors[cluster[0]];
  order[1] = neighbors[cluster[1]];  
  order[2] = neighbors[cluster[2]];

  //VEC3F& orth = (vecgrpList[distances[0].first.second].first);

  //cout << " Ortho: (" << orth[0] << "," << orth[1] << "," << orth[2] << ")" << endl;

  //float AdotO = avg[0]*orth[0] + avg[1]*orth[1] + avg[2]*orth[2];
  //cout << " Avg dot Ortho: " << AdotO << endl << endl;

  // Third 'ortho' vector is the neighbor that was closest to 90 degrees from the avg vector
  order[3] = neighbors[distances[0].first.second]; // */
}


//////////////////////////////////////////////////////////////////////
// Generate n1 and n2 frame textures for the rest pose of the mesh.
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH::genRestFrame()
{
  // Temporary texture to hold rest positions of all vertices (same as _tetVertexTexID in the INSTANCE 
  // class, only this one is not deformed). Needed as input to the next pass.
  GLuint restPosTex;
  glGenTextures(1, &restPosTex);
  if (restPosTex == 0)
  {
    cout << "Could not allocate temporary memory for rest pos texture object" << endl;
    return;
  }
  glBindTexture(GL_TEXTURE_2D, restPosTex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, _tetVertexTexWidth,
               _tetVertexTexHeight, 0, GL_RGBA, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

  // Create texture objects to hold rest frames
  glGenTextures(1, &_n1RestTexID);
  if (_n1RestTexID == 0)
  {
    cout << "Could not allocate memory for n1Rest texture object" << endl;
    return;
  }
  glBindTexture(GL_TEXTURE_2D, _n1RestTexID);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, _tetVertexTexWidth,
               _tetVertexTexHeight, 0, GL_RGBA, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

  glGenTextures(1, &_n2RestTexID);
  if (_n2RestTexID == 0)
  {
    cout << "Could not allocate memory for n2Rest texture object" << endl;
    return;
  }  
  glBindTexture(GL_TEXTURE_2D, _n2RestTexID);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, _tetVertexTexWidth,
               _tetVertexTexHeight, 0, GL_RGBA, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  
  glBindTexture(GL_TEXTURE_2D, 0);

  // Make a temporary framebuffer for the texture generation
  GLuint restFrameFbo;
  glGenFramebuffersEXT(1, &restFrameFbo);
  if (restFrameFbo == 0)
  {
    cout << " Could not generate tet vertex framebuffer object." << endl;
    return;
  }

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, restFrameFbo);
  GLenum buf[2] = { GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT };
  glDrawBuffers(2, buf);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
                            GL_TEXTURE_2D, restPosTex, 0);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1_EXT,
                            GL_TEXTURE_2D, _n2RestTexID, 0);
 
  GLenum fboStatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  if (fboStatus != GL_FRAMEBUFFER_COMPLETE_EXT)
  {
    cout << " Framebuffer status is not complete; an error occurred at ";
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
    return; 
  }

  // Make shader to create temporary rest position texture
  GLSL_SHADER *genRestPosShader = new GLSL_SHADER();
  genRestPosShader->attachVert("src/rendering/glsl/tetrestpos.vert");
  genRestPosShader->attachFrag("src/rendering/glsl/tetrestpos.frag");
  genRestPosShader->compile();

  // Set up drawing to texture  
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);

  GLint windowViewport[4];
  glGetIntegerv(GL_VIEWPORT, windowViewport);

  glMatrixMode(GL_MODELVIEW);
  glViewport(0, 0, _tetVertexTexWidth, _tetVertexTexHeight);

  glClear(GL_COLOR_BUFFER_BIT);

  glPointSize(1.0);
  
  glBindBuffer(GL_ARRAY_BUFFER, _tetMeshVboID);
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE0);
  glBindBuffer(GL_ARRAY_BUFFER, _tetMeshInfoVboID);
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(3, GL_FLOAT, 0, 0);

  glUseProgram(genRestPosShader->getHandle());

  // Create temporary rest pose texture
  glDrawArrays(GL_POINTS, 0, _tetMeshSize);

  delete genRestPosShader;

  // Now that we have the rest pose texture, we can use it to generate the rest frame textures
  // that are used later.

  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
                            GL_TEXTURE_2D, _n1RestTexID, 0);

  GLSL_SHADER *genRestRotShader = new GLSL_SHADER();
  genRestRotShader->attachVert("src/rendering/glsl/tetrestrot.vert");
  genRestRotShader->attachFrag("src/rendering/glsl/tetrestrot.frag");
  genRestRotShader->compile();

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, restPosTex);

  glClear(GL_COLOR_BUFFER_BIT);

  glBindBuffer(GL_ARRAY_BUFFER, _surfaceVertVboID);
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(2, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE0);
  glBindBuffer(GL_ARRAY_BUFFER, _surfaceVertNeighbor0VboID);
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(4, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE1);
  glBindBuffer(GL_ARRAY_BUFFER, _surfaceVertNeighbor1VboID);
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(4, GL_FLOAT, 0, 0);

  glUseProgram(genRestRotShader->getHandle());
  glUniform1i(genRestRotShader->uniformLoc("tetMeshTex"), 0);
  glUniform2f(genRestRotShader->uniformLoc("offset"),
    (0.5f / (GLfloat)_tetVertexTexWidth), (0.5f / (GLfloat)_tetVertexTexHeight));

  glDrawArrays(GL_POINTS, 0, _surfaceVertSize);

  // Clean up all the temporary objects we created.
  glUseProgram(0);

  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glClientActiveTexture(GL_TEXTURE0);
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);

  glBindTexture(GL_TEXTURE_2D, 0);

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
  glDeleteFramebuffersEXT(1, &restFrameFbo);
  glDeleteTextures(1, &restPosTex);

  delete genRestRotShader;

  glViewport(windowViewport[0],windowViewport[1],windowViewport[2],windowViewport[3]);

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

//////////////////////////////////////////////////////////////////////
// Copy tetVertexTex width and height to the given pointers
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH::getVertexTexDims(GLsizei *width, GLsizei *height,
                                     GLsizei *embeddedwidth,
                                     GLsizei *embeddedheight)
{
  (*width) = _tetVertexTexWidth;
  (*height) = _tetVertexTexHeight;
  (*embeddedwidth) = _embeddedVertexTexWidth;
  (*embeddedheight) = _embeddedVertexTexHeight;
}

//////////////////////////////////////////////////////////////////////
// Simple accessors, with asserts
//////////////////////////////////////////////////////////////////////
SUBSPACE_TET_MESH *GLSL_TET_MESH::tetMesh()      
{ 
  assert(_tetMesh != NULL);
  return _tetMesh;    
}

GLSL_SHADER *GLSL_TET_MESH::tetVertexShader()
{
  assert(_tetVertexShader != NULL);
  return _tetVertexShader;  
}

GLSL_SHADER *GLSL_TET_MESH::genRotationShader()  
{ 
  assert(_genRotationShader != NULL);
  return _genRotationShader;  
}
GLSL_SHADER *GLSL_TET_MESH::drawBasesShader()
{
  assert(_drawBasesShader != NULL);
  return _drawBasesShader;  
}
GLSL_SHADER *GLSL_TET_MESH::drawFacesShader()
{
  assert(_drawFacesShader != NULL);
  return _drawFacesShader;  
}
GLSL_SHADER *GLSL_TET_MESH::drawNeighborsShader()
{
  assert(_drawNeighborsShader != NULL);
  return _drawNeighborsShader;  
}
GLSL_SHADER *GLSL_TET_MESH::drawAOVShader()
{
  assert(_drawAOVShader != NULL);
  return _drawAOVShader;  
}
GLSL_SHADER *GLSL_TET_MESH::drawAOV2Shader()
{
  assert(_drawAOV2Shader != NULL);
  return _drawAOV2Shader;  
}
GLSL_SHADER *GLSL_TET_MESH::drawAOV3Shader()
{
  assert(_drawAOV3Shader != NULL);
  return _drawAOV3Shader;  
}
GLSL_SHADER *GLSL_TET_MESH::renderFromTextureTest()
{
  assert(_renderFromTextureTest != NULL);
  return _renderFromTextureTest;  
}
GLuint GLSL_TET_MESH::tetMeshVboID()         
{
  assert(_tetMeshVboID != 0);
  return _tetMeshVboID; 
}
GLuint GLSL_TET_MESH::tetMeshInfoVboID()     
{ 
  assert(_tetMeshInfoVboID != 0);
  return _tetMeshInfoVboID;  
}
GLsizei GLSL_TET_MESH::tetMeshSize()         
{ 
  return _tetMeshSize;  
}

GLuint GLSL_TET_MESH::tetFaceInfo0VboID()    
{
  assert(_tetFaceInfo0VboID != 0);
  return _tetFaceInfo0VboID; 
}
GLuint GLSL_TET_MESH::tetFaceInfo1VboID()  
{ 
  assert(_tetFaceInfo1VboID != 0);
  return _tetFaceInfo1VboID;
}
GLsizei GLSL_TET_MESH::tetFaceSize()        
{
  return _tetFaceSize; 
}

GLuint GLSL_TET_MESH::embeddedMeshVboID()     
{
  assert(_embeddedMeshVboID != 0);
  return _embeddedMeshVboID;
}
GLuint GLSL_TET_MESH::embeddedNormalVboID()  
{ 
  assert(_embeddedNormalVboID != 0);
  return _embeddedNormalVboID; 
}
GLuint GLSL_TET_MESH::embeddedVertLoc0VboID()  
{
  assert(_embeddedVertLoc0VboID != 0);
  return _embeddedVertLoc0VboID; 
}
GLuint GLSL_TET_MESH::embeddedVertLoc1VboID() 
{
  assert(_embeddedVertLoc1VboID != 0);
  return _embeddedVertLoc1VboID;
}
GLuint GLSL_TET_MESH::embeddedVertLoc2VboID() 
{
  assert(_embeddedVertLoc2VboID != 0);
  return _embeddedVertLoc2VboID;
}

GLuint GLSL_TET_MESH::embeddedMeshIboID() 
{
  assert(_embeddedMeshIboID != 0);
  return _embeddedMeshIboID;
}
GLsizei GLSL_TET_MESH::embeddedMeshIboSize()  
{
  return _embeddedMeshIboSize;
}
GLsizei GLSL_TET_MESH::embeddedVertSize()  
{
  return _embeddedVertSize;
}

GLuint GLSL_TET_MESH::surfaceVertVboID()      
{
  assert(_surfaceVertVboID != 0);
  return _surfaceVertVboID; 
}
GLuint GLSL_TET_MESH::surfaceVertNeighbor0VboID() 
{
  assert(_surfaceVertNeighbor0VboID != 0);
  return _surfaceVertNeighbor0VboID; 
}
GLuint GLSL_TET_MESH::surfaceVertNeighbor1VboID()  
{ 
  assert(_surfaceVertNeighbor1VboID != 0);
  return _surfaceVertNeighbor1VboID; 
}
GLsizei GLSL_TET_MESH::surfaceVertSize()        
{
  return _surfaceVertSize; 
}

GLuint GLSL_TET_MESH::n1RestTexID()    
{
  assert(_n1RestTexID != 0);
  return _n1RestTexID; 
}
GLuint GLSL_TET_MESH::n2RestTexID()  
{ 
  assert(_n2RestTexID != 0);
  return _n2RestTexID;
}

GLuint GLSL_TET_MESH::UcoordVboID()     
{ 
  assert(_UcoordVboID != 0);
  return _UcoordVboID;
}
GLuint GLSL_TET_MESH::UbasisTexID()     
{
  assert(_UbasisTexID != 0);
  return _UbasisTexID;
}

