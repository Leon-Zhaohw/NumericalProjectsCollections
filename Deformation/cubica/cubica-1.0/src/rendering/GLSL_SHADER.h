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
// GLSL_SHADER.h: Interface for OpenGL GLSL wrapper class
//
/////////////////////////////////////

#ifndef GLSL_SHADER_H
#define GLSL_SHADER_H

#if _WIN32
#include <gl/glew.h>
#include <gl/glut.h>
#elif USING_OSX
#include <GL/glew.h>
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#endif

#include <iostream>
#include <map>
#include <vector>
#include <cstdio>
#include <string>

using namespace std;

//////////////////////////////////////////////////////////////////////
// OpenGL GLSL Shader Program wrapper class 
//////////////////////////////////////////////////////////////////////
class GLSL_SHADER {

public:
  // Constructor for empty shader
  GLSL_SHADER();
  // Constructor for shader program given a vertex and fragment shader filename
  GLSL_SHADER(const char* vertFilename, const char* fragFilename);
  // Destructor
  virtual ~GLSL_SHADER();

  // Attach a vertex shader to the shader program
  bool attachVert(const char* filename);   
  // Attach a fragment shader to the shader program
  bool attachFrag(const char* filename);
  // Attach a geometry shader to the shader program
  bool attachGeom(const char* filename);

  // Get the uniform location in the compiled shader
  GLint uniformLoc(const char* uniformName);
  // Set a program constant (#define 'constName' 'constValue'), in vert=1, frag=2, vert&frag=3.
  void setProgramConst(const char* constName, char* constValue, unsigned short whichShader=3);
  void setProgramConst(const char* constName, int constValue, unsigned short whichShader=3);
  void setProgramConst(const char* constName, float constValue, unsigned short whichShader=3);

  // Compile and link the shader program using the previously attached shaders
  bool compile();
  // Return the OpenGL handle to the shader program
  GLuint getHandle();

  void setGeomInput(GLenum input);
  void setGeomOutput(GLenum output);
  void setGeomVerticesOut(GLint max);

  // Simple function to load a text file, returned as a char array
  static char* loadShader(const char* filename);

protected:
  // The GLSL shader program handle
  GLuint _handle;
  // Flag to keep track of whether or not the shader has been compiled
  bool _isCompiled;
  // List of shader handles used to construct shader program
  vector<GLuint> _shaderHandles;

  // Geometry shader constants
  GLenum _geomInput;
  GLenum _geomOutput;
  GLint _geomVerticesOut;

  // Map of uniform-value pairs
  map<string, GLint> _uniformMap;
  // Map of constant-value pairs for the vertex shaders
  map<string, string> _constVertMap;
  // Map of constant-value pairs for the fragment shaders
  map<string, string> _constFragMap;
    
  // Vector of vertex shader strings
  map<string,string> _vertexShaders;
  // Vector of fragment shader strings
  map<string,string> _fragmentShaders;

  // Map of constant-value pairs for the geometry shaders
  map<string, string> _constGeomMap;
  // Vector of fragment shader strings
  map<string,string> _geometryShaders;

  // cache the shader filenames
  string _vertFilename;
  string _fragFilename;
};

#endif

