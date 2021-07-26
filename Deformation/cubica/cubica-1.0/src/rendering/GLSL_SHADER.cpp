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
// GLSL_SHADER.cpp
//
//////////////////////////////////////////////////////////////////////

#include <GLSL_SHADER.h>

#include <cassert>
#include <cstdio>

//////////////////////////////////////////////////////////////////////
// Constructor for empty shader
//////////////////////////////////////////////////////////////////////
GLSL_SHADER::GLSL_SHADER() 
  : _handle(0), _isCompiled(false),
  _vertFilename(""),
  _fragFilename("")
{
  _handle = glCreateProgram();
  if (_handle == 0)
  {
    cout << " An error occurred when trying to create a GLSL Shader Program." << endl;
  }
  else 
  {
    _geomInput = GL_POINTS;
    _geomOutput = GL_LINE_STRIP;
    _geomVerticesOut = 0;

    /*_uniformMap.clear();
    _constVertMap.clear();
    _constFragMap.clear();
    _constGeomMap.clear();
    _vertexShaders.clear();
    _fragmentShaders.clear();
    _geometryShaders.clear();
    _shaderHandles.clear();*/
  } 
}

//////////////////////////////////////////////////////////////////////
// Constructor for shader program given a vertex and fragment shader filename
//////////////////////////////////////////////////////////////////////
GLSL_SHADER::GLSL_SHADER(const char* vertFilename, const char* fragFilename)
  : _handle(0), _isCompiled(false),
  _vertFilename(vertFilename),
  _fragFilename(fragFilename)
{
  _handle = glCreateProgram();
  if (_handle == 0)
  {
    cout << " An error occurred when trying to create a GLSL Shader Program." << endl;
  }
  attachVert(vertFilename);
  attachFrag(fragFilename);
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
GLSL_SHADER::~GLSL_SHADER()
{
  if (_shaderHandles.size() > 0)
  {
    for (int x = 0; x < (int)_shaderHandles.size(); x++)
      glDeleteShader(_shaderHandles[x]);
    _shaderHandles.clear();
  }

  if (_handle != 0)
    glDeleteProgram(_handle);

  _uniformMap.clear();
  _constVertMap.clear();
  _constFragMap.clear();
  _vertexShaders.clear();
  _fragmentShaders.clear();
  _geometryShaders.clear();
} 

//////////////////////////////////////////////////////////////////////
// Attach a vertex shader to the shader program
//////////////////////////////////////////////////////////////////////
bool GLSL_SHADER::attachVert(const char* filename)
{
  _vertFilename = string(filename);
  char *shaderSource = loadShader(filename);
  if (shaderSource != NULL)
  {
    if (_vertexShaders.find(string(filename)) == _vertexShaders.end())
    {
      _vertexShaders[string(filename)] = string(shaderSource);
      return true;
    }
  }
  return false;
}

//////////////////////////////////////////////////////////////////////
// Attach a fragment shader to the shader program
//////////////////////////////////////////////////////////////////////
bool GLSL_SHADER::attachFrag(const char* filename)
{
  _fragFilename = string(filename);
  char *shaderSource = loadShader(filename);
  if (shaderSource != NULL)
  {
    if (_fragmentShaders.find(string(filename)) == _fragmentShaders.end())
    {
      _fragmentShaders[string(filename)] = string(shaderSource);
      return true;
    }
  }
  return false;
}

//////////////////////////////////////////////////////////////////////
// Attach a geometry shader to the shader program
//////////////////////////////////////////////////////////////////////
bool GLSL_SHADER::attachGeom(const char* filename)
{
  char *shaderSource = loadShader(filename);
  if (shaderSource != NULL)
  {
    if (_geometryShaders.find(string(filename)) == _geometryShaders.end())
    {
      _geometryShaders[string(filename)] = string(shaderSource);
      return true;
    }
  }
  return false;
}// */

//////////////////////////////////////////////////////////////////////
// Get the uniform location in the compiled shader
//////////////////////////////////////////////////////////////////////
GLint GLSL_SHADER::uniformLoc(const char* uniformName)
{
  if (!_isCompiled)
  {
    cout << " Shader program must be compiled before a uniform location can be extracted." << endl;
    return -1;
  }

  if (_uniformMap.find(string(uniformName)) == _uniformMap.end())
  {
    GLint uniformLoc = glGetUniformLocation(_handle, uniformName);
    if (uniformLoc < 0)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Shader program returned an invalid uniform location (" << uniformLoc << ") for \'" << string(uniformName) << "\'" << endl;
      cout << " Vertex program: " << _vertFilename.c_str() << endl;
      cout << " Fragment program: " << _fragFilename.c_str() << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    }
    _uniformMap[string(uniformName)] = uniformLoc;
  }
  return _uniformMap[string(uniformName)];
}

//////////////////////////////////////////////////////////////////////
// Set a program constant (#define 'constName' 'constValue'), in vert=1, frag=2, vert&frag=3.
//////////////////////////////////////////////////////////////////////
void GLSL_SHADER::setProgramConst(const char* constName, char* constValue, unsigned short whichShader)
{
  // whichShader is a bit vector: 1 is vertex, 2 is frag, 4 is geom.
  if ((whichShader & 1) == 1)   // vertex
  {
    if (_constVertMap.find(string(constName)) == _constVertMap.end())
    {
      _constVertMap[string(constName)] = string(constValue);
    }
  }
  if ((whichShader & 2) == 2)   // fragment
  {
    if (_constFragMap.find(string(constName)) == _constFragMap.end())
    {
      _constFragMap[string(constName)] = string(constValue);
    }
  }
  if ((whichShader & 4) == 4)   // geometry
  {
    if (_constGeomMap.find(string(constName)) == _constGeomMap.end())
    {
      _constGeomMap[string(constName)] = string(constValue);
    }
  }// */

  // Constants are inserted into the code, so changing them requires a recompile
  if (_isCompiled)
  {
    cout << " Constant changed on compiled shader program; recompiling." << endl;
    _isCompiled = false;
    _uniformMap.clear();
    compile();
  }
}

//////////////////////////////////////////////////////////////////////
// Set a program constant (#define 'constName' 'constValue'), in vert=1, frag=2, vert&frag=3.
//////////////////////////////////////////////////////////////////////
void GLSL_SHADER::setProgramConst(const char* constName, int constValue, unsigned short whichShader)
{
  const size_t strSize = 16;
  char strValue[strSize];
  int result = 0;
  result = snprintf(strValue, strSize, "%d", constValue);
  if (result < 0)
  {
    cout << " An error occurred encoding int to string during setProgramConst" << endl;
  }
  if ((size_t)result >= strSize)
  {
    cout << " Shader program constant " << string(constName) << " had its value truncated to " << string(strValue) << endl;
  }
  setProgramConst(constName, strValue, whichShader);
}

//////////////////////////////////////////////////////////////////////
// Set a program constant (#define 'constName' 'constValue'), in vert=1, frag=2, vert&frag=3.
//////////////////////////////////////////////////////////////////////
void GLSL_SHADER::setProgramConst(const char* constName, float constValue, unsigned short whichShader)
{
  const size_t strSize = 16;
  char strValue[strSize];
  int result = 0;
  if (constValue > 10.0f)
  {
    result = snprintf(strValue, strSize, "%.1f", constValue);
  } else {
    result = snprintf(strValue, strSize, "%.5f", constValue);
  }
    if (result < 0)
  {
    cout << " An error occurred encoding float to string during setProgramConst" << endl;
  }
  if ((size_t)result >= strSize)
  {
    cout << " Shader program constant " << string(constName) << " had its value truncated to " << string(strValue) << endl;
  }
  setProgramConst(constName, strValue, whichShader);
}


//////////////////////////////////////////////////////////////////////
// Compile and link the shader program using the previously attached shaders
//////////////////////////////////////////////////////////////////////
bool GLSL_SHADER::compile()
{
  GLchar infoLog[1024];
  GLsizei infoLogLength;

  GLint success;
  GLuint shaderHandle = 0;
  string shaderSource = "";
  map<string, string>::iterator constIter;
  map<string, string>::iterator shaderIter;

  // If shaderHandles already stored, delete them and clear
  if (_shaderHandles.size() > 0)
  {
    for (int x = 0; x < (int)_shaderHandles.size(); x++)
      glDeleteShader(_shaderHandles[x]);
    _shaderHandles.clear();
  }
  
  ////////////////////////////////////////////////////////////////////
  // Compile and attach VERTEX SHADERS
  for (shaderIter = _vertexShaders.begin(); shaderIter != _vertexShaders.end(); shaderIter++)
  {
    // Clear the source from last time through loop (if any)
    shaderSource.clear();

    // Loop through all constant in the vertex constant map
    for (constIter = _constVertMap.begin(); constIter != _constVertMap.end(); constIter++)
    {
      shaderSource.append("#define ");
      shaderSource.append(constIter->first);
      shaderSource.append(" ");
      shaderSource.append(constIter->second);
      shaderSource.append("\n");
    }

    // Append actual shader source after constants
    shaderSource.append(shaderIter->second);

    // Create new vertex shader object
    shaderHandle = glCreateShader(GL_VERTEX_SHADER);	
    if (shaderHandle == 0)
    {
      cout << " Could not create a vertex shader during compile." << endl;
    }

    const char *srcloc = shaderSource.c_str();

    //fprintf(stdout, "\n===== Vertex source: =====\n%s\n", srcloc);

    // Copy shader source to GPU, and compile
    glShaderSource(shaderHandle, 1, (const GLchar**) &srcloc, NULL);
    assert(shaderHandle != 0);   
    glCompileShader(shaderHandle);

    // Check status of compilation
    glGetShaderiv(shaderHandle, GL_COMPILE_STATUS, &success);
    if (success)
    {
      assert(_handle != 0);
      glAttachShader(_handle, shaderHandle);
      _shaderHandles.push_back(shaderHandle);
    } 
    else 
    {
      cout << " Vertex shader " << shaderIter->first << " didn't compile!" << endl;
      cout << endl << " Vertex Shader Source: " << endl;
      cout << shaderSource << endl;
      glGetShaderInfoLog(shaderHandle, 1023, &infoLogLength, infoLog);
      cout << " Info log: " << infoLog << endl;
      glDeleteShader(shaderHandle);
    }
  }

  ////////////////////////////////////////////////////////////////////
  // Compile and attach FRAGMENT SHADERS
  for (shaderIter = _fragmentShaders.begin(); shaderIter != _fragmentShaders.end(); shaderIter++)
  {
    // Clear the source from last time through loop (if any)
    shaderSource.clear();

    // Loop through all constant in the fragment constant map
    for (constIter = _constFragMap.begin(); constIter != _constFragMap.end(); constIter++)
    {
      shaderSource.append("#define ");
      shaderSource.append(constIter->first);
      shaderSource.append(" ");
      shaderSource.append(constIter->second);
      shaderSource.append("\n");
    }

    // Append actual shader source after constants
    shaderSource.append(shaderIter->second);

    // Create new vertex shader object
    shaderHandle = glCreateShader(GL_FRAGMENT_SHADER);	
    if (shaderHandle == 0)
    {
      cout << " Could not create a fragment shader during compile." << endl;
    }

    const char *srcloc = shaderSource.c_str();

    //fprintf(stdout, "\n===== Fragment source: =====\n%s\n", srcloc);

    // Copy shader source to GPU, and compile
    glShaderSource(shaderHandle, 1, (const GLchar**) &srcloc, NULL);
    assert(shaderHandle != 0);
    glCompileShader(shaderHandle);

    // Check status of compilation
    glGetShaderiv(shaderHandle, GL_COMPILE_STATUS, &success);
    if (success)
    {
      assert(_handle != 0);
      glAttachShader(_handle, shaderHandle);
      _shaderHandles.push_back(shaderHandle);
    } 
    else 
    {
      cout << " Fragment shader " << shaderIter->first << " didn't compile!" << endl;
      cout << endl << " Fragment Shader Source: " << endl;
      cout << shaderSource << endl;

      glGetShaderInfoLog(shaderHandle, 1023, &infoLogLength, infoLog);
      cout << " Info log: " << infoLog << endl;
      glDeleteShader(shaderHandle);
    }
  }
  
  ////////////////////////////////////////////////////////////////////
  // Compile and attach GEOMETRY SHADERS
  for (shaderIter = _geometryShaders.begin(); shaderIter != _geometryShaders.end(); shaderIter++)
  {
    // Clear the source from last time through loop (if any)
    shaderSource.clear();

    // Loop through all constant in the geometry constant map
    for (constIter = _constGeomMap.begin(); constIter != _constGeomMap.end(); constIter++)
    {
      shaderSource.append("#define ");
      shaderSource.append(constIter->first);
      shaderSource.append(" ");
      shaderSource.append(constIter->second);
      shaderSource.append("\n");
    }

    // Append actual shader source after constants
    shaderSource.append(shaderIter->second);

    // Create new vertex shader object
    shaderHandle = glCreateShader(GL_GEOMETRY_SHADER_EXT);	
    if (shaderHandle == 0)
    {
      cout << " Could not create a geometry shader during compile." << endl;
    }

    const char *srcloc = shaderSource.c_str();

    //fprintf(stdout, "\n===== Geometry source: =====\n%s\n", srcloc);

    // Copy shader source to GPU, and compile
    glShaderSource(shaderHandle, 1, (const GLchar**) &srcloc, NULL);	 
    glCompileShader(shaderHandle);

    // Check status of compilation
    glGetShaderiv(shaderHandle, GL_COMPILE_STATUS, &success);
    if (success)
    {
      glAttachShader(_handle, shaderHandle);
      _shaderHandles.push_back(shaderHandle);

      glProgramParameteriEXT(_handle, GL_GEOMETRY_INPUT_TYPE_EXT, _geomInput);
  		glProgramParameteriEXT(_handle, GL_GEOMETRY_OUTPUT_TYPE_EXT, _geomOutput);
      glProgramParameteriEXT(_handle, GL_GEOMETRY_VERTICES_OUT_EXT, _geomVerticesOut);      
    } 
    else 
    {
      cout << " Geometry shader " << shaderIter->first << " didn't compile!" << endl;
      cout << endl << " Geometry Shader Source: " << endl;
      cout << shaderSource << endl;

      glGetShaderInfoLog(shaderHandle, 1023, &infoLogLength, infoLog);
      cout << " Info log: " << infoLog << endl;
      glDeleteShader(shaderHandle);
    }
  }// */

  ////////////////////////////////////////////////////////////////////
  // LINK shaders to program object
  assert(_handle != 0);
  glLinkProgram(_handle);
  glGetProgramiv(_handle, GL_LINK_STATUS, &success);
  if (!success)
  {
    cout << " Linking shader program failed!" << endl;
    glGetProgramInfoLog(_handle, 1023, &infoLogLength, infoLog);
    cout << " Info log: " << infoLog << endl;
  }

  if (success)
  {
    _isCompiled = true;
    return true;
  }
  
  return false;
}

//////////////////////////////////////////////////////////////////////
// Return the OpenGL handle to the shader program
//////////////////////////////////////////////////////////////////////
GLuint GLSL_SHADER::getHandle()
{
  if (_isCompiled)
    return _handle;
  cout << " Shader program is not compiled!" << endl;
  return 0;
}

//////////////////////////////////////////////////////////////////////
// Set the program parameter for the geometry shaders input
//////////////////////////////////////////////////////////////////////
void GLSL_SHADER::setGeomInput(GLenum input)
{
  if ( input == GL_POINTS
    || input == GL_LINES 
    || input == GL_TRIANGLES
    || input == GL_LINES_ADJACENCY_EXT 
    || input == GL_TRIANGLES_ADJACENCY_EXT )
  {
    _geomInput = input;

    if (_isCompiled)
    {
      cout << " Geometry shader constant changed on compiled shader program; recompiling." << endl;
      _isCompiled = false;
      _uniformMap.clear();
      compile();
    }
  }
  else
  {
    cout << " Error! Unknown geometry shader input type." << endl;
  }
} // */
 
//////////////////////////////////////////////////////////////////////
// Set the program parameter for the geometry shaders output
//////////////////////////////////////////////////////////////////////
void GLSL_SHADER::setGeomOutput(GLenum output)
{
  if ( output == GL_POINTS
    || output == GL_LINE_STRIP 
    || output == GL_TRIANGLE_STRIP )
  {
    _geomOutput = output;

    if (_isCompiled)
    {
      cout << " Geometry shader constant changed on compiled shader program; recompiling." << endl;
      _isCompiled = false;
      _uniformMap.clear();
      compile();
    }

  }
  else
  {
    cout << " Error! Unknown geometry shader output type." << endl;
  }
}// */

//////////////////////////////////////////////////////////////////////
// Set the program parameter for the geometry shaders maximum number of output vertices
//////////////////////////////////////////////////////////////////////
void GLSL_SHADER::setGeomVerticesOut(GLint max)
{
  int temp;
	glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES_EXT,&temp);

  if (max > temp)
  {
    cout << " setGeomVerticesOut() was given " << max;
    cout << ", but cannot me higher than " << temp << endl;
    _geomVerticesOut = temp;
  } 
  else
  {
    _geomVerticesOut = max;
  }

  if (_isCompiled)
  {
    cout << " Geometry shader constant changed on compiled shader program; recompiling." << endl;
    _isCompiled = false;
    _uniformMap.clear();
    compile();
  }
}//  */


//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Simple function to load a text file, returned as a char array
//////////////////////////////////////////////////////////////////////
char* GLSL_SHADER::loadShader(const char* filename)
{
  //char tempChar = NULL;
  long charCount = 0;
  FILE *filePtr;


  filePtr = fopen(filename, "rb");				// Open shader file
  if (filePtr == NULL)
  {
    cout << " Could not open file " << filename << "!" << endl;
    return NULL;
  }

  fseek(filePtr, 0, SEEK_END);					// goto end of file...
  charCount = ftell(filePtr);					// and ask how big it is
  char *shaderSource = new char[charCount+1];			// make the shader buffer

  fseek(filePtr, 0, SEEK_SET);					// go back to the beginning of the file...
  fread(shaderSource, sizeof(char), charCount, filePtr);	// read in the entire file...
  fclose(filePtr);						// then close, cause we're done with it
  shaderSource[charCount] = '\0';				// add the NULL termination to the end

  return shaderSource;						// return the shader
}


