///////////////////////////////////////////////////////////////////////////////////////////////////
// I. LICENSE CONDITIONS
//
// Copyright (c) 2018 by Disney-Pixar
//
// Permission is hereby granted to use this software solely for non-commercial applications
// and purposes including academic or industrial research, evaluation and not-for-profit media
// production. All other rights are retained by Pixar. For use for or in connection with
// commercial applications and purposes, including without limitation in or in connection with
// software products offered for sale or for-profit media production, please contact Pixar at
// tech-licensing@pixar.com.
//
// THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
// NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PIXAR OR ITS AFFILIATES BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef CUBE_H
#define CUBE_H

#include "Settings.h"

#include <array>

namespace CubeSim
{

class Material;

// The vertex indices are:          
//                                  
//               4               6  
//                 o---------------o
//                /|              /|
//               / |             / |
//              /  |            /  |
//           5 /   |         7 /   |
//            o---------------o    |
//            |  0 |          |  2 |
//            |    o----------|----o
//            |   /           |   / 
//            |  /            |  /  
//            | /             | /   
//            |/              |/    
//            o---------------o     
//           1               3      
//                                  
//                                  
//     +z                           
//      |                           
//      o--- +y                     
//     /                            
//   +x                             

class Cube final
{

public:

    Cube() = default;
    Cube(const bool cubeFixed, const std::array<int,8>& vertexIndices,
         const std::vector<Vector3>& verts, const std::vector<Vector3>& restVerts);

    Cube(Cube&) = delete;
    Cube& operator=(Cube&) = delete;

    Cube(Cube&&) = default;
    Cube& operator=(Cube&&) = default;

    // True if all corners of the cube are fixed.
    bool AllCornersFixed() const;

    // Computes the deformation gradient at any of the cube's eight Gauss points.
    Matrix3 ComputeF(const std::vector<Vector3>& verts, const int gpIdx) const;

    // Integrates the strain energy density over the cube.
    Scalar ComputeInternalEnergy(const Material& material, const std::vector<Vector3>& verts) const;

    // Computes the forces on each corner of the cube. 3 DoFs per vertex.
    Vector24 ComputeForces(const Material& material, const std::vector<Vector3>& verts) const;

    // Computes the isotropic material force Jacobian.
    // TODO: rename the fullJacobian parameter...
    Matrix2424 ComputeForceJacobian(const Material& material, const std::vector<Vector3>& verts) const;

    // Given a vertex (corner) number in [0,8), returns the global index of that vertex.
    int VertexIndex(const int vertexNumber) const;

private:

    // Computes the isotropic material force at a given Gauss point.
    Vector24 _ComputeForces(const Material& material, const std::vector<Vector3>& verts, const int gpIdx) const;

    // Computes the force Jacobian at a given Gauss point.
    Matrix2424 _ComputeForceJacobian(const Material& material,
                                     const std::vector<Vector3>& verts,
                                     const int gaussPoint) const;

    // True if all corners of the cube are fixed.
    bool _cubeFixed;

    // Vertex indices in the *global* DoFs.
    std::array<int,8> _vertexIndices;

    // Rest pose (material) matrix.
    Matrix38 _Dm;

    // Rest volume of the cube.
    Scalar _restVolume;

    // Cache some per-gauss point quantities.
    std::array<Matrix38,8> _Bmg;
    std::array<Matrix83,8> _HDmHInv;
    std::array<Matrix924,8> _pFpuMatrix;

};

}

#endif
