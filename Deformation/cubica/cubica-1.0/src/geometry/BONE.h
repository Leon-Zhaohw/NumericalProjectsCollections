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
// BONE.h: interface for the BONE class.
//
//////////////////////////////////////////////////////////////////////

#ifndef BONE_H
#define BONE_H

#include <VEC3.h>
#include <MATRIX3.h>
#include <vector>
#include <QUATERNION.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Bone class for a skinned skeleton
//////////////////////////////////////////////////////////////////////
class BONE {
public:
  BONE(BONE* parent, VEC3F& translation, MATRIX3& rotation) {
    _translation = translation;
    _rotation = QUATERNION(rotation);
    _parent = parent;
    _globalTransform.resizeAndWipe(4,4);
    if (parent != NULL)
      _parent->children().push_back(this);

    // compute the total translation
    _globalRestPosition = _translation;
    BONE* ancestor = parent;
    while (ancestor != NULL)
    {
      _globalRestPosition += ancestor->_translation;
      ancestor = ancestor->_parent;
    }

    _translationOriginal = _translation;
    _rotationOriginal = _rotation;
    _globalTransformOriginal = _globalTransform;
    _globalRestPositionOriginal = _globalRestPosition;
  };
  virtual ~BONE() {};
 
  VEC3F& translation()   { return _translation; };
  VEC3F& translationOriginal()   { return _translationOriginal; };
  VEC3F& globalRestPosition()   { return _globalRestPosition; };
  MATRIX3 rotation()            { return _rotation.toExplicitMatrix3x3(); };
  MATRIX3 rotationOriginal()    { return _rotationOriginal.toExplicitMatrix3x3(); };
  QUATERNION& quaternion()            { return _rotation; };
  QUATERNION& quaternionOriginal()    { return _rotationOriginal; };
  VEC3F& globalPosition()     { return _globalPosition; };
  vector<BONE*>& children()   { return _children; };
  BONE*& parent()             { return _parent; };
  MATRIX& globalTransform()   { return _globalTransform; };
 
  void reset() {
    _translation = _translationOriginal;
    _rotation = _rotationOriginal;
    _globalTransform = _globalTransformOriginal;
    _globalRestPosition = _globalRestPositionOriginal;
  };

private:
  VEC3F   _translation;
  QUATERNION _rotation;
  VEC3F   _globalRestPosition;
  VEC3F   _globalPosition;
  MATRIX  _globalTransform;

  // reset versions (hacky, but assume we don't have a lot of bones)
  VEC3F   _translationOriginal;
  QUATERNION _rotationOriginal;
  VEC3F   _globalRestPositionOriginal;
  VEC3F   _globalPositionOriginal;
  MATRIX  _globalTransformOriginal;

  vector<BONE*> _children;
  BONE* _parent;
};

#endif
