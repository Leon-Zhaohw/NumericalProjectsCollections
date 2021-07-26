/*************************************************************************
 *                                                                       *
 * Open Dynamics Engine, Copyright (C) 2001,2002 Russell L. Smith.       *
 * All rights reserved.  Email: russ@q12.org   Web: www.q12.org          *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of EITHER:                                  *
 *   (1) The GNU Lesser General Public License as published by the Free  *
 *       Software Foundation; either version 2.1 of the License, or (at  *
 *       your option) any later version. The text of the GNU Lesser      *
 *       General Public License is included with this library in the     *
 *       file LICENSE.TXT.                                               *
 *   (2) The BSD-style license that is included with this library in     *
 *       the file LICENSE-BSD.TXT.                                       *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the files    *
 * LICENSE.TXT and LICENSE-BSD.TXT for more details.                     *
 *                                                                       *
 *************************************************************************/

#include <ode/ode.h>
#include <iostream>
#include "VEC3.h"
#include "MATRIX3.h"
#include <vector>
#include <map>

#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#include <glvu.hpp>

using namespace std;
GLVU glvu;
bool animate = false;
bool singleStep = false;

////////////////////////////////////////////////////////////////////////
// The ragdoll code is adapted from the Python code at:
//
// http://monsterden.net/software/ragdoll-pyode-tutorial
////////////////////////////////////////////////////////////////////////
VEC3F norm3(VEC3F v)
{
  Real l = norm(v);
  if (l > 0.0) 
    return VEC3F(v[0] / l, v[1] / l, v[2] / l);
  return VEC3F(0.0, 0.0, 0.0);
}

Real acosdot3(VEC3F a, VEC3F b)
{
  Real x = (a * b);
  if (x < -1.0) 
    return M_PI;
  else if (x > 1.0) 
    return 0.0;
  return acos(x);
}

VEC3F rotate3(MATRIX3 m, VEC3F v)
{
  return VEC3F(v[0] * m(0,0) + v[1] * m(0,1) + v[2] * m(0,2),
    v[0] * m(1,0) + v[1] * m(1,1) + v[2] * m(1,2),
    v[0] * m(2,0) + v[1] * m(2,1) + v[2] * m(2,2));
}

// repack an ODE matrix into MATRIX3
MATRIX3 repack(const dReal* data)
{
  return MATRIX3(
      VEC3F(data[0], data[1], data[2]),
      VEC3F(data[4], data[5], data[6]),
      VEC3F(data[8], data[9], data[10]));
}

void unpack(const MATRIX3& input, dReal* data)
{
  data[0] = input(0,0);
  data[1] = input(0,1);
  data[2] = input(0,2);
  data[3] = 0;

  data[4] = input(1,0);
  data[5] = input(1,1);
  data[6] = input(1,2);
  data[7] = 0;

  data[8] = input(2,0);
  data[9] = input(2,1);
  data[10] = input(2,2);
  data[11] = 0;
}

MATRIX3 calcRotMatrix(VEC3F axis, Real angle)
{
  Real cosTheta = cos(angle);
  Real sinTheta = sin(angle);
  Real t = 1.0 - cosTheta;
  MATRIX3 final;
  final(0,0) = t * axis[0] * axis[0] + cosTheta;
  final(0,1) = t * axis[0] * axis[1] - sinTheta * axis[2];
  final(0,2) = t * axis[0] * axis[2] + sinTheta * axis[1];
  final(1,0) = t * axis[0] * axis[1] + sinTheta * axis[2];
  final(1,1) = t * axis[1] * axis[1] + cosTheta;
  final(1,2) = t * axis[1] * axis[2] - sinTheta * axis[0];
  final(2,0) = t * axis[0] * axis[2] - sinTheta * axis[1];
  final(2,1) = t * axis[1] * axis[2] + sinTheta * axis[0];
  final(2,2) = t * axis[2] * axis[2] + cosTheta;
  return final;
}

void makeOpenGLMatrix(const dReal* r, const dReal* p, float* data)
{
  data[0]  = r[0]; data[1]  = r[4]; data[2]  = r[8]; data[3]  = 0;
  data[4]  = r[1]; data[5]  = r[5]; data[6]  = r[9]; data[7]  = 0;
  data[8]  = r[2]; data[9]  = r[6]; data[10] = r[10]; data[11] = 0;
  data[12] = p[0]; data[13] = p[1]; data[14] = p[2]; data[15] = 1;
}

VEC3F getBodyRelVec(dBodyID b, VEC3F v)
{
  const dReal* bRotation = dBodyGetRotation(b);
  MATRIX3 rotation;
  rotation(0,0) = bRotation[0];
  rotation(0,1) = bRotation[1];
  rotation(0,2) = bRotation[2];
  rotation(1,0) = bRotation[4];
  rotation(1,1) = bRotation[5];
  rotation(1,2) = bRotation[6];
  rotation(2,0) = bRotation[8];
  rotation(2,1) = bRotation[9];
  rotation(2,2) = bRotation[10];

  return rotate3(rotation.transpose(), v);
}

VEC3F angularVelocity(dBodyID body)
{
  const dReal* velocity = dBodyGetAngularVel(body);
  return VEC3F(velocity[0], velocity[1], velocity[2]);
}

////////////////////////////////////////////////////////////////////////
// Rag doll body parameters
////////////////////////////////////////////////////////////////////////
MATRIX3 rightRot(VEC3F(0.0, 0.0, -1.0), VEC3F(0.0, 1.0, 0.0), VEC3F(1.0, 0.0, 0.0));
MATRIX3 leftRot(VEC3F(0.0, 0.0, 1.0), VEC3F(0.0, 1.0, 0.0), VEC3F(-1.0, 0.0, 0.0));
MATRIX3 upRot(VEC3F(1.0, 0.0, 0.0), VEC3F(0.0, 0.0, -1.0), VEC3F(0.0, 1.0, 0.0));
MATRIX3 downRot(VEC3F(1.0, 0.0, 0.0), VEC3F(0.0, 0.0, -1.0), VEC3F(0.0, 1.0, 0.0));
MATRIX3 bkwdRot(VEC3F(1.0, 0.0, 0.0), VEC3F(0.0, 1.0, 0.0), VEC3F(0.0, 0.0, 1.0));

VEC3F rightAxis(1.0, 0.0, 0.0);
VEC3F leftAxis(-1.0, 0.0, 0.0);
VEC3F upAxis(0.0, 1.0, 0.0);
VEC3F downAxis(0.0, -1.0, 0.0);
VEC3F bkwdAxis(0.0, 0.0, 1.0);
VEC3F fwdAxis(0.0, 0.0, -1.0);

Real UPPER_ARM_LEN = 0.30;
Real FORE_ARM_LEN = 0.25;
Real HAND_LEN = 0.13;
Real FOOT_LEN = 0.18;
Real HEEL_LEN = 0.05;

Real BROW_H = 1.68;
Real MOUTH_H = 1.53;
Real NECK_H = 1.50;
Real SHOULDER_H = 1.37;
Real CHEST_H = 1.35;
Real HIP_H = 0.86;
Real KNEE_H = 0.48;
Real ANKLE_H = 0.08;

Real SHOULDER_W = 0.41;
Real CHEST_W = 0.36;
Real LEG_W = 0.28;
Real PELVIS_W = 0.25;

VEC3F R_SHOULDER_POS(-SHOULDER_W * 0.5, SHOULDER_H, 0.0);
VEC3F L_SHOULDER_POS(SHOULDER_W * 0.5, SHOULDER_H, 0.0);
VEC3F R_ELBOW_POS = R_SHOULDER_POS - VEC3F(UPPER_ARM_LEN, 0.0, 0.0);
VEC3F L_ELBOW_POS = L_SHOULDER_POS + VEC3F(UPPER_ARM_LEN, 0.0, 0.0);
VEC3F R_WRIST_POS = R_ELBOW_POS - VEC3F(FORE_ARM_LEN, 0.0, 0.0);
VEC3F L_WRIST_POS = L_ELBOW_POS + VEC3F(FORE_ARM_LEN, 0.0, 0.0);
VEC3F R_FINGERS_POS = R_WRIST_POS - VEC3F(HAND_LEN, 0.0, 0.0);
VEC3F L_FINGERS_POS = L_WRIST_POS + VEC3F(HAND_LEN, 0.0, 0.0);

VEC3F R_HIP_POS(-LEG_W * 0.5, HIP_H, 0.0);
VEC3F L_HIP_POS(LEG_W * 0.5, HIP_H, 0.0);
VEC3F R_KNEE_POS(-LEG_W * 0.5, KNEE_H, 0.0);
VEC3F L_KNEE_POS(LEG_W * 0.5, KNEE_H, 0.0);
VEC3F R_ANKLE_POS(-LEG_W * 0.5, ANKLE_H, 0.0);
VEC3F L_ANKLE_POS(LEG_W * 0.5, ANKLE_H, 0.0);
VEC3F R_HEEL_POS = R_ANKLE_POS - VEC3F(0.0, 0.0, HEEL_LEN);
VEC3F L_HEEL_POS = L_ANKLE_POS - VEC3F(0.0, 0.0, HEEL_LEN);
VEC3F R_TOES_POS = R_ANKLE_POS + VEC3F(0.0, 0.0, FOOT_LEN);
VEC3F L_TOES_POS = L_ANKLE_POS + VEC3F(0.0, 0.0, FOOT_LEN);

////////////////////////////////////////////////////////////////////////
// RAGDOLL class
////////////////////////////////////////////////////////////////////////
class RAG_DOLL {

private:
  struct JOINT 
  {
    dJointID jointID;
    VEC3F baseAxis;
    VEC3F baseTwistUp;
    VEC3F baseTwistUp2;
    Real flexLimit;
    Real twistLimit;
    Real flexForce;
    Real twistForce;
    string style;
  };

public:
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  RAG_DOLL(dWorldID world,
           dSpaceID space,
           Real density,
           VEC3F offset = VEC3F(0,0,0))
  {
    _world = world;
    _space = space;
    _density = density;
    _totalMass = 0.0;
    _offset = offset;
    _jointGroup = dJointGroupCreate(0);
    _contactGroup = dJointGroupCreate(0);

    // armadillo settings
    VEC3F chestLeft  = VEC3F(0.323926, 0.298774, 0.404543);
    VEC3F chestRight = VEC3F(-0.255926, 0.301226, 0.417457);
    _chest = addBody(chestLeft ,chestRight, 0.3);
    VEC3F bellyLeft  = VEC3F(0.032, 0.170256, 0.473662);
    VEC3F bellyRight = VEC3F(0.032, -0.248256, 0.438338);
    _belly = addBody(bellyLeft ,bellyRight, 0.2);
    VEC3F pelvisLeft  = VEC3F(-0.224915, -0.186, 0.443342);
    VEC3F pelvisRight = VEC3F(0.294915, -0.186, 0.456658);
    _pelvis = addBody(pelvisLeft ,pelvisRight, 0.12);

    VEC3F headLeft  = VEC3F(0.049, 0.585756, -0.165875);
    VEC3F headRight = VEC3F(0.049, 0.300244, 0.315875);

    _head = addSphere(headLeft ,headRight, 0.28);

    VEC3F leftUpperArmLeft  = VEC3F(-0.220985, 0.234347, 0.344812);
    VEC3F leftUpperArmRight = VEC3F(-0.455015, 0.259653, 0.193188);
    _leftUpperArm = addBody(leftUpperArmLeft ,leftUpperArmRight, 0.09);
    VEC3F leftForeArmLeft  = VEC3F(-0.421589, 0.258004, 0.219155);
    VEC3F leftForeArmRight = VEC3F(-0.530411, 0.347996, 0.000844961);
    _leftForeArm = addBody(leftForeArmLeft ,leftForeArmRight, 0.09);
    VEC3F leftHandLeft  = VEC3F(-0.523, 0.343871, -0.00408064);
    VEC3F leftHandRight = VEC3F(-0.523, 0.380129, -0.159919);
    _leftHand = addBody(leftHandLeft ,leftHandRight, 0.06);
    VEC3F rightUpperArmLeft  = VEC3F(0.497623, 0.221, 0.233973);
    VEC3F rightUpperArmRight = VEC3F(0.288377, 0.221, 0.420027);
    _rightUpperArm = addBody(rightUpperArmLeft ,rightUpperArmRight, 0.1);
    VEC3F rightForeArmLeft  = VEC3F(0.624517, 0.407923, 0.0889622);
    VEC3F rightForeArmRight = VEC3F(0.501483, 0.246077, 0.251038);
    _rightForeArm = addBody(rightForeArmLeft ,rightForeArmRight, 0.11);
    VEC3F rightHandLeft  = VEC3F(0.609, 0.532853, -0.0555617);
    VEC3F rightHandRight = VEC3F(0.609, 0.429147, 0.0915617);
    _rightHand = addBody(rightHandLeft ,rightHandRight, 0.06);
    VEC3F rightUpperLegLeft  = VEC3F(0.359602, -0.473623, 0.375354);
    VEC3F rightUpperLegRight = VEC3F(0.210398, -0.226377, 0.456646);
    _rightUpperLeg = addBody(rightUpperLegLeft ,rightUpperLegRight, 0.11);
    VEC3F rightLowerLegLeft  = VEC3F(0.241704, -0.728029, 0.547673);
    VEC3F rightLowerLegRight = VEC3F(0.326296, -0.467971, 0.424327);
    _rightLowerLeg = addBody(rightLowerLegLeft ,rightLowerLegRight, 0.1);
    VEC3F rightFootLeft  = VEC3F(0.399786, -0.81094, 0.376629);
    VEC3F rightFootRight = VEC3F(0.194214, -0.79906, 0.535371);
    _rightFoot = addBody(rightFootLeft ,rightFootRight, 0.06);
    VEC3F leftUpperLegLeft  = VEC3F(-0.304147, -0.437885, 0.40273);
    VEC3F leftUpperLegRight = VEC3F(-0.139853, -0.214115, 0.43927);
    _leftUpperLeg = addBody(leftUpperLegLeft ,leftUpperLegRight, 0.11);
    VEC3F leftLowerLegLeft  = VEC3F(-0.268254, -0.748709, 0.50418);
    VEC3F leftLowerLegRight = VEC3F(-0.265746, -0.463291, 0.41182);
    _leftLowerLeg = addBody(leftLowerLegLeft ,leftLowerLegRight, 0.1);
    VEC3F leftFootLeft  = VEC3F(-0.3828, -0.805, 0.351165);
    VEC3F leftFootRight = VEC3F(-0.2332, -0.805, 0.538835);
    _leftFoot = addBody(leftFootLeft ,leftFootRight, 0.06);

    VEC3F diff, diff2, crossed, mean;

    _midSpine = addFixedJoint(_chest, _belly);
    _lowSpine = addFixedJoint(_belly, _pelvis);
    diff = headRight - headLeft;
    diff.normalize();

    diff2 = chestRight - chestLeft;
    diff2.normalize();

    crossed = cross(diff, diff2);
    crossed.normalize();
    _neck = addBallJoint(_chest, _head, headRight, diff, crossed, M_PI * 0.25, M_PI* 0.25, 150.0, 100.0);

    diff = pelvisRight - pelvisLeft;
    diff2 = rightUpperLegRight - rightUpperLegLeft;
    diff2.normalize();
    crossed = cross(diff, diff2);
    crossed.normalize();

    mean = (pelvisRight + rightUpperLegLeft);
    mean *= 0.5;
    _rightHip = addHingeJoint(_pelvis, _rightUpperLeg, rightUpperLegRight, -crossed, -M_PI * 0.5, 0);

    diff = pelvisRight - pelvisLeft;
    diff2 = leftUpperLegRight - leftUpperLegLeft;
    diff2.normalize();
    crossed = cross(diff, diff2);
    crossed.normalize();

    mean = (pelvisRight  + leftUpperLegLeft);
    mean *= 0.5;
    _leftHip = addHingeJoint(_pelvis, _leftUpperLeg, leftUpperLegRight, crossed, -M_PI * 0.5, 0);

    mean = (rightUpperLegRight + rightLowerLegLeft);
    mean *= 0.5;
    _rightKnee = addHingeJoint(_rightUpperLeg, _rightLowerLeg, mean, rightAxis, 0.0, M_PI * 0.5);

    mean = (leftUpperLegRight + leftLowerLegLeft);
    mean *= 0.5;
    _leftKnee = addHingeJoint(_leftUpperLeg, _leftLowerLeg, mean, rightAxis, 0.0, M_PI * 0.5);

    mean = (rightLowerLegRight + rightFootLeft);
    mean *= 0.5;
    _rightAnkle = addHingeJoint(_rightLowerLeg, _rightFoot, mean, leftAxis, -0.1 * M_PI, 0.05 * M_PI);

    mean = (leftLowerLegRight + leftFootLeft);
    mean *= 0.5;
    _leftAnkle = addHingeJoint(_leftLowerLeg, _leftFoot, mean, leftAxis, -0.1 * M_PI, 0.05 * M_PI);

    diff = rightUpperArmRight - rightUpperArmLeft;
    diff.normalize();
    diff2 = chestRight - chestLeft;
    diff2.normalize();
    crossed = cross(diff, diff2);
    crossed.normalize();
    Real limit = M_PI * 0.25;
    _rightShoulder = addBallJoint(_chest, _rightUpperArm, rightUpperArmRight, diff, crossed, limit, limit, 500.0, 500.0);

    diff = leftUpperArmRight - leftUpperArmLeft;
    diff.normalize();
    crossed = cross(diff, diff2);
    crossed.normalize();
    _leftShoulder = addBallJoint(_chest, _leftUpperArm, leftUpperArmLeft, diff, crossed, limit, limit, 500.0, 500.0);

    diff = (rightUpperArmRight - rightUpperArmLeft);
    diff2 = (rightForeArmLeft - rightForeArmRight);
    crossed = cross(diff, diff2);
    crossed.normalize();

    mean = rightUpperArmLeft + rightForeArmRight;
    mean *= 0.5;
    _rightElbow = addHingeJoint(_rightUpperArm, _rightForeArm, mean, crossed, 0.0, 0.6 * M_PI);
  
    diff = (leftUpperArmLeft - leftUpperArmRight);
    diff2 = (leftForeArmRight - leftForeArmLeft);
    crossed = cross(diff, diff2);
    crossed.normalize();

    mean = leftUpperArmRight + leftForeArmLeft;
    mean *= 0.5;
    _leftElbow = addHingeJoint(_leftUpperArm, _leftForeArm, mean, crossed, 0.0, 0.6 * M_PI);

    mean = rightForeArmLeft + rightHandRight;
    mean *= 0.5;

    diff = (rightForeArmRight - rightForeArmLeft);
    diff2 = (rightHandLeft - rightHandRight);
    crossed = cross(diff, diff2);
    crossed.normalize();

    _rightWrist = addHingeJoint(_rightForeArm, _rightHand, mean, crossed, -0.1 * M_PI, 0.2 * M_PI);

    mean = leftForeArmRight + leftHandLeft;
    mean *= 0.5;

    diff = (leftForeArmLeft - leftForeArmRight);
    diff2 = (leftHandRight - leftHandLeft);
    crossed = cross(diff, diff2);
    crossed.normalize();

    _leftWrist = addHingeJoint(_leftForeArm, _leftHand, mean, bkwdAxis, -0.1 * M_PI, 0.2 * M_PI);
  }


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  dBodyID addSphere(VEC3F p1, VEC3F p2, Real radius)
  {
    p1 = p1 + _offset;
    p2 = p2 + _offset;

    Real cyllen = norm(p1 - p2) - radius;

    dBodyID body = dBodyCreate(_world);
    dMass m;

    dMassSetSphere(&m, _density, 0.11);
    dBodySetMass (body, &m);

    dGeomID geom;
    geom = dCreateSphere(_space, radius);
    dGeomSetBody(geom, body);

    VEC3F za = norm3(p2 - p1);
    VEC3F xa(0.0, 1.0, 0.0);
    if (fabs((za * VEC3F(1.0, 0.0, 0.0))) < 0.7) 
      xa = VEC3F(1.0, 0.0, 0.0);
    VEC3F ya = cross(za, xa);
    xa = norm3(cross(ya, za));
    ya = cross(za, xa);
    MATRIX3 rot(VEC3F(xa[0], ya[0], za[0]), 
                VEC3F(xa[1], ya[1], za[1]), 
                VEC3F(xa[2], ya[2], za[2]));

    VEC3F bodyPosition = p1 + p2;
    bodyPosition *= 0.5;
    dBodySetPosition(body, bodyPosition[0], bodyPosition[1], bodyPosition[2]);

    dMatrix3 bodyRotation;
    bodyRotation[0] = xa[0];
    bodyRotation[1] = ya[0];
    bodyRotation[2] = za[0];
    bodyRotation[3] = 0;

    bodyRotation[4] = xa[1];
    bodyRotation[5] = ya[1];
    bodyRotation[6] = za[1];
    bodyRotation[7] = 0;

    bodyRotation[8] = xa[2];
    bodyRotation[9] = ya[2];
    bodyRotation[10] = za[2];
    bodyRotation[11] = 0;

    dBodySetRotation(body, bodyRotation);

    _bodies.push_back(body);
    _radii.push_back(radius);
    _lengths.push_back(cyllen);
    _geoms.push_back(geom);

    // add to the hash
    _bodyIndices[body] = _bodies.size() - 1;

    _totalMass += m.mass;
    return body;
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  dBodyID addBody(VEC3F p1, VEC3F p2, Real radius)
  {
    p1 = p1 + _offset;
    p2 = p2 + _offset;

    Real cyllen = norm(p1 - p2) - radius;

    dBodyID body = dBodyCreate(_world);
    dMass m;

    dMassSetCapsule(&m, _density, 3, radius, cyllen);
    dBodySetMass (body, &m);

    dGeomID geom = dCreateCapsule(_space, radius, cyllen);
    dGeomSetBody(geom, body);

    VEC3F za = norm3(p2 - p1);
    VEC3F xa(0.0, 1.0, 0.0);
    if (fabs((za * VEC3F(1.0, 0.0, 0.0))) < 0.7) 
      xa = VEC3F(1.0, 0.0, 0.0);
    VEC3F ya = cross(za, xa);
    xa = norm3(cross(ya, za));
    ya = cross(za, xa);
    MATRIX3 rot(VEC3F(xa[0], ya[0], za[0]), 
                VEC3F(xa[1], ya[1], za[1]), 
                VEC3F(xa[2], ya[2], za[2]));

    VEC3F bodyPosition = p1 + p2;
    bodyPosition *= 0.5;
    dBodySetPosition(body, bodyPosition[0], bodyPosition[1], bodyPosition[2]);

    dMatrix3 bodyRotation;
    bodyRotation[0] = xa[0];
    bodyRotation[1] = ya[0];
    bodyRotation[2] = za[0];
    bodyRotation[3] = 0;

    bodyRotation[4] = xa[1];
    bodyRotation[5] = ya[1];
    bodyRotation[6] = za[1];
    bodyRotation[7] = 0;

    bodyRotation[8] = xa[2];
    bodyRotation[9] = ya[2];
    bodyRotation[10] = za[2];
    bodyRotation[11] = 0;

    dBodySetRotation(body, bodyRotation);

    _bodies.push_back(body);
    _radii.push_back(radius);
    _lengths.push_back(cyllen);
    _geoms.push_back(geom);

    // add to the hash
    _bodyIndices[body] = _bodies.size() - 1;

    _totalMass += m.mass;
    return body;
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  JOINT addFixedJoint(dBodyID body1, dBodyID body2)
  {
    dJointID jointID = dJointCreateFixed(_world, _jointGroup);
    dJointAttach(jointID, body1, body2);
    dJointSetFixed(jointID);

    JOINT joint;
    joint.jointID = jointID;
    joint.style = string("fixed");

    _joints.push_back(joint);
    _parentBodies[_bodyIndices[body2]] = _bodyIndices[body1];
    return joint;
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  JOINT addBallJoint(dBodyID body1, dBodyID body2, VEC3F anchor, VEC3F baseAxis, VEC3F baseTwistUp,
    Real flexLimit = M_PI, Real twistLimit = M_PI, Real flexForce = 0.0, Real twistForce = 0.0)
  {
    anchor = anchor + _offset;

    dJointID jointID = dJointCreateBall(_world, _jointGroup); 
    dJointAttach(jointID, body1, body2);
    dJointSetBallAnchor(jointID, anchor[0], anchor[1], anchor[2]);

    JOINT joint;
    joint.jointID = jointID;

    joint.baseAxis = getBodyRelVec(body1, baseAxis);
    VEC3F tempTwistUp = getBodyRelVec(body1, baseTwistUp);
    VEC3F baseSide = norm3(cross(tempTwistUp, joint.baseAxis));
    joint.baseTwistUp = norm3(cross(joint.baseAxis, baseSide));
    joint.baseTwistUp2 = getBodyRelVec(body2, baseTwistUp);

    joint.flexLimit = flexLimit;
    joint.twistLimit = twistLimit;
    joint.flexForce = flexForce;
    joint.twistForce = twistForce;

    joint.style = string("ball");
    _joints.push_back(joint);
    _parentBodies[_bodyIndices[body2]] = _bodyIndices[body1];

    return joint;
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  JOINT addUniversalJoint(dBodyID body1, dBodyID body2, VEC3F anchor, VEC3F axis1, VEC3F axis2, 
      Real loStop1 = -dInfinity, Real hiStop1 = dInfinity,
      Real loStop2 = -dInfinity, Real hiStop2 = dInfinity)
  {
    anchor = anchor + _offset;
    dJointID jointID = dJointCreateUniversal(_world, _jointGroup);

    dJointSetUniversalAnchor(jointID, anchor[0], anchor[1], anchor[2]);
    dJointSetUniversalAxis1(jointID, axis1[0], axis1[1], axis1[2]);
    dJointSetUniversalAxis2(jointID, axis2[0], axis2[1], axis2[2]);
    dJointAttach(jointID, body1, body2);

    dJointSetUniversalParam(jointID, dParamLoStop, loStop1);
    dJointSetUniversalParam(jointID, dParamHiStop, hiStop1);
    dJointSetUniversalParam(jointID, dParamLoStop2, loStop2);
    dJointSetUniversalParam(jointID, dParamHiStop2, hiStop2);

    JOINT joint;
    joint.jointID = jointID;
    joint.style = string("univ");
    _joints.push_back(joint);
    _parentBodies[_bodyIndices[body2]] = _bodyIndices[body1];

    return joint;
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  JOINT addHingeJoint(dBodyID body1, dBodyID body2, VEC3F anchor, VEC3F axis, 
      Real loStop = dInfinity,
      Real hiStop = dInfinity)
  {
    anchor = anchor + _offset;

    dJointID jointID = dJointCreateHinge(_world, _jointGroup);
    dJointAttach(jointID, body1, body2);
    dJointSetHingeAnchor(jointID, anchor[0], anchor[1], anchor[2]);
    dJointSetHingeAxis(jointID, axis[0], axis[1], axis[2]);
    dJointSetHingeParam(jointID, dParamLoStop, loStop);
    dJointSetHingeParam(jointID, dParamHiStop, hiStop);

    JOINT joint;
    joint.jointID = jointID;
    joint.style = string("hinge");

    _joints.push_back(joint);
    _parentBodies[_bodyIndices[body2]] = _bodyIndices[body1];

    return joint;
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  void update()
  {
    for (unsigned int x = 0; x < _joints.size(); x++)
    {
      if (_joints[x].style.compare("ball") != 0) continue;

      JOINT& j = _joints[x];
      dJointID jointID = j.jointID;
     
      dBodyID body0 = dJointGetBody(jointID, 0);
      dBodyID body1 = dJointGetBody(jointID, 1);

      const dReal* rotation0 = dBodyGetRotation(body0);
      const dReal* rotation1 = dBodyGetRotation(body1);

      VEC3F baseAxis = rotate3(repack(rotation0), j.baseAxis);
      MATRIX3 repacked = repack(rotation1);
      VEC3F currAxis(repacked(0,2), repacked(1,2), repacked(2,2));
      
      // get angular velocity of attached body relative to fixed body
      VEC3F relAngVel = angularVelocity(body1) - angularVelocity(body0);
      VEC3F twistAngVel = relAngVel * (norm3(relAngVel) * currAxis);
      VEC3F flexAngVel = relAngVel - twistAngVel;

      //restrict limbs rotating too far from base axis
      Real angle = acosdot3(currAxis, baseAxis);
      if (angle > j.flexLimit)
      {
        // add torque to push body back towards base axis
        VEC3F torque =  norm3(cross(currAxis, baseAxis)) *
          (angle - j.flexLimit) * j.flexForce;
        dBodyAddTorque(body1, torque[0], torque[1], torque[2]);

        // dampen flex to prevent bounceback
        torque = flexAngVel * (-0.01f * j.flexForce);
        dBodyAddTorque(body1, torque[0], torque[1], torque[2]);
      }

      // determine the base twist up vector for the current attached
      // body by applying the current joint flex to the fixed body's
      // base twist up vector
      VEC3F baseTwistUp = rotate3(repack(rotation0), j.baseTwistUp);
      MATRIX3 base2current = calcRotMatrix(norm3(cross(baseAxis, currAxis)), acosdot3(baseAxis, currAxis));
      VEC3F projBaseTwistUp = rotate3(base2current, baseTwistUp);

      // determine the current twist up vector from the attached body
      VEC3F actualTwistUp = rotate3(repack(rotation1), j.baseTwistUp2);

      // restrict limbs twisting
      angle = acosdot3(actualTwistUp, projBaseTwistUp);
      if (angle > j.twistLimit)
      {
        // add torque to rotate body back towards base angle
        VEC3F torque = 
          norm3(cross(actualTwistUp, projBaseTwistUp)) *
          (angle - j.twistLimit) * j.twistForce;
        dBodyAddTorque(body1, torque[0], torque[1], torque[2]);

        // dampen twisting
        torque = twistAngVel * (-0.01f * j.twistForce);
        dBodyAddTorque(body1, torque[0], torque[1], torque[2]);
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  dJointGroupID& contactGroup() { return _contactGroup; };
  dSpaceID& space()             { return _space; };
  vector<dBodyID>& bodies()     { return _bodies; };
  vector<Real>& radii()         { return _radii; };
  vector<Real>& lengths()       { return _lengths; };

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  void printRigids()
  {
    for (unsigned int x = 0; x < _bodies.size(); x++)
    {
      const dReal* translation = dBodyGetPosition(_bodies[x]);
      const dReal* rotation = dBodyGetQuaternion(_bodies[x]);

      cout << " Body " << x << endl;
      cout << " Translation: " << translation[0] << " " 
                               << translation[1] << " "
                               << translation[2] << endl;
      cout << " Rotation:    " << rotation[1] << " " 
                               << rotation[2] << " "
                               << rotation[3] << " "
                               << rotation[0] << endl;
    }
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  void writeSkeleton(const char* filename)
  {
    FILE* file = fopen(filename, "w");
    if (file == NULL)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Could not open skeleton file " << filename << "!!!" << endl;
      exit(0);
    }
    cout << " Writing state file: " << filename << endl;

    // write out the number of total bones
    int totalBones = _bodies.size();
    //fprintf(file, "%i\n", totalBones);

    // get the parent, and the rigid transform
    for (int x = 0; x < totalBones; x++)
    {
      const dReal* translation = dBodyGetPosition(_bodies[x]);
      const dReal* rotation = dBodyGetQuaternion(_bodies[x]);
      //int parent = _parentBodies[x];
      int parent = -1;
      if (_parentBodies.find(x) != _parentBodies.end())
        parent = _parentBodies[x];

      fprintf(file, "%i %f %f %f %i\n", x, translation[0], translation[1], translation[2], parent);

      // use cubica quaternion ordering
      fprintf(file, "%f %f %f %f\n", rotation[1], rotation[2], rotation[3], rotation[0]);

      // dump out the capsule dimensions as well
      fprintf(file, "%f %f\n", _lengths[x], _radii[x]);
    }

    fclose(file);
  }

  ////////////////////////////////////////////////////////////////////////
  // get the translation of a specific body
  ////////////////////////////////////////////////////////////////////////
  VEC3F bodyTranslation(int body)
  {
    if (body > _bodies.size() - 1)
      return VEC3F();
    const dReal* translation = dBodyGetPosition(_bodies[body]);

    VEC3F final;
    final[0] = translation[0];
    final[1] = translation[1];
    final[2] = translation[2];

    return final;
  }

private:
  dWorldID _world;
  dSpaceID _space;
  Real _density;
  Real _totalMass;
  VEC3F _offset;

  vector<dBodyID> _bodies;
  vector<dGeomID> _geoms;
  vector<JOINT> _joints;
  vector<Real> _radii;
  vector<Real> _lengths;

  dBodyID _chest;
  dBodyID _belly;
  dBodyID _pelvis;
  dBodyID _head;
  dBodyID _rightUpperLeg;
  dBodyID _leftUpperLeg;
  dBodyID _rightLowerLeg;
  dBodyID _leftLowerLeg;
  dBodyID _rightFoot;
  dBodyID _leftFoot;
  dBodyID _rightUpperArm;
  dBodyID _leftUpperArm;
  dBodyID _rightForeArm;
  dBodyID _leftForeArm;
  dBodyID _rightHand;
  dBodyID _leftHand;
  JOINT _midSpine;
  JOINT _lowSpine;
  JOINT _neck;
  JOINT _rightHip;
  JOINT _leftHip;
  JOINT _rightKnee;
  JOINT _leftKnee;
  JOINT _rightAnkle;
  JOINT _leftAnkle;
  JOINT _rightShoulder;
  JOINT _leftShoulder;
  JOINT _rightElbow;
  JOINT _leftElbow;
  JOINT _rightWrist;
  JOINT _leftWrist;

  dJointGroupID _jointGroup;
  dJointGroupID _contactGroup;

  // keep a hash of all the bones for output
  map<dBodyID, int> _bodyIndices;

  // keep a hash of all the parent bones, i.e.
  // _parentBodies[body] gives the ID of the parent of 'body'
  map<int, int> _parentBodies;
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
RAG_DOLL* ragdoll;
dWorldID world;

vector<VEC3F> stairSizes;
vector<VEC3F> stairPositions;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void draw_body(dBodyID body, Real length, Real radius)
{
  // polygon resolution for capsule bodies
  int CAPSULE_SLICES = 16;
  int CAPSULE_STACKS = 12;

  const dReal* rotationd = dBodyGetRotation(body);
  const dReal* positiond = dBodyGetPosition(body);
  float rot[16];
  makeOpenGLMatrix(rotationd, positiond, &rot[0]);
  glPushMatrix();
    glMultMatrixf(rot);
    Real cylHalfHeight = length / 2.0;
    glBegin(GL_QUAD_STRIP);
    for (int i = 0; i < CAPSULE_SLICES + 1; i++)
    {
      Real angle = i / float(CAPSULE_SLICES) * 2.0 * M_PI;
      Real ca = cos(angle);
      Real sa = sin(angle);
      glNormal3f(ca, sa, 0);
      glVertex3f(radius * ca, radius * sa, cylHalfHeight);
      glVertex3f(radius * ca, radius * sa, -cylHalfHeight);
    }
    glEnd();
    glTranslated(0, 0, cylHalfHeight);
    glutSolidSphere(radius, CAPSULE_SLICES, CAPSULE_STACKS);
    glTranslated(0, 0, -2.0 * cylHalfHeight);
    glutSolidSphere(radius, CAPSULE_SLICES, CAPSULE_STACKS);
  glPopMatrix();
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void drawStairs()
{
  for (unsigned int x = 0; x < stairSizes.size(); x++)
  {
    glPushMatrix();
      glTranslatef(stairPositions[x][0], stairPositions[x][1], stairPositions[x][2]);
      glScalef(stairSizes[x][0], stairSizes[x][1], stairSizes[x][2]);
      glutSolidCube(1);
    glPopMatrix();
  }
}

void createStairs()
{
  float downSpacing = 1.65;
  float leftSpacing = downSpacing * 0.75;

  for (int y = 0; y < 2; y++)
    for (int x = 0; x < 10; x++)
    {
      VEC3F size(0.55,0.125,5);
     
      float jitter = (y % 2) ? 0.5 * downSpacing : 0; 

      VEC3F position(0.5 -y * leftSpacing,-1.2 - x * downSpacing - jitter, 1);

      stairSizes.push_back(size);
      stairPositions.push_back(position);

      dGeomID stair = dCreateBox(ragdoll->space(), size[0], size[1], size[2]);
      dGeomSetPosition(stair, position[0], position[1], position[2]);
    }

  {
  // left wall 
  VEC3F size(0.1, 20, 5);
  VEC3F position(0.85 ,-10, 1);

  stairSizes.push_back(size);
  stairPositions.push_back(position);
  dGeomID stair = dCreateBox(ragdoll->space(), size[0], size[1], size[2]);
  dGeomSetPosition(stair, position[0], position[1], position[2]);
  }

  {
  // right wall 
  VEC3F size(0.1, 20, 5);
  VEC3F position(-1.1 ,-10, 1);

  stairSizes.push_back(size);
  stairPositions.push_back(position);
  dGeomID stair = dCreateBox(ragdoll->space(), size[0], size[1], size[2]);
  dGeomSetPosition(stair, position[0], position[1], position[2]);
  }

  // bottom wall
  {
  VEC3F size(1.5, 1, 5);
  VEC3F position((-1 + 0.75) * 0.5 ,-20, 1);

  stairSizes.push_back(size);
  stairPositions.push_back(position);
  dGeomID stair = dCreateBox(ragdoll->space(), size[0], size[1], size[2]);
  dGeomSetPosition(stair, position[0], position[1], position[2]);
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void near_callback(void* args, dGeomID geom1, dGeomID geom2)
{
  dBodyID body1 = dGeomGetBody(geom1);
  dBodyID body2 = dGeomGetBody(geom2);

  // first make sure neither is 0, i.e. the static environment
  if (body1 && body2)
    if (dAreConnected(body1, body2))
      return;

  // check if the objects collide
  dContact contact;

  contact.surface.mode = dContactBounce | dContactSoftCFM;
  // friction parameter
  contact.surface.mu = 0;
  // bounce is the amount of "bouncyness".
  contact.surface.bounce = 0.3;
  // bounce_vel is the minimum incoming velocity to cause a bounce
  contact.surface.bounce_vel = 0.1;
  // constraint force mixing parameter
  //
  // setting to zero prevents all obstacle penetration
  contact.surface.soft_cfm = 0.0;

  int numc = dCollide (geom1,geom2,1,&contact.geom,sizeof(dContact));
  if (numc > 0)
  {
    dJointID c = dJointCreateContact (world,ragdoll->contactGroup(), &contact);
    dJointAttach (c,body1,body2);
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void onIdle()
{
  if (animate)
  {
    static int frame = 0;

    // Detect collisions and create contact joints
    dSpaceCollide(ragdoll->space(), world, near_callback);

    // Simulation step (with slo motion)
    dWorldStep(world, 1.0 / 120.0);

    // apply internal ragdoll forces
    ragdoll->update();

    // Remove all contact joints
    dJointGroupEmpty(ragdoll->contactGroup());

    // dump out skeleton state
    char buffer[256];
    sprintf(buffer, "%04i", frame);

    string filename("./data/ode.motion.");
    filename = filename + string(buffer);
    filename = filename + string(".skeleton");

    // dump out the skeleton data
    ragdoll->writeSkeleton(filename.c_str());

    // always point camera at the middle of the torso
    Camera* camera = glvu.GetCurrentCam();
    Vec3f eye;
    Vec3f lookat;
    Vec3f up;
    camera->GetLookAtParams(&eye, &lookat, &up);

    Vec3f direction = eye - lookat;
    direction *= 4;

    VEC3F center = ragdoll->bodyTranslation(0);
    lookat[0] = center[0];
    lookat[1] = center[1];
    lookat[2] = center[2];

    eye = lookat + direction;

    camera->LookAt(eye, lookat, up);

    if (singleStep)
    {
      singleStep = false;
      animate = false;
    }

    frame++;

    if (frame % 100 == 0)
    {
      cout << " Frame " << frame << endl;
      if (frame == 2000)
        exit(0);
    }
  }
	
  glutPostRedisplay();
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void onKey(unsigned char c, int x, int y)
{
  if (c == 'q' || c == 'Q')
    exit(0);

  if (c == 'a')
    animate = !animate;

  if (c == 's')
  {
    singleStep = true;
    animate = true;
  }

  if (c == 'v')
  {
    Camera* camera = glvu.GetCurrentCam();
    Vec3f eye;
    Vec3f lookat;
    Vec3f up;
    camera->GetLookAtParams(&eye, &lookat, &up);
    cout << " Eye(" << eye[0] << ", " << eye[1] << ", " << eye[2] << "), " ;
    cout << " LookAtCntr(" << lookat[0] << ", " << lookat[1] << ", " << lookat[2] << "), " ;
    cout << " Up(" << up[0] << ", " << up[1] << ", " << up[2] << ");" << endl;
  }

  glvu.Keyboard(c,x,y);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void onDraw()
{
  glvu.BeginFrame();
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_MULTISAMPLE);

  vector<dBodyID>& bodies = ragdoll->bodies();
  vector<Real> radii = ragdoll->radii();
  vector<Real> lengths = ragdoll->lengths();

  for (unsigned int x = 0; x < bodies.size(); x++)
    draw_body(bodies[x], lengths[x], radii[x]); 

  drawStairs();

  glvu.EndFrame();
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void createCapsule(dWorldID world, dSpaceID space, Real density, Real length, Real radius, dBodyID& body, dGeomID& geom)
{
	body = dBodyCreate(world);
  dMass m;
  dMassSetCapsule(&m, density, 3, radius, length);
  dBodySetMass (body, &m);

  geom = dCreateCapsule(space, radius, length);
  dGeomSetBody(geom, body);
}

//////////////////////////////////////////////////////////////////////////////
// open the GLVU window
//////////////////////////////////////////////////////////////////////////////
int glvuWindow()
{
  glvu.Init("Ragdoll Viewer",
            GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH,
            //0, 0, 640, 400);
            0, 0, 800, 600);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  float position[] = {0,0,1,0};
  Real white = 0.75;
  float diffuse[] = {white, white, white, 1.0};
	glLightfv(GL_LIGHT0,GL_POSITION,position);
	glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuse);
	glLightfv(GL_LIGHT0,GL_SPECULAR,diffuse);
	glEnable(GL_LIGHT0);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
  glEnable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);
	glClearColor(0.8, 0.8, 0.9, 0.0);

  glutDisplayFunc(onDraw);
  glutKeyboardFunc(onKey);
  glutIdleFunc(onIdle);

  Vec3f ModelMin(-10,-10,-10), ModelMax(10,10,10), 
        Eye(0.0297215, -0.221213, -3.88846),  LookAtCntr(-0.00337817, -0.221785, -2.88901),  Up(-0.0577071, 0.998334, -0.00134061);

  float Yfov = 45;
  float Aspect = 1;
  float Near = 0.001f;
  float Far = 10.0f;
  glvu.SetAllCams(ModelMin, ModelMax, Eye,LookAtCntr,Up, Yfov,Aspect, Near,Far);

  // reset center
  Vec3f center(0.5f, 0.5f, 0.5f);
  glvu.SetWorldCenter(center);
  
  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
int main (int argc, char **argv)
{
  //create an ODE world object
  dInitODE2(0);
  world = dWorldCreate();

  dWorldSetGravity(world, 0, -9.81,0); 
  dWorldSetERP(world, 0.1);
  dWorldSetCFM(world, 1e-4);

  // create an ODE space object
  dSpaceID space = dHashSpaceCreate(0);

  // create the rag doll
  ragdoll = new RAG_DOLL(world, space, 500);

  // create stairs to fall down
  createStairs();

  glutInit(&argc, argv);
  glvuWindow();

  dWorldDestroy (world);
  dCloseODE();
  return 0;
}
