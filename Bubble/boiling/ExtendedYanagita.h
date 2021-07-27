/* (C)Copyright 2007                                               */
/* International Business Machines Corporation,                    */
/* All Rights Reserved.                                            */
/*                                                                 */
/* This program is made available under the terms of the           */
/* Common Public License v1.0 which accompanies this distribution. */
/* --------------------------------------------------------------- */
/* PROLOG END TAG zYx                                              */

////////////////////////////////////////////////////////////////////////////
// This project is a 2D implementation of the boiling module from:
//
//   T. Kim and M. Carlson, A Simple Boiling Module. Proceedings of
//    ACM SIGGRAPH / Eurographics Symposium on Computer Animation 2007
//
// This file contains the actual implementation of the module.
//
// Please direct any questions or comments to Ted Kim: twkim@us.ibm.com
////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Calculate local curvature
//////////////////////////////////////////////////////////////////////
void calcCurve(int i, int offset, float& curve, int& seen) {
  if (crossings[i + offset]) {
    float diffx = (offset > 0) ? normalX[i] - normalX[i + offset] : 
                                 normalX[i + offset] - normalX[i];
    float diffy = (offset > 0) ? normalY[i] - normalY[i + offset] : 
                                 normalY[i + offset] - normalY[i];
    float magnitude = diffx * diffx + diffy * diffy;

    if (offset == 1) diffx -= 1.0f;
    else diffy -= 1.0f;
    
    float sign = (diffx * diffx + diffy * diffy < 1.0f) ? -1.0f : 1.0f;
    curve += (crossings[i]) ? sign * sqrtf(magnitude) : 0.0f;
    seen++;
  }
}

//////////////////////////////////////////////////////////////////////
// Extended Yanagita model designed in the paper
//////////////////////////////////////////////////////////////////////
void stepExtendedYanagita() {
  // step heat diffusion
  float beta = 1.0f / (1.0f + 4.0f * diffConst);
  for (int z = 0; z < 20; z++)
    for (int y = 1, i = res + 1; y < res - 1; y++, i += 2)
      for (int x = 1; x < res - 1; x++, i++)
        heat1[i] = (heat[i] + (heat1[i - 1]   + heat1[i + 1] + 
                               heat1[i - res] + heat1[i + res]) * diffConst) * beta;

  // step buoyancy
  for (int y = 1, i = res + 1; y < res - 1; y++, i += 2)
    for (int x = 1; x < res - 1; x++, i++)
      heat2[i] = heat1[i] * (1.0f - (tanh(alpha * (heat1[i + res] - Tc)) - 
                                     tanh(alpha * (heat1[i - res] - Tc))) * buoyConst);
 
  // step latent heat
  for (int y = 1, i = res + 1; y < res - 1; y++, i += 2)
    for (int x = 1; x < res - 1; x++, i++) {
      float delta = 0.0f;
      delta = (heat2[i] < Tc && heat[i] > Tc) ? eta : delta;
      heat2[i] += delta * latentConst;
    }
  // swap arrays
  float* temp = heat; heat = heat2; heat2 = temp;

  // calc phase
  for (int x = 0; x < res * res; x++)
    phase[x] = (tanh(alpha * (heat[x] - bottom)) + 0.5f) * 2.0f;

  // calc crossings and normals
  for (int y = 1, i = res + 1; y < res - 1; y++, i += 2)
    for (int x = 1; x < res - 1; x++, i++) {
      crossings[i] = (phase[i] * phase[i + 1] <= 0.0f || 
                      phase[i] * phase[i - 1] <= 0.0f || 
                      phase[i] * phase[i + res] <= 0.0f || 
                      phase[i] * phase[i - res] <= 0.0f) ? true : false;
      if (crossings[i]) {
        float dx = phase[i + 1]   - phase[i - 1];
        float dy = phase[i + res] - phase[i - res];
        float magnitude = sqrtf(dx * dx + dy * dy);
        magnitude = (magnitude > 0.0f) ? -1.0f / magnitude : 0.0f;
        normalX[i] = dx * magnitude;
        normalY[i] = dy * magnitude;
      }
    }

  // calc curvature (store it in heat2)
  for (int y = 1, i = res + 1; y < res - 1; y++, i += 2)
    for (int x = 1; x < res - 1; x++, i++) {
      float curve = 0.0f;
      int seen = 0;
      calcCurve(i, 1,    curve, seen);
      calcCurve(i, res,  curve, seen);
      calcCurve(i, -1,   curve, seen);
      calcCurve(i, -res, curve, seen);
      heat2[i] = (seen != 0) ? curve / seen : 0.0f;
    }

  // add surface tension
  for (int y = 1, i = res + 1; y < res - 1; y++, i += 2)
    for (int x = 1; x < res - 1; x++, i++)
      heat[i] += (crossings[i]) ? -heat2[i] * tensionConst : 0.0f;
}
