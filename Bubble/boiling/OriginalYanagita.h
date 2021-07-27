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
// For comparison purposes, this file contains an implementation of:
//
//   Yanagita, T. Phenomenology of boiling: A coupled map lattice model. 
//    Chaos, 2 3. 343-350. 1992.
//
// Please direct any questions or comments to Ted Kim: twkim@us.ibm.com
////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Original Yanagita model
//////////////////////////////////////////////////////////////////////
void stepOriginalYanagita() {
  // set bottom to hot
  //for (int x = 0; x < factor * res; x++) 
  //  heat[x] = heat1[x] = heat2[x] = bottom;

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
    for (int x = 1; x < res - 1; x++, i++)  {
      float delta = 0.0f;      
      delta += (heat2[i + 1] < Tc && heat[i + 1] > Tc) ? eta : 0.0;
      delta += (heat2[i - 1] < Tc && heat[i - 1] > Tc) ? eta : 0.0;
      delta += (heat2[i + res] < Tc && heat[i + res] > Tc) ? eta : 0.0;
      delta += (heat2[i - res] < Tc && heat[i - res] > Tc) ? eta : 0.0;
      heat2[i] += delta * latentConst;
    }

  // swap arrays
  float* temp = heat; heat = heat2; heat2 = temp;

  // calc phase
  for (int x = 0; x < res * res; x++)
    phase[x] = (tanh(alpha * (heat[x] - bottom)) + 0.5f) * 2.0f;

  // calc crossings and normals -- just so they can render too
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
}
