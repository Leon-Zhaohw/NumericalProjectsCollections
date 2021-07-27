#include <cstdlib>
#include <cassert>

#include "helper_cuda.h"
#include "NONLINEAR_SLICE_CUDA.h"
#include "TIMER.h"

// stop g++ from choking when Thrust is included
#undef _GLIBCXX_ATOMIC_BUILTINS
#undef _GLIBCXX_USE_INT128

#include <nppi.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>

__constant__ double4 topRoots[MAX_COEFFS];
__constant__ double4 bottomRoots[MAX_COEFFS];
__constant__ double topPowers[MAX_COEFFS];
__constant__ double bottomPowers[MAX_COEFFS];

__constant__ unsigned int totalTopRoots;
__constant__ unsigned int totalBottomRoots;

__constant__ int xRes;
__constant__ int yRes;
__constant__ int zRes;
__constant__ Real xLength;
__constant__ Real yLength;
__constant__ Real zLength;
__constant__ Real dx;
__constant__ Real dy;
__constant__ Real dz;
__constant__ Real xCenter;
__constant__ Real yCenter;
__constant__ Real zCenter;
__constant__ Real escape;
__constant__ Real expScaling;
__constant__ Real quaternionSlice;
__constant__ Real isosurface;
__constant__ int maxIterations;
__constant__ int binaryBandwidth;
__constant__ int curvatureBandwidth;

__device__ unsigned int blockCounter;   // global counter, initialized to zero before kernel launch

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void setConsts(const int xResHost, const int yResHost, const int zResHost,
               const Real dxHost, const Real dyHost, const Real dzHost,
               const Real xCenterHost, const Real yCenterHost, const Real zCenterHost,
               const Real xLengthHost, const Real yLengthHost, const Real zLengthHost,
               const Real escapeHost, const int maxIterationsHost, const Real expScalingHost,
               const Real quaternionSliceHost, const Real isosurfaceHost)
{
  cudaMemcpyToSymbol(xRes, &xResHost, sizeof(int));
  getLastCudaError("xRes copy failed");
  cudaMemcpyToSymbol(yRes, &yResHost, sizeof(int));
  getLastCudaError("yRes copy failed");
  cudaMemcpyToSymbol(zRes, &zResHost, sizeof(int));
  getLastCudaError("zRes copy failed");

  cudaMemcpyToSymbol(dx, &dxHost, sizeof(Real));
  getLastCudaError("dx copy failed");
  cudaMemcpyToSymbol(dy, &dyHost, sizeof(Real));
  getLastCudaError("dy copy failed");
  cudaMemcpyToSymbol(dz, &dzHost, sizeof(Real));
  getLastCudaError("dz copy failed");

  cudaMemcpyToSymbol(xCenter, &xCenterHost, sizeof(Real));
  getLastCudaError("xCenter copy failed");
  cudaMemcpyToSymbol(yCenter, &yCenterHost, sizeof(Real));
  getLastCudaError("yCenter copy failed");
  cudaMemcpyToSymbol(zCenter, &zCenterHost, sizeof(Real));
  getLastCudaError("zCenter copy failed");

  cudaMemcpyToSymbol(xLength, &xLengthHost, sizeof(Real));
  getLastCudaError("xLength copy failed");
  cudaMemcpyToSymbol(yLength, &yLengthHost, sizeof(Real));
  getLastCudaError("yLength copy failed");
  cudaMemcpyToSymbol(zLength, &zLengthHost, sizeof(Real));
  getLastCudaError("zLength copy failed");

  cudaMemcpyToSymbol(escape, &escapeHost, sizeof(Real));
  getLastCudaError("escape copy failed");

  cudaMemcpyToSymbol(maxIterations, &maxIterationsHost, sizeof(int));
  getLastCudaError("maxIterations copy failed");
  
  cudaMemcpyToSymbol(expScaling, &expScalingHost, sizeof(Real));
  getLastCudaError("expScaling copy failed");
  
  cudaMemcpyToSymbol(quaternionSlice, &quaternionSliceHost, sizeof(Real));
  getLastCudaError("quaternionSlice copy failed");
  
  cudaMemcpyToSymbol(isosurface, &isosurfaceHost, sizeof(Real));
  getLastCudaError("isosurface copy failed");
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
cudaError_t setTopRoots(double4* roots, unsigned int total)
{
  assert(total < MAX_COEFFS);
  cudaMemcpyToSymbol(totalTopRoots, &total,
                     sizeof(unsigned int));

  return cudaMemcpyToSymbol(topRoots, roots,
                            sizeof(double4) * total);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
cudaError_t setTopPowers(double* powers, unsigned int total)
{
  assert(total < MAX_COEFFS);

  return cudaMemcpyToSymbol(topPowers, powers,
                            sizeof(double) * total);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
cudaError_t setBottomRoots(double4* roots, unsigned int total)
{
  assert(total < MAX_COEFFS);
  cudaMemcpyToSymbol(totalBottomRoots, &total,
                     sizeof(unsigned int));

  return cudaMemcpyToSymbol(bottomRoots, roots,
                            sizeof(double4) * total);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
cudaError_t setBottomPowers(double* powers, unsigned int total)
{
  assert(total < MAX_COEFFS);

  return cudaMemcpyToSymbol(bottomPowers, powers,
                            sizeof(double) * total);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
__device__ inline double4 make_quaternion(const double& w,
                                          const double& x,
                                          const double& y,
                                          const double& z)
{
  return make_double4(x,y,z,w);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
__device__ inline double magnitudeSq4(const double4& data)
{
  return data.x * data.x + data.y * data.y + data.z * data.z + data.w * data.w;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
__device__ inline double magnitude4(const double4& data)
{
  return sqrt(data.x * data.x + data.y * data.y + data.z * data.z + data.w * data.w);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
__device__ inline double4 conjugateScaled(const double4& data, const double& scale)
{
  return make_quaternion(scale * data.w, 
                        -scale * data.x, 
                        -scale * data.y, 
                        -scale * data.z);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
__device__ inline double4 conjugate(const double4& data)
{
  return make_quaternion(data.w, -data.x, -data.y, -data.z);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
__device__ inline void multiply(const double4& left, const double4& right, double4& result)
{
  result.x = left.y * right.z - left.z * right.y + right.w * left.x + left.w * right.x;
  result.y = left.z * right.x - left.x * right.z + right.w * left.y + left.w * right.y;
  result.z = left.x * right.y - left.y * right.x + right.w * left.z + left.w * right.z;
  result.w = left.w * right.w - left.x * right.x - right.y * left.y - left.z * right.z;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
__device__ inline double3 add(const double3& left, const double3& right)
{
  double3 result;
  result.x = left.x + right.x;
  result.y = left.y + right.y;
  result.z = left.z + right.z;
  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
__device__ inline void subtract(const double4& left, const double4& right, double4& result)
{
  result.x = left.x - right.x;
  result.y = left.y - right.y;
  result.z = left.z - right.z;
  result.w = left.w - right.w;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
__device__ inline void scale(double4& left, const double& factor)
{
  left.x *= factor;
  left.y *= factor;
  left.z *= factor;
  left.w *= factor;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
__device__ inline void scale(double3& left, const double& factor)
{
  left.x *= factor;
  left.y *= factor;
  left.z *= factor;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
__device__ double4 pow(const double4& q, const double& exponent)
{
  const double partial = q.x * q.x + q.y * q.y + q.z * q.z;
  const double qMagnitude = sqrt(partial + q.w * q.w);
  const double vMagnitude = sqrt(partial);
  const double vMagnitudeInv = (vMagnitude > 0.0) ? 1.0 / vMagnitude : 0.0;

  const double scale = exponent * acos(q.w / qMagnitude) * vMagnitudeInv;

  const double magnitude = scale * vMagnitude;
  const double magnitudeInv = (magnitude > 0.0) ? 1.0 / magnitude : 0.0;

  const double exps = std::exp(exponent * std::log(qMagnitude));

  double sMag,cMag;
  sincos(magnitude, &sMag, &cMag);
  const double scale2 = scale * exps * magnitudeInv * sMag;
  return make_quaternion(exps * cMag, 
                         scale2 * q.x, 
                         scale2 * q.y, 
                         scale2 * q.z);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
__device__ inline double4 evaluateBottom(const double4& iterate)
{
  double4 result, temp, temp2;

  result = make_quaternion(1.0, 0.0, 0.0, 0.0);

  for (int x = 0; x < totalBottomRoots; x++)
  {
    const double4& root = bottomRoots[x];
    temp.x = iterate.x - root.x;
    temp.y = iterate.y - root.y;
    temp.z = iterate.z - root.z;
    temp.w = iterate.w - root.w;
    
    // raise temp to a power
    temp2 = pow(temp, bottomPowers[x]);

    multiply(result, temp2, temp);
    result = temp;
  }

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
__device__ inline double4 evaluateTop(const double4& iterate)
{
  double4 result, temp, temp2;
  subtract(iterate, topRoots[0], temp);

  // raise temp to a power
  result = pow(temp, topPowers[0]);

  for (int x = 1; x < totalTopRoots; x++)
  {
    const double4& root = topRoots[x];
    temp.x = iterate.x - root.x;
    temp.y = iterate.y - root.y;
    temp.z = iterate.z - root.z;
    temp.w = iterate.w - root.w;
    
    // raise temp to a power
    temp2 = pow(temp, topPowers[x]);

    multiply(result, temp2, temp);
    result = temp;
  }

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
__device__ inline double3 cellCenter(const int x, const int y, const int z)
{
  return make_double3(xCenter - 0.5 * xLength + x * dx + dx * 0.5,
                      yCenter - 0.5 * yLength + y * dy + dy * 0.5,
                      zCenter - 0.5 * zLength + z * dz + dz * 0.5);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void __global__ computeCompactedFlags(const int* flags, 
                                      const int* indexTranslation, 
                                      int* compactedFlags)
{
  const int x = blockDim.x * blockIdx.x + threadIdx.x;
  const int y = blockDim.y * blockIdx.y + threadIdx.y;

  // why doesn't this fix non-power-of-two problems?
  if (x >= xRes || y >= yRes) return;
  const int index = x + y * xRes;
  if (flags[index] == -1) return;

  compactedFlags[indexTranslation[index]] = flags[index];
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void __global__ computeFlags(const double* slab0, 
                             const double* slab1, 
                             unsigned char* flags,
                             int* occupancies)
{
  const int x = blockDim.x * blockIdx.x + threadIdx.x;
  const int y = blockDim.y * blockIdx.y + threadIdx.y;

  // why doesn't this fix non-power-of-two problems?
  //if (x >= xRes || y >= yRes) return;
  if (x >= xRes - 1 || y >= yRes - 1) return;

  const int index = x + y * xRes;
  //const double NNN = slab0(x,y);
  const double NNN = slab0[index];
  //const double NNP = slab1(x,y);
  const double NNP = slab1[index];
  //const double NPN = slab0(x,y + 1);
  const double NPN = slab0[index + xRes];
  //const double NPP = slab1(x,y + 1);
  const double NPP = slab1[index + xRes];
  //const double PNN = slab0(x + 1,y);
  const double PNN = slab0[index + 1];
  //const double PNP = slab1(x + 1,y);
  const double PNP = slab1[index + 1];
  //const double PPN = slab0(x + 1,y + 1);
  const double PPN = slab0[index + 1 + xRes];
  //const double PPP = slab1(x + 1,y + 1);
  const double PPP = slab1[index + 1 + xRes];
        
  unsigned char flag =    ((NNN > 0) + 2 *   (NNP > 0) + 4  * (NPN > 0) +
                       8 * (NPP > 0) + 16 *  (PNN > 0) + 32 * (PNP > 0) +
                       64 *(PPN > 0) + 128 * (PPP > 0));

  flags[index] = flag;

  int occupied = (flag == 0 || flag == 255) ? 0 : 1;
  occupancies[index] = occupied;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void __global__ computeFlags(const double* slab0, const double* slab1, unsigned char* flags)
{
  const int x = blockDim.x * blockIdx.x + threadIdx.x;
  const int y = blockDim.y * blockIdx.y + threadIdx.y;

  // why doesn't this fix non-power-of-two problems?
  //if (x >= xRes || y >= yRes) return;
  if (x >= xRes - 1 || y >= yRes - 1) return;

  const int index = x + y * xRes;
  //const double NNN = slab0(x,y);
  const double NNN = slab0[index];
  //const double NNP = slab1(x,y);
  const double NNP = slab1[index];
  //const double NPN = slab0(x,y + 1);
  const double NPN = slab0[index + xRes];
  //const double NPP = slab1(x,y + 1);
  const double NPP = slab1[index + xRes];
  //const double PNN = slab0(x + 1,y);
  const double PNN = slab0[index + 1];
  //const double PNP = slab1(x + 1,y);
  const double PNP = slab1[index + 1];
  //const double PPN = slab0(x + 1,y + 1);
  const double PPN = slab0[index + 1 + xRes];
  //const double PPP = slab1(x + 1,y + 1);
  const double PPP = slab1[index + 1 + xRes];
        
  unsigned char flag =    ((NNN > 0) + 2 *   (NNP > 0) + 4  * (NPN > 0) +
                       8 * (NPP > 0) + 16 *  (PNN > 0) + 32 * (PNP > 0) +
                       64 *(PPN > 0) + 128 * (PPP > 0));

  flags[index] = flag;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
__device__ inline double nonlinearValue(const double3& center)
{
  double4 iterate = make_quaternion(center.x, center.y, center.z, 
                                    quaternionSlice);
  double magnitude = magnitude4(iterate);
  int totalIterations = 0;

  double bottomPowerSum = bottomPowers[0];
  for (int x = 1; x < totalBottomRoots; x++)
    bottomPowerSum += bottomPowers[x];

  const double4 bail = make_double4(DBL_MAX, DBL_MAX, 
                                    DBL_MAX, DBL_MAX);
  bool bailed = false;
  while (magnitude < escape && totalIterations < maxIterations)
  {
    const double4 topEval = evaluateTop(iterate);
    const double4 bottomEval = evaluateBottom(iterate);
    const double bottomMagnitudeSq = magnitudeSq4(bottomEval);
    const double bottomMagnitude = sqrt(bottomMagnitudeSq);

    // compute guard value from Eqn. 2
    const double rhs = bottomPowerSum + 
                       log(bottomMagnitude) / log(10.0);
    const double topLimit = 308.0 - 1.05 * rhs;

    // if the division is tiny, bail
    if (bottomMagnitude < pow(10.0, topLimit))
    {
      const double4 bottomInv = 
        conjugateScaled(bottomEval, 1.0 / bottomMagnitudeSq);
      multiply(topEval, bottomInv, iterate);
    }
    else
    {
      iterate = bail;
      bailed = true;
    }

    const double scaling = (!bailed) ? expScaling : 1.0;
    scale(iterate, scaling);
    magnitude = (!bailed) ? magnitude4(iterate) : DBL_MAX;
    totalIterations++;
  }

  return log(magnitude);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void __global__ computeSlice(const int z, double* field)
{
  const int x = blockDim.x * blockIdx.x + threadIdx.x;
  const int y = blockDim.y * blockIdx.y + threadIdx.y;

  if (x >= xRes || y >= yRes) return;

  const int index = x + y * xRes;
  const double3 center = cellCenter(x,y,z);

  field[index] = nonlinearValue(center);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
double3 __device__ midpointSearchForLoop(const double3& positiveVertex, 
                                         const double& positiveValue, 
                                         const double3& negativeVertex, 
                                         const double& negativeValue)
{
  double3 pVert = positiveVertex;
  double3 nVert = negativeVertex;

  double3 midpointVertex;
  double midpointValue;

  for (int x = 0; x <= 6; x++)
  {
      midpointVertex = add(pVert, nVert);
      scale(midpointVertex, 0.5);

      midpointValue = nonlinearValue(midpointVertex);
      if (fabs(midpointValue) < 1e-8) return midpointVertex;

      if (midpointValue < 0)
        nVert = midpointVertex;
      else
        pVert = midpointVertex;
  }
  return midpointVertex;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void __global__ computeEdges(const int3* firstVertices, 
                             const int3* secondVertices,
                             const int size,
                             double3* finalVertices)
{
  const int index  = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= size) return;
    
  const int3& firstXYZ = firstVertices[index];
  const int3& secondXYZ = secondVertices[index];
  double3 firstVertex = cellCenter(firstXYZ.x, firstXYZ.y, firstXYZ.z);
  double3 secondVertex = cellCenter(secondXYZ.x, secondXYZ.y, secondXYZ.z);

  double firstValue = nonlinearValue(firstVertex);
  double secondValue = nonlinearValue(secondVertex);

  double3 positiveVertex = firstVertex;
  double positiveValue = firstValue;
  double3 negativeVertex = secondVertex;
  double negativeValue = secondValue;

  if (firstValue < 0)
  {
    positiveVertex = secondVertex;
    positiveValue = secondValue;
    negativeVertex = firstVertex;
    negativeValue = firstValue;
  }

  // this turns the midpoint search on and off. If you want to compare to just traditional
  // marching cubes with linear interpolation, set this to 0.
  double3 finalVertex = midpointSearchForLoop(positiveVertex, positiveValue, negativeVertex, negativeValue);
  finalVertices[index] = finalVertex;
}

#define BLOCK_SIZE 16
#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16
void computeSliceOnGPU(const int xRes, 
                       const int yRes,
                       const int z, 
                       double* field)
{
    TIMER functionTimer(__FUNCTION__);
    dim3 gridSize((xRes + BLOCK_SIZE_X - 1) / BLOCK_SIZE_X, (yRes + BLOCK_SIZE_Y - 1) / BLOCK_SIZE_Y);
    dim3 blockSize(BLOCK_SIZE_X, BLOCK_SIZE_Y);
    computeSlice<<<gridSize, blockSize>>>(z, field);

    cudaDeviceSynchronize();
}

void computeFlagsOnGPU(const int xRes, 
                       const int yRes,
                       const double* slab0,
                       const double* slab1,
                       unsigned char* flags)
{
    dim3 gridSize((xRes + BLOCK_SIZE_X - 1) / BLOCK_SIZE_X, (yRes + BLOCK_SIZE_Y - 1) / BLOCK_SIZE_Y);
    dim3 blockSize(BLOCK_SIZE_X, BLOCK_SIZE_Y);
    computeFlags<<<gridSize, blockSize>>>(slab0, slab1, flags);
}

void computeEdgesOnGPU(const int3* firstVertices,
                       const int3* secondVertices,
                       const int size,
                       double3* finalVertices)
{
    TIMER functionTimer(__FUNCTION__);
    if (size == 0) return;

    computeEdges<<<size / 64 + 1,64>>>(firstVertices, secondVertices, size, finalVertices);
    cudaDeviceSynchronize();
}

int countFlagsOnGPU(const int xRes,
                    const int yRes,
                    const double* slab0,
                    const double* slab1,
                    unsigned char* flags,
                    int* occupancies,
                    int* indexTranslation,
                    int* compactedFlags)
{
    TIMER functionTimer(__FUNCTION__);
    dim3 gridSize((xRes + BLOCK_SIZE_X - 1) / BLOCK_SIZE_X, (yRes + BLOCK_SIZE_Y - 1) / BLOCK_SIZE_Y);
    dim3 blockSize(BLOCK_SIZE_X, BLOCK_SIZE_Y);
    computeFlags<<<gridSize, blockSize>>>(slab0, slab1, flags, occupancies);

    TIMER reduceTimer("Thrust::reduce");
    thrust::device_ptr<int> oPtr(occupancies);
    int finalCount = thrust::reduce(oPtr, oPtr + xRes * yRes);
    reduceTimer.stop();

    return finalCount;
}
