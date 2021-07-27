#ifndef NONLINEAR_SLICE_CUDA_H 
#define NONLINEAR_SLICE_CUDA_H

#include "SETTINGS.h"

// pull in CUDA types
#include "vector_types.h"

// CUDA can't handle this many
#define MAX_COEFFS 200

void checkForCudaErrors();
cudaError_t setTopRoots(double4* coeffs, unsigned int total);
cudaError_t setBottomRoots(double4* coeffs, unsigned int total);
cudaError_t setTopPowers(double* powers, unsigned int total);
cudaError_t setBottomPowers(double* powers, unsigned int total);

void setConsts(const int xResHost, const int yResHost, const int zResHost, 
               const Real dxHost, const Real dyHost, const Real dzHost,
               const Real xCenterHost, const Real yCenterHost, const Real zCenterHost,
               const Real xLengthHost, const Real yLengthHost, const Real zLengthHost,
               const Real escapeHost, const int maxIterationsHost, 
               const Real expScalingHost, const Real quaternionSliceHost,
               const Real isosurface);

void computeSliceOnGPU(const int xRes,
                       const int yRes,
                       const int z,
                       double* field);

void computeEdgesOnGPU(const int3* firstVertices,
                       const int3* secondVertices,
                       const int size,
                       double3* finalVertices);

int countFlagsOnGPU(const int xRes,
                    const int yRes,
                    const double* slab0,
                    const double* slab1,
                    unsigned char* flags,
                    int* occupancies,
                    int* indexTranslation,
                    int* compactedFlags);

void computeFlagsOnGPU(const int xRes, 
                       const int yRes,
                       const double* slab0,
                       const double* slab1,
                       unsigned char* flags);

#endif
