import numpy as np
from Cube import *
resolution = 4

Nx = resolution + 1
Nvertex = Nx**3
Ncube = resolution**3
mu = 1
lam = 10

assert mu > 0
assert lam > 0
assert (mu / lam) > 0

dx = 2 / resolution
cubes = [Cube() for i in range(Ncube)]
vertexPos = np.zeros((Nvertex,3))
vertexf = np.zeros((Nvertex,3))
fixedvertex = np.zeros((Nvertex),dtype = bool)
dofcnt = 0
for k in range(Nx):
    for j in range(Nx):
        for i in range(Nx):
            idx = int(k*Nx*Nx + j*Nx + i)
            vertexPos[idx,0] = i * dx - 1
            vertexPos[idx,1] = j * dx - 1
            vertexPos[idx,2] = k * dx - 1
            if j == 0 or j == Nx - 1:
                fixedvertex[idx] = True
            else:
                fixedvertex[idx] = False
vertexPosRest = vertexPos.copy()

for k in range(resolution):
    for j in range(resolution):
        for i in range(resolution):
            cidx = k*resolution*resolution + j*resolution + i
            pidx = k*Nx*Nx + j*Nx + i
            cubes[cidx].generate(pidx, Nx, dx)
            for p in range(8):
                idx = cubes[cidx].vertexIndex[p]
                cubes[cidx].vertexPos[p,:] = vertexPos[idx,:]
            cubes[cidx].ComputeValue()
            
Ndofs = 3 * (Nx - 2) * Nx * Nx
fvec = np.zeros((Ndofs))
Kmat = np.zeros((Ndofs,Ndofs))

def dotmat(mat0,mat1):
    n = mat0.shape[0]
    m = mat0.shape[1]
    res = 0
    for i0 in range(n):
        for j0 in range(m):
            res += mat0[i0,j0] * mat1[i0,j0]
    return res

ite = 0
ite_max = 100
dt = 0.1
while(ite < ite_max):
    ite += 1
    vertexf[:,:] = 0
    # UpdateBoundary
    for k in range(Nx):
        for i in range(Nx):
           j = 0
           idx = int(k*Nx*Nx + j*Nx + i)
           vertexPos[idx,1] = - 1 - ite * dt
           j = Nx - 1
           idx = int(k*Nx*Nx + j*Nx + i)
           vertexPos[idx,1] = 1 + ite * dt
    # 妈的，这玩意怎么传到Cube里面去？
    for cidx in range(Ncube):
        for i in range(8):
            idx = cubes[cidx].vertexIndex[i]
            cubes[cidx].vertexPos[i,:] = vertexPos[idx,:]
        cubes[cidx].ComputeForce()
        for i in range(8):
            vidx = cubes[cidx].vertexIndex[i]
            if fixedvertex[vidx] == False:
                vertexf[vidx,:] += cubes[cidx].elementForce[i*3:(i+1)*3]
                
    residual = dotmat(vertexf, vertexf)
    if residual < 1e-6:
        break
    
    # ComputeStiffnessMatrixSparse 之中的Local完成了
    for cidx in range(Ncube):
        cubes[cidx].ComputeForceJacobian()
            
            