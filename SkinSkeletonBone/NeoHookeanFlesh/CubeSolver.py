import numpy as np
from Cube import *

class CubeSolver:
    
    def __init__(self,dofs,Nvertex,Ncube):
        self.DOF = dofs
        self.u = np.zeros((dofs))
        self.uOld = np.zeros((dofs))
        self.uDelta = np.zeros((dofs))
        self.f = np.zeros((dofs))
        self.K = np.zeros((dofs,dofs))
        self.vertexToSimulatedDof = np.zeros((Nvertex))
        self.elementEnergy = np.zeros((Ncube))
        self.elementForce = np.zeros((Ncube,24))
        self.localIK = np.zeros((Ncube,24,24))
        self.cubes = [Cube(index, pos)] * 1000
        
    def InitValue(self):
        curDof = 0
        Nx = 11
        for k in range(Nx):
            for j in range(Nx):
                for i in range(Nx):
                    idx = k * Nx * Nx + j * Nx + i
                    if j == 0 or j == Nx - 1:
                        vertexToSimulatedDof[idx] = - 1
                    else:
                        vertexToSimulatedDof[idx] = curDof
                    curDof += 3

    def ComputeForces(f):
        
        for cidx in range(self.Ncube):
            self.elementForce[cidx] = self.cubes.ComputeForces()
            for y in range(8):
                vidx = self.cubes[cidx].Vertex[y]
                if self.vertexToSimulatedDof[vidx] == -1:
                    continue
                f[vidx,:] = elementEnergy[3*y:3*(y+1)]
        
        