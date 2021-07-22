import numpy as np
import scipy.io as scio
import math
"""
Elastoviscoplastic

D:\FluidSim\FluidSim\FEMNEW\2002-CC_KR-Elastoviscoplastic_FE_Analysis_in_Matlab\Software3\example4

完成状态：第115行有bug.n
"""
coordinates = scio.loadmat('coordinates.mat')['coordinates']
elements = scio.loadmat('elements.mat')['elements'] - 1
neumann = scio.loadmat('neumann.mat')['neumann'] - 1
dirichlet = scio.loadmat('dirichlet.mat')['dirichlet'] - 1
lam=1.107438169066076e+05
mu=8.019379844961240e+04
C1=lam+2*mu/3
sigma_y=450
th = 1
dt = 0.05
nu = 0
sigma = np.zeros((32,9))
for t in range(1):
    C2=nu/(nu/(2*mu)+th*dt)
    C3=th*dt*sigma_y/(nu/(2*mu)+th*dt);
    trsigma = sigma[:,0] + sigma[:,4] + sigma[:,8]
    devsigma = np.zeros((32,9))
    devarray = np.array([1,0,0,0,1,0,0,0,1])
    e0 = np.zeros((32,9))
    for i in range(32):
        for j in range(9):
            devsigma[i,j] = sigma[i,j] - trsigma[i] * devarray[j]
            e0[i,j] = trsigma[i] * devarray[j] / (9 * lam + 6 * mu) + devsigma[i,j] / (2 * mu)
            
    N = coordinates.shape[0]
    DF = np.zeros((3 * N,3 * N))
    Q = np.zeros((3 * N))
    P = np.zeros((3 * N))
    numOfElem = elements.shape[0]
    
    u0 = np.zeros((12))
    u1 = np.zeros((12))
    
    for k in range(numOfElem):
        x0 = coordinates[elements[k,0],0]
        y0 = coordinates[elements[k,0],1]
        z0 = coordinates[elements[k,0],2]
        x1 = coordinates[elements[k,1],0]
        y1 = coordinates[elements[k,1],1]
        z1 = coordinates[elements[k,1],2]
        x2 = coordinates[elements[k,2],0]
        y2 = coordinates[elements[k,2],1]
        z2 = coordinates[elements[k,2],2]
        x3 = coordinates[elements[k,3],0]
        y3 = coordinates[elements[k,3],1]
        z3 = coordinates[elements[k,3],2]
        
        G = np.array([[1,1,1,1],[x0,x1,x2,x3],[y0,y1,y2,y3],[z0,z1,z2,z3]])
        temp0 = np.linalg.inv(G)
        PhiGrad = np.dot(temp0,np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]))
        
        Deta1 = np.zeros((12,9))
        idx0 = np.array([0,3,6,9],dtype = int)
        idx1 = np.array([0,1,2],dtype = int)
        for i in range(4):
            for j in range(3):
                Deta1[idx0[i],idx1[j]] = PhiGrad[i,j]
                
        idx0 = np.array([1,4,7,10],dtype = int)
        idx1 = np.array([3,4,5],dtype = int)
        for i in range(4):
            for j in range(3):
                Deta1[idx0[i],idx1[j]] = PhiGrad[i,j]
                
        idx0 = np.array([2,5,8,11],dtype = int)
        idx1 = np.array([6,7,8],dtype = int)
        for i in range(4):
            for j in range(3):
                Deta1[idx0[i],idx1[j]] = PhiGrad[i,j]
                
        Deta2 = np.zeros((12,9))
        idx0 = np.array([0,3,6,9],dtype = int)
        idx1 = np.array([0,3,6],dtype = int)
        for i in range(4):
            for j in range(3):
                Deta2[idx0[i],idx1[j]] = PhiGrad[i,j]
                
        idx0 = np.array([1,4,7,10],dtype = int)
        idx1 = np.array([1,4,7],dtype = int)
        for i in range(4):
            for j in range(3):
                Deta2[idx0[i],idx1[j]] = PhiGrad[i,j]
                
        idx0 = np.array([2,5,8,11],dtype = int)
        idx1 = np.array([2,5,8],dtype = int)
        for i in range(4):
            for j in range(3):
                Deta2[idx0[i],idx1[j]] = PhiGrad[i,j]
                
        eps = (Deta1 + Deta2)/2
        T = np.linalg.det(G) / 6
        devv = np.zeros((9))
        if np.linalg.norm(devv) - sigma_y / (2 * mu) > 0:
            C5 = C2 + C3 / np.linalg.norm(devv)
            C6 = C3 / np.linalg.norm(devv)**3 # 后面的没写
        else:
            C5 = 2 * mu
            C6 = np.zeros((12))
        
        treps = eps[:,0] + eps[:,4] + eps[:,8]
        deveps = np.zeros((eps.shape[0],9))
        for i in range(eps.shape[0]):
            for j in range(9):
                deveps[i,j] = eps[i,j] - treps[i] * devarray[j]
                

        dM0 = (C1 * np.dot(treps,np.transpose(treps))
               + C5 * np.dot(deveps,np.transpose(eps)) 
               - C6 * np.dot(devv,np.transpose(eps)))
        
        idx = np.zeros((12),dtype = int)
        idx[0] = elements[k,0] * 3
        idx[1] = elements[k,0] * 3 + 1
        idx[2] = elements[k,0] * 3 + 2
        idx[3] = elements[k,0] * 3
        idx[4] = elements[k,0] * 3 + 1
        idx[5] = elements[k,0] * 3 + 2
        idx[6] = elements[k,0] * 3
        idx[7] = elements[k,0] * 3 + 1
        idx[8] = elements[k,0] * 3 + 2
        idx[9] = elements[k,0] * 3
        idx[10] = elements[k,0] * 3 + 1
        idx[11] = elements[k,0] * 3 + 2
        for i in range(12):
            for j in range(12):
                idxi = idx[i]
                idxj = idx[j]
                DF[idxi,idxj] += dM0[i,j] * T
                
    