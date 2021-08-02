import numpy as np

class Cube:
    restGaussPoints = np.zeros((8,3))

    mu = 1
    lam = 10
    Nx = 5
    dx = 0.1
    restVolume = 0
    
    def __init__(self):
        invSqrt3 = 1 / np.sqrt(3)
        self.restGaussPoints[0,:] = np.array([-1,-1,-1]) * invSqrt3
        self.restGaussPoints[1,:] = np.array([1,-1,-1]) * invSqrt3
        self.restGaussPoints[2,:] = np.array([-1,1,-1]) * invSqrt3
        self.restGaussPoints[3,:] = np.array([1,1,-1]) * invSqrt3
        self.restGaussPoints[4,:] = np.array([-1,-1,1]) * invSqrt3
        self.restGaussPoints[5,:] = np.array([1,-1,1]) * invSqrt3
        self.restGaussPoints[6,:] = np.array([-1,1,1]) * invSqrt3
        self.restGaussPoints[7,:] = np.array([1,1,1]) * invSqrt3
        self.vertexIndex = np.zeros((8),dtype = int)
        self.vertexPos = np.zeros((8,3))
        self.vertexPosRest = np.zeros((8,3))
        self.elementForce = np.zeros((24))
        self.localIKs = np.zeros((24,24))
        
    def generate(self,idx,Nx,dx):
        self.Nx = Nx
        self.dx = dx
        self.restVolume = dx * dx * dx
        self.vertexIndex[0] = idx
        self.vertexIndex[1] = idx + 1
        self.vertexIndex[2] = idx + Nx
        self.vertexIndex[3] = idx + 1 + Nx
        self.vertexIndex[4] = idx + Nx * Nx
        self.vertexIndex[5] = idx + Nx * Nx + 1
        self.vertexIndex[6] = idx + Nx * Nx + Nx
        self.vertexIndex[7] = idx + Nx * Nx + 1 + Nx
        
    def mulMatMat(self,mat0,mat1):
        n0 = mat0.shape[0]
        m0 = mat0.shape[1]
        n1 = mat1.shape[0]
        m1 = mat1.shape[1]
        assert m0 == n1
        res = np.zeros((n0,m1))
        for i in range(n0):
            for j in range(m1):
                for k in range(m0):
                    res[i,j] += mat0[i,k] * mat1[k,j]
        return res
    
    def ComputeH(self,idx):
        H = np.zeros((8,3))
        x1m = 1 - self.restGaussPoints[idx,0]
        x1p = 1 + self.restGaussPoints[idx,0]
        x2m = 1 - self.restGaussPoints[idx,1]
        x2p = 1 + self.restGaussPoints[idx,1]
        x3m = 1 - self.restGaussPoints[idx,2]
        x3p = 1 + self.restGaussPoints[idx,2]
            
        H[0,0] = - x2m * x3m
        H[0,1] = - x1m * x3m
        H[0,2] = - x1m * x2m
        H[1,0] = x2m * x3m
        H[1,1] = - x1p * x3m
        H[1,2] = - x1p * x2m
        
        H[2,0] = - x2p * x3m
        H[2,1] = x1m * x3m
        H[2,2] = - x1m * x2p
        H[3,0] = x2p * x3m
        H[3,1] = x1p * x3m
        H[3,2] = - x1p * x2p
        
        H[4,0] = - x2m * x3p
        H[4,1] = - x1m * x3p
        H[4,2] = x1m * x2m
        H[5,0] = x2m * x3p
        H[5,1] = - x1p * x3p
        H[5,2] = x1p * x2m
        
        H[6,0] = - x2p * x3p
        H[6,1] = x1m * x3p
        H[6,2] = x1m * x2p
        H[7,0] = x2p * x3p
        H[7,1] = x1p * x3p
        H[7,2] = x1p * x2p
        
        H /= 8
        return H
        
    # Dm 是现在顶点的位置，
    def ComputeValue(self):
        self.Bmg = np.zeros((8,3,8))
        self.HDmHinv = np.zeros((8,8,3))
        self.pFpuMatrix = np.zeros((8,9,24))
        
        
        for i in range(8):
            H = self.ComputeH(i)
            DmHg = np.dot(np.transpose(self.vertexPos),H)
            DmHginv = np.linalg.inv(DmHg)
            self.Bmg[i,:,:] = np.dot(np.transpose(DmHginv),np.transpose(H)) * self.restVolume
            self.HDmHinv[i,:,:] = np.dot(H,DmHginv)
            temp = self.HDmHinv[i,:,:]
            
        for i in range(8):
           for j in range(24):
               
               H = self.ComputeH(i)
               DmH = np.dot(np.transpose(self.vertexPos),H)
               A = np.dot(H,np.linalg.inv(DmH)) # 8x3 matrix
               delta = np.zeros((3,8))
               delta[int(j%3),int(j/3)] = 1
               pFpu = np.dot(delta,A) # 3x3 matrix
               self.pFpuMatrix[i,0:3,j] = pFpu[:,0]
               self.pFpuMatrix[i,3:6,j] = pFpu[:,1]
               self.pFpuMatrix[i,6:9,j] = pFpu[:,2]
               
    def cross(self,u,v):
        w = np.zeros((3))
        w[0] = u[1]*v[2] - u[2]*v[1]
        w[1] = u[2]*v[0] - u[0]*v[2]
        w[2] = u[0]*v[1] - u[1]*v[0]
        return w

    def PartialJpartialF(self,F):
        res = np.zeros((3,3))
        res[0,:] = self.cross(F[1,:], F[2,:])
        res[1,:] = self.cross(F[2,:], F[0,:])
        res[2,:] = self.cross(F[0,:], F[1,:])
        return res
    
    def Flatten(self,mat):
        n = mat.shape[0]
        m = mat.shape[1]
        res = np.zeros((n * m))
        for j in range(m):
            for i in range(n):
                idx = j * n + i
                res[idx] = mat[i,j]
        return res
                
               
    def ComputeForce(self):
        self.elementForce[:] = 0
        forcevector = np.zeros((24))
        for gpidx in range(8):
            # F 是 3x3 的力矩阵
            sv = self.vertexPos
            sd = self.HDmHinv[gpidx,:,:]
            F = np.dot(np.transpose(self.vertexPos),self.HDmHinv[gpidx,:,:])
            pJpF = self.PartialJpartialF(F)
            ratio = 0.1
            Jminus = np.linalg.det(F) - 1 - ratio
            Pk1 =  self.mu * F + self.lam * Jminus * pJpF
            G = - np.dot(Pk1,self.Bmg[gpidx,:,:]) # 3 x 8
            forcevector += self.Flatten(G)
        self.elementForce[:] = forcevector.copy()           
            
    def RotationVariantSVD(self,F):
        U,sigma,Vt = np.linalg.svd(F)
        V = np.transpose(Vt)
        if np.linalg.det(U) < 0:
            U[:,0] = - U[:,0]
            sigma[0] = - sigma[0]
        if np.linalg.det(V) < 0:
            V[:,0] = - V[:,0]
            sigma[0] = - sigma[0]
        return U,sigma,V
        
    def ProjectHessian(self,F,U,V,S):
        eigenvalues = np.zeros((9))
        eigenvectors = np.zeros((9,9))
        J = np.linalg.det(F)
        lam = 10
        mu = 1
        
        # Build Twist And Flip Eigenvectors
        eigenvalues[0:3] = S[:]
        eigenvalues[3:6] = -S[:]
        evScale = lam * (J - 1) - mu
        eigenvalues[0:6] = eigenvalues[0:6]*evScale + mu
        
        scale = 1 / np.sqrt(2)
        sV = np.transpose(V * scale)
        U = np.transpose(U)
        
        A = np.array([sV[0,2]*U[0,1] , sV[0,2]*U[1,1] , sV[0,2]*U[2,1] ,
                      sV[1,2]*U[0,1] ,sV[1,2]*U[1,1] ,sV[1,2]*U[2,1] ,
                      sV[2,2]*U[0,1] ,sV[2,2]*U[1,1] ,sV[2,2]*U[2,1]])
        
        # 硬编码，真有你的哦
        B = np.array([sV[0,1]*U[0,2] , sV[0,1]*U[1,2] , sV[0,1]*U[2,2] ,
                       sV[1,1]*U[0,2] , sV[1,1]*U[1,2] , sV[1,1]*U[2,2] ,
                        sV[2,1]*U[0,2], sV[2,1]*U[1,2], sV[2,1]*U[2,2]])
        
        C = np.array([sV[0,2]*U[0,0] , sV[0,2]*U[1,0] , sV[0,2]*U[2,0] ,
                      sV[1,2]*U[0,0] , sV[1,2]*U[1,0] , sV[1,2]*U[2,0] ,
                      sV[2,2]*U[0,0], sV[2,2]*U[1,0], sV[2,2]*U[2,0]])
        
        D = np.array([sV[0,0]*U[0,2] , sV[0,0]*U[1,2] , sV[0,0]*U[2,2] ,
                      sV[1,0]*U[0,2] , sV[1,0]*U[1,2] , sV[1,0]*U[2,2] ,
                      sV[2,0]*U[0,2], sV[2,0]*U[1,2], sV[2,0]*U[2,2]])
        
        E = np.array([sV[0,1]*U[0,0] , sV[0,1]*U[1,0] , sV[0,1]*U[2,0] ,
                      sV[1,1]*U[0,0] , sV[1,1]*U[1,0] , sV[1,1]*U[2,0] ,
                      sV[2,1]*U[0,0] , sV[2,1]*U[1,0], sV[2,1]*U[2,0]])
        
        F = np.array([sV[0,0]*U[0,1] , sV[0,0]*U[1,1] , sV[0,0]*U[2,1] ,
                      sV[1,0]*U[0,1] , sV[1,0]*U[1,1] , sV[1,0]*U[2,1] ,
                      sV[2,0]*U[0,1], sV[2,0]*U[1,1], sV[2,0]*U[2,1]])
        
        eigenvectors[:,0] = B - A
        eigenvectors[:,1] = D - C
        eigenvectors[:,2] = F - E
        
        eigenvectors[:,3] = A + B
        eigenvectors[:,4] = C + D
        eigenvectors[:,5] = E + F
        
        Amat = np.zeros((3,3))
        s0s0 = S[0]*S[0]
        s1s1 = S[1]*S[1]
        s2s2 = S[2]*S[2]
        
        Amat[0,0] = mu + lam * s1s1 * s2s2
        Amat[1,1] = mu + lam * s0s0 * s2s2
        Amat[2,2] = mu + lam * s0s0 * s1s1
        evScale = lam * (2 * J - 1) - mu
        Amat[1,0] = Amat[0,1] = evScale * S[2]
        Amat[2,0] = Amat[0,2] = evScale * S[1]
        Amat[1,2] = Amat[2,1] = evScale * S[0]
        value,vector = np.linalg.eig(Amat)
        # numpy 库 和 Eigen 库 各种意义上是反的 ...
        eigenvalues[6] = value[2]
        eigenvalues[7] = value[1]
        eigenvalues[8] = value[0]
        U = np.transpose(U)
        diag = np.zeros((3,3))
        diag[0,0] = vector[0,2]
        diag[1,1] = vector[1,2]
        diag[2,2] = vector[2,2]
        eigenvectors[:,6] = self.Flatten(np.dot(np.dot(U,diag),np.transpose(V)))
        diag[0,0] = vector[0,1]
        diag[1,1] = vector[1,1]
        diag[2,2] = vector[2,1]
        eigenvectors[:,7] = self.Flatten(np.dot(np.dot(U,diag),np.transpose(V)))
        diag[0,0] = vector[0,0]
        diag[1,1] = vector[1,0]
        diag[2,2] = vector[2,0]
        eigenvectors[:,8] = self.Flatten(np.dot(np.dot(U,diag),np.transpose(V)))
        
        diag = np.zeros((9,9))
        for i0 in range(9):
            if eigenvalues[i0] < 0:
                eigenvalues[i0] = 0
            diag[i0,i0] = eigenvalues[i0]

        return np.dot(np.dot(eigenvectors,diag),np.transpose(eigenvectors))                

        
        
        
    def ComputeForceJacobian(self):
        K = np.zeros((24,24))
        for i in range(8):
            F = np.dot(np.transpose(self.vertexPos),self.HDmHinv[i,:,:])
            U,sigma,V = self.RotationVariantSVD(F)
            
            matrix1 = np.transpose(self.pFpuMatrix[i,:,:]) # 24 x 9
            matrix2 = self.ProjectHessian(F,U.copy(),V,sigma) # 9 x 9 
            matrix3 = self.pFpuMatrix[i,:,:] # 9 x 24
            res = np.dot(np.dot(matrix1,matrix2),matrix3) * (-self.restVolume)
            K += res
        self.localIKs = K.copy()
            
               
               
               
               
               
               
        
        
        