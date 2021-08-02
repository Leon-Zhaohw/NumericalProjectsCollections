import numpy as np

def axisAngle(axis,theta):
    R = np.zeros((3,3))
    norm = np.sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
    axis /= norm
    cosTheta = np.cos(theta)
    sinTheta = np.sin(theta)
    versine = 1 - cosTheta
    
    R[0,0] = axis[0] * axis[0] * versine + cosTheta
    R[0,1] = axis[0] * axis[1] * versine - axis[2] * sinTheta
    R[0,2] = axis[0] * axis[2] * versine + axis[1] * sinTheta
    
    R[1,0] = axis[1] * axis[0] * versine + axis[2] * sinTheta
    R[1,1] = axis[1] * axis[1] * versine + cosTheta
    R[1,2] = axis[1] * axis[2] * versine - axis[0] * sinTheta
    
    R[2,0] = axis[2] * axis[0] * versine - axis[1] * sinTheta
    R[2,1] = axis[2] * axis[1] * versine + axis[0] * sinTheta
    R[2,2] = axis[2] * axis[2] * versine + cosTheta
    
    return R

def mulMatMat(mat0,mat1):
    n = mat0.shape[0]
    res = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            
            for k in range(n):
                res[i,j] += mat0[i,k] * mat1[k,j]
    return res

def mulTranspose(vec0,vec1):
    n = len(vec0)
    res = np.zeros((n,n))
    for i in range(n):
        res[i,:] = vec0[i] * vec1[:]
    return res

# Rotation-variant SVD which pulls any reflections out of
# U and V and pushes them to Sigma
def svd_rv(mat):
    U,sigma,Vt = np.linalg.svd(mat)
    V = np.transpose(Vt)
    L = np.array([[1,0,0],[0,1,0],[0,0,1]],dtype = float)
    L[2,2] = np.linalg.det(np.dot(U,Vt))
    
    detU = np.linalg.det(U)
    detV = np.linalg.det(V)
    
    if detU < 0 and detV > 0:
        U = np.dot(U,L)
    elif detU > 0 and detV < 0:
        V = np.dot(V,L)
    elif detU < 0 and detV < 0:
        L = np.array([[1,0,0],[0,1,0],[0,0,-1]],dtype = float)
        U = np.dot(U,L)
        V = np.dot(V,L)
        L = np.array([[1,0,0],[0,1,0],[0,0,1]],dtype = float)
    si = np.array([[sigma[0],0,0],[0,sigma[1],0],[0,0,sigma[2]]])
    si = np.dot(si,L)
    return U,si,V
        
    

# Random Sigma
sigma = np.array([[2,0,0],[0,5,0],[0,0,11]],dtype = float)
a = np.array([1,2,3],dtype = float)
norma = np.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
a /= norma

# Random rotation and axisAngle

angleU = 95.1413
U = axisAngle(np.array([1,0.5,0.25]), angleU)

angleV = 212.1818
V = axisAngle(np.array([0.5,2.125,3.25]), angleV)

F = np.dot(np.dot(U,sigma), np.transpose(V))

PK1direct = 2 * np.dot(F, mulTranspose(a, a))

PK1numerical = np.zeros((3,3))

eps = 1e-2
previousdiff = 100000
converged = True

for epoch in range(5):
    F0 = F.copy()
    # a 3 x 1 matirx
    # P = a' * F' * F * a
    
    atFt = np.zeros((3))
    atFt = a[0]*F0[:,0] + a[1]*F0[:,1] + a[2]*F0[:,2]
    atFtF = np.zeros((3))
    atFtF = atFt[0]*F0[0,:] + atFt[1]*F0[1,:] + atFt[2]*F0[2,:]
    Psi0 = atFtF[0]*a[0] + atFtF[1]*a[1] + atFtF[2]*a[2]
    
    for y in range(3):
        for x in range(3):
            F1 = F.copy()
            F1[x,y] = F1[x,y] + eps
            U1,sigma1,V1 = svd_rv(F1)
            
            Ftemp = np.dot(np.dot(U1,sigma1), np.transpose(V1))
            atFt = np.zeros((3))
            atFt = a[0]*Ftemp[:,0] + a[1]*Ftemp[:,1] + a[2]*Ftemp[:,2]
            atFtF = np.zeros((3))
            atFtF = atFt[0]*Ftemp[0,:] + atFt[1]*Ftemp[1,:] + atFt[2]*Ftemp[2,:]
            Psi1 = atFtF[0]*a[0] + atFtF[1]*a[1] + atFtF[2]*a[2]
            
            PK1numerical[x,y] = (Psi1 - Psi0) / eps
            
    PK1diff = PK1direct - PK1numerical
    diff = np.linalg.norm(PK1diff) / 9
    eps *= 0.1
    if diff >= previousdiff / 2:
        converged = False
    previousdiff = diff
    print('eps : ',eps,'diff : ',diff,' pass ? :',converged)
