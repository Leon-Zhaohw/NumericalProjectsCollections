import numpy as np

Nx = 4
Ny = 6
Nz = 4

dt = 0.1
buoyancy = 0.1
heatDiffusion = 1e-3
vorticityEps = 2

maxRes = max(max(Nx,Ny),Nz)
dx = 1 / maxRes
scaling = 64 / maxRes
if scaling < 1:
    scaling = 1
vorticityEps /= scaling

divergence = np.zeros((Nx,Ny,Nz))
vorticity = np.zeros((Nx,Ny,Nz))
xVel = np.zeros((Nx,Ny,Nz))
yVel = np.zeros((Nx,Ny,Nz))
zVel = np.zeros((Nx,Ny,Nz))
xForce = np.zeros((Nx,Ny,Nz))
yForce = np.zeros((Nx,Ny,Nz))
zForce = np.zeros((Nx,Ny,Nz))
xVort = np.zeros((Nx,Ny,Nz))
yVort = np.zeros((Nx,Ny,Nz))
zVort = np.zeros((Nx,Ny,Nz))
density = np.zeros((Nx,Ny,Nz))
pressure = np.zeros((Nx,Ny,Nz))
heat = np.zeros((Nx,Ny,Nz))
obstacle = np.zeros((Nx,Ny,Nz),dtype = int)

# 前后左右有围墙，而上下没有
obstacle[0,:,:] = 1
obstacle[Nx-1,:,:] = 1
obstacle[:,:,0] = 1
obstacle[:,:,Nx-1] = 1

sphereCenter = [0,375,0.5,0.375]
sphereRadixVels = 0.1
for k in range(Nz):
    for j in range(Ny):
        for i in range(Nz):
            diffx = i * dx - sphereCenter[0]
            diffy = j * dx - sphereCenter[1]
            diffz = k * dx - sphereCenter[2]
            mag = diffx**2 + diffy**2 + diffz**2
            if mag*mag < sphereRadixVels*sphereRadixVels:
                obstacle[i,j,k] = 1
                
def addSmoke(field):
    heightMin = 0.05
    heightMax = 0.1
    yMin = int(heightMin * Ny)
    yMax = int(heightMax * Ny)
    for k in range(Nz):
        for j in range(yMin,yMax+1):
            for i in range(Nx):
                xLength = (i - 0.4 * Nx) * dx
                zLength = (k - 0.5 * Nz) * dx
                rad = np.sqrt(xLength**2 + zLength**2)
                if rad < 0.075 * dx * Nx:
                    field[i,j,k] = 1
    
def ZeroBoundary(field):
    field[0,:,:] = 0
    field[Nx-1,:,:] = 0
    field[:,0,:] = 0
    field[:,Ny-1,:] = 0
    field[:,:,0] = 0
    field[:,:,Nz-1] = 0
                
def WipeBoundary():
    ZeroBoundary(xVel)
    ZeroBoundary(yVel)
    ZeroBoundary(zVel)
    ZeroBoundary(density)
    
def addVorticity():
    for k in range(1,Nz-1):
        for j in range(1,Ny-1):
            for i in range(1,Nx-1):
                
                right = i + (obstacle[i+1,j,k] == 0)
                left = i - (obstacle[i-1,j,k] == 0)
                up = j + (obstacle[i,j+1,k] == 0)
                down = j - (obstacle[i,j-1,k] == 0)
                front = k + (obstacle[i,j,k+1] == 0)
                back = k - (obstacle[i,j,k-1] == 0)
                
                diffx = max(right - left,1)
                diffy = max(up - down,1)
                diffz = max(front - back,1)
                
                xVort[i,j,k] = (zVel[i,up,k]-zVel[i,down,k])/diffy + (-yVel[i,j,front]+yVel[i,j,back])/diffz
                yVort[i,j,k] = (xVel[i,j,front]-xVel[i,j,back])/diffz + (-zVel[right,j,k]+zVel[left,j,k])/diffx
                zVort[i,j,k] = (yVel[right,j,k]-yVel[left,j,k])/diffx + (-xVel[i,up,k]+xVel[i,down,k])/diffy   
                vorticity[i,j,k] = np.sqrt(xVort[i,j,k]**2 + yVort[i,j,k]**2 + zVort[i,j,k]**2)
             
    for k in range(1,Nz-1):
        for j in range(1,Ny-1):
            for i in range(1,Nx-1):
                if obstacle[i,j,k] == 1:
                    continue
                
                right = i + (obstacle[i+1,j,k] == 0)
                left = i - (obstacle[i-1,j,k] == 0)
                up = j + (obstacle[i,j+1,k] == 0)
                down = j - (obstacle[i,j-1,k] == 0)
                front = k + (obstacle[i,j,k+1] == 0)
                back = k - (obstacle[i,j,k-1] == 0)
                
                diffx = max(right - left,1)
                diffy = max(up - down,1)
                diffz = max(front - back,1)
                
                n0 = (vorticity[right,j,k] - vorticity[left,j,k])/diffx
                n1 = (vorticity[i,up,k] - vorticity[i,down,k])/diffy
                n2 = (vorticity[i,j,front] - vorticity[i,j,back])/diffz
                
                mag = np.sqrt(n0**2 + n1**2 + n2**2)
                if mag > 0:
                    n0 /= mag
                    n1 /= mag
                    n2 /= mag
                    xForce[i,j,k] = (n1*zVort[i,j,k] - n2*yVort[i,j,k])
                    yForce[i,j,k] = -(n0*zVort[i,j,k] - n2*xVort[i,j,k])
                    zForce[i,j,k] = (n0*yVort[i,j,k] - n1*xVort[i,j,k])
                else:
                    xForce[i,j,k] = yForce[i,j,k] = zForce[i,j,k]
                
def addBuoyancy():
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                yForce[i,j,k] += buoyancy * heat[i,j,k]
    
def applyForce():
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                xVel[i,j,k] += dt * xForce[i,j,k]
                yVel[i,j,k] += dt * yForce[i,j,k]
                zVel[i,j,k] += dt * zForce[i,j,k]
    
def pressureAverage():
    for k in range(1,Nz-1):
        for j in range(1,Ny-1):
            for i in range(1,Nx-1):
                if obstacle[i,j,k] == 1:
                    # 只求流体能够分一点压力给自己
                    pressure[i,j,k] = 0
                    xVel[i,j,k] = yVel[i,j,k] = zVel[i,j,k]
                    cnt = 0
                    if obstacle[i+1,j,k] == 0:
                        pressure[i,j,k] += pressure[i+1,j,k]
                        cnt += 1
                    if obstacle[i-1,j,k] == 0:
                        pressure[i,j,k] += pressure[i-1,j,k]
                        cnt += 1
                    if obstacle[i,j+1,k] == 0:
                        pressure[i,j,k] += pressure[i,j+1,k]
                        cnt += 1
                    if obstacle[i,j-1,k] == 0:
                        pressure[i,j,k] += pressure[i,j-1,k]
                        cnt += 1
                    if obstacle[i,j,k+1] == 0:
                        pressure[i,j,k] += pressure[i,j,k+1]
                        cnt += 1
                    if obstacle[i,j,k-1] == 0:
                        pressure[i,j,k] += pressure[i,j,k-1]
                        cnt += 1
                    if cnt > 0:
                        pressure[i,j,k] /= cnt
    
def computeDiv():
    for k in range(1,Nz-1):
        for j in range(1,Ny-1):
            for i in range(1,Nx-1):
                
                right = i + (obstacle[i+1,j,k] == 0)
                left = i - (obstacle[i-1,j,k] == 0)
                up = j + (obstacle[i,j+1,k] == 0)
                down = j - (obstacle[i,j-1,k] == 0)
                front = k + (obstacle[i,j,k+1] == 0)
                back = k - (obstacle[i,j,k-1] == 0)
                
                diffx = max(right - left,1)
                diffy = max(up - down,1)
                diffz = max(front - back,1)
                
                term0 = (xVel[right,j,k] - xVel[left,j,k]) / diffx
                term1 = (yVel[i,up,k] - yVel[i,down,k]) / diffy
                term2 = (zVel[i,j,front] - zVel[i,j,back]) / diffz
                
                divergence[i,j,k] = - (term0 + term1 + term2) * dx
                
# 给大家表演一下 如何不存储那个超大的二维矩阵A
def ConjudateGradient(result,source,Ascale):
    residual = np.zeros((Nx,Ny,Nz))
    for k in range(1,Nz-1):
        for j in range(1,Ny-1):
            for i in range(1,Nx-1):
                A = Ascale
                totalp = 0
                if obstacle[i,j,k] == 0:
                    if obstacle[i+1,j,k] == 0:
                        A += 1
                        totalp += pressure[i+1,j,k]
                    if obstacle[i-1,j,k] == 0:
                        totalp += pressure[i-1,j,k]
                        A += 1
                    if obstacle[i,j+1,k] == 0:
                        totalp += pressure[i,j+1,k]
                        A += 1
                    if obstacle[i,j-1,k] == 0:
                        totalp += pressure[i,j-1,k]
                        A += 1
                    if obstacle[i,j,k+1] == 0:
                        totalp += pressure[i,j,k+1]
                        A += 1
                    if obstacle[i,j,k-1] == 0:
                        totalp += pressure[i,j,k-1]
                        A += 1
                    residual[i,j,k] = source[i,j,k] - (A * pressure[i,j,k] - totalp)
                    
    direction = residual.copy()
    deltaNew = np.dot(residual,residual)
    ite = 0
    ite_max = 100
    eps = 100
    eps_min = 0.01
    q = np.zeros((Nx,Ny,Nz))
    while(ite < ite_max and eps < eps_min):
        for k in range(1,Nz-1):
            for j in range(1,Ny-1):
                for i in range(1,Nx-1):
                    A = 0
                    totalp = 0
                    if obstacle[i,j,k] == 0:
                        if obstacle[i+1,j,k] == 0:
                            A += 1
                            totalp += pressure[i+1,j,k]
                        if obstacle[i-1,j,k] == 0:
                            totalp += pressure[i-1,j,k]
                            A += 1
                        if obstacle[i,j+1,k] == 0:
                            totalp += pressure[i,j+1,k]
                            A += 1
                        if obstacle[i,j-1,k] == 0:
                            totalp += pressure[i,j-1,k]
                            A += 1
                        if obstacle[i,j,k+1] == 0:
                            totalp += pressure[i,j,k+1]
                            A += 1
                        if obstacle[i,j,k-1] == 0:
                            totalp += pressure[i,j,k-1]
                            A += 1
                        q[i,j,k] = A * pressure[i,j,k] - totalp
                    else:
                        q[i,j,k] = 0
                        
        alpha = np.dot(direction,q)
        if alpha > 0:
            alpha = deltaNew / alpha
        
        result += alpha * direction
        residual -= alpha * q
        
        deltaOld = deltaNew
        deltaNew = np.dot(residual,residual)
        beta = deltaNew / deltaOld
        
        direction = residual + beta * direction
        eps = max(residual)
        
def project():
    for k in range(1,Nz-1):
        for j in range(1,Ny-1):
            for i in range(1,Nx-1):
                xVel[i,j,k] += (pressure[i-1,j,k] - pressure[i+1,j,k]) / 2 / dx
                yVel[i,j,k] += (pressure[i,j-1,k] - pressure[i,j+1,k]) / 2 / dx
                zVel[i,j,k] += (pressure[i,j,k-1] - pressure[i,j,k+1]) / 2 / dx
                
def solveHeat():
    heat[0,:,:] = heat[1,:,:]
    heat[Nx-1,:,:] = heat[Nx-2,:,:]
    heat[:,0,:] = heat[:,1,:]
    heat[:,Ny-1,:] = heat[:,Ny-2,:]
    heat[:,:,0] = heat[:,:,1]
    heat[:,:,Nz-1] = heat[:,:,Nz-2]
    
    heatConst = dt * heatDiffusion / dx / dx
    heatOld = heat.copy()
    # 好奇怪，热源竟是我自己？
    ConjudateGradient(heat, heatOld, 1 / heatConst)
    
    for k in range(1,Nz-1):
        for j in range(1,Ny-1):
            for i in range(1,Nx-1):
                if obstacle[i,j,k] == 1:
                    heat[i,j,k] = 0
                    
def minD(a0,a1,a2):
    t0 = (a1 - a0) * dx
    t1 = (a2 - a1) * dx
    t2 = (a2 - a0) * dx / 2
    return min(t0,min(t1,t2))
                   
# 此函数细节未必正确 
def computeLU(A):
    m = A.shape[0]
    n = A.shape[1]
    piv = np.zeros((m))
    for i in range(m):
        piv[i] = i
    pivsign = 1
    LUcolj = np.zeros((m))
    LU = A.copy()
    for j in range(n):
        # Make a copy of the j-th column to localize references.
        for i in range(m):
            LUcolj[i] = LU[i][j]
        # Apply previous transformations.
        for i in range(m):
            LUrowi = LU[i,:]
            # Most of the time is spent in the following dot product.
            kmax = min(i,j)
            s = 0
            for k in range(kmax):
                s += LUrowi[k] * LUcolj[k]
            LUcolj[i] -= s
            LUrowi[j] = LUcolj[i]
        # Find pivot and exchange if necessary.
        p = j
        for i in range(j+1,m):
            if abs(LUcolj[i]) > abs(LUcolj[p]):
                p = i
        if p != j:
            k = 0
            for k in range(n):
                t = LU[p,k]
                LU[p,k] = LU[j,k]
                LU[j,k] = t
            k = piv[p]
            piv[p] = piv[j]
            piv[j] = k
            pivsign = - pivsign
            
        # compute multipliers
        if j < m and LU[j,j] != 0 :
            for i in range(j+1,m):
                LU[i,j] /= LU[j,j]
    return LU
            
def isNonSingular(LU):
    n = LU.shape[0]
    for i in range(n):
        if LU[i,i] == 0:
            return False
    return True
    
def permuteCopy(B,piv):
    n = len(piv)
    x = np.zeros((n))
    for i in range(n):
        x[i] = B[piv[i]]
    return x

# solve L * U * x = b[piv,:]
def solveLUX(B):
    n = len(B)
    LU = np.zeros((n,n))
    piv = np.zeros((n))
    x = permuteCopy(B, piv)
    
    # solve L * Y = B
    for k in range(n):
        for i in range(k+1,n):
            for j in range(n):
                x[i,j] -= x[k,j]*LU[i,k]
                
    # solve U * X = Y
    for k in range(n-1,-1,-1):
        for j in range(n):
            x[k,j] /= LU[k,k]
        for i in range(k):
            for j in range(n):
                x[i,j] -= x[k,j] * LU[i,k]
    
        
# eig 3x1 matrix A 3x3 matrix
def EigenValue(eig,A):
    n = A.shape[0]
    V = np.zeros((n,n))
    d = np.zeros((n))
    e = np.zeros((n))
    
    isSymmetric = True
    for j in range(n):
        for i in range(n):
            if A[i,j] != A[j,i]:
                isSymmetric = False
                break
    
    if isSymmetric == True:
        for i in range(n):
            for j in range(n):
                V[i,j] = A[i,j]
    

def computeEigenValues():
    maxeig = -1
    mineig = 10
    for k in range(1,Nz-1):
        for j in range(1,Ny-1):
            for i in range(1,Nx-1):
                minxVelDx = minD(xVel[i-1,j,k],xVel[i,j,k],xVel[i+1,j,k])
                minyVelDx = minD(yVel[i-1,j,k],yVel[i,j,k],yVel[i+1,j,k])
                minzVelDx = minD(zVel[i-1,j,k],zVel[i,j,k],zVel[i+1,j,k])
                minxVelDy = minD(xVel[i,j-1,k],xVel[i,j,k],xVel[i,j+1,k])
                minyVelDy = minD(yVel[i,j-1,k],yVel[i,j,k],yVel[i,j+1,k])
                minzVelDy = minD(zVel[i,j-1,k],zVel[i,j,k],zVel[i,j+1,k])
                minxVelDz = minD(xVel[i,j,k-1],xVel[i,j,k],xVel[i,j,k+1])
                minyVelDz = minD(yVel[i,j,k-1],yVel[i,j,k],yVel[i,j,k+1])
                minzVelDz = minD(zVel[i,j,k-1],zVel[i,j,k],zVel[i,j,k+1])
                jac = np.array([[minxVelDx,minyVelDx,minzVelDx],
                                [minxVelDy,minyVelDy,minzVelDy],
                                [minxVelDz,minyVelDz,minzVelDz]])
                
                LU = computeLU(jac)
                if isNonSingular(LU):
                    eigenvalues = [1,1,1]
                    
                    
                    
                
time = 0
timeFinal = 1
while(time < timeFinal):
    time = 100
    addSmoke(density)
    addSmoke(heat)
    # 还有个zVelayVeleTxVelrbxVellence不知道在搞什么
    WipeBoundary()
    addVorticity()
    addBuoyancy()
    applyForce()
    pressureAverage()
    computeDiv()
    ConjudateGradient(pressure, divergence, 0)
    project()
    solveHeat()
    
    
    
    
    
