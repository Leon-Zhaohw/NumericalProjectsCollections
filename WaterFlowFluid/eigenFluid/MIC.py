import numpy as np
# Builds the Modified Incomplete Cholesky Preconditioner per Bridson
# 成功，果然是总散度的关系
heat = False # 是否解的是热力学方程
dt = 0.1
heatDiffusion = 1
dx = 0.1
Nx = 16
Ny = 16
Adiag = np.zeros((Nx,Ny))
Aplusi = np.zeros((Nx,Ny))
Aplusj = np.zeros((Nx,Ny))
precon = np.zeros((Nx,Ny))
obstacle = np.zeros((Nx,Ny),dtype = bool)
def buildA():
    global Adiag
    global Aplusi
    global Aplusj
    Aoff = 1
    if heat == True:
        Aoff = dt * heatDiffusion / dx / dx
    for j in range(1,Ny-1):
        for i in range(1,Nx-1):
            if obstacle[i,j] == False:
                Adiag[i,j] = 0
                if heat == True:
                    Adiag[i,j] = 1
                    
                # 记住，Aplusi最后一列是没用的
                # Aplusi[i,j]的位置，实际上是在[i+0.5,j]
                if obstacle[i+1,j] == False:
                    Adiag[i,j] += Aoff
                    Aplusi[i,j] = -Aoff
                if obstacle[i-1,j] == False:
                    Adiag[i,j] += Aoff
                    Aplusi[i-1,j] = -Aoff
                if obstacle[i,j+1] == False:
                    Adiag[i,j] += Aoff
                    Aplusj[i,j] = -Aoff
                if obstacle[i,j-1] == False:
                    Adiag[i,j] += Aoff
                    Aplusj[i,j-1] = -Aoff
                    
def buildMICPreconitioner():
    global precon
    global Adiag
    global Aplusi
    global Aplusj
    tau = 0.97 # tuning constant
    sig = 0.25 # safty constant
    for j in range(1,Ny-1):
        for i in range(1,Nx-1):
            pmini = precon[i-1,j]**2
            pminj = precon[i,j-1]**2
            e0 = Aplusi[i-1,j]**2*pmini + Aplusj[i,j-1]**2*pminj
            e1 = tau * (Aplusi[i-1,j] * Aplusj[i-1,j] * pmini
                + Aplusj[i,j-1] * Aplusi[i,j-1] * pminj)
            e = Adiag[i,j] - e0 - e1
            # if we fall under safey
            if e < sig * Adiag[i,j]:
                e = Adiag[i,j] 
            if abs(Adiag[i,j]) < 1e-6: # 如果什么都没有
                precon[i,j] = 0
            else:
                precon[i,j] = 1 / np.sqrt(e)
        
def zeroBoundary(val):
    global obstacle
    for j in range(Ny):
        for i in range(Nx):
            if obstacle[i,j] == True:
                val[i,j] = 0
    
def multiA(xf,bf):
    global Aplusi
    global Aplusj
    global Adiag
    for j in range(1,Ny-1):
        for i in range(1,Nx-1):
            bf[i,j] = Adiag[i,j] * xf[i,j] + (
                xf[i-1,j] * Aplusi[i-1,j] + 
                xf[i+1,j] * Aplusi[i,j] + 
                xf[i,j-1] * Aplusj[i,j-1] + 
                xf[i,j+1] * Aplusj[i,j])
            
def applyMICPreconditioner(zf,r,q):
    global precon
    global Aplusi
    global Aplusj
    # L q = r
    for j in range(1,Ny-1):
        for i in range(1,Nx-1):
            t = r[i,j] - Aplusi[i-1,j] * precon[i-1,j] * q[i-1,j]
            t -= Aplusj[i,j-1] * precon[i,j-1] * q[i,j-1]
            q[i,j] = t * precon[i,j]
            
    # L^T z = q
    for j in range(Ny-2,0,-1):
        for i in range(Nx-2,0,-1):
            pre = precon[i,j]
            t = q[i,j] - Aplusi[i,j] * pre * zf[i+1,j]
            t -= Aplusj[i,j] * pre * zf[i,j+1]
            zf[i,j] = t * pre
    
def scaleDot(a,b):
    res = 0
    for j in range(Ny):
        for i in range(Nx):
            res += a[i,j]*b[i,j]
    return res

obstacle[0,:] = True
obstacle[Nx-1,:] = True
obstacle[:,0] = True
obstacle[:,Ny-1] = True
buildA()
buildMICPreconitioner()
residual = np.zeros((Nx,Ny))
b = np.zeros((Nx,Ny)) # 散度
b[5,5] = 1
b[5,6] = -1
x = np.zeros((Nx,Ny))
z = np.zeros((Nx,Ny))
q = np.zeros((Nx,Ny))
multiA(x,residual)
for j in range(1,Ny-1):
    for i in range(1,Nx-1):
        residual[i,j] = b[i,j] - residual[i,j]
zeroBoundary(residual)
applyMICPreconditioner(z, residual, q)
direction = z.copy()
deltaNew = scaleDot(residual,z)
eps_min = 1e-10
eps = 100
ite = 0
ite_max = 100
while(ite < ite_max and eps > eps_min):
    multiA(direction,z)
    alpha = scaleDot(direction,z)
    if abs(alpha) > 0:
        alpha = deltaNew / alpha
    else:
        print("alpha broken down.solver failed")
        break
    x = x + alpha * direction
    residual = residual - alpha * z
    applyMICPreconditioner(z, residual, q)
    zeroBoundary(q)
    
    # residual[0,:] = 0
    # residual[Nx-1,:] = 0
    # residual[:,0] = 0
    # residual[:,Ny-1] = 0
    
    eps = np.max(abs(residual))
    deltaOld = deltaNew
    deltaNew = scaleDot(residual,z)
    beta = deltaNew / deltaOld
    direction = z + beta * direction
    ite += 1