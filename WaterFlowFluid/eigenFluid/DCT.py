import numpy as np
# Use DCT to enforce the divergence free condition. 
# Note only Dirichlet BC are supported.
# 看不懂意思，结果也不对
Nx = 128
Ny = 128
u = np.zeros((Nx,Ny))
v = np.zeros((Nx,Ny))
u[6:12,8:10] = 100
uFreq = np.fft.fft2(u).real
vFreq = np.fft.fft2(v).real
for j in range(Ny):
    for i in range(Nx):
        K = np.array([i / Nx,j / Ny]) # waveNumber
        W = np.zeros((2))
        if i > 0:
            W[0] = uFreq[i-1,j]
        if j > 0:
            W[1] = vFreq[i,j-1]
        if i > 0 or j > 0:
            Ks = np.sqrt(K[0]**2 + K[1]**2)
            Kdw = W[0]*K[0] + W[1]*K[1]
            W = W - 1 / Ks * Kdw * K
        if i > 0:
            uFreq[i-1,j] = W[0]
        if j > 0:
            vFreq[i,j-1] = W[1]
u = np.fft.ifft2(uFreq).real
v = np.fft.ifft2(vFreq).real
div = np.zeros((Nx,Ny))
for j in range(1,Ny-1):
    for i in range(1,Nx-1):
        div[i,j] = (u[i+1,j] - u[i-1,j] + v[i,j+1] - v[i,j-1])/2