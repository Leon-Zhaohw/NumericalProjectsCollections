# 将表面的三角形存储起来
surfaceTriangle = np.zeros((resolution*resolution*2*6,3))
cnt = 0
# -y
for k in range(resolution):
    j = 0
    for i in range(resolution):
        c0 = int(k*Ny*Nx + j*Nx + i)
        c1 = c0 + 1
        c4 = c0 + Nx * Ny
        c5 = c1 + Nx * Ny
        surfaceTriangle[cnt,:] = np.array([c0,c1,c4])
        cnt += 1
        surfaceTriangle[cnt,:] = np.array([c4,c1,c5])
        cnt += 1
# +y
for k in range(resolution):
    j = resolution - 1
    for i in range(resolution):
        c0 = int(k*Ny*Nx + j*Nx + i)
        c1 = c0 + 1
        c2 = c0 + Nx
        c3 = c1 + Nx
        c6 = c2 + Nx * Ny
        c7 = c3 + Nx * Ny
        surfaceTriangle[cnt,:] = np.array([c2,c6,c3])
        cnt += 1
        surfaceTriangle[cnt,:] = np.array([c3,c6,c7])
        cnt += 1
# -z   
k = 0
for j in range(resolution):
    for i in range(resolution):
        c0 = int(k*Ny*Nx + j*Nx + i)
        c1 = c0 + 1
        c2 = c0 + Nx
        c3 = c1 + Nx
        surfaceTriangle[cnt,:] = np.array([c0,c2,c1])
        cnt += 1
        surfaceTriangle[cnt,:] = np.array([c1,c2,c3])
        cnt += 1
# +z
k = resolution - 1
for j in range(resolution):
    for i in range(resolution):
        c0 = int(k*Ny*Nx + j*Nx + i)
        c1 = c0 + 1
        c2 = c0 + Nx
        c3 = c1 + Nx
        c4 = c0 + Nx * Ny
        c5 = c1 + Nx * Ny
        c6 = c2 + Nx * Ny
        c7 = c3 + Nx * Ny
        surfaceTriangle[cnt,:] = np.array([c4,c5,c6])
        cnt += 1
        surfaceTriangle[cnt,:] = np.array([c5,c7,c6])
        cnt += 1
# -x
i = 0
for k in range(resolution):
    for j in range(resolution):
        c0 = int(k*Ny*Nx + j*Nx + i)
        c2 = c0 + Nx
        c4 = c0 + Nx * Ny
        c6 = c2 + Nx * Ny
        surfaceTriangle[cnt,:] = np.array([c0,c6,c2])
        cnt += 1
        surfaceTriangle[cnt,:] = np.array([c0,c4,c6])
        cnt += 1
# +x
i = resolution - 1
for k in range(resolution):
    for j in range(resolution):
        c0 = int(k*Ny*Nx + j*Nx + i)
        c1 = c0 + 1
        c3 = c1 + Nx
        c5 = c1 + Nx * Ny
        c7 = c3 + Nx * Ny
        surfaceTriangle[cnt,:] = np.array([c1,c3,c7])
        cnt += 1
        surfaceTriangle[cnt,:] = np.array([c1,c7,c5])
        cnt += 1