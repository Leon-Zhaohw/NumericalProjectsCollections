import numpy as np

Nx = 1
Ny = 1
Nz = 1

vortx = np.zeros((Nx,Ny,Nz))
vorty = np.zeros((Nx,Ny,Nz))
vortz = np.zeros((Nx,Ny,Nz))
gradrho = np.zeros((Nx,Ny,Nz))
vortpos = np.zeros((Nx,Ny,Nz))
'''
作用：多项式插值
参数：r : 距离，浮点数
返回: 权重
'''
def poly_shape(r):
    
    rr = abs(r)
    # 第一种实现方法
    if rr < 1:
        return 1 - rr
    return 0
    # 第二种实现方法
    # if rr < 1:
    #     return (1 - rr)**4
    # return 0
    # 第三种实现方法
    # res = 0
    # if rr < 2:
    #     res = ((5 - 2*rr) - np.sqrt(-4*rr*rr + 12*rr - 7))
    # if rr < 1:
    #     res = ((3 - 2*rr) + np.sqrt(1 + 4*rr - 4*rr*rr))
    # return res / 8

def poly_shape2(r):
    rr = abs(r)
    if rr < 1:
        return (1 - rr)**4
    return 0

'''
作用：计算网格位置与粒子位置之间的权重
参数：grid : 网格位置，三维向量
     particle: 粒子位置，三维向量
返回: 权重
'''
def inperpolate_weights(grid,particle):
    dirx = poly_shape(grid[0] - particle[0])
    diry = poly_shape(grid[1] - particle[1])
    dirz = poly_shape(grid[2] - particle[2])
    return np.array([dirx,diry,dirz])

def inperpolate_weights2(grid,particle):
    dirx = poly_shape(grid[0] - particle[0])
    diry = poly_shape(grid[1] - particle[1])
    dirz = poly_shape(grid[2] - particle[2])
    return poly_shape2(dirx**2 + diry**2 + dirz**2)

def cross_uxv(u,v):
    term0 = u[1]*v[2] - u[2]*v[1] # UyVz - UzVy
    term1 = u[2]*v[0] - u[0]*v[2] # UzVx - UxVz
    term2 = u[0]*v[1] - u[1]*v[0] # UxVy - UyVx
    return np.array([term0,term1,term2])

dx = 1 / Nx
dy = 1 / Ny
dz = 1 / Nz
ParticleNum = 1
ParticlePos = np.zeros((ParticleNum,3))
dt = 0
'''
作用：计算粒子斜压
参数：ppos: 涡粒子位置
返回: 权重
'''
def computeBaroclinic(pidx):
    px = ParticlePos[pidx,0] / dx
    py = ParticlePos[pidx,1] / dy
    pz = ParticlePos[pidx,2] / dz
    bx = int(px)
    by = int(py)
    bz = int(pz)
    summ = np.zeros((3))
    for k in range(bz-2,bz+3):
        for j in range(by-2,by+3):
            for i in range(bx-2,bx+3):
                pijk = np.array([i+0.5,j+0.5,k+0.5])
                w = inperpolate_weights(pijk*dx,ParticlePos[pidx,0])
                summ += gradrho * w
    baroclinic = cross_uxv(summ,[0,10,0])
    vortx[pidx] += dt * 0.1 * baroclinic[0]
    vorty[pidx] += dt * 0.1 * baroclinic[1]
    vortz[pidx] += dt * 0.1 * baroclinic[2]
      
def computeBoundaryVortex():
    Vel = np.zeros((3))
    Normal = np.zeros((3))
    normal_component = Vel[0]*Normal[0] + Vel[1]*Normal[1] + Vel[2]*Normal[2]
    tangent = np.zeros((3))
    tangent[0] = Vel[0] - normal_component * Normal[0]
    tangent[1] = Vel[1] - normal_component * Normal[1]
    tangent[2] = Vel[2] - normal_component * Normal[2]
    normal = np.array([])
    intensity = cross_uxv(Normal,tangent)
    dt_c = 0
    area = 0
    vort_x = dt_c * area * intensity[0]
    vort_y = dt_c * area * intensity[1]
    vort_z = dt_c * area * intensity[2]
    