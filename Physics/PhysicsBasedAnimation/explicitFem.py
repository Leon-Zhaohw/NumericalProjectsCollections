import numpy as np
import scipy.linalg
# 已完成，数据完全一致，稍微有点精度问题

dt = 0.00033333
lam = 1e4
mu = 1e7
damping = 0.01
density = 1e4

'''    读取文件  开始    '''

f = open("code/inputs/input.mesh","r")   #设置文件对象
line = f.readline()
num_elem = 0
num_particle = 1
num_tri = 0
while line:   
   line = f.readline()
   st = line.split(' ')
   if st[0] == 'p':
       num_particle += 1
   elif st[0] == 'e':
       num_elem += 1
   elif st[0] == 't':
       num_tri += 1
       
f = open("code/inputs/input.mesh","r")   #设置文件对象
pos = np.zeros((num_particle,3))
force = np.zeros((num_particle,3))
vel = np.zeros((num_particle,3))
mass = np.zeros((num_particle))
normal = np.zeros((num_elem,4,3))
elements = np.zeros((num_elem,4),dtype = int)
basis = np.zeros((num_elem,3,3))
tri = np.zeros((num_tri,3),dtype = int)
pcnt = 0
ecnt = 0
tcnt = 0
line = f.readline()
st = line.split(' ')
pos[pcnt,0] = st[1]
pos[pcnt,1] = st[2]
pos[pcnt,2] = st[3]
pcnt += 1
while line:   
   line = f.readline()
   st = line.split(' ')
   if st[0] == 'p':
       pos[pcnt,0] = st[1]
       pos[pcnt,1] = st[2]
       pos[pcnt,2] = st[3]
       pcnt += 1
   elif st[0] == 'e':
       elements[ecnt,0] = st[1]
       elements[ecnt,1] = st[2]
       elements[ecnt,2] = st[3]
       elements[ecnt,3] = st[4]
       ecnt += 1
   elif st[0] == 't':
       tri[tcnt,0] = st[1]
       tri[tcnt,1] = st[2]
       tri[tcnt,2] = st[3]
       tcnt += 1
       
f.close() #关闭文件
'''    读取文件  结束    '''

def cross(u,v):
    t0 = u[1]*v[2] - u[2]*v[1]
    t1 = u[2]*v[0] - u[0]*v[2]
    t2 = u[0]*v[1] - u[1]*v[0]
    return np.array([t0,t1,t2])

def scaledot(mat,vec):
    res = np.zeros((3))
    for i in range(3):
        for j in range(3):
            res[i] += mat[i,j] * vec[j]
    return res

def init():
    for ie in range(elements.shape[0]):
        x0 = pos[elements[ie,0],:]
        x1 = pos[elements[ie,1],:]
        x2 = pos[elements[ie,2],:]
        x3 = pos[elements[ie,3],:]
    
        X = np.zeros((3,3))
        X[0,:] = x1 - x0
        X[1,:] = x2 - x0
        X[2,:] = x3 - x0
        
        normal[ie,0,:] = cross(x2 - x1,x3 - x1)
        normal[ie,1,:] = cross(x3 - x0,x2 - x0)
        normal[ie,2,:] = cross(x1 - x0,x3 - x0)
        normal[ie,3,:] = cross(x2 - x0,x1 - x0)
        
        det = np.linalg.det(X)
        m = density * det / 24
        basis[ie,:,:] = np.linalg.inv(X)
        mass[elements[ie,0]] += m
        mass[elements[ie,1]] += m
        mass[elements[ie,2]] += m
        mass[elements[ie,3]] += m
        
init()
time = 0
timeFinal = 8


while(time < timeFinal):
    
    for p in range(num_particle):
        force[p,:] = 0
        force[p,2] = -9.8 * mass[p]
    
    for ie in range(elements.shape[0]):
        x0 = pos[elements[ie,0],:]
        x1 = pos[elements[ie,1],:]
        x2 = pos[elements[ie,2],:]
        x3 = pos[elements[ie,3],:]
    
        v0 = vel[elements[ie,0],:]
        v1 = vel[elements[ie,1],:]
        v2 = vel[elements[ie,2],:]
        v3 = vel[elements[ie,3],:]
        X = np.zeros((3,3))
        X[0,:] = x1 - x0
        X[1,:] = x2 - x0
        X[2,:] = x3 - x0
        
        # 注意Eigen库和numpy算乘法的方式不一样，不过怎么会不一样呢？
        F = np.dot(basis[ie],X)
        
        # scipy 的 svd 解法的精度问题
        for i in range(3):
            for j in range(3):
                if F[i,j] < 1e-10:
                    F[i,j] = 0
        Q = np.zeros((3,3)) # 旋转矩阵
        # 注意scipy 算出来的是V的转置，所以是V_transpose
        U,sigma,Vt = scipy.linalg.svd(F)
        
        Q = np.dot(Vt,U)
        Ftilde = np.dot(F,np.transpose(Q))
        
        iden = np.zeros((3,3))
        iden[0,0] = iden[1,1] = iden[2,2] = 1
        strain = (Ftilde + np.transpose(Ftilde)) / 2 - iden
        tr = strain[0,0] + strain[1,1] + strain[2,2]
        stress = lam * tr * iden + 2 * mu * strain
        
        Xdot = np.zeros((3,3))
        Xdot[0,:] = v1 - v0
        Xdot[1,:] = v2 - v0
        Xdot[2,:] = v3 - v0
        
        Fdot = np.dot(basis[ie,:,:],Xdot)
        Fdottilde = np.dot(Fdot,np.transpose(Q))
        
        strainrate = (Fdottilde + np.transpose(Fdottilde)) / 2 - iden
        tr1 = strainrate[0,0] + strainrate[1,1] + strainrate[2,2]
        stressrate = damping * (lam * tr1 * iden + 2 * mu * strainrate)
        
        qs = np.dot((stress + stressrate),Q)
        
        force[elements[ie,0],:] += scaledot(qs,normal[ie,0,:]) / 6
        force[elements[ie,1],:] += scaledot(qs,normal[ie,1,:]) / 6
        force[elements[ie,2],:] += scaledot(qs,normal[ie,2,:]) / 6
        force[elements[ie,3],:] += scaledot(qs,normal[ie,3,:]) / 6
        
    for p in range(num_particle):
        vel[p,:] += dt *(force[p,:] / mass[p])
        pos[p,:] += dt * vel[p,:]
        
        if(pos[p,2] < 0):
            pos[p,2] = 0
            vel[p,2] = 0
        