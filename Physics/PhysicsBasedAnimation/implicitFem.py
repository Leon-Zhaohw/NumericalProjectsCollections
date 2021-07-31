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
elementQ = np.zeros((num_elem,3,3))
stiffness = np.zeros((num_elem,16,9))
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

def vecdot(vec0,vec1):
    res = 0
    for i in range(len(vec0)):
        res += vec0[i] * vec1[i]
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
        
        for i in range(4):
            for j in range(4):
                ni = normal[ie,i,:]
                nj = normal[ie,j,:]
                nnt = np.zeros((3,3))
                for i0 in range(3):
                    for j0 in range(3):
                        nnt[i0,j0] = ni[i0] * nj[j0]
                mat = np.zeros((3,3))
                idx = int(i * 4 + j)
                mat = lam * nnt + mu * np.transpose(nnt)
                diag = mu * vecdot(ni,nj)
                mat[0,0] += diag
                mat[1,1] += diag
                mat[2,2] += diag
                mat /= (6 * det)
                
                for j0 in range(3):
                    for i0 in range(3):
                        idx0 = int(j0 * 3 + i0)
                        stiffness[ie,idx,idx0] = mat[i0,j0]
             
def dotMat(mat0,mat1):
    n = mat0.shape[0]
    m = mat0.shape[1]
    res = 0
    for i in range(n):
        for j in range(m):
           res += mat0[i,j] * mat1[i,j] 
    return res

def solve(b,x):
    # 又是共轭梯度，我都快背下来了
    tol = 1e-4
    ite = 0
    ite_max = 100
    r = np.zeros((num_particle,3))
    d = np.zeros((num_particle,3))
    q = np.zeros((num_particle,3))
    
    mvmul(x,q)
    
    r = b - q
    d = r.copy()
    
    sigma = dotMat(r,r)
    tol *= sigma
    
    while(ite < ite_max):
        ite += 1
        mvmul(d, q)
        alpha = sigma / dotMat(d,q)
        x += alpha * d
        r -= alpha * q
        sigma_old = sigma
        sigma = dotMat(r,r)
        if sigma < tol or alpha < 1e-20:
            break
        beta = sigma / sigma_old
        q = r.copy()
        q += beta * d
        d = q.copy()
        
        
    
def mvmul(x,y):
    global mass
    for p in range(num_particle):
        y[p,:] = mass[p] * x[p,:]
    scale = dt * dt + dt * damping
    for ie in range(num_elem):
        for i in range(4):
            for j in range(4):
                
                mat = np.zeros((3,3))
                idx = int(4 * i + j)
                for j0 in range(3):
                    for i0 in range(3):
                        idx0 = j0 * 3 + i0
                        mat[i0,j0] = stiffness[ie,idx,idx0]
                        
                term = np.dot(np.dot(elementQ[ie,:,:],mat),np.transpose(elementQ[ie,:,:]))
                y[elements[ie,i],:] += scale * scaledot(term[:,:],x[elements[ie,j],:])
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

        qs = np.dot(stress,Q)
        
        elementQ[ie,:,:] = Q[:,:]
        
        force[elements[ie,0],:] += scaledot(qs,normal[ie,0,:]) / 6
        force[elements[ie,1],:] += scaledot(qs,normal[ie,1,:]) / 6
        force[elements[ie,2],:] += scaledot(qs,normal[ie,2,:]) / 6
        force[elements[ie,3],:] += scaledot(qs,normal[ie,3,:]) / 6
        
    x = np.zeros((num_particle,3))
    b = np.zeros((num_particle,3))
    
    for p in range(num_particle):
        b[p,:] = mass[p] * vel[p,:] + dt * force[p,:]
        x[p,:] = vel[p,:]
        
    solve(b,x)
        
    for p in range(num_particle):
        vel[p,:] = x[p,:]
        pos[p,:] += dt * vel[p,:]
        
        if(pos[p,2] < 0):
            pos[p,2] = 0
            vel[p,2] = 0
        