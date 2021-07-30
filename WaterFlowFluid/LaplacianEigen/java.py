import numpy as np
Nx = 32
Ny = 32
N = 16

Velx = np.zeros((Nx+1,Ny+1))
Vely = np.zeros((Nx+1,Ny+1))
density = np.zeros((Nx+1,Ny+1))

coef = np.zeros((N))
forces_dw = np.zeros((N))

N_sqrt = int(np.sqrt(N))
lookup_table = np.zeros((N,2))
rlookup_table = np.zeros((N_sqrt+1,N_sqrt+1))
vel_basis_x = np.zeros((N,Nx+1,Ny+1))
vel_basis_y = np.zeros((N,Nx+1,Ny+1))
Ck = np.zeros((N,N,N)) # 应该用稀疏矩阵来存。稀疏矩阵是个很简单的数据结构，两个变长数组，一个存索引，一个存值

amp = 1

def coefdensity(a1,b1,a2,b2,c,tt):
    if tt == 0:
        # SS x SS
        if c == 0:
            return - (a1*b2 - a2*b1) / 4
        if c == 1:
            return (a1*b2 + a2*b1) / 4
        if c == 2:
            return -(a1*b2 + a2*b1) / 4
        if c == 3:
            return (a1*b2 - a2*b1) / 4
    elif tt == 1:
        # SC x SS
        if c == 0:
            return - (a1*b2 - a2*b1) / 4
        if c == 1:
            return - (a1*b2 + a2*b1) / 4
        if c == 2:
            return (a1*b2 + a2*b1) / 4
        if c == 3:
            return (a1*b2 - a2*b1) / 4
    elif tt == 2:
        # CS x SS
        if c == 0:
            return - (a1*b2 - a2*b1) / 4
        if c == 1:
            return - (a1*b2 + a2*b1) / 4
        if c == 2:
            return (a1*b2 + a2*b1) / 4
        if c == 3:
            return (a1*b2 - a2*b1) / 4
    elif tt == 3:
        # CS x SS
        if c == 0:
            return - (a1*b2 - a2*b1) / 4
        if c == 1:
            return - (a1*b2 + a2*b1) / 4
        if c == 2:
            return (a1*b2 + a2*b1) / 4
        if c == 3:
            return (a1*b2 - a2*b1) / 4
        
def cur_energy():
    energy = 0
    global eigs_inv
    global coef
    for i in range(N):
        energy += eigs_inv[i]*coef[i]*coef[i]
    return energy

# fill_lookup_table

rlookup_table[:,:] = -1
idx = 0
for k1 in range(N_sqrt+1):
    for k2 in range(N_sqrt+1):
        if k1 > N_sqrt or k1 < 1 or k2 > N_sqrt or k2 < 1:
            continue
        lookup_table[idx,0] = k1
        lookup_table[idx,1] = k2
        rlookup_table[k1,k2] = idx
        idx += 1

# precompute_basis_fields

for k in range(N):
    k1 = lookup_table[k,0]
    k2 = lookup_table[k,1]
    
    xfact = 1
    if k1 != 0:
        xfact = -1 / (k1**2 + k2**2)
    yfact = 1
    if k2 != 0:
        yfact = -1 / (k1**2 + k2**2)
    
    deltax = np.pi / Nx
    deltay = np.pi / Ny
    
    for j in range(Ny+1):
        for i in range(Nx+1):
            vel_basis_x[k,i,j] = amp * xfact * -k2 * np.sin(k1*i) * np.cos(k2*(j+deltay/2))
            vel_basis_y[k,i,j] = amp * yfact * k1 * np.sin(k1*(i+deltax/2)) * np.cos(k2*j)
            
# precompute_dynamics
eigs = np.zeros((N))
eigs_inv = np.zeros((N))
eigs_inv_root = np.zeros((N))

for k in range(N):
    k1 = lookup_table[k,0]
    k2 = lookup_table[k,1]
    eigs[k] = (k1*k1 + k2*k2)
    eigs_inv = 1 / (k1*k1 + k2*k2)
    eigs_inv_root = 1 / np.sqrt(k1*k1 + k2*k2)
    
for i in range(N):
    a1 = lookup_table[i,0]
    a2 = lookup_table[i,1]
    lambda_a = -(a1*a1 + a2*a2)
    inv_lambda_a = - 1 / (a1*a1 + a2*a2)
    for j in range(N):
        b1 = lookup_table[j,0]
        b2 = lookup_table[j,1]
        lambda_b = -(b1*b1 + b2*b2)
        inv_lambda_b = - 1 / (b1*b1 + b2*b2)
        
        k1 = rlookup_table[a1,a2]
        k2 = rlookup_table[b1,b2]
        
        antiparis = np.zeros((4,2),dtype = int)
        antiparis[0,0] = a1 - b1
        antiparis[0,1] = a2 - b2
        antiparis[1,0] = a1 - b1
        antiparis[1,1] = a2 + b2
        antiparis[2,0] = a1 + b1
        antiparis[2,1] = a2 - b2
        antiparis[3,0] = a1 + b1
        antiparis[3,1] = a2 + b2
        
        for c in range(4):
            i = antiparis[c,0]
            j = antiparis[c,1]
            idx = rlookup_table[i,j]
            
            if idx != -1:
                coe = coefdensity(a1, b1, a2, b2, c, 0) * inv_lambda_b
                Ck[idx,k1,k2] = coe
                Ck[idx,k2,k1] = coe * - lambda_b / lambda_a
                
def MatMulVec(mat,vec):
    n = mat.shape[0]
    res = np.zeros((n))
    for i in range(n):
        for j in range(n):
            res[i] += mat[i,j]*vec[j]
    return res

def scaleDot(vec0,vec1):
    n = len(vec0)
    res = 0
    for i in range(n):
        res += vec0[i]*vec1[i]
        
def set_energy(desired_e):
    global coef
    cur_e = cur_energy()
    fact = np.sqrt(desired_e) / np.sqrt(cur_e)
    coef *= fact
                
def expand_basis():
    global Velx
    global Vely
    Velx[:,:] = 0
    Vely[:,:] = 0
    for k in range(N):
        for j in range(Ny+1):
            for i in range(Nx+1):
                Velx[i,j] += coef[k] * vel_basis_x[k,i,j]
                Vely[i,j] += coef[k] * vel_basis_y[k,i,j]
    
time = 0
timeFinal = 100
dt = 0.1
visc = 1
while(time < timeFinal):
    dw = np.zeros((N))
    pre_e = cur_energy()
    dwt = np.zeros((4,N))
    qn = np.zeros((4,N))
    for k in range(N):
        dwt[0,k] = scaleDot(qn[0,:],MatMulVec(Ck[k,:,:], qn[0,:]))
        qn[1,k] = qn[0,k] + dt * dwt[0,k] / 2
    for k in range(N):
        dwt[1,k] = scaleDot(qn[1,:],MatMulVec(Ck[k,:,:], qn[1,:]))
        qn[2,k] = qn[0,k] + dt * dwt[1,k] / 2
    for k in range(N):
        dwt[2,k] = scaleDot(qn[2,:],MatMulVec(Ck[k,:,:], qn[2,:]))
        qn[3,k] = qn[0,k] + dt * dwt[2,k]
    for k in range(N):
        dwt[3,k] = scaleDot(qn[3,:],MatMulVec(Ck[k,:,:], qn[3,:]))
        dw[k] = (dwt[0,k] + 2 * dwt[1,k] + 2 * dwt[2,k] + dwt[3,k]) / 6
        
    coef += dw * dt
    set_energy(pre_e)
    coef = coef * np.exp(-1 * eigs * dt * visc) + forces_dw
    forces_dw[:] = 0
    expand_basis()
        