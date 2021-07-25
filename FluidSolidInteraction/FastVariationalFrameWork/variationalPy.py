import numpy as np
import random
# A Variational FrameWork , 
# 照着抄，抄了8个小时...不得不说把c++代码改编成python代码还是需要一点时间的

Nx = 20
Ny = 20
dx = 1/ Nx
u = np.zeros((Nx+1,Ny))
v = np.zeros((Nx,Ny+1))

u_weights = np.zeros((Nx+1,Ny))
v_weights = np.zeros((Nx,Ny+1))


valid = np.zeros((Nx+1,Ny+1))
old_valid = np.zeros((Nx+1,Ny+1))

cen0 = [0.5,0.5]
cen1 = [0.7,0.5]
cen2 = [0.3,0.35]
cen3 = [0.5,0.7]
rad0 = 0.4
rad1 = 0.1
rad2 = 0.1
rad3 = 0.1




def circlePhi(pos,cen,rad):
    return np.sqrt((pos[0]-cen[0])**2 + (pos[1]-cen[1])**2) - rad

def boundaryPhi(pos):
    phi0 = - circlePhi(pos,cen0,rad0)
    phi1 = circlePhi(pos,cen1,rad1)
    phi2 = circlePhi(pos,cen2,rad2)
    phi3 = circlePhi(pos,cen3,rad3)
    return min(min(phi0,phi1),min(phi2,phi3))

def set_boundary():
    nodal_solid_phi = np.zeros((Nx+1,Ny+1))
    for j in range(Ny+1):
        for i in range(Nx+1):
            nodal_solid_phi[i,j] = boundaryPhi([i*dx,j*dx])
    return nodal_solid_phi
            
def randhashf(seed):
    return np.sin(seed / 12.653)
            
def init_particle():
    num = 0
    for i in range(Nx*Nx):
        x = random.random()
        y = random.random()
        x = randhashf(i*2)
        y = randhashf(i*2+1)
        if boundaryPhi([x,y]) > 0:
            num += 1
    idx = 0
    ParticlePos = np.zeros((num,2))
    for i in range(Nx*Nx):
        x = randhashf(i*2)
        y = randhashf(i*2+1)
        if boundaryPhi([x,y]) > 0:
            ParticlePos[idx,:] = [x,y]
            idx += 1
    return ParticlePos
            
def interpolate_value(basef,val):
    
    basex = int(basef[0])
    basey = int(basef[1])
    basef[0] = basef[0] - basex
    basef[1] = basef[1] - basey
    
    if basex < 0:
        basex = 0
        basef[0] = 0
    elif basex > val.shape[0] - 2:
        basex = val.shape[0] - 2
        basef[0] = 1
    
    if basey < 0:
        basey = 0
        basef[1] = 0
    elif basey > val.shape[1] - 2:
        basey = val.shape[1] - 2
        basef[1] = 1
    
    term0 = val[basex,basey]*(1 - basef[0]) + val[basex+1,basey]*basef[0]
    term1 = val[basex,basey+1]*(1 - basef[0]) + val[basex+1,basey+1]*basef[0]
    return term0 * (1 - basef[1]) + term1 * basef[1]
    
def interpolate_gradient(basef,val):
    
    basex = int(basef[0])
    basey = int(basef[1])
    basef[0] = basef[0] - basex
    basef[1] = basef[1] - basey
    
    if basex < 0:
        basex = 0
        basef[0] = 0
    elif basex > val.shape[0] - 2:
        basex = val.shape[0] - 2
        basef[0] = 1
    
    if basey < 0:
        basey = 0
        basef[1] = 0
    elif basey > val.shape[1] - 2:
        basey = val.shape[1] - 2
        basef[1] = 1
        
    dx0 = val[basex+1,basey] - val[basex,basey]
    dx1 = val[basex+1,basey+1] - val[basex,basey+1]
    dy0 = val[basex,basey+1] - val[basex,basey]
    dy1 = val[basex+1,basey+1] - val[basex+1,basey]
    
    gradx = dx0 * (1 - basef[1]) + dx1 * basef[1]
    grady = dy0 * (1 - basef[0]) + dy1 * basef[0]
    
    return np.array([gradx,grady])

def get_velocity(pos):
    global u
    global v
    u_value = interpolate_value(pos / dx - [0,0.5], u)
    v_value = interpolate_value(pos / dx - [0.5,0], v)
    return np.array([u_value,v_value])
    
def trace_rk2(pos,dt0):
    vel = get_velocity(pos)
    vel2 = get_velocity(pos + dt0 * vel / 2)
    return pos + dt0 * vel2
        
def advect_particles():
    global ParticlePos
    global nodal_solid_phi
    ParticleNum = ParticlePos.shape[0]
    for p in range(ParticleNum):
        pos = trace_rk2(ParticlePos[p],0.005)
        phi_value = interpolate_value(pos/dx,nodal_solid_phi)
        # 如果 phi 小于零，也就是碰到固体了,那么把它弹回去
        if phi_value < 0:
            norm = interpolate_gradient(pos / dx, nodal_solid_phi)
            norm /= np.sqrt(norm[0]**2 + norm[1]**2)
            pos -= phi_value * norm
        ParticlePos[p] = pos
            
def advect():
    global u
    global v
    u_temp = u.copy()
    v_temp = v.copy()
    for j in range(Ny):
        for i in range(Nx+1):
            pos = np.array([i*dx,(j+0.5)*dx])
            pos = trace_rk2(pos,-dt)
            u_temp[i,j] = get_velocity(pos)[0]
    for j in range(Ny+1):
        for i in range(Nx):
            pos = np.array([(i+0.5)*dx,j*dx])
            pos = trace_rk2(pos,-dt)
            v_temp[i,j] = get_velocity(pos)[1]
    u = u_temp.copy()
    v = v_temp.copy()
            
        
def add_force():
    global v
    mid = Nx // 2
    v[mid:mid+2,mid] = 10

def fraction_inside(phi_left,phi_right):
    if phi_left < 0 and phi_right < 0:
        return 1
    if phi_left < 0 and phi_right >= 0:
        return phi_left / (phi_left - phi_right)
    if phi_left >= 0 and phi_right < 0:
        return phi_right / (phi_right - phi_left)
    return 0

def compute_weights():
    global u_weights
    global v_weights
    for j in range(Ny):
        for i in range(Nx+1):
            u_weights[i,j] = 1 - fraction_inside(nodal_solid_phi[i,j+1],nodal_solid_phi[i,j])
    for j in range(Ny+1):
        for i in range(Nx):
            v_weights[i,j] = 1 - fraction_inside(nodal_solid_phi[i+1,j],nodal_solid_phi[i,j])
            
            
def pcgsolve(lhs,rhs,toleranace_factor,marker):
    residual_out = max(abs(rhs))
    if residual_out == 0:
        return  True
    tol = toleranace_factor * residual_out
    
    invdiag,colstart,rowindex,value = form_preconditioner(lhs,rhs,marker)
    z = apply_preconditioner(rhs, invdiag, colstart, rowindex, value)
    rho = np.dot(z,rhs)
    if rho == 0:
        print("Matrix Error")
        return False
    s = z.copy()
    result = np.zeros((Nx*Ny))
    matrixValue = np.zeros((Nx*Ny*5))
    matrixColIndex = np.zeros((Nx*Ny*5),dtype = int)
    rowStart = np.zeros((Nx*Ny + 1),dtype = int)
    idx = 0
    for i in range(Nx*Ny):
        for j in range(Nx*Ny):
            if marker[i,j] == 1:
                matrixColIndex[idx] = j
                matrixValue[idx] = lhs[i,j]
                idx += 1
            rowStart[i + 1] = idx
    ite = 0
    ite_max = 100
    
    # 终于到主循环了
    while(ite < ite_max):
        ite += 1
        
        for i in range(Nx*Ny):
            z[i] = 0
            for j in range(rowStart[i],rowStart[i+1]):
                idx = matrixColIndex[j]
                z[i] += matrixValue[j] * s[idx]
        
        alpha = rho / np.dot(s,z)
        result += alpha * s
        rhs -= alpha * z
        residual_out = max(abs(rhs))
        if residual_out <= tol:
            break
        z = apply_preconditioner(rhs, invdiag, colstart, rowindex, value)
        rho_new = np.dot(z,rhs)
        beta = rho_new / rho
        z += beta * s
        temp = z.copy()
        z = s.copy()
        s = temp.copy()
        rho = rho_new
    return result
        
        
        
    
def form_preconditioner(lhs,rhs,marker):
    
    colstart = np.zeros((Nx*Ny),dtype = int)
    summ = 0
    adiag = np.zeros((Nx*Ny))
    invdiag = np.zeros((Nx*Ny))
    value = np.zeros((Nx*Ny*2))
    rowindex = np.zeros((Nx*Ny*2),dtype = int)
    matrixindex = np.ones((Nx*Ny,5)) * (-1)
    idx = 0
    for i in range(Nx*Ny):
        adiag[i] = lhs[i,i]
        invdiag[i] = lhs[i,i]
        idx0 = 0
        for j in range(Nx*Ny):
            if marker[i,j] == 1:
                matrixindex[i,idx0] = j
                idx0 += 1
                if j > i:
                    value[idx] = lhs[i,j]
                    rowindex[idx] = j
                    idx += 1
        colstart[i] = idx
    min_diagonal_ratio = 0.25
    for k in range(Nx*Ny):
        if adiag[k] == 0:
            continue
        if invdiag[k] < min_diagonal_ratio * adiag[k]:
            invdiag[k] = 1 / np.sqrt(adiag[k])
        else:
            invdiag[k] = 1 / np.sqrt(invdiag[k])
        for p in range(colstart[k-1],colstart[k]):
            value[p] *= invdiag[k]
        for p in range(colstart[k-1],colstart[k]):
            j = rowindex[p]
            multiplier = value[p]
            missing = 0
            a = colstart[k-1]
            while(a < colstart[k] and rowindex[a] < j):
                idx0 = 0
                while(matrixindex[j,idx0] != -1):
                    if matrixindex[j,idx0] < rowindex[a]:
                        idx0 += 1
                    elif matrixindex[j,idx0] == rowindex[a]:
                        break
                    else:
                        missing += value[a]
                        break
                a += 1
            
            if(a < colstart[k] and rowindex[a] == j):
                invdiag[j] -= multiplier * value[a]
                
            a += 1
            b = colstart[j-1]
            while(a < colstart[k] and b < colstart[j]):
                if rowindex[b] < rowindex[a]:
                    b += 1
                elif rowindex[b] == rowindex[a]:
                    value[b] -= multiplier * value[a]
                    a += 1
                    b += 1
                else:
                    missing += value[a]
                    a += 1
            while(a < colstart[k]):
                missing += value[a]
                a += 1
            
            modifiction = 0.97
            invdiag[j] -= modifiction * multiplier * missing
                
    return invdiag,colstart,rowindex,value

def apply_preconditioner(rhs,invdiag,colstart,rowindex,value):
    result = rhs.copy()
    for i in range(Nx*Ny):
        result[i] *= invdiag[i]
        if i == 0:
            continue
        for j in range(colstart[i-1],colstart[i]):
            result[int(rowindex[j])] -= value[j] * result[i]
            
    i = Nx*Ny - 1
    while(i > 0):
        for j in range(colstart[i-1],colstart[i]):
            result[i] -= value[j] * result[int(rowindex[j])]
        result[i] *= invdiag[i]
        i -= 1
    return result


        
def solve_pressure():
    global u_weights
    global v_weights
    global u
    global v
    rhs = np.zeros((Nx*Ny))
    pressure = np.zeros((Nx*Ny))
    lhs = np.zeros((Nx*Ny,Nx*Ny))
    
    marker = np.zeros((Nx*Ny,Nx*Ny))
    
    for j in range(1,Ny-1):
        for i in range(1,Nx-1):
            idx = j * Nx + i
            
            term = u_weights[i+1,j] * dt / dx / dx
            lhs[idx,idx] += term
            lhs[idx,idx+1] -= term
            rhs[idx] -= u_weights[i+1,j]*u[i+1,j] / dx
            
            term = u_weights[i,j] * dt / dx / dx
            lhs[idx,idx] += term
            lhs[idx,idx-1] -= term
            rhs[idx] += u_weights[i,j]*u[i,j] / dx
            
            term = v_weights[i,j+1] * dt / dx / dx
            lhs[idx,idx] += term
            lhs[idx,idx+Nx] -= term
            rhs[idx] -= v_weights[i,j+1]*v[i,j+1] / dx
            
            term = v_weights[i,j] * dt / dx / dx
            lhs[idx,idx] += term
            lhs[idx,idx-Nx] -= term
            rhs[idx] += v_weights[i,j]*v[i,j] / dx
            
            marker[idx,idx] = 1
            marker[idx,idx+1] = 1
            marker[idx,idx-1] = 1
            marker[idx,idx+Nx] = 1
            marker[idx,idx-Nx] = 1
            
    # 效率低，建议换成preconditioned conjugate gradient
    pressure = pcgsolve(lhs, rhs,1e-5,marker)
    
    for j in range(Ny):
        for i in range(Nx+1):
            idx = j * Nx + i
            if u_weights[i,j] > 0:
                u[i,j] -= dt * (pressure[idx] - pressure[idx-1]) / dx
            else:
                u[i,j] = 0
    
    for j in range(Ny+1):
        for i in range(Nx):
            idx = j * Nx + i
            if v_weights[i,j] > 0:
                v[i,j] -= dt * (pressure[idx] - pressure[idx-Nx]) / dx
            else:
                v[i,j] = 0    
    
    return pressure
            
def extrapolate(val,weights):
    nx = weights.shape[0]
    ny = weights.shape[1]
    valid = np.zeros((nx,ny),dtype = bool)
    valid1 = np.zeros((nx,ny),dtype = bool)
    for j in range(ny):
        for i in range(nx):
            valid1[i,j] = weights[i,j] > 0
    for ite in range(5):
        valid = valid1.copy()
        for j in range(1,ny-1):
            for i in range(1,nx-1):
                summ = 0
                count = 0
                if(valid[i,j] == False):
                    if(valid[i+1,j] == True):
                        summ += val[i+1,j]
                        count += 1
                    if(valid[i-1,j] == True):
                        summ += val[i-1,j]
                        count += 1
                    if(valid[i,j+1] == True):
                        summ += val[i,j+1]
                        count += 1
                    if(valid[i,j-1] == True):
                        summ += val[i,j-1]
                        count += 1
                    if count > 0:
                        val[i,j] = summ / count
                        valid1[i,j] = True
            
def constrain_veloicty():
    global u
    global v
    u_temp = u.copy()
    v_temp = v.copy()
    normal = np.zeros((2,Nx+1,Ny))
    for j in range(Ny):
        for i in range(Nx+1):
            if u_weights[i,j] == 0:
                pos = np.array([i*dx,(j+0.5)*dx])
                vel = get_velocity(pos)
                
                norm = interpolate_gradient(pos/dx,nodal_solid_phi)
                norm /= np.sqrt(norm[0]**2 + norm[1]**2)
                
                prep = vel[0]*norm[0] + vel[1]*norm[1]
                normal[:,i,j] = vel
                vel -= prep * norm
                
                u_temp[i,j] = vel[0]
    for j in range(Ny+1):
        for i in range(Nx):
            if v_weights[i,j] == 0:
                pos = np.array([(i+0.5)*dx,j*dx])
                vel = get_velocity(pos)
                norm = interpolate_gradient(pos/dx,nodal_solid_phi)
                norm /= np.sqrt(norm[0]**2 + norm[1]**2)
                prep = vel[0]*norm[0] + vel[1]*norm[1]
                vel -= prep * norm
                v_temp[i,j] = vel[1]    
    u = u_temp.copy()
    v = v_temp.copy()
    return normal
            
nodal_solid_phi = set_boundary()
ParticlePos = init_particle()

time = 0
timeFianl = 1
dt = 0.005
while(time < timeFianl):
    time += 0.5
    advect_particles()
    advect()
    add_force()
    compute_weights()
    solve_pressure()
    extrapolate(u,u_weights)
    extrapolate(v,v_weights)
    normal = constrain_veloicty()