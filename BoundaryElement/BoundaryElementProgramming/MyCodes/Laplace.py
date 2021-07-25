import numpy as np

'''
Fundamental solution for scalar wave equation Pressure

args:
   r - Distance between source and field point
   k - wave number
   dim - Cartesian dimension (2-D,3-D)
return:
   
'''
def UW(r,k,dim):
    if dim == 2:
        C0 = 0.25j
    elif dim == 3:
        C0 = k*r*1j
        uw = 1/(4*k*r)*np.exp(C0)
        
'''
Fundamental solution for Potential problems Temperature/Potential

args:
   r - Distance between source and field point
   k - wave number
   dim - Cartesian dimension (2-D,3-D)
return:
   
'''

def U(r,k,dim):
    res = 0
    if dim == 2:
        res = 1 / (2 * np.pi * k) * np.log(1 / r)
    elif dim == 3:
        res = 1 / (4 * np.pi * r * k)
    return res
        
'''
Fundamental solution for Potential problems Normal gradient

args:
   r - Distance between source and field point
   dxr - rx/r , ry/r , rz/r
   vnorm - Normal vector
   dim - Cartesian dimension (2-D,3-D)
return:
   
'''
def T(r,dxr,vnorm,dim):
    res = 0
    if dim == 2:
        res = - np.dot(vnorm,dxr) / (2 * np.pi * r)
    elif dim == 3:
        res = - np.dot(vnorm,dxr) / (4 * np.pi * r * r)
    return res
        
'''
Derivatives of Fundamental solution for Potential problems 
Temperature/Potential

args:
   r - Distance between source and field point
   dxr - Distances in Cartesian directions divided by r
   dim - Cartesian dimension (2-D,3-D)
return:
   
'''
def dU(r,dxr,dim):
    res = np.zeros((dim))
    if dim == 2:
        C = 1 / (2 * np.pi * r)
        res[0] = -C * dxr[0]
        res[1] = -C * dxr[1]
    elif dim == 3:
        C = 1 / (4 * np.pi * r**2)
        res[0] = C * dxr[0]
        res[1] = C * dxr[1]
        res[2] = C * dxr[2]
    return res
    
'''
derivatives of Fundamental solution for Potential problems 
Normal gradient

args:
   r - Distance between source and field point
   dxr - rx/r , ry/r , rz/r
   vnorm - Normal vector
   dim - Cartesian dimension (2-D,3-D)
return:
   
'''
def dT(r,dxr,vnorm,dim):
    res = np.zeros((dim))
    costh = np.dot(vnorm,dxr)
    if dim == 2:
        C = 1 / (2 * np.pi * r**2)
        res[0] = C * (2 * dxr[0] * costh - vnorm[0])
        res[1] = C * (2 * dxr[1] * costh - vnorm[1])
    elif dim == 3:
        C = 3 / (4 * np.pi * r**3)
        res[0] = C * costh * dxr[0]
        res[1] = C * costh * dxr[1]
        res[2] = C * costh * dxr[2]
        