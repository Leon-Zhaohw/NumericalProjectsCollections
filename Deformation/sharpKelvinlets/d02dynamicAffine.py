import numpy as np

axis_element = 50
num_element = axis_element * axis_element
axis_nodes = axis_element + 1
num_nodes = axis_nodes * axis_nodes
L = 10
dx = L / axis_element

nodePos = np.zeros((num_nodes,3)) # 假的，实际上是二维的
for i in range(num_nodes):
    x = int(i % axis_nodes)
    y = int(i / axis_nodes)
    nodePos[i,0] = x * dx - L / 2
    nodePos[i,1] = y * dx - L / 2
    
eps = 1.5
nu = 0.45
mu = 5
_beta  = np.sqrt(mu)
_alpha = _beta * np.sqrt(1 + 1 / (1 - 2 * nu))
if nu == 0.5:
    _alpha = np.inf
force = np.zeros((3,3))
action = 0 # 0 twist, 1 uniform scale , 2 pinch
if action == 0:
    axisAngle = np.pi / 4
    force[0,1] = -axisAngle
    force[1,0] = axisAngle
elif action == 1:
    force[0,0] = force[1,1] = force[2,2] = 1
elif action == 2:
    force[0,0] = 1
    force[1,1] = -1

e4 = eps * eps * eps * eps
q = np.zeros((3))
q[0] = (force[2,1] - force[1,2])/2
q[1] = (force[0,2] - force[2,0])/2
q[2] = (force[1,0] - force[0,1])/2

s = (force[0,0] + force[1,1] + force[2,2]) / 3
iden = np.array([[1,0,0],[0,1,0],[0,0,1]])
pinch = (force + np.transpose(force)) / 2 - s * iden

tFactor = - 10 * _beta * e4
q *= tFactor / 5
if _alpha == np.inf:
    s = 0
    pinch *= tFactor / 2
else:
    sFactor = -10 * _alpha * e4
    s *= sFactor / 5
    pinch /= (2/tFactor + 3/sFactor)




def AssembleSkewSymMatrix(A):
    res = np.zeros((3,3))
    res[0,1] = -A[2]
    res[0,2] = A[1]
    res[1,0] = A[2]
    res[1,2] = -A[0]
    res[2,0] = -A[1]
    res[2,1] = A[0]
    return res

force = s * iden
force += AssembleSkewSymMatrix(q)
force += pinch

def EvalDen(r,a):
    return 16 * np.pi * a * pow(r,3)

def EvalU(r,t,a):
    if a == np.inf:
        return 0
    if t == np.inf:
        return 0
    if r < 1e-4:
        return EvalU0(t,a)
    
    fP = EvalW(r,r + a * t)
    fN = EvalW(r,r - a * t)
    return (fP - fN) / EvalDen(r,a)

def EvalGradU(r,t,a):
    if a == np.inf:
        return 0
    if t == np.inf:
        return 0
    if r < 1e-4:
        return 0
    
    fP = EvalW(r,r + a * t)
    fN = EvalW(r,r - a * t)
    gP = EvalGradW(r,r + a * t)
    gN = EvalGradW(r,r - a * t)
    return ((gP - gN) - 3*(fP - fN) / r) / EvalDen(r, a)

def EvalHessU(r,t,a):
    if a == np.inf:
        return 0
    if t == np.inf:
        return 0
    if r < 1e-4:
        return 0
    r2 = r * r
    fP = EvalW(r,r + a * t)
    fN = EvalW(r,r - a * t)
    gP = EvalGradW(r,r + a * t)
    gN = EvalGradW(r,r - a * t)
    hP = EvalHessW(r,r + a * t)
    hN = EvalHessW(r,r - a * t)
    return ((hP - hN) - 6*(gP - gN)/r + 12 * (fP - fN) / r2) / EvalDen(r, a)

def EvalU0(t,a):
    a2t2 =pow(a*t,2)
    e2 = pow(eps,2)
    num = 5 * t * e2 * e2
    den = 8 * np.pi * np.sqrt(pow(a2t2+e2,7))
    return num / den

def EvalW(r,R):
    e2 = pow(eps,2)
    R2 = R*R
    Re2 = R2 + e2
    Re = np.sqrt(Re2)
    
    num = e2 + 2 * R2 - r * R * (3 - R2 / Re2)
    den = Re
    return num / den

def EvalGradW(r,R):
    e2 = pow(eps,2)
    R2 = R*R
    Re2 = R2 + e2
    Re = np.sqrt(Re2)
    
    num = -3 * e2 * e2 * r
    den = Re2 * Re2 * Re
    return num / den

def EvalHessW(r,R):
    e2 = pow(eps,2)
    R2 = R*R
    Re2 = R2 + e2
    Re = np.sqrt(Re2)
    
    num = -3 * e2 * e2 * (Re2 - 5 * r * R)
    den = Re2 * Re2 * Re2 * Re
    return num / den
    
def EvalA(cacheua,cacheub,r):
    return cacheua[0] + 2 * cacheub[0] + r * cacheub[1]

def EvalB(cacheua,cacheub,r):
    if r < 1e-4:
        return 0
    return (cacheua[1] - cacheub[1]) / r

def EvalGradA(cacheua,cacheub,r):
    return cacheua[1] + 3 * cacheub[1] + r * cacheub[2]

def EvalGradB(cacheua,cacheub,r):
    B = EvalB(cacheua,cacheub,r)
    return (cacheua[2] - cacheub[2] - B) / r

def MulMatVec(mat,vec):
    n = mat.shape[0]
    res = np.zeros((3))
    for i in range(3):
        for j in range(3):
            res[i] += mat[i,j] * vec[j]
    return  res

def EvalDisp(x,t):
    if t <= 0:
        return np.zeros((3))
    r = np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
    if r < 1e-4:
        return np.zeros((3))
    
    CacheUa = np.zeros((3))
    CacheUb = np.zeros((3))
    
    CacheUa[1] = EvalGradU(r,t,_alpha)
    CacheUa[2] = EvalHessU(r,t,_alpha)
    CacheUb[1] = EvalGradU(r,t,_beta)
    CacheUb[2] = EvalHessU(r,t,_beta)
    
    Fx = MulMatVec(force, x)
    xdotFx = x[0]*Fx[0] + x[1]*Fx[1] + x[2]*Fx[2] 
    result = np.zeros((3))
    result += EvalGradA(CacheUa, CacheUb, r) * Fx / r
    result += EvalGradB(CacheUa, CacheUb, r) * xdotFx * x / r
    trace = force[0,0] + force[1,1] + force[2,2]
    result += EvalB(CacheUa, CacheUb, r) * (MulMatVec(np.transpose(force), x) + trace * x)
    return result
    
    
    
    

time = 0
timeFinal = 4
dt = 1
nSteps = int(timeFinal / dt)
nodePost = np.zeros((num_nodes,3,nSteps))
while(time < timeFinal):
    for i in range(num_nodes):
        temp = np.zeros((3))
        v0 = EvalDisp(nodePos[i,:],time)
        v1 = EvalDisp(nodePos[i,:] + 0.5 * v0,time)
        v2 = EvalDisp(nodePos[i,:] + 0.5 * v1,time)
        v3 = EvalDisp(nodePos[i,:] + v2,time)
        # 它的时间居然是纯线性的，不过按照代码来看，也是纯线性的
        # nodePos[i,:] += (v0 + 2 * v1 + 2 * v2 + v3) / 6
        nodePost[i,:,int(time/dt)] = nodePos[i,:] + (v0 + 2 * v1 + 2 * v2 + v3) / 6
    time += dt

