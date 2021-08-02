import numpy as np

axis_element = 100
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
    
eps = 4
nu = 0.45
mu = 5
force = np.zeros((3))
_a = 1 / 4 / np.pi / mu
_b = _a / (1 - nu) / 4

def EvalA(r):
    R = r / eps
    if R > 5:
        return 0
    eps5 = pow(eps,5)
    sqR = R * R
    R1 = np.sqrt(sqR + 1)
    R1_3 = pow(R1,3)
    R1_9 = pow(R1_3,3)
    R1_10 = R1_9 * R1
    
    iterm0 = 512 * sqR + 2304
    iterm0 = iterm0*sqR + 4032
    iterm0 = iterm0*sqR + 3360
    iterm0 = iterm0*sqR + 1260
    iterm0 = iterm0*sqR + 105
    iterm0 = iterm0 * R1 - 512 * R * R1_10
    iterm0 = iterm0 * _a
    
    iterm1 = 448 + 128 * sqR
    iterm1 = iterm1 * sqR + 560
    iterm1 = iterm1 * sqR + 280
    iterm1 = iterm1 * sqR + 35
    iterm1 = 128 * R * R1_10 - iterm1 * R1_3
    iterm1 = iterm1 * 2 * _b
    
    return (iterm0 + iterm1) * 9 / (eps5 * R1_10)


def EvalB(r):
    R = r / eps
    if R > 5:
        return 0
    eps5 = pow(eps,5)
    eps7 = pow(eps,7)
    
    sqR = R * R
    R1 = np.sqrt(sqR + 1)
    R1_3 = pow(R1,3)
    R1_9 = pow(R1_3,3)
    
    Rterm = 576 + 128 * sqR
    Rterm = Rterm * sqR + 1008
    Rterm = Rterm * sqR + 840
    Rterm = Rterm * sqR + 315
    Rterm = 128 * R1_9 - Rterm * R
    return Rterm * (18 * _b) / (eps7 * R * R1_9)
    

force[1] = 1 / EvalA(0)

for i in range(num_nodes):
    temp = np.zeros((3))
    x = nodePos[i,:]
    r = np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
    u = EvalA(r) * force
    if r < 1e-4:
        temp = u.copy()
    else:
        dotres = force[0]*x[0] + force[1]*x[1] + force[2]*x[2]
        u += EvalB(r) * dotres * x
        temp = u.copy()
    nodePos[i,:] += temp

