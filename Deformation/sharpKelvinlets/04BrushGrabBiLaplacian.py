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
    
    eps2 = eps * eps
    eps4 = eps2 * eps2
    r2 = r * r
    re2 = r2 + eps2
    re = np.sqrt(re2)
    re11 = pow(re,11)
    
    term = 52.5 * eps4 * (3 * _a * (eps2 - 2 * r2) - 2 * _b * re2)
    return term / re11
    


def EvalB(r):
    
    eps2 = eps * eps
    eps4 = eps2 * eps2
    r2 = r * r
    re2 = r2 + eps2
    re = np.sqrt(re2)
    re11 = pow(re,11)
    
    return 945 * _b * eps4 / re11
    

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

