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

def EvalRadial(r,e,a,b):
    e2 = e * e
    re2 = r * r + e2
    re = np.sqrt(re2)
    return (a - b) / re + a * e2 / re / re2 / 2

def EvalBulge(r,e,b):
    e2 = e * e
    re2 = r * r + e2
    re = np.sqrt(re2)
    return b / re / re2

def EvalA(r):
    A0 = EvalRadial(r,eps,_a,_b)
    A1 = EvalRadial(r,1.1*eps,_a,_b)
    return A0 - A1

def EvalB(r):
    B0 = EvalBulge(r,eps,_b)
    B1 = EvalBulge(r,1.1*eps,_b)
    return B0 - B1

force[1] = 1 / EvalA(0)
for i in range(num_nodes):
    temp = np.zeros((3))
    x = nodePos[i,:]
    r = np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
    u = EvalA(r) * force
    if r < 1e-4:
        temp = u
    dotres = force[0]*x[0] + force[1]*x[1] + force[2]*x[2]
    u += EvalB(r) * dotres * x
    temp = u.copy()
    nodePos[i,:] += temp

