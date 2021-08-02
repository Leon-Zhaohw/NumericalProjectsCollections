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
force = np.zeros((3,3))
force[0,0] = 0.75
force[1,1] = -0.75
_a = 1 / 4 / np.pi / mu
_b = _a / (1 - nu) / 4

def EvalA(r):
    eps2 = eps * eps
    eps4 = eps2 * eps2
    r2 = r * r
    re2 = r2 + eps2
    re = np.sqrt(re2)
    re7 = pow(re,7)
    return (_b * re2 * (5 * eps2 + 2 * r2) - 7.5 * _a * eps4)/re7


def EvalB(r):
    eps2 = eps * eps
    r2 = r * r
    re2 = r2 + eps2
    re = np.sqrt(re2)
    re7 = pow(re,7)
    return -3 * _b * (7 * eps2 + 2 * r2) / re7

def scaledot(vec0,vec1):
    n = len(vec0)
    res = np.zeros((n))
    for i in range(n):
        res[i] = vec0[i] * vec1[i]
    return res

def scalerDot(vec0,vec1):
    n = len(vec0)
    res = 0
    for i in range(n):
        res = vec0[i] * vec1[i]
    return res

def mulMatVec(mat,vec):
    n = len(mat)
    res = np.zeros((n))
    for i in range(n):
        for j in range(n):
            res[i] += mat[i,j] * vec[j]
    return res

def AssembleSkewSymMatrix(A):
    res = np.zeros((3,3))
    res[0,1] = -A[2]
    res[0,2] = A[1]
    res[1,0] = A[2]
    res[1,2] = -A[0]
    res[2,0] = -A[1]
    res[2,1] = A[0]
    return res

e3 = eps * eps * eps
q = np.zeros((3))
q[0] = (force[2,1] - force[1,2])/2
q[1] = (force[0,2] - force[2,0])/2
q[2] = (force[1,0] - force[0,1])/2

s = (force[0,0] + force[1,1] + force[2,2]) / 3
iden = np.array([[1,0,0],[0,1,0],[0,0,1]])
pinch = (force + np.transpose(force)) / 2 - s * iden

q = q * -0.4 * e3 / _a
if nu == 0.5:
    s = 0
else:
    s = s * 0.4 * e3 / (2 * _b - _a)
pinch = pinch * (2 * e3 / (4 * _b - 5 * _a))

force = s * iden
force += AssembleSkewSymMatrix(q)
force += pinch
        

for i in range(num_nodes):
    temp = np.zeros((3))
    x = nodePos[i,:]
    r = np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
    if r < 1e-4:
        continue
    re2 = r*r + eps*eps
    re = np.sqrt(re2)
    re3= re2*re
    summ = np.zeros((3))
    mulres = mulMatVec(force, x)
    summ -= _a / re3 * (1 + 1.5*eps*eps/re2) * mulres
    summ += _b / re3 * mulMatVec((force + np.transpose(force)), x)
    forceTrace = force[0,0] + force[1,1] + force[2,2]
    term = scalerDot(x, mulres)
    summ += _b / re3 * (forceTrace - 3 / re2 * term) * x
    temp = summ.copy()
    nodePos[i,:] += temp

