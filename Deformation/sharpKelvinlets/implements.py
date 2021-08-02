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
_a = 0
_b = 0
force = np.zeros((3))
action = 0
# << "ACTION CAN BE\n"
# << "...0 - BrushGrab\n"
# << "...1 - BrushGrabBiScale\n"
# << "...2 - BrushGrabTriScale\n"
# << "...3 - BrushGrabLaplacian\n"
# << "...4 - BrushGrabBiLaplacian\n"
# << "...5 - BrushGrabCusp\n"
# << "...6 - BrushGrabCuspLaplacian\n"
# << "...7 - BrushGrabCuspBiLaplacian\n"
# << "...8 - BrushAffine (Twist)\n"
# << "...9 - BrushAffine (Scale)\n"
# << "..10 - BrushAffine (Pinch)\n";

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

def BrushGrabInit():
    global _a
    global _b
    global force
    _a = (np.pi * mu) / 4
    _b = _a / (1 - nu) / 4
    te = _a
    force[1] = 1 / EvalRadial(0, eps, _a, _b)

def BrushGrabEval(x):
    r = np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
    
    u = EvalRadial(r, eps, _a, _b) * force
    if r < 1e-5:
        return u
    dotres = force[0]*x[0] + force[1]*x[1] + force[2]*x[2]
    u += EvalBulge(r,eps,_b) * x * dotres
    return u
def BrushGrabBiScaleInit():
    return 0
def BrushGrabBiScaleEval(x):
    return 0
def BrushGrabTriScale():
    return 0
def BrushGrabLaplacian():
    return 0
def BrushGrabBiLaplacian():
    return 0
def BrushGrabCusp():
    return 0
def BrushGrabCuspLaplacian():
    return 0
def BrushGrabCuspBiLaplacian():
    return 0
def BrushAffineTwist():
    return 0
def BrushAffineScale():
    return 0
def BrushAffinePinch():
    return 0

if action == 0:
    BrushGrabInit()








for i in range(num_nodes):
    temp = np.zeros((3))
    if action == 0:
        temp = BrushGrabEval(nodePos[i,:])
    elif action == 1:
        temp = BrushGrabBiScale(nodePos[i,:])
    elif action == 2:
        temp = BrushGrabTriScale(nodePos[i,:])
    elif action == 3:
        temp = BrushGrabLaplacian(nodePos[i,:])
    elif action == 4:
        temp = BrushGrabBiLaplacian(nodePos[i,:])
    elif action == 5:
        temp = BrushGrabCusp(nodePos[i,:])
    elif action == 6:
        temp = BrushGrabCuspLaplacian(nodePos[i,:])
    elif action == 7:
        temp = BrushGrabCuspBiLaplacian(nodePos[i,:])
    elif action == 8:
        temp = BrushAffineTwist(nodePos[i,:])
    elif action == 9:
        temp = BrushAffineScale(nodePos[i,:])
    elif action == 10:
        temp = BrushAffinePinch(nodePos[i,:])
    nodePos[i,:] += temp

