import numpy as np
# 有限元法与MATLAB程序设计 郭吉坦 薛齐文
# 未完成，搞不懂那些面力之类的
E = 2.06e5 # 弹性模量
nu = 0.3 # 泊松比
gravity = 7.7e-05 # 重力
t0 = 2 # 厚度

num_nodes = 8 # 节点数
coord = np.array([[20,0],
               [40,0],
               [60,0],
               [0,20],
               [20,20],
               [40,20],
               [0,40],
               [20,40]])

num_elem = 7 # 单元总数
element = np.array([[6,3,4],
                    [4,7,6],
                    [7,4,5],
                    [0,4,3],
                    [4,0,1],
                    [1,5,4],
                    [5,1,2]])

num_constraint = 5 # 位移约束数
boundary = np.array([[6,0],
                     [3,0],
                     [0,1],
                     [1,1],
                     [2,1]])

num_force = 2
force = np.array([[2,1000,0],[7,0,2000]])

num_volumeForce = 0

mx = 2 # 作用表面力线段的个数
qmx = np.array([[2,5,0],[5,7,0]])

md = 3 # 作用表面里结点的个数
qmd = np.array([[3,1,40],[6,1,30],[8,1,20]])

D = np.zeros([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])*E/(1-nu*nu) # 平面应力
D = np.zeros([[1-nu,nu,0],[nu,1-nu,0],[0,0,(1-2*nu)/2]])*E/(1+nu)/(1-2*nu) # 平面应变

Kmat = np.zeros((2*num_nodes,2*num_nodes))
bmat = np.zeros((2*num_nodes))

# 计算总体刚度矩阵
for ie in range(num_elem):
    p0 = element[ie,0]
    p1 = element[ie,1]
    p2 = element[ie,2]
    
    b1 = coord[p1,1] - coord[p2,1]
    b2 = coord[p2,1] - coord[p0,1]
    b3 = coord[p0,1] - coord[p1,1]
    
    c1 = coord[p2,0] - coord[p1,0]
    c2 = coord[p1,0] - coord[p0,0]
    c3 = coord[p0,0] - coord[p2,0]
    
    A = (b2 * c3 - b3 * c2) / 2 # 面积
    B = np.array([[b1,0,b2,0,b3,0],
                  [0,c1,0,c2,0,c3],
                  [c1,b1,c2,b2,c3,b3]]) / (2 * A)
    
    Ke = np.dot(np.transpose(B),D)
    Ke = np.dot(Ke,B) * A * t0
    
    idx = np.zeros((6),dtype = int)
    idx[0] = p0 * 2
    idx[1] = p0 * 2 + 1
    idx[2] = p1 * 2
    idx[3] = p1 * 2 + 1
    idx[4] = p2 * 2
    idx[5] = p2 * 2 + 1
    for i in range(6):
        for j in range(6):
            Kmat[idx[i],idx[j]] += Ke[i,j]
    
# 计算体力的等效结点载荷
for ie in range(num_elem):
    p0 = element[ie,0]
    p1 = element[ie,1]
    p2 = element[ie,2]
    
    x0 = coord[p0,0]
    x1 = coord[p1,0]
    x2 = coord[p2,0]
    y0 = coord[p0,1]
    y1 = coord[p1,2]
    y2 = coord[p2,3]
    
    A = (x0*(y1 - y2) + x1*(y2 - y0) + x2*(x0 - y2))/2
    
    # 重力的影响
    ga = - A * t0 * gravity / 3
    bmat[p0*2+1] += ga
    bmat[p1*2+1] += ga
    bmat[p2*2+1] += ga
    
# 组装节点集中力
for i in range(num_force):
    idx = int(force[i,0])
    bmat[idx,0] += force[i,1]
    bmat[idx,1] += force[i,2]
    
    