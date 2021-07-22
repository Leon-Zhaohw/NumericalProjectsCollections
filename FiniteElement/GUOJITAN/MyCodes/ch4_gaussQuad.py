import numpy as np
# 有限元法与MATLAB程序设计 郭吉坦 薛齐文
# 未完成，搞不懂那些面力之类的
E = 2.06e5 # 弹性模量
nu = 0.3 # 泊松比
gravity = 7.7e-05 # 重力
t0 = 2 # 厚度

num_nodes = 12 # 节点数
coord = np.array([[0,20],
               [0,0],
               [20,20],
               [20,0],
               [40,0],
               [60,20],
               [80,20],
               [80,0],
               [100,20],
               [100,0]])

num_elem = 5 # 单元总数
element = np.array([[0,1,3,2],
                    [2,3,5,4],
                    [4,5,7,6],
                    [6,7,9,8],
                    [8,9,10,11]])

num_constraint = 4 # 位移约束数
boundary = np.array([[0,0],
                     [1,0],
                     [0,1],
                     [1,1]])

num_force = 2 # 集中力数据
force = np.array([[10,1000,0],[11,-1000,0]])


D = np.zeros([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])*E/(1-nu*nu) # 平面应力
D = np.zeros([[1-nu,nu,0],[nu,1-nu,0],[0,0,(1-2*nu)/2]])*E/(1+nu)/(1-2*nu) # 平面应变

Kmat = np.zeros((2*num_nodes,2*num_nodes))
bmat = np.zeros((2*num_nodes))

# 计算总体刚度矩阵
for ie in range(num_elem):
    p0 = element[ie,0]
    p1 = element[ie,1]
    p2 = element[ie,2]
    p3 = element[ie,3]
    
    X0 = [coord[p0,0],coord[p1,0],coord[p2,0],coord[p3,0]]
    Y0 = [coord[p0,1],coord[p1,1],coord[p2,1],coord[p3,1]]
    
    ke = np.zeros((8,8))
    gx = [-0.57735,0.57735]
    w = [1,1]
    for i in range(2):
        for j in range(2):
            xi = gx[i]
            eta = gx[j]
            N_xi = [-(1-eta),(1-eta),(1+eta),-(1+eta)]/4
            N_eta = [-(1-xi),-(1+xi),(1+xi),(1-xi)]/4
            
            J11 = np.dot(N_xi,X0)
            J12 = np.dot(N_xi,Y0)
            J21 = np.dot(N_eta,X0)
            J22 = np.dot(N_eta,Y0)
            
            J = np.array([[J11,J12],[J21,J22]])
            det = np.linalg.det(J)
            Jinv = np.linalg.inv(J)
            
            N_x = Jinv[0,0] * N_xi + Jinv[0,1] * N_eta
            N_y = Jinv[1,0] * N_xi + Jinv[1,1] * N_eta
            
            B = np.zeros([[N_x[0],0,N_x[1],0,N_x[2],0,N_x[3],0],
                          [0,N_y[0],0,N_y[1],0,N_y[2],0,N_y[3]],
                          [N_y[0],N_x[0],N_y[1],N_x[1],N_y[2],N_x[2],N_y[3],N_x[3]]])
            
            ke = ke + t0 * np.dot(np.dot(np.transpose(B),D),B) * det * w[i] * w[j]
            

    