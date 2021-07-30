import numpy as np

# 幂方法求主特征向量
A = np.array([[0,11,-5],[-2,17,-7],[-4,26,-10]])
eigvec = np.array([1,1,1],dtype = float) # 猜一个主特征向量
eig = 0 # 猜一个主特征值
diff = 100
eps = 1e-8
while(diff > eps):
    eigold = eig
    eigvec = np.dot(A,eigvec)
    eig = -1 # 求绝对值最大的分量，值
    for k in range(3):
        if eig < abs(eigvec[k]):
            eig = eigvec[k]
    for k in range(3):
        eigvec[k] /= eig
    diff = abs(eig - eigold)
        