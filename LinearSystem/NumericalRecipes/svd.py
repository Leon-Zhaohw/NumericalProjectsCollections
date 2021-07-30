import numpy as np
# numerical recipes in c source code second edition 
# sigular value decomposition
# Singular-value decomposition
from scipy.linalg import svd
# define a matrix
A = np.array([[1, 2], [3, 4], [5, 6]])
# SVD
sci_U, sci_s, sci_Vt = svd(A)

m = A.shape[0]
n = A.shape[1]
U = A.copy()
w = np.zeros((n))
V = np.zeros((n,n)) # 未转置，未转置，未转置

rvl = np.zeros((n))
g = 0
scale = 0
anorm = 0 # Householder reduction to bidiagonal form
for i in range(n):
    l = i + 1
    rvl[i] = scale * g
    g = s = scale = 0
    if i < m:
        for k in range(i,m):
            scale += abs(U[k,i])
        if scale > 0:
            for k in range(i,m):
                U[k,i] /= scale
                s += U[k,i] * U[k,i]
            
            f = U[i,i]
            g = 0 # g = - SIGN(sqrt(s),f) 本行代码是错的
            h = f * g - s
            U[i,i] = f - g
            for j in range(n):
                s = 0
                for k in range(i,m):
                    s += U[k,i] * U[k,j]
                f = s / h
                for k in range(i,m):
                    U[k,j] += f * U[k,i]
            for k in range(i,m):
                U[k,i] *= scale
                
    w[i] = scale * g
    g = s = scale = 0
    if i < m and i < n - 1:
        for k in range(l,n):
            scale += abs(U[i,k])
        if scale > 0:
            for k in range(l,n):
                U[i,k] /= scale
                s += U[i,k] * U[i,k]
        f = U[i,l]
        g = 0 # g = - SIGN(sqrt(s),f) 本行代码是错的
        h = f * g - s
        U[i,l] = f - g
        for k in range(l,n):
            rvl[k] = U[i,k] / h
        for j in range(l,m):
            s = 0
            for k in range(l,n):
                s += U[j,k] * U[i,k]
            for k in range(l,n):
                U[j,k] += s * rvl[k]
        for k in range(l,n):
            U[i,k] *= scale
    # 接下来还有两页的算法，看不太懂，停止写了
                
        