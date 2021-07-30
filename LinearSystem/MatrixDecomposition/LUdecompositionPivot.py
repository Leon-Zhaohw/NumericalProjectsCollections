# -*- coding: utf-8 -*-
"""
Created on Tue May  9 11:44:53 2017
@author: Leandro
"""

import numpy as np

def pivotear(A, indexs, j):
    coluna = abs(np.transpose(A)[j])
    maiorj = j
    for i in range(j+1,len(coluna)):
        if(coluna[i] > coluna[j]):
            coluna[j], coluna[i] = coluna[i], coluna[j]
            maiorj = i
    if(maiorj != indexs[j]):
        indexs[j], indexs[maiorj] = indexs[maiorj], indexs[j]
        for k in range(len(indexs)):
            A[indexs[j]][k], A[indexs[maiorj]][k] = A[indexs[maiorj]][k], A[indexs[j]][k]

def LU_decomp(A):
    LU = np.copy(A)
    n  = len(A)
    indexs = np.array([i for i in range(n)])
    for j in range(n-1):
        pivotear(LU, indexs, j)
        for i in range(j+1,n):
            LU[i][j] /= LU[j][j]
            for k in range(j+1,n):
                LU[i][k] -= LU[i][j]*LU[j][k]
    return LU, indexs

def LU_forwardsub(L,b): #triangular inferior com diagonal composta apenas por uns
    n = len(L)
    y = np.copy(b)
    for i in range(1,n):
        for j in range(i):
            y[i] -= L[i][j]*y[j]
    return y

def LU_backwardsub(U,y): #triangular superior
    n = len(U)
    x = np.copy(y)
    for i in range(n-1,-1,-1): #de n-1 até 0
        for j in range(i+1,n):
            x[i] -= U[i][j]*x[j]
        x[i] /= U[i][i]
    return x

def LU_solve(LU,b):
    y = LU_forwardsub(LU, b) #Ly = b
    x = LU_backwardsub(LU, y) #Ux = y
    return x

def LU_det(U): #det(LU)=det(L)*det(U)=1*det(U)
    det = 1
    for i in range(len(U)):
        det *= U[i][i]
    return det

A = np.array([[1.,-1., 1.,-1.], \
              [1., 0., 0., 0.], \
              [1., 1., 1., 1.], \
              [1., 2., 4., 8.]])
b = np.array([-2., 0., 0., 6.])
LU, indexs = LU_decomp(A)
x =  LU_solve(LU, b[indexs])
print("Solução:\n x = {}".format(x.round(2)))

A2 = np.array([[1., 2.],\
               [3., 4.]])
b2 = np.array([1., 2.])
LU, indexs = LU_decomp(A2)
x =  LU_solve(LU, b2[indexs])
print("Solução:\n x = {}".format(x.round(2)))