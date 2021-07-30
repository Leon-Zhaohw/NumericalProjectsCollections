import numpy as np
# numerical recipes in c source code second edition lu decomposition
n = 4
A = np.array([[2, 5, 8, 7], [5, 2, 2, 8], [7, 5, 6, 6], [5, 4, 4, 8]])
b = np.array([22,17,24,21])
LU = np.zeros((n,n))
lower = np.zeros((n,n))
upper = np.zeros((n,n))
for j in range(n):
    lower[j,j] = 1
    for i in range(j+1):
        summ = 0
        for k in range(i):
            summ += lower[i,k] * upper[k,j]
        upper[i,j] = A[i,j] - summ
    for i in range(j+1,n):
        summ = 0
        for k in range(j):
            summ += lower[i,k] * upper[k,j]
        lower[i,j] = (A[i,j] - summ) / upper[j,j]
        
from scipy.linalg import lu

p0, l0, u0 = lu(A)