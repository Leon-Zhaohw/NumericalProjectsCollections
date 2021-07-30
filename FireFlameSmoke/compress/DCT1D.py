import numpy as np
import scipy.fft
# scipy 的库函数结果似乎有问题
N = 8
arr = np.array([0, 0.707, 1, 0.707, 0, -0.707, -1, -0.707])
dctres = np.zeros((N))
idctres = np.zeros((N))
C = np.zeros((N,N))
for k in range(N):
    for n in range(N):
        if k == 0:
            C[k,n] = np.sqrt(1/N)
        else:
            C[k,n] = np.sqrt(2/N)*np.cos((np.pi*k*(1/2+n))/N)
for i in range(N):
    for j in range(N):
        dctres[i] += C[i,j]*arr[j]
for i in range(N):
    for j in range(N):
        idctres[i] += C[j,i]*dctres[j]