https://zhuanlan.zhihu.com/p/114626779

离散余弦变换
$$
[X(0),...,X(N-1)] = [x(0),...,x(N-1)]\begin{bmatrix}C_{0,0} & ... & C_{0,N-1} \\ \\ \\ C_{N-1,0} & ... & C_{N-1,N-1} \end{bmatrix}
$$
其中C是系数矩阵，如下
$$
C = \sqrt{\frac{2}{N}}\begin{bmatrix} \sqrt{1/2}\cos(\pi/2N) & ... & \cos ((N-1)\pi/2N) \\ \sqrt{1/2}\cos(3\pi/2N) & ... & \cos ((N-1)3\pi/2N) \\ \\ \sqrt{1/2}\cos((2N-1)\pi/2N) & ... & \cos ((N-1)(2N-1)\pi/2N)\end{bmatrix}
$$
