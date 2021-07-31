物体的形变力，由旋转，缩放，反射，位移组成。

要移去位移，那么
$$
\frac{\partial \phi(x)}{\partial x} = \frac{\partial (Fx+t)}{\partial x} = F
$$
第一个尝试是Dirichlet energy 计算如下
$$
\psi_{Dirichlet} = ||\bold F||_F^2 = \sum_{i=0}^2 \sum_{i=0}^2 f_{i,j}^2
$$
当F是零的时候，Dirichlet能量当然是零。但是当F = I，也没有形变，因此Dirichlet 能量也应当等于零，但是计算出来不是这样。

接下来尝试Neo-Hookean
$$
\psi_{neo}= ||\bold F||_F^2 - 3
$$
然而其它模式也返回零，也可能返回负数。



我们希望某个函数的计算结果，当F是纯旋转的时候，这个函数的计算结果为零。这就是St.Venant Kirichhoff
$$
\psi = \frac{1}{2}||\bold F^T\bold F - \bold I||_F^2 \\
= ||\bold F^T\bold F||_F^2 + tr\bold I - 2 tr(\bold F^T\bold F)
$$
同样还有几个别有深意的名字。首先是The right Cauchy-Green tensor。
$$
C = F^TF
$$
以及GreenStrain
$$
E =\frac{1}{2}(F^TF - I)
$$


少一点非线性版本的
$$
\psi_{ARAP} = ||\bold F - \bold R||^2_F = ||\bold F||_F^2 + ||\bold R||^2_F - tr(\bold F^T\bold R) \\= ||\bold F||_F^2 + 3 - 2tr(\bold S)
$$
the deformations measures directly generates forces.

计算力
$$
\frac{\partial \psi}{\partial x} = \frac{\partial \bold F}{\partial \bold x} : \frac{\partial \psi}{\partial \bold F}
$$
也叫first Piola-Kirchhoff stress tensor，对于Dirichlet力来说
$$
\frac{\partial \psi_D}{\partial \bold F} = \bold P(\bold F) = 2\bold F
$$
对于St.Kirchhoff-stretch来说
$$
\psi = \frac{1}{4}||\bold F^T\bold F - \bold I||^2_F \\
= \frac{1}{4}(||\bold F^T\bold F||_F^2 + tr\bold I - 2tr(\bold F^T\bold F))
$$
那么PK1力推导如下
$$
\bold P(\bold F) = \frac{1}{4}(4\bold F \bold F^T\bold F - 4\bold F) = \bold F(\bold F^T\bold F - \bold I) = \bold F\bold E
$$
计算Force Gradient
$$
\frac{\partial^2 \psi}{\partial x^2} = \frac{\partial \bold F^T}{\partial \bold x}\frac{\partial^2 \psi}{\partial \bold F^2}\frac{\partial \bold F}{\partial \bold x}
$$
那么Cauchy-Green不变量

如何把一个矩阵变成一个标量？可以使用柯西格林不变量。
$$
I_C = ||\bold F||^2_F \qquad II_C = ||\bold F^T\bold F||_F^2 \qquad III_C = det(\bold F^T\bold F)
$$
如果我们能算出三个不变量的梯度和Hessian矩阵，那么计算能量的力梯度将变得容易。剩下的只需简单的运算即可得出。这就像在写“快乐暑假”的作业的时候，正在吐槽作业最后怎么没参考答案的时候，一份详细的标准答案突然出现在面前，最后只需要简单的复制粘贴就可以了。

These are usually specified using the Lamé parameters, µ and λ. The µ parameter controls
length preservation, and is also sometimes called the shear modulus or Lamé’s second
parameter. The λ parameter controls volume preservation, and is sometimes called Lamé’s
first parameter.6 If you want lots of volume preservation, you set λ to be much larger than
µ. Conversely, if you want more length preservation, you push up the value of µ  

就最后这个我看懂了

F 可以分成极化分解，首先使用SV分解
$$
F = U\Sigma V^T \qquad R = UV^T \qquad S = V\Sigma V^T
$$
U和V都是酉矩阵unitary Matrix
$$
RS = UV^TV\Sigma V^T = UI \Sigma V^T = U\Sigma V^T =F
$$
