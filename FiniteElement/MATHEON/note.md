https://www2.mathematik.hu-berlin.de/~cc/cc_homepage/software/software.shtml

第一篇是*Remarks around 50 lines of Matlab: short finite element implementation*。这篇文章及代码



最后要解的形式如下
$$
\bold A \bold u = b
$$
刚度矩阵



要由计算外力所造成的位移，需要解三个方程。首先第一个重要的方程是平衡方程
$$
\frac{\partial \sigma_x}{\partial x} + \frac{\partial \tau_{xy}}{\partial y} + F_x = 0 \qquad \frac{\partial \tau_{xy}}{\partial x} + \frac{\partial \sigma_{y}}{\partial y} + F_y = 0 
$$
其中sigma是应力分类，F单位体力。也就是由外力计算出应力。写成矩阵形式如下
$$
\begin{bmatrix} \partial /\partial x & 0 & \partial /\partial y \\ 0 & \partial /\partial y & \partial /\partial x\end{bmatrix}\begin{bmatrix} \sigma_x \\ \sigma_y \\ \tau_{xy}\end{bmatrix} = -\begin{bmatrix} F_x \\ F_y\end{bmatrix}
$$
简写如下
$$
[\bold A]^T[\bold \sigma] = \{\bold f\}
$$
第二个重要的方程是本构关系，也就是平面应力。
$$
\begin{bmatrix} \sigma_x \\ \sigma_y \\ \tau_{xy}\end{bmatrix} = \frac{E}{1-v^2}\begin{bmatrix} 1 & v & 0 \\ v & 1 & 0 \\ 0 & 0 & \frac{1-v}{2}\end{bmatrix}\begin{bmatrix} \varepsilon_x \\ \varepsilon_y \\ \gamma_{xy}\end{bmatrix}
$$
其中E为弹性模量，nu为泊松比，varepsilon和gamma都是应变分量。也就是由上一步算出的应力，计算出应变。也可简写如下
$$
\{\sigma\} = [D]\{\varepsilon\}
$$
第三个重要的方程是应变-位移关系
$$
\begin{bmatrix} \varepsilon_x \\ \varepsilon_y \\ \gamma_{xy}\end{bmatrix} = \begin{bmatrix} \partial /\partial x & 0 \\ 0 & \partial /\partial y \\ \partial  /\partial y & \partial /\partial x\end{bmatrix}\begin{bmatrix} u \\ v\end{bmatrix}
$$
其中u和v为位移分量。也就是由上一步计算出的应变，算出最终的位移。简写如下
$$
\{\varepsilon\} = [A]\{\bold u\}
$$
不过注意A矩阵，在上面的数学表达式中，它是三行两列的。但实际上我们更经常用的是形函数来表达这个关系，对于线性三角形来说，需要三个形函数，那么它就是三行六列的。


$$
\begin{bmatrix} \varepsilon_x \\ \varepsilon_y \\ \gamma_{xy}\end{bmatrix} = \begin{bmatrix} \partial N_1/\partial x & 0 & \partial N_2/\partial x & 0 & \partial N_3/\partial x & 0 \\ 0 & \partial N_1/\partial y & 0 & \partial N_2/\partial y & 0 & \partial N_3/\partial y \\ \partial N_1 /\partial y & \partial N_1/\partial x & \partial N_2 /\partial y & \partial N_2/\partial x  & \partial N_3 /\partial y & \partial N_3/\partial x \end{bmatrix}\begin{bmatrix} u_1 \\ v_1 \\ u_2 \\ v_2 \\ u_3 \\ v_3\end{bmatrix}
$$
右边那个u1v1u2v2u3v3就是三角形三个顶点。而中间
$$
N_1 = \xi  \qquad N_2 = \eta \qquad N_3 = 1 - \xi - \eta
$$

$$
\frac{\partial N_1}{\partial x} = \frac{\partial N_1}{\partial \xi}\frac{\partial \xi}{\partial x}  + \frac{\partial N_1}{\partial \eta}\frac{\partial \eta}{\partial x} = \frac{1}{2|Area|}(y_2-y_3)\\
\frac{\partial N_1}{\partial y} = \frac{\partial N_1}{\partial \xi}\frac{\partial \xi}{\partial y}  + \frac{\partial N_1}{\partial \eta}\frac{\partial \eta}{\partial y} = \frac{1}{2|Area|}(x_3-x_2)
$$
等式右边那坨奇怪的式子，来源如下。首先
$$
\begin{bmatrix} dN/d\xi \\dN/d\eta \end{bmatrix} = \begin{bmatrix} dx/d\xi & dy/d\xi\\dx/d\eta & dx/d\eta \end{bmatrix}\begin{bmatrix} dN/dx \\dN/dy \end{bmatrix}
$$
其中雅可比矩阵为
$$
J = \frac{\partial(x,y)}{\partial (\xi,\eta)} = \begin{bmatrix} x2 - x1 & y2 - y1 \\x3 - x1 & y3 - y1\end{bmatrix} \qquad |J| = 2|Area|
$$
那么上式反过来变换就是
$$
\begin{bmatrix} dN/dx \\dN/dy \end{bmatrix} = \frac{1}{2|Area|}\begin{bmatrix} y_3 - y_1 & y_1 - y_2\\x_1 - x_3 & x_2 - x_1 \end{bmatrix}\begin{bmatrix} dN/d\xi \\dN/d\eta \end{bmatrix}
$$
然后
$$
\frac{\partial x}{\partial \xi} = x_1 - x_3 \qquad \frac{\partial y}{\partial \xi} = y_1 - y_3 \\ \frac{\partial x}{\partial \eta} = x_2 - x_3 \qquad \frac{\partial y}{\partial \eta} = y_2 - y_3
$$

$$
M_{ij} = \rho\int_\Omega \nabla N_i \cdot \nabla N_j d\Omega
$$
问题一，上面这个和下面这个有什么不同？
$$
\bold K = \iint \bold B^T \bold D \bold B dx \qquad \bold B = \frac{d \bold N}{dx}
$$
右边项
$$
b_j = \sum_{T \in \mathcal{T}}f\eta_j dx + \sum_{E \in \Gamma_N}\int_E g\eta_j ds - \sum_{k=1}^NU_k \sum_{T\in \mathcal{T}}\int_T \nabla \eta_j \cdot \nabla \eta_kdx
$$
可以使用最小势能原理给出弹性平面单元刚度方程式。单位厚度上单元的应变能为
$$
U = \iint \frac{1}{2}\{\bold \sigma\}^T\{\bold \epsilon\}dxdy = \\
\frac{1}{2}\{\bold u\}^T \iint([\bold A][\bold S]^T)\bold D([\bold A][\bold S])dxdy\{\bold u\}\\
=\frac{1}{2}\{\bold u\}^T \iint [\bold B]^T[\bold D][\bold B]dxdy\{\bold u\}
$$
其中
$$
[\bold A] = \begin{bmatrix} \partial /\partial x & 0 \\ 0 & \partial /\partial y \\ \partial/\partial y & \partial /\partial x\end{bmatrix}
$$


问题二：这个的逆究竟代表什么

```python
G = np.array([[1,1,1],[x0,x1,x2],[y0,y1,y2]])
G0 = np.array([[0,0],[1,0],[0,1]])
PhiGrad = np.dot(np.linalg.inv(G),G0)
```

$$
G_0 = \begin{bmatrix} 1 & 1 & 1 \\ x_1 & x_2 & x_3 \\ y_1 & y_2 & y_3\end{bmatrix}
$$

首先求det(G0)，也就是
$$
det(G_0) = x_2y_3 - y_2x_3 + x_3y_1 - x_1y_3 + x_1y_2 - x_2y_1
$$
https://zh.wikihow.com/%E6%B1%823X3%E7%9F%A9%E9%98%B5%E7%9A%84%E9%80%86%E7%9F%A9%E9%98%B5

首先转置
$$
G_0^T = \begin{bmatrix} 1 & x_1 & y_1 \\ 1 & x_2 & y_2 \\ 1 & x_3 & y_3\end{bmatrix}
$$
求它的伴随矩阵，从左到右，从上到下依次是
$$
\begin{bmatrix}(x_2y_3 - x_3y_2) & -(y_3-y_2)& (x_3 - x_2) \\ -(x_1y_3 - x_3y_1) & (y_3 - y_1) & -(x_3 - x_1) \\ (x_1y_2 - x_2y_1) & -(y_2 - y_1) & (x_2 - x_1)\end{bmatrix}
$$
真正的G如下
$$
G = \begin{bmatrix} y_2 - y_3& x_3 - x_2 \\  y_3 - y_1 & x_1 - x_3 \\  y_1 - y_2 & x_2 - x_1\end{bmatrix}
$$
这就和pdf上的那个公式对应上了
$$
M_{jk} = \int_T \nabla \eta_j (\nabla \eta_k)^T dx = \frac{|T|}{(2|T|)^2}\begin{bmatrix} y_{j+1} - y_{j+2},x_{j+2} - x_{j+1}\end{bmatrix}\begin{bmatrix} y_{k+1} - y_{k+2} \\ x_{k+2} - x_{k+1}\end{bmatrix}
$$
不过它写了
$$
\begin{bmatrix} \nabla N_1 \\ \nabla N_2 \\ \nabla N_3\end{bmatrix} = \frac{1}{2|T|}\begin{bmatrix} y_2 - y_3& x_3 - x_2 \\  y_3 - y_1 & x_1 - x_3 \\  y_1 - y_2 & x_2 - x_1\end{bmatrix} = \begin{bmatrix} 1 & 1 & 1 \\ x_1 & x_2 & x_3 \\ y_1 & y_2 & y_3\end{bmatrix}^{-1}\begin{bmatrix} 0 & 0 \\ 1 & 0 \\ 0 & 1\end{bmatrix}
$$
