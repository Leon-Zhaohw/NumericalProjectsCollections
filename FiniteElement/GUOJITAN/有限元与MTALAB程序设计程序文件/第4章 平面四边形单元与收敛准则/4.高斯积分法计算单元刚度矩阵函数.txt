

 	function KE =Gauss_Stiff_Matrix( ie,D )   
 	%   用高斯积分法求解平面等参数单元的刚度矩阵
 	%  输入参数：   ie——单元号
 	%  返回量：  KE ——单元刚度矩阵
 	global  pm  E  nv  t0  nd  ne  XY  EL     %用全局变量定义结点坐标、单元信息等
 	KE = zeros( 8,8) ;
 	  gx = [-0.577350269189626,  0.577350269189626] ;       % 2×2  高斯积分点和权系数
 	  w = [1, 1] ;                      %将gs、w输入不同阶次数据，适用相应阶次高斯积分
 	  n0=length(w);     %确定积分阶数，可根据w的元素个数，自行调整积分阶次  
 	for i=1: n0
 	  for j=1: n0
 	     [B, dt] = Quad_B4_Matrix( ie, gx(i), gx(j) );   %调用函数计算几何矩阵、雅克比行列式 
 	        KE =KE +t0*transpose(B)*D*B*dt*w(i)*w(j);            %一个高斯点的单刚 
 	   end
 	end
 	return














