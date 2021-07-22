

 	function  QE=Equivalent_Nodal_Force_Body（ie，P0）  
 	%高斯积分法计算平面四点等参元体力等效结点载荷向量；
 	%输入量： ie——单元号，P0—体力分量或体力与自重合力
 	%返回量：体力的等效载荷向量QE
 	global  nd  ne  nj  md  mx   XY  EL  
 	   %全局变量： t0  lou ----实数： 厚度、密度；
 	  QE = zeros(8,1 ) ;      % 结构的总载荷列阵       
 	   gx = [-0.577350269189626,  0.577350269189626] ;       % 2×2  高斯积分点和权系数
 	   w = [1, 1] ;                     %将gs、w输入不同阶次数据，适用相应阶次高斯积分
 	   n0=length(w);     %确定积分阶数，可根据w的元素个数，自行调整积分阶次
 	     X0(1:4,1)=XY(EL(ie,1:4),1);         %将单元ie的4个结点x坐标值，赋予列向量X0(4)         
 	     Y0(1:4,1)=XY(EL(ie, 1:4),2);         % y坐标值，赋予列向量Y0(4)
 	for i=1: n0
 	   cxi= gx(i);
 	  for j=1: n0
 	    eta= gx(j);
 	     N=[(1-cxi) *(1-eta),  (1+cxi) *(1-eta),  (1+cxi) *(1+eta),  (1-cxi)*(1+eta)]/4;
 	    N_cxi = [-(1-eta),  (1-eta),  (1+eta),  -(1+eta)]/4;      %形函数对无量纲局部坐标偏导数
 	    N_eta =[-(1-cxi),  -(1+cxi),  (1+cxi),  (1-cxi)]/4;
 	        JA =  [ N_cxi*X0,  N_cxi*Y0           %  计算雅克比矩阵
 	               N_eta*X0,   N_eta*Y0 ];
 	       dt=det(JA);         %  计算雅克比行列式  
 	     for r=1:1:4
 	       QE(2*r-1,1)= QE(2*r-1,1) + t0*P0(1)*N(r)*dt *w(i)*w(j);       %一个高斯点的体力等效结点载荷
 	       QE(2*r,1)= QE(2*r,1) + t0*P0(2) *N(r)*dt *w(i)*w(j);      
 	     end
 	   end
 	end
 	return
