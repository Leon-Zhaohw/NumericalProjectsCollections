function D =Elastic_Matrix (pm, E, nv)           
 	% 根据问题属性生成平面问题的弹性矩阵
 	%输入量:，E-弹性模量，nv—泊松比 
 	%pm-问题属性: 平面应力问题—1, 平面应变问题—2, 轴对称问题—3, 空间问题—4。
 	%返回值：弹性矩阵D
 	switch pm
 	 case 1                                    %pm=1，平面应力问题
 	     D = [1  nv  0 ; 
 	         nv  1  0 ; 
 	         0  0 (1-nv)/2] *E/(1 - nv*nv);
 	 case 2                                    %pm=2，平面应变问题
 	     D = [1-nv  nv  0 ; 
 	          nv  1-nv  0 ; 
 	          0    0 (1-2*nv)/2] *E/(1+nv)/(1-2*nv);
 	 case 3                                     %pm=3，轴对称问题
 	     D = [1-nv  nv  nv   0 ; 
 	          nv  1-nv  nv  0 ; 
 	          nv  nv   1-nv  0 ; 
 	           0  0 (1-2*nv)/2] *E/(1+nv)/(1-2*nv);
 	  case 4                                     %pm=4，空间问题
 	     D11 = [1-nv  nv  nv ;
 	            nv  1-nv  nv ; 
 	            nv  nv   1-nv ] *E/(1+nv)/(1-2*nv);
 	     D12=zeros(3);
 	     D22=eye(3)*E/(1+nv)/2;
 	     D = [D11, D12; 
 	         D12, D22];
 	%disp('  输入的问题属性错误，（平面应力问题-1，平面应变问题-2）请修改！')
 	end
 	return



