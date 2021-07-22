

 	function  Sig12 = Main_Strees(a1,a2,a3) 
 	%计算主应力、Mises应力
 	%输入量：3个应力分量---sigmx、sigmy、tauxy
 	%输出量：第一主应力、第二主应力、Mises应力
 	   b=sqrt((a1- a2)^2 + 4* a3^2 )/2;  
 	   Sig12(1) = (a1+a2)/2+b;                      %第一主应力
 	   Sig12(2) = (a1+a2)/2-b;                        %第二主应力
 	   Sig12(3) = sqrt(a1^2 + a2^2 -a1*a2 +3*a3^2);    % Mises应力，
 	  return










