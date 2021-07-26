% convert axis-angle representation into a 3x3 rotation matrix
function [R] = axisAngle (axis, theta)

R = zeros(3,3);
axis = axis / norm(axis, 2);
cosTheta = cos(theta); 
sinTheta = sin(theta); 
versine = (1.0 - cosTheta);

R(1,1) = axis(1) * axis(1) * versine + cosTheta;
R(1,2) = axis(1) * axis(2) * versine - axis(3) * sinTheta;  
R(1,3) = axis(1) * axis(3) * versine + axis(2) * sinTheta;  

R(2,1) = axis(2) * axis(1) * versine + axis(3) * sinTheta;  
R(2,2) = axis(2) * axis(2) * versine + cosTheta;
R(2,3) = axis(2) * axis(3) * versine - axis(1) * sinTheta;  

R(3,1) = axis(3) * axis(1) * versine - axis(2) * sinTheta;  
R(3,2) = axis(3) * axis(2) * versine + axis(1) * sinTheta;  
R(3,3) = axis(3) * axis(3) * versine + cosTheta; 

end
