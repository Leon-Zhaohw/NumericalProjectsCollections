function [R] = axisAngle2D(theta)

cosTheta = cos(theta); 
sinTheta = sin(theta); 

R = [cosTheta -sinTheta;
     sinTheta cosTheta];

end
