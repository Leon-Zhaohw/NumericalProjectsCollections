function [ newVelocityField ] = applyForces( velocityField, forceField, deltaTime )
%APPLYFORCES Apply force to the given velocityField
%   return newVelocityField: matrix(yNodes, xNodes, dim)
%   velocity: matrix(yNodes, xNodes, dim)
%   forceField: matrix(yNodes, xNodes, dimension), force applied to
%       fluid in m/s^2
%   deltaTime: float

    newVelocityField = velocityField + forceField * deltaTime;
    
end

