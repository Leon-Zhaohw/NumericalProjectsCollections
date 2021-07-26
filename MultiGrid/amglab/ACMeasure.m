function u = ACMeasure(M,B,gamma) 

% Divya Bansal
% Department of Computational Sciences
% University of Kentucky


u = (log(normest(B) * M) - log10((normest(B) * M) - 1)) * (1 - gamma);