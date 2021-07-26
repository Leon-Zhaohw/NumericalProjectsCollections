function  CC = crossCorrAbs(f,g)

% A function for calculating the cross correlation function for functions f and g.
% The functions are describes in N points, with sampling frequency fs.
% See details in a chapter in the textbook
% The Physics of Oscillations and Waves, Arnt Inge Vistnes, 2018.

N = length(f);
sf2 = sum(f.*f);
sg2 = sum(g.*g);
normConstant = sqrt(sf2*sg2);
CC = abs(sum((f.*g))/normConstant);

end

