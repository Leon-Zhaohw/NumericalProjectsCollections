function DDCheck(at_R)
%Plot the discretization
y=[-1 0 1];
for i=1:length(at_R)-1
    x=ones(1,3)*at_R(i);
    plot(x,y,'k')
    hold on
end