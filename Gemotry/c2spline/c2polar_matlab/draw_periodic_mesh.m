% Helper function to display polar configuration.

function draw_periodic_mesh(n, polcfg_x, polcfg_y, polcfg_z)

CIRC=[1:n,1];
for i=1:size(polcfg_x,1)
    plot3(polcfg_x(i,CIRC), polcfg_y(i,CIRC), polcfg_z(i,CIRC));
end
for i=1:n
    plot3(polcfg_x(:,i), polcfg_y(:,i), polcfg_z(:,i));
end

