n = 12;
polcfg_R3 = [
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
-3.0 0.0 -1.8
-2.59807621135332 1.5 -1.8
-1.5 2.59807621135332 -1.8
-1.83690953073357e-16 3.0 -1.8
1.5 2.59807621135332 -1.8
2.59807621135332 1.5 -1.8
3.0 3.67381906146713e-16 -1.8
2.59807621135332 -1.5 -1.8
1.5 -2.59807621135332 -1.8
5.5107285922007e-16 -3.0 -1.8
-1.5 -2.59807621135332 -1.8
-2.59807621135332 -1.5 -1.8
-5.0 0.0 -5.0
-4.33012701892219 2.5 -5.0
-2.5 4.33012701892219 -5.0
-3.06151588455594e-16 5.0 -5.0
2.5 4.33012701892219 -5.0
4.33012701892219 2.5 -5.0
5.0 6.12303176911189e-16 -5.0
4.33012701892219 -2.5 -5.0
2.5 -4.33012701892219 -5.0
9.18454765366783e-16 -5.0 -5.0
-2.5 -4.33012701892219 -5.0
-4.33012701892219 -2.5 -5.0
-7.0 0.0 -9.8
-6.06217782649107 3.5 -9.8
-3.5 6.06217782649107 -9.8
-4.28612223837832e-16 7.0 -9.8
3.5 6.06217782649107 -9.8
6.06217782649107 3.5 -9.8
7.0 8.57224447675664e-16 -9.8
6.06217782649107 -3.5 -9.8
3.5 -6.06217782649107 -9.8
1.2858366715135e-15 -7.0 -9.8
-3.5 -6.06217782649107 -9.8
-6.06217782649107 -3.5 -9.8
];

polcfg_x = polcfg_R3(:,1); polcfg_x = polcfg_x([1:n; n+1:2*n; 2*n+1:3*n; 3*n+1:4*n]);
polcfg_y = polcfg_R3(:,2); polcfg_y = polcfg_y([1:n; n+1:2*n; 2*n+1:3*n; 3*n+1:4*n]);
polcfg_z = polcfg_R3(:,3); polcfg_z = polcfg_z([1:n; n+1:2*n; 2*n+1:3*n; 3*n+1:4*n]);
b33x = c1_polar_spline(n, polcfg_x);
b33y = c1_polar_spline(n, polcfg_y);
b33z = c1_polar_spline(n, polcfg_z);
b36x = c2_polar_spline(n, polcfg_x);
b36y = c2_polar_spline(n, polcfg_y);
b36z = c2_polar_spline(n, polcfg_z);

figure(1); hold on;
draw_periodic_mesh(n, polcfg_x, polcfg_y, polcfg_z); % polar configuration
%draw_periodic_mesh(n, b33x, b33y, b33z); % b-spline control net
draw_polar_spline33(n, b33x, b33y, b33z);
axis equal;
rotate3d on;
hold off;

figure(2); hold on;
draw_periodic_mesh(n, polcfg_x, polcfg_y, polcfg_z); % polar configuration
%draw_periodic_mesh(4*n, b36x, b36y, b36z); % b-spline control net
draw_polar_spline36(n, b36x, b36y, b36z);
axis equal;
rotate3d on;
hold off;

