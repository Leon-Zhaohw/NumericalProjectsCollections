function [ u ] = Advec1D ( u, FinalTime )

%*****************************************************************************80
%
%% ADVEC1D integrates the 1D advection problem until FinalTime.
%
%  Licensing:
%
%    Permission to use this software for noncommercial
%    research and educational purposes is hereby granted
%    without fee.  Redistribution, sale, or incorporation
%    of this software into a commercial product is prohibited.
%
%    THE AUTHORS OR PUBLISHER DISCLAIMS ANY AND ALL WARRANTIES
%    WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
%    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
%    PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHORS OR
%    THE PUBLISHER BE LIABLE FOR ANY SPECIAL, INDIRECT OR
%    CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
%    RESULTING FROM LOSS OF USE, DATA OR PROFITS.
%
%  Modified:
%
%    23 September 2018
%
%  Author:
%
%    Original version by Jan Hesthaven, Tim Warburton.
%    Some modifications by John Burkardt.
%
%  Reference:
%
%    Jan Hesthaven, Tim Warburton,
%    Nodal Discontinuous Galerkin Methods: 
%    Algorithms, Analysis, and Applications,
%    Springer, 2007,
%    ISBN: 978-0387720654.
%
  Globals1D;

  time = 0.0;

  Elements = 10;
%
%  Runge-Kutta residual storage.
%
  resu = zeros ( Np, K );
%
%  Compute time step size.
%
  xmin = min ( abs ( x(1,:) - x(2,:) ) );
  CFL = 0.75;
  dt = CFL / ( 2 * pi ) * xmin;
  dt = 0.5 * dt;
  Nsteps = ceil ( FinalTime / dt );
  dt = FinalTime / Nsteps;
%
%  Advection speed.
%
  a = 2.0 * pi;
%
%  Outer time step loop.
%
  figure ( 1 );
  shapestr = { '-o','-x' };

  for tstep = 1 : Nsteps
    for INTRK = 1 : 5
      timelocal = time + rk4c ( INTRK ) * dt;
      [rhsu] = AdvecRHS1D ( u, timelocal, a );
      resu = rk4a ( INTRK ) * resu + dt * rhsu;
      u = u + rk4b ( INTRK ) * resu;
    end;
%
%  Increment time.
%
    time = time + dt;
%
%  Display solution on every 10th time step.
%
    if ( rem ( tstep, 10 ) == 0 )
      for i = 1 : Elements
        plot ( x(:,i), u(:,i), shapestr{1+rem(i,2)}, ...
          'Markersize', 8, 'LineWidth', 2 );
        hold all
      end
      grid ( 'on' );
      axis ( [ 0.0, 2.0, -3.0, 3.0 ] );
      drawnow;
      hold off
    end

  end

  return
end
