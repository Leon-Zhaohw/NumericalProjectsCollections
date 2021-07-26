function [ E, H ] = Maxwell1D ( E, H, eps, mu, FinalTime );

%*****************************************************************************80
%
%% MAXWELL1D integrates 1D Maxwell’s equations.
%
%  Discussion:
%
%    Integration is from time 0 until FinalTime, starting with
%    conditions (E(t=0),H(t=0)) and materials (eps,mu).
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
%    24 September 2018
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
%
%  Runge-Kutta residual storage.
%
  resE = zeros(Np,K);
  resH = zeros(Np,K);
%
%  Compute time step size.
% 
  xmin = min ( abs ( x(1,:) - x(2,:) ) );
  CFL = 1.0;
  dt = CFL * xmin;
  Nsteps = ceil ( FinalTime / dt );
  dt = FinalTime / Nsteps;
%
%  Outer time step loop.
%
  for tstep = 1 : Nsteps
    for INTRK = 1 : 5
      [ rhsE, rhsH ] = MaxwellRHS1D ( E, H, eps, mu );
      resE = rk4a(INTRK) * resE + dt * rhsE;
      resH = rk4a(INTRK) * resH + dt * rhsH;
      E = E + rk4b(INTRK) * resE;
      H = H + rk4b(INTRK) * resH;
    end
    time = time + dt;
  end

  return
end
