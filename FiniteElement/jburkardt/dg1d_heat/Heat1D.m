function [ u, time ] = Heat1D ( u, FinalTime )

%*****************************************************************************80
%
%% HEAT1D integrates the 1D heat equation up to a given final time.
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
%    17 September 2018
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
%  Runge-Kutta residual storage  
%
  resu = zeros(Np, K); 
%
%  Compute time step size.
%
  xmin = min ( abs(x(1,:)-x(2,:)) );
  CFL = 0.25;
  dt = CFL * (xmin)^2;
  Nsteps = ceil ( FinalTime / dt );
  dt = FinalTime / Nsteps;
%
%  Outer time step loop 
%
  for tstep=1:Nsteps
    for INTRK = 1:5

      timelocal = time + rk4c(INTRK)*dt;        
%
%  Compute right hand side of 1D advection equations.
%
      [rhsu] = HeatCRHS1D(u,timelocal);
%
%  Initiate and increment Runge-Kutta residuals
%
      resu = rk4a(INTRK) * resu + dt * rhsu;  
%    
%  Update fields
%
      u = u + rk4b(INTRK) * resu;

    end;
%
%  Increment time.
%
    time = time+dt;
  end

  return
end
