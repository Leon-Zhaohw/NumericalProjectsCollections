import numpy as np
import matplotlib.pyplot as plt
import h5py
from dedalus import public as de
from dedalus.extras import flow_tools
import time
from IPython import display
import logging
root = logging.root
for h in root.handlers:
    h.setLevel("INFO")
    
logger = logging.getLogger(__name__)

#Aspect ratio 2
Lx, Ly = (1., 1.)
nx, ny = (40, 40)

# Create bases and domain
x_basis = de.SinCos('x', nx, interval=(0, Lx), dealias=3/2)
y_basis = de.Chebyshev('y', ny, interval=(0, Ly), dealias=3/2)
domain = de.Domain([x_basis, y_basis], grid_dtype=np.float64)

problem = de.IVP(domain, variables=['p','u','v','w'])

Reynolds = 1e8
Schmidt = 1.
Diff = 0.00


problem.parameters['Re'] = Reynolds
problem.parameters['Sc'] = Schmidt
problem.parameters['Diff'] = Diff

problem.meta['u']['x']['parity'] = -1
problem.meta['v']['x']['parity'] = 1
problem.meta['p']['x']['parity'] = 1
problem.meta['w']['x']['parity'] = -1
#N-S equation 
problem.add_equation("dt(u) + dx(p) = w*v")
problem.add_equation("dt(v) + dy(p) = -w*u")
problem.add_equation("w = dx(v) - dy(u)")
#Incompressible constrain
problem.add_equation("dx(u) + dy(v) = 0")
#problem.add_equation("uy - dy(u) = 0")
#problem.add_equation("vy - dy(v) = 0")
'''
problem.add_bc("left(u) = 0", )
problem.add_bc("right(u) = 0")
problem.add_bc("left(vy) = 0", condition="(nx != 0)")
problem.add_bc("right(vy) = 0", condition="(nx != 0)")
#problem.add_bc("right(v) = 0", condition="(nx == 0)")
#problem.add_bc("right(p) = 0", condition="(nx == 0)")
problem.add_bc("integ(p,'y') = 0", condition="(nx == 0)")
problem.add_bc("integ(v,'y') = 0", condition="(nx == 0)")
'''
problem.add_bc("left(dy(v)) = 0", condition="(nx != 0)")
problem.add_bc("right(dy(v)) = 0", condition="(nx != 0)")
#problem.add_bc("integ(p,'y') = 0", condition="(nx == 0)")
#problem.add_bc("integ(v,'y') = 0", condition="(nx == 0)")
problem.add_bc("left(v) = 0", condition="(nx == 0)")
problem.add_bc("left(p) = 0", condition="(nx == 0)")

#problem.add_bc("left(u) = 0")
#problem.add_bc("right(u) = 0")
#problem.add_bc("left(v) = 0")
#problem.add_bc("right(v) = 0", condition="(nx != 0)")
#problem.add_bc("integ(p,'y') = 0", condition="(nx == 0)")

ts = de.timesteppers.RK443
solver =  problem.build_solver(ts)

x = domain.grid(0)
y = domain.grid(1)
u = solver.state['u']
v = solver.state['v']
#vy = solver.state['vy']
#uy = solver.state['uy']

#Output the init u and v
v['g'] = 0.02*np.cos(x*np.pi)*np.cos(y*np.pi)
#Dirichlet direction, sine, sine
u['g'] = 0.02*np.sin(x*np.pi)*np.sin(y*np.pi)
#v['g'] += 0.02*np.cos(4*x*np.pi)*np.cos(4*y*np.pi)
#u['g'] += 0.02*np.sin(4*x*np.pi)*np.sin(4*y*np.pi)
#aa = np.zeros(v['g'].shape)
#aa[30:40,30:40] = 1
#v['g'] = 0.2*aa
#u.differentiate('y',out=uy)
#v.differentiate('y',out=vy)

solver.stop_sim_time = 10.00
solver.stop_wall_time = np.inf
solver.stop_iteration = np.inf

initial_dt = 0.2*Lx/nx
cfl = flow_tools.CFL(solver,initial_dt,safety=0.8)
cfl.add_velocities(('u','v'))

analysis = solver.evaluator.add_file_handler('analysis_tasks', sim_dt=0.1, max_writes=50)
#analysis.add_task('s')
analysis.add_task('u')

# Make plot of scalar field
x = domain.grid(0,scales=domain.dealias)
y = domain.grid(1,scales=domain.dealias)
xm, ym = np.meshgrid(x,y)
fig, axis = plt.subplots(figsize=(6,6))
dt = 1e-6
solver.step(dt)
p = axis.pcolormesh(xm, ym, v['g'].T, cmap='RdBu_r');
axis.set_xlim([0,1.])
axis.set_ylim([0,1.])

logger.info('Starting loop')
start_time = time.time()
time_array = []

while solver.iteration <= 500:
    p.set_array(np.ravel(v['g'][:-1,:-1].T))
    fname = "./SineCosOutput2/%05d.png" % (solver.iteration - 1)
    plt.savefig(fname)
    dt = cfl.compute_dt()
    #dt = 1e-2
    solver.step(dt)
    time_array.append(dt)
    logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))

end_time = time.time()
time_array = np.asarray(time_array)
p.set_array(np.ravel(v['g'][:-1,:-1].T))
display.clear_output()
# Print statistics
logger.info('Run time: %f' %(end_time-start_time))
logger.info('Iterations: %i' %solver.iteration)