DDF Fluid solver with Turbulence extensions
---------------------------------------------

Sample implementation for the paper 
 "Synthetic turbulence using Artificial Boundary layers"
 T.Pfaff, N.Thuerey, A.Selle, M.Gross
 ACM SIGGRAPH Asia 2009


 License
------------

The code is provided under the General Public License (GPL v2).
Please see the included file COPYRIGHT for more information.


 Installing
------------

You will need the GLUT and ZLIB libraries in order to build the project.
If they aren't installed in the normal include, lib paths, you'll need to
manually add their install path to the makefile.
For rendering, you also need to have the (open source) PBRT raytracer 
installed (see www.pbrt.org for details).

By default, the code uses OpenMP for some operations. If your compiler
doesn't handle OpenMP, disable it as indicated in the Makefile.


 Running
------------

The makefiles will build two executables, the command line tool turbCmd 
and the GUI interface turbGui.
They include two hard-coded scenes, a static one, similar to the fuel tank
scene in the paper, and a dynamic one, similar to the car scene in the
paper. Both executables take the following command line arguments:

turbCmd/turbGui [mode]
  where mode is
    scene       : generate scene files. has to be run first.
    pre-static  : precompute static scene
    pre-dynamic : precompute dynamic scene
    run-static  : run static scene
    run-dynamic : run dynamic scene


 GUI
-------------

The GUI can be used to visualize vortex particles, scalar and velocity grids.
It shows a 2D slice of the 3D volume grid.
Additional Information (which slice is displayed, display ranges etc. 
are shown on the console)

Basic Navigation:
- Press p to start/stop the simulation, l performs a single step
- Move camera position with WASD, zoom in using Q/E.
- Hold left mouse button to rotate camera, right mouse button to translate
- Select displayed slide with +/-
- Select x/y/z slice with *
- Press ESC to exit

Grid display:
- Vector grids (velocity etc.) are hidden by default. Press V to enable display.
- Shift+V toggles component-wise display and centered vector display of vector grids.
- Use [ ] keys to adjust the viewing range(scaling) of scalar grids, { } keys 
  for vector grids.
- Use X to select the next scalar grid, Z to select the next vector grid.
- The simulations include different "solvers", e.g. to handle the normal and the
  upscaled grids. Use ~ to switch between solvers.

Additional options:
- By default, all vortex particles are displayed. The Shift+G key toggles the display
  of vortex particles in the current slice only.


 Rendering
--------------

The program saves gzipped pbrt or df3 (povray) smoke density files, which can be used to 
render the smoke volumes. An example script for rendering using PBRT is included 
("render.sh" shell script). You need to edit the PBRT path in the script.
Use "./render.sh sta" or "./render.sh dyn" to render the static/dynamic scene.

Note that using  PBRT rendering of volume data is fairly slow, thus the render
process will take a while. Also, the rendering of animated meshes is not easily possible,
therefore the car mesh is omitted in the dynamic scene.


