Dynamic Kelvinlets: Secondary Motions based on Fundamental Solutions of Elastodynamics
ACM SIGGRAPH 2018

++++++++++++++++
+ DEPENDENCIES +
++++++++++++++++

* EIGEN -- eigen.tuxfamily.org
* [OPTIONAL] TBB -- www.threadingbuildingblocks.org

+++++++++
+ BUILD +
+++++++++

mkdir build; 
cd build; 
cmake -DEIGEN3_INCLUDE_DIR=<PATH_TO_EIGEN> ..; 
make;

# With TBB
cmake -DEIGEN3_INCLUDE_DIR=<PATH_TO_EIGEN> -DUSE_TBB=ON ..;

+++++++
+ RUN +
+++++++

# Help
./kelvinlet -h

# Default: 60 steps of Impulse-Grab applied to a grid.
./kelvinlet

# Example: 10 steps of Push-Twist applied to a grid.
./kelvinlet -a 1 -t 1 -s 0.1 -p

+++++++++++++++++
+ VISUALIZATION +
+++++++++++++++++

The output is a list of OBJs of a deforming grid.
You can visualize these files using the Houdini file viewer.hipnc.
Set the prefix path to the output files in the "currentFrame" node.
Download Houdini from: www.sidefx.com/products/houdini-apprentice

++++++++++++++++++
+ CODE STRUCTURE +
++++++++++++++++++

main.cpp -- generates a grid, deforms it using Dynamic Kelvinlets, 
            and writes each frame to an OBJ file.

kelvinlet/ 
* types           -- Eigen typedefs
* dynaBase        -- abstract class with Dynamic Kelvinlet API.
* dynaPulseBase   -- abstract class with Impulse Potentials.
* dynaPulseGrab   -- Dynamic Kelvinlet with Impulse-Grab deformation.
* dynaPulseAffine -- Dynamic Kelvinlet with Impulse-Affine deformation.
* dynaPushBase    -- abstract class with Push Potentials.
* dynaPushGrab    -- Dynamic Kelvinlet with Push-Grab deformation.
* dynaPushAffine  -- Dynamic Kelvinlet with Push-Affine deformation.

Fernando de Goes and Doug L. James
Apr/2018
