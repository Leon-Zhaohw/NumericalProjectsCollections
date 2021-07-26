Sharp Kelvinlets: Elastic Deformations with Cusps and Localized Falloffs
Fernando de Goes and Doug L. James
DigiPro 2019

++++++++++++++++
+ DEPENDENCIES +
++++++++++++++++

* EIGEN -- eigen.tuxfamily.org
* [OPTIONAL] TBB -- www.threadingbuildingblocks.org

+++++++++
+ BUILD +
+++++++++

cd exampleBrush; 
[or] 
cd exampleDyna;

mkdir build; 
cd build; 

cmake -DEIGEN3_INCLUDE_DIR=<PATH_TO_EIGEN> ..;  
[or] 
cmake -DEIGEN3_INCLUDE_DIR=<PATH_TO_EIGEN> -DUSE_TBB=ON ..;

make;

+++++++
+ RUN +
+++++++

# Help
./kelvinlet -h

++++++++++++++++
+ exampleBrush +
++++++++++++++++

# Default: Grab applied to a grid.
./kelvinlet

# Example: Grab-Cusp applied to a grid.
./kelvinlet -a 5 

+++++++++++++++
+ exampleDyna +
+++++++++++++++

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

Fernando de Goes and Doug L. James
June/2019
