#LDFLAGS_COMMON = -framework Accelerate -framework GLUT -framework OpenGL -L/usr/local/lib/x86_64 -L/opt/local/lib -lstdc++ -lpng -lz -ljpeg -fopenmp -msse2 -lmpfr -lgmp
LDFLAGS_COMMON = -L/opt/local/lib -lstdc++ -lpng -lz -ljpeg
MY_INCLUDES = -I../../src/integrators -I../../src/algebra -I../../src/geometry -I../../src/util -I./ -I/opt/local/include -I../../ -I./eigen/Eigen -I../../src/solvers/
#CFLAGS_COMMON = -c -Wall -I../../src/solvers -I../../src/TOUPEE -I../../src/integrators -I../../src/algebra -I../../src/geometry -I../../src/util -I../../src/solvers -I./ -I/opt/local/include -I../../ -I../../src/Eigen/ -O3 -msse2 -fopenmp
CFLAGS_COMMON = -c -Wall -I../../src/solvers -I../../src/TOUPEE -I../../src/integrators -I../../src/algebra -I../../src/geometry -I../../src/util -I../../src/solvers -I./ -I/opt/local/include -I../../ -I./eigen/ -O3 -msse2 -std=c++11
