LDFLAGS_COMMON = -framework Accelerate -framework GLUT -framework OpenGL -L/usr/local/lib -lstdc++ -lpng -lz -ljpeg -msse2 -lfftw3 -lm
CFLAGS_COMMON = -v -c -Wall -Wno-deprecated -I../../src/integrators -I../../src/linearalgebra -I../../src/geometry -I../../src/util -I../../src/glvu -DDO_PNG_OUT=0 -I./ -I/usr/local/include -I../../ -I../../src/Eigen/ -O3 -DNO_FFT -fopenmp -msse2 -lstdc++ -I/opt/local/include/

