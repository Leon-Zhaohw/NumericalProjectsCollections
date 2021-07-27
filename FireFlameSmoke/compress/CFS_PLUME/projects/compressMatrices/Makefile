include ../include.mk

# calls:
CC         = g++
CFLAGS     = ${CFLAGS_COMMON} -I../../src/cubature -I../../src/Eigen/ -I../../src/compression
LDFLAGS    = ${LDFLAGS_COMMON}
EXECUTABLE = ../../bin/compressMatrices

SOURCES    = compressMatrices.cpp \
						 ../../src/glvu/glvu.cpp \
						 ../../src/util/TIMER.cpp \
						 ../../src/util/SIMPLE_PARSER.cpp \
						 ../../src/linearalgebra/SPARSE_MATRIX.cpp \
						 ../../src/linearalgebra/SPARSE_MATRIX_ARRAY.cpp \
						 ../../src/linearalgebra/VECTOR.cpp \
						 ../../src/linearalgebra/MATRIX.cpp \
						 ../../src/linearalgebra/BIG_MATRIX.cpp \
						 ../../src/linearalgebra/MATRIX3.cpp \
						 ../../src/linearalgebra/NNLS.cpp \
						 ../../src/linearalgebra/EIGEN.cpp \
						 ../../src/cubature/CUBATURE_GENERATOR_EIGEN.cpp \
						 ../../src/geometry/FIELD_3D.cpp \
						 ../../src/geometry/VECTOR3_FIELD_3D.cpp \
						 ../../src/geometry/BOX.cpp \
						 ../../src/integrators/SUBSPACE_FLUID_3D_EIGEN.cpp \
						 ../../src/integrators/FLUID_3D_MIC.cpp \
						 ../../src/geometry/FIELD_2D.cpp \
						 ../../src/geometry/VECTOR3_FIELD_2D.cpp \
						 ../../src/util/MIN_HEAP.cpp \
						 ../../src/compression/COMPRESSION.cpp \
						 ../../src/compression/COMPRESSION_DATA.cpp \
						 ../../src/compression/MATRIX_COMPRESSION_DATA.cpp \
						 ../../src/geometry/INTEGER_FIELD_3D.cpp

OBJECTS    = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
	rm -f ../../bin/compressMatrices
