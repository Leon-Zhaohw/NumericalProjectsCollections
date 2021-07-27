include ../include.mk

# calls:
CC         = g++
CFLAGS     = ${CFLAGS_COMMON}
LDFLAGS    = ${LDFLAGS_COMMON}
EXECUTABLE = ../../bin/buildProducts

SOURCES    = buildProducts.cpp \
						 ../../src/glvu/glvu.cpp \
						 ../../src/util/TIMER.cpp \
						 ../../src/util/SIMPLE_PARSER.cpp \
						 ../../src/linearalgebra/SPARSE_MATRIX.cpp \
						 ../../src/linearalgebra/SPARSE_MATRIX_ARRAY.cpp \
						 ../../src/linearalgebra/VECTOR.cpp \
						 ../../src/linearalgebra/MATRIX.cpp \
						 ../../src/linearalgebra/MATRIX3.cpp \
						 ../../src/linearalgebra/EIGEN.cpp \
						 ../../src/geometry/FIELD_3D.cpp \
						 ../../src/geometry/VECTOR3_FIELD_3D.cpp \
						 ../../src/integrators/SUBSPACE_FLUID_3D_EIGEN.cpp \
						 ../../src/integrators/FLUID_3D_MIC.cpp \
						 ../../src/geometry/BOX.cpp \
						 ../../src/geometry/FIELD_2D.cpp \
						 ../../src/geometry/VECTOR3_FIELD_2D.cpp \
						 ../../src/util/MIN_HEAP.cpp \

OBJECTS    = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
	rm -f ../../bin/buildProducts
