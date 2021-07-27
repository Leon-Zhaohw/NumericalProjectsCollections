#include ../include.mk

# calls:
CC         = g++
LDFLAGS = -framework GLUT -framework OpenGL -L/usr/local/lib/x86_64 -fopenmp -lz
CFLAGS = -c -Wall -I./src -I./src/glvu -I./ -fopenmp -O3 -DNDEBUG
HOUDINI = iWaveHoudini
PHYSBAM = iWavePhysBAM

SOURCES    = ./src/IWAVE_3D.cpp \
						 ./src/FIELD_2D.cpp \
						 ./src/FIELD_3D.cpp \
						 ./src/VECTOR3_FIELD_3D.cpp \
						 ./src/BOX.cpp \
						 ./src/TRIANGLE.cpp \
						 ./src/MATRIX3.cpp \
						 ./src/glvu/glvu.cpp \
						 ./src/TIMER.cpp \
						 ./src/TRIANGLE_MESH.cpp
HOUDINI_OBJECTS    = iWaveHoudini.o $(SOURCES:.cpp=.o)
PHYSBAM_OBJECTS    = iWavePhysBAM.o $(SOURCES:.cpp=.o)

all: $(SOURCES) $(HOUDINI) $(PHYSBAM)
	
$(HOUDINI): $(HOUDINI_OBJECTS) 
	$(CC) $(HOUDINI_OBJECTS) $(LDFLAGS) -o $@

$(PHYSBAM): $(PHYSBAM_OBJECTS) 
	$(CC) $(PHYSBAM_OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o;rm -f ./src/*.o;rm -f ./src/glvu/*.o;rm $(HOUDINI); rm $(PHYSBAM)
