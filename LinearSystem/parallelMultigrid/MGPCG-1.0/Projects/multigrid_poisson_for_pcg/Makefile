CC=g++


CFLAGS= -c
#if you want to use unoptimized version, uncomment this
#CFLAGS+=-DMGPCG_UNOPTIMIZED
#if you want to write opengl output, uncomment this
#CFLAGS+=-DMGPCG_WRITE_OPENGL_OUTPUT
LINK_LIBS=
LIBPATHS=

PhysBAM_build=../../PhysBAM-build
PhysBAM_libpath=$(PhysBAM_build)/lib

CFLAGS+=-I$(PhysBAM_build)/include/physbam
LIBPATHS+=-L$(PhysBAM_libpath) 
LINK_LIBS+=-lPhysBAM

SOURCES=main.cpp MULTIGRID_POISSON.cpp MULTIGRID_POISSON_REFINEMENT.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mgpcg_test

all:$(SOURCES) $(EXECUTABLE)

$(EXECUTABLE):$(OBJECTS)
	$(CC) $(LIBPATHS) $(LINK_LIBS) -Wl,-rpath=$(PhysBAM_libpath) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(OBJECTS) $(EXECUTABLE)