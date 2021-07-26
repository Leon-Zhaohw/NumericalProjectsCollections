C     = gcc
CC    = g++
FLAGS = -g -Wall 

# DEPENDENCIES
TNT     = ../TNT
TNT_INC = -I$(TNT)

INC = $(TNT_INC)
LIB = -lm 

# FILES
HEADERS = rply.h heap.h tml.h r3.h io.h edge.h facet.h vertex.h trimesh.h 

SOURCES = rply.c heap.cc quaternion.cc              \
          basic_io.cc ply_io.cc ifs_io.cc obj_io.cc \
          edge.cc facet.cc vertex.cc 

OBJECTS = rply.o heap.o quaternion.o            \
          basic_io.o ply_io.o ifs_io.o obj_io.o \
          edge.o facet.o vertex.o 

# OUTPUT
LIBNAME = libtml.a
LIB_DIR = ../lib
INC_DIR = ../include

# TO COMPILE

.c.o: $<
	$(C) -c $< -o $@ $(FLAGS) $(INC)

.cc.o: $<
	$(CC) -c $< -o $@ $(FLAGS) $(INC) 

$(LIBNAME): $(HEADERS) $(OBJECTS) 
	ar ruv $(LIBNAME) $(OBJECTS)
	ranlib $(LIBNAME)

install: $(LIBNAME)
	mv $(LIBNAME) $(LIB_DIR)
	cp $(HEADERS) $(INC_DIR)

clean:
	rm -f *.o *~ $(LIBNAME)
