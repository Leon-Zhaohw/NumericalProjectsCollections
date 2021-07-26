CC    = g++
FLAGS = -g -Wall # -DCHECKER

# DEPENDENCIES
TNT     = ../TML/TNT
TNT_INC = -I$(TNT)

TML     = ../TML
TML_INC = -I$(TML)/include $(TNT_INC)
TML_LIB = -L$(TML)/lib -ltml -lm

INC = $(TML_INC)
LIB = $(TML_LIB)

# FILES
HEADERS = adm.h edge.h facet.h vertex.h oracle.h criteria.h mesh.h

SOURCES = a48elmts.cc stochastic.cc adapt.cc front.cc stellar.cc triquad.cc mesh.cc color.cc

OBJECTS = a48elmts.o stochastic.o adapt.o front.o stellar.o triquad.o mesh.o color.o

# OUTPUT
LIBNAME = libadm.a
LIB_DIR = ../lib
INC_DIR = ../include

# TO COMPILE

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
