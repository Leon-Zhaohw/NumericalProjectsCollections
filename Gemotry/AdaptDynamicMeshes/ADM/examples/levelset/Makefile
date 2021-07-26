include ../ADMmake.defs

DIR  = levelset

CPLS3D = ../../CPLS
CPLS3D_INC = -I$(CPLS3D)/include/3D
CPLS3D_LIB = -L$(CPLS3D)/lib -lCPLS3D

INC = $(VIEWER) $(SURFACES) $(CPLS3D_INC) $(ADM_INC) $(FLTK_INC) $(MAGICK_INC) 
LIB = $(CPLS3D_LIB) $(ADM_LIB) $(FLTK_LIB) $(MAGICK_LIB)

HDRS = 
SRCS = transform.cc jiao.cc enright.cc twospheres.cc
OBJS = transform.o  jiao.o  enright.o  twospheres.o
BIN  = transform    jiao    enright    twospheres

all: $(BIN)

transform: 
	$(CC) -c transform.cc -o transform.o $(FLAGS) $(TML_INC)
	$(CC) transform.o -o $@ $(FLAGS) $(TML_INC) $(TML_LIB)

jiao: jiao.o
	$(CC) $< -o $@ $(FLAGS) $(INC) $(LIB)

enright: enright.o
	$(CC) $< -o $@ $(FLAGS) $(INC) $(LIB)

twospheres: twospheres.o
	$(CC) $< -o $@ $(FLAGS) $(INC) $(LIB)

clean:
	rm -f *~ *.o $(BIN) ../*~
