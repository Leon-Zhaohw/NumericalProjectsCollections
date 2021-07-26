include ../ADMmake.defs

DIR  = extrudelp

INC = $(VIEWER) $(SURFACES) $(ADM_INC) $(FLTK_INC) $(MAGICK_INC)
LIB = $(ADM_LIB) $(FLTK_LIB) $(MAGICK_LIB)

HDRS = extrudelp.h
SRCS = main.cc
OBJS = main.o
BIN  = extrudelp

all: $(OBJS)
	$(CC) $(OBJS) -o $(BIN) $(FLAGS) $(INC) $(LIB)

clean:
	rm -f *~ *.o $(BIN) ../*~
