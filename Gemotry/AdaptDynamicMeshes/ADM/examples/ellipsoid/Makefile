include ../ADMmake.defs

DIR  = ellipsoid

INC = $(VIEWER) $(SURFACES) $(ADM_INC) $(FLTK_INC) $(MAGICK_INC)
LIB = $(ADM_LIB) $(FLTK_LIB) $(MAGICK_LIB)

HDRS = ellipsoid.h
SRCS = main.cc
OBJS = main.o
BIN  = ellipsoid

all: $(OBJS)
	$(CC) $(OBJS) -o $(BIN) $(FLAGS) $(INC) $(LIB)

clean:
	rm -f *~ *.o $(BIN) ../*~
