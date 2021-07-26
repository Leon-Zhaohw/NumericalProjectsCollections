include ../ADMmake.defs

DIR  = sinc

INC = $(VIEWER) $(SURFACES) $(ADM_INC) $(FLTK_INC) $(MAGICK_INC)
LIB = $(ADM_LIB) $(FLTK_LIB) $(MAGICK_LIB)

HDRS = sinc.h
SRCS = bump.cc height.cc
OBJS = bump.o  height.o
BIN  = bump height

all: $(BIN)

bump: bump.o
	$(CC) $< -o $@ $(FLAGS) $(INC) $(LIB)

height: height.o
	$(CC) $< -o $@ $(FLAGS) $(INC) $(LIB)

clean:
	rm -f *~ *.o $(BIN) ../*~
