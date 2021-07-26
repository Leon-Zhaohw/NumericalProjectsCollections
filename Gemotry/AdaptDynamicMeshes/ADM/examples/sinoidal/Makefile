include ../ADMmake.defs

DIR  = sinoidal

INC = $(VIEWER) $(SURFACES) $(ADM_INC) $(FLTK_INC) $(MAGICK_INC) 
LIB = $(ADM_LIB) $(FLTK_LIB) $(MAGICK_LIB)

HDRS = sinoidal.h
SRCS = main.cc
OBJS = main.o
BIN  = sin

all: $(OBJS)
	$(CC) $(OBJS) -o $(BIN) $(FLAGS) $(INC) $(LIB)

clean:
	rm -f *~ *.o $(BIN) ../*~
