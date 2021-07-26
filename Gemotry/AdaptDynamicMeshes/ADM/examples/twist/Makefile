include ../ADMmake.defs

DIR  = twist

INC = $(VIEWER) $(SURFACES) $(ADM_INC) $(FLTK_INC) $(MAGICK_INC)
LIB = $(ADM_LIB) $(FLTK_LIB) $(MAGICK_LIB)

HDRS = twist.h
SRCS = twist.cc main.cc
OBJS = twist.o  main.o
BIN  = twist

all: $(OBJS)
	$(CC) $(OBJS) -o $(BIN) $(FLAGS) $(INC) $(LIB)

clean:
	rm -f *~ *.o $(BIN) ../*~
