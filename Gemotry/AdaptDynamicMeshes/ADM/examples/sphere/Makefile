include ../ADMmake.defs

DIR  = sphere

INC = $(VIEWER) $(SURFACES) $(ADM_INC) $(FLTK_INC) $(MAGICK_INC)
LIB = $(ADM_LIB) $(FLTK_LIB) $(MAGICK_LIB)

HDRS = sphere.h multi_criteria.h
SRCS = multi.cc unit.cc
OBJS = multi.o  unit.o
BIN  = multi unit

all: $(BIN)

multi: multi.o
	$(CC) $< -o $@ $(FLAGS) $(INC) $(LIB)

unit: unit.o
	$(CC) $< -o $@ $(FLAGS) $(INC) $(LIB)

clean:
	rm -f *~ *.o $(BIN) ../*~
