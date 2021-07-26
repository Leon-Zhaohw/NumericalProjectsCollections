include ../ADMmake.defs

DIR  = cubelp

INC = $(VIEWER) $(SURFACES) $(ADM_INC) $(FLTK_INC) #$(MAGICK_INC)
LIB = $(ADM_LIB) $(FLTK_LIB) #$(MAGICK_LIB)

HDRS = lp3d.h multi_criteria.h
SRCS = main.cc multi.cc
OBJS = main.o multi.o
BIN  = cubelp multi

all: $(BIN)

cubelp: main.o
	$(CC) main.o -o cubelp $(FLAGS) $(INC) $(LIB)

multi: multi.o
	$(CC) $< -o $@ $(FLAGS) $(INC) $(LIB)

clean:
	rm -f *~ *.o $(BIN) ../*~
