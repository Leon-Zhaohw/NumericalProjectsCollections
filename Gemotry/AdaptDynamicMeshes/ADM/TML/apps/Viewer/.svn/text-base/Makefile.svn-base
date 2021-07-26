include ../TMLmake.defs

MAGICK_INC = `Magick++-config --cxxflags --cppflags`
MAGICK_LIB = `Magick++-config --ldflags --libs`

FLTK_INC = `fltk-config --use-forms --use-gl --cflags --cxxflags`
FLTK_LIB = `fltk-config --use-forms --use-gl --ldflags --libs`

INC   = $(TML_INC) $(FLTK_INC) #$(MAGICK_INC)
LIB   = $(TML_LIB) $(FLTK_LIB) #$(MAGICK_LIB)

# FILES
OBJS = 
BIN  = viewer viewer2

# TO COMPILE
CC    = g++
FLAGS = -g -Wall

.cc.o: $<
	$(CC) -c $< -o $@ $(FLAGS) $(INC)

all: $(BIN)

#-----------------------------------------------------#

OBJS_V1 = viewer.o
viewer: $(OBJS_V1)
	$(CC) -o $@ $(OBJS_V1) $(FLAGS) $(INC) $(LIB)

OBJS_V2 = TMViewer.o Drawing.o viewer2.o
viewer2: $(OBJS_V2)
	$(CC) -o $@ $(OBJS_V2) $(FLAGS) $(INC) $(LIB)

#-----------------------------------------------------#

clean:
	rm -rf *.o *~ $(BIN)
