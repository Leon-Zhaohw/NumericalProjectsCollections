
# User provide definitions
#
include make.inc


TARGETS := spline2d
all: $(TARGETS)


OBJ_DIR = build/

VPATH := ./src

SOURCES  = $(wildcard src/*.cpp)
OBJECTS  = $(addprefix $(OBJ_DIR), $(addsuffix .o, $(basename $(notdir $(SOURCES)))))

LDFLAGS := $(GL_LIB)

CPPFLAGS := -I$(ATB_INC) -I$(EIGEN_INC) -I$(IAGM_INC) -Isrc -g -O3 -Wall -Woverloaded-virtual -Wuninitialized -Wimplicit
LDLIBS   := $(ATB_LIB) $(IAGM_LIB) 

$(OBJ_DIR)%.o: %.cpp
	$(shell mkdir -p $(OBJ_DIR))
	$(CXX) $(CPPFLAGS) -c -o $@ $<

spline2d: $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@ $(CPPFLAGS) $(LDLIBS)

clean: 
	-rm -f $(TARGETS)
	-rm -rf $(OBJECTS)
	-rm -f .depend

depend:
	$(CXX) -MM $(CPPFLAGS) $(SOURCES) | sed 's/.*\.o/build\/&/' > .depend

rebuild: clean depend all

ifeq (.depend, $(wildcard .depend))
  include .depend
endif

