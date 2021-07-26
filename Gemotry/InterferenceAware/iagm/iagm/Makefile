
# Things for users to define
#
include make.inc


IAGM_LIB := libiagm.$(LIB_EXT)
all: $(IAGM_LIB)

CXXFLAGS := -Wall -Wextra -pipe -fPIC -shared -O3

INCS := -Iinc -Isrc -I$(EIGEN_INC)
LIBS := 

VPATH := ./src
SRCS := $(wildcard ./src/*.cpp)

BUILD_DIR := build/
OBJS := $(addprefix $(BUILD_DIR), $(addsuffix .o, $(basename $(notdir $(SRCS)))))


$(IAGM_LIB): $(OBJS)
	mkdir -p lib
	$(CXX) -dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-compatibility_version,1.0,-current_version,1.0,-install_name,./lib/$(IAGM_LIB) -o lib/$(IAGM_LIB) $(OBJS) $(LIBS)

$(BUILD_DIR)%.o: %.cpp
	$(shell mkdir -p $(BUILD_DIR))
	$(CXX) $(CXXFLAGS) $(INCS) -c -o $@ $<

depend:
	$(CXX) $(INCS) -MM $(SRCS) | sed 's/.*\.o/build\/&/' > .depend

ifeq (.depend, $(wildcard .depend))
  include .depend
endif

clean: 
	-rm -f lib/$(IAGM_LIB) $(BUILD_DIR)*.o .depend

rebuild: clean depend all

