#
#	Makefile
#	
#	Created by Ryoichi Ando
#	Last Updated on 1st May 2013
#	Email: and@verygood.aid.design.kyushu-u.ac.jp
#

RELEASE := 1
OPT := -Wall
ifeq ($(RELEASE),1)
OPT += -Os
else
OPT += -g
endif

NAME := flip2d
DEBUG := 0
BINDIR := ./bin
SRCDIR := ./src

SRC := $(shell find $(SRCDIR) -name "*.cpp")
HEADER := $(shell find $(SRCDIR) -name "*.h")
INC	:= $(addprefix -I , $(shell find $(SRCDIR) -name ".svn" -prune -o -type d))
ARCH := $(shell uname)
SUM	= $(shell find . -name "*.h" -print0 -o -name "*.cpp" -print0 | xargs -0 wc -l | grep total )

ifeq ($(ARCH),Darwin)
	LIB := -framework OpenGL -framework GLUT
endif

ifeq ($(ARCH),Linux)
	LIB := -lglut -lGLEW
endif

all:
	@make depend
	@make $(BINDIR)/$(NAME) DEP=1
	@echo "\nLine sum: ${SUM}"

$(BINDIR)/$(NAME) : $(SRC:.cpp=.o)
	@mkdir -p $(BINDIR)
	@$(CXX) -o $@ $(SRC:.cpp=.o) $(LIB) $(OPT)
	@echo "Linking $@..."

run: all
	$(BINDIR)/$(NAME)

depend : $(SRC:.cpp=.d)
%.d : %.cpp
	@echo "$@ $(dir $@)$$(gcc $(INC) -MM -MG $<)" > $@

%.o : %.cpp
	@echo "Compiling $< ..."
	@$(CXX) $(OPT) $(INC) -o $@ -c $<

.PHONY: clean
clean:
	rm -rf $(shell find . -name "*.o" -o -name "*.d")
	rm -rf $(BINDIR)
	rm -rf GPATH GRTAGS GSYMS GTAGS HTML

ifeq ($(DEP), 1)
-include $(SRC:.cpp=.d)
endif
