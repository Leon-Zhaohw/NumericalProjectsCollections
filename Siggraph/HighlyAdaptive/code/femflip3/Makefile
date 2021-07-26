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

NAME := flip3d
DEBUG := 0
BINDIR := ./bin
SRCDIR := ./src
EXTHEAD := /opt/local/include
EXTLIB := /opt/local/lib

SRC := $(shell find $(SRCDIR) -name "*.cpp")
HEADER := $(shell find $(SRCDIR) -name "*.h")
INC	:= $(addprefix -I , $(shell find $(SRCDIR) -name ".svn" -prune -o -type d)) -I$(EXTHEAD)
ARCH := $(shell uname)
SUM	= $(shell find . -name "*.h" -print0 -o -name "*.cpp" -print0 | xargs -0 wc -l | grep total )

ifeq ($(ARCH),Darwin)
	LIB := -framework OpenGL -framework GLUT -L$(EXTLIB) $(SRCDIR)/tetgen/libtet.a -lgsl -lgslcblas -lz
endif

ifeq ($(ARCH),Linux)
	LIB := -lglut -lGLEW -fopenmp $(SRCDIR)/tetgen/libtet.a -lgsl -lgslcblas -lz
endif

all:
	@make depend
	@make $(BINDIR)/$(NAME) DEP=1
	@echo "\nLine sum: ${SUM}"

error:
	@make 2>&1 | grep 'error'

$(BINDIR)/$(NAME) : $(SRC:.cpp=.o) $(SRCDIR)/tetgen/libtet.a $(SRCDIR)/rply/rply.o
	@mkdir -p $(BINDIR)
	@$(CXX) -o $@ $(SRC:.cpp=.o) $(SRCDIR)/rply/rply.o $(LIB) $(OPT)
	@echo "Linking $@..."

$(SRCDIR)/tetgen/libtet.a: 
	@cd $(SRCDIR)/tetgen; make tetlib

run: all
	$(BINDIR)/$(NAME)

depend : $(SRC:.cpp=.d)
%.d : %.cpp
	@echo "$@ $(dir $@)$$(gcc $(INC) -MM -MG $<)" > $@

%.o : %.cpp
	@echo "Compiling $< ..."
	@$(CXX) $(OPT) $(INC) -o $@ -c $<

%.o : %.c
	@echo "Compiling $< ..."
	@$(CC) $(OPT) $(INC) -o $@ -c $<

.PHONY: clean
clean:
	rm -rf render
	rm -rf $(SRCDIR)/tetgen/libtet.a
	rm -rf $(shell find . -name "*.o" -o -name "*.d")
	rm -rf $(BINDIR)
	rm -rf GPATH GRTAGS GSYMS GTAGS HTML

ifeq ($(DEP), 1)
-include $(SRC:.cpp=.d)
endif
