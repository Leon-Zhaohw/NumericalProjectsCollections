CXX := g++
CXXFLAGS += -O3 -Wall -g -Wunused-variable
LDFLAGS = -g -fopenmp

TARGETS = targetrdf
SOURCES = $(wildcard *.cc)

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
CXXFLAGS += -I/opt/local/include
LDFLAGS += -L/opt/local/lib
endif

ifeq ($(UNAME), Linux)
CXXFLAGS += -frounding-math
endif

OBJS = $(patsubst %.cc, %.o, $(SOURCES))
DEPS = $(patsubst %.cc, .deps/%.P, $(SOURCES))

.PHONY: all clean install
all: $(TARGETS)

clean:
	rm -f $(TARGETS) $(OBJS) $(DEPS)

targetrdf: targetrdf.o common.o param.o
	$(CXX) $(LDFLAGS) -o $@ $^

%.o : %.cc
	@mkdir -p $(dir .deps/$*)
	$(CXX) -MD -Wp,-MD,.deps/$*.P -MT '$@' -MP -c $(CXXFLAGS) -o $@ $<

-include $(DEPS)
