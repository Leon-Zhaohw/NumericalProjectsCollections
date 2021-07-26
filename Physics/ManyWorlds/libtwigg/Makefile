#
#  Makefile for fltk applications
#

#CXXFLAGS=-Wall  -g -O0 -fno-inline -pthread -DUSE_GSL
CXXFLAGS=-Wall  -fno-math-errno -fno-trapping-math -fno-signaling-nans -Os -march=pentium4 -pthread -DNDEBUG -DUSE_GSL
CXXFLAGS_DEBUG=-Wall  -O0 -fno-inline -pthread -DUSE_GSL
CFLAGS=-Wall  -O2 -march=pentium2 -msse -mfpmath=sse -DNDEBUG
#CXXFLAGS=-Wall  -g -O2 -march=pentium4 -msse2 -msse -mfpmath=sse -pthread -DNDEBUG `wx-config --cxxflags`
#CXXFLAGS=-Wnon-virtual-dtor  -g -O2 -march=pentium4 -msse2 -msse -mfpmath=sse -pthread -DNDEBUG `wx-config --cxxflags`
LDFLAGS=-pthread `wx-config --libs`

INCLUDE = -I/nfs/hn00/cdtwigg/libtwigg -I/nfs/hn00/cdtwigg/cllux07ddl/include  -I/opt/intel/mkl61/include -I/nfs/hn00/cdtwigg/MagicSoftware/WildMagic2/Include -pthread -I/usr/include/freetype2 -I../ -I/nfs/hn00/cdtwigg/lib/wx/include/gtk2u-2.4 -I/nfs/hn00/cdtwigg/qslim-2.1/mixkit/src/ -I/nfs/hn00/cdtwigg/tristripper-1.1.0-beta-4/include
LIBDIR = -L/nfs/hn00/cdtwigg/lib -L/opt/intel/mkl61/lib/32 -L/usr/X11R6/lib -L/nfs/hn00/cdtwigg/tristripper-1.1.0-beta-4/lib

DEPDIR = .deps
df = $(DEPDIR)/$(*F)

#CFLAGS = -Wall -g -fno-inline -march=pentium4 -msse2 -msse -mfpmath=sse -DSAMPLE_SOLUTION -DBSP_ACCEL -Isrc/fileio/ -Isrc/parser -Isrc/scene -Isrc/SceneObjects/ -Isrc/ui -Isrc/vecmath -Isrc/ -I/afs/cs/academic/class/15864-s04/include

CXX = g++
CC = gcc


.SUFFIXES: .o .cpp

CXX_SRCS = blockmat.cpp boundingbox.cpp camera.cpp exception.cpp gslWrappers.cpp imagelib.cpp keyable.cpp linalg.cpp nagWrappers.cpp objfile.cpp pathUtils.cpp prt.cpp renderUtils.cpp stdafx.cpp TRAFile.cpp random.cpp ioutil.cpp interpolation.cpp util.cpp util.cpp volumeIntegrals.cpp
C_SRCS = tr.c EulerAngles.c
OBJS = $(CXX_SRCS:%.cpp=%.o) $(C_SRCS:%.c=%.o)

all: libtwigg.a libtwiggd.a

libtwiggd.a: $(CXX_SRCS:%.cpp=debug/%.o) $(C_SRCS:%.c=%.o)
	ar cru libtwiggd.a $(CXX_SRCS:%.cpp=debug/%.o) $(C_SRCS:%.c=%.o)

libtwigg.a: $(CXX_SRCS:%.cpp=release/%.o) $(C_SRCS:%.c=%.o)
	ar cru libtwigg.a $(CXX_SRCS:%.cpp=release/%.o) $(C_SRCS:%.c=%.o)

install: libtwigg.a libtwiggd.a
	cp -f libtwigg.a /nfs/hn00/cdtwigg/lib
	cp -f libtwiggd.a /nfs/hn00/cdtwigg/lib

tr.o: tr.c
	$(CC) -c $(CFLAGS) $(INCLUDE) -o $@ $<

EulerAngles.o: EulerAngles.c
	$(CC) -c $(CFLAGS) $(INCLUDE) -o $@ $<

release/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -MD -o $@ $<
	@cp release/$*.d $*.P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		-e '/^$$/ d' -e 's/$$/ :/' < release/$*.d >> $*.P; \
	rm -f release/$*.d

debug/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS_DEBUG) $(INCLUDE) -o $@ $<

-include $(CXX_SRCS:%.cpp=%.P) $(C_SRCS:%.c=%.P)

clean:
	rm -f $(C_SRCS:.c=.o) $(CXX_SRCS:%.cpp=release/%.o) $(CXX_SRCS:.cpp=.P) $(CXX_SRCS:%.cpp=debug/%.o)


