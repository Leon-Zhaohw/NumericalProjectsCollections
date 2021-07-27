# Minimal LBM demonstration program by Nils Thuerey 
# See http://www.ntoken.com and http://www.matthiasmueller.info/realtimephysics/ for further info
# The code is released under the GPL, see LICENSE.txt for details

EXE=synthturb
SOURCES = $(wildcard source/*.cpp util/*.cpp)
HEADERS = $(wildcard source/*.h util/*.h)
OBJCMD = $(patsubst %.cpp,build_cmd/%.o,$(SOURCES))
OBJGUI = $(patsubst %.cpp,build_gui/%.o,$(SOURCES))

CXX=g++
CFLAGS     = -c -O3 -Isource -Iutil -DDDF_OPENMP=0
LDFLAGS    = -lz

# comment out the following two lines for disabling OPENMP
LDFLAGS    = -lz -fopenmp -lgomp
CFLAGS     = -c -fopenmp -O3 -Isource -Iutil -DDDF_OPENMP=1
# uncomment this line for macs
#LDFLAGS    = -lz -framework OpenGL -framework GLUT -lgomp

DEFGUI     = -DDDF_GLUTGUI=1
DEFCMD     = -DDDF_GLUTGUI=0

default: all

all: cmd gui

cmd: $(OBJCMD)	 
	$(CXX) $(LDFLAGS) -o turbCmd $(OBJCMD)

gui: $(OBJGUI)
	$(CXX) $(LDFLAGS) -lglut -o turbGui $(OBJGUI)

build_cmd/%.o : %.cpp $(HEADERS)
	$(CXX) $(CFLAGS) $(DEFCMD) -o $@ $<

build_gui/%.o : %.cpp $(HEADERS)
	$(CXX) $(CFLAGS) $(DEFGUI) -o $@ $<

clean:
	 rm -f *.o build_gui/source/*.o build_gui/util/*.o build_cmd/source/*.o build_gui/util/*.o turbCmd turbGui

