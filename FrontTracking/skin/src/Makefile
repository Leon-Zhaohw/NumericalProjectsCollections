#!smake 

#-----------------------------------------
#OpenVDB Directory -----------------------
openvdb = yes
OpenVDB_DIR = /opt
ifeq (yes,$(strip $(openvdb)))
	OpenVDB_FLAGS = -DOPENVDB
else
	OpenVDB_FLAGS =
endif
#-----------------------------------------
#Basic Stuff -----------------------------
CC          = `which c++` -std=c++0x 
cc          = `which c++` -std=c++0x

#-----------------------------------------
#Misc Flags ------------------------------
MISC_FLAGS = -Wall -Wno-c++11-extensions 

#-----------------------------------------
#Optimization ----------------------------
OPT   = -O3  -g

#-----------------------------------------
# Debug mode -----------------------------
DEBUG = 
LIB_MOD = 

TARGETS = particleskinner
#-----------------------------------------
#-----------------------------------------
OBJECTS = slVector.o slUtil.o main.o particleIO.o

ifeq (yes,$(strip $(openvdb)))
	OBJECTS += smoothingFilter.o
else	
	OBJECTS += smoothingGrid.o marchingCube.o kdTree.o slMatrix.o
endif

#-----------------------------------------

LIBS =  
INCS =  -I/usr/local/include/ -I/usr/include 

ifeq (yes,$(strip $(openvdb)))
	LIBS += -L/usr/local/lib -L/usr/lib/ -L$(OpenVDB_DIR)/lib -lopenvdb -lHalf -ltbb 
	INCS += -I$(OpenVDB_DIR)/include
endif

CCOPTS = $(OPT) $(MISC_FLAGS) $(DEBUG) $(INCS) $(OpenVDB_FLAGS)
LDOPTS = $(OPT) $(MISC_FLAGS) $(DEBUG) $(LIBS)

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS)

clean:
	/bin/rm -rf *.o particleskinner

#-----------------------------------------
particleskinner: $(OBJECTS)
	$(CC) $(LDOPTS) $(OBJECTS) $(LIBS) -o particleskinner

#-----------------------------------------
.C.o: 
	$(CC) $(CCOPTS) -c $< 

.cpp.o: 
	$(CC) $(CCOPTS) -c $< 

.c.o: 
	$(cc) $(CCOPTS) -c $< 

.o: $(OBJECTS)
	$(CC)  $(LDOPTS) $(OBJS) $(OBJECTS) $< $(LIBS) -o $@ 

.C: $(OBJECTS)
	$(CC) $(LDOPTS)   $(OBJS) $(OBJECTS) $< $(LIBS) $(FOR_LIB) -o $@

#-----------------------------------------
#-----------------------------------------














