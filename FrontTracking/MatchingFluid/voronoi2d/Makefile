
# --------------------
#
# Makefile for Talpa, and example application using El Topo.
#
# Create a file called Makefile.local_defs before running make.  See Makefile.example_defs for suggestions.
# You must also build the El Topo library before building Talpa.
#
# --------------------

# This is for GNU make; other versions of make may not run correctly.

# Local machine settings
include Makefile.local_defs

# --------------------
# Executable
# --------------------

MAIN_PROGRAM = voronoi2d

INCLUDE_PATH = -I. -Icommon -Inewsparse -Isparse -Itriangle

RELEASE_FLAGS += -DTRILIBRARY -DANSI_DECLARATORS
DEBUG_FLAGS += -DTRILIBRARY -DANSI_DECLARATORS

SRC = accelerationgrid.cpp eltopo.cpp pressure_vor2d.cpp \
      broadphasegrid.cpp   framestepper.cpp     sampleseeder.cpp \
      dist_funcs.cpp       geom_routines.cpp    scenes.cpp \
      dualfluidsim.cpp     interpolator.cpp     surftrack.cpp \
      dynamicsurface.cpp   main.cpp             triangle_wrapper.cpp \
      edgemesh.cpp         meshlevelset.cpp     trimesh2d.cpp \
      common/collisionqueries.cpp   common/openglutils.cpp \
      common/cubic_ccd_wrapper.cpp  common/predicates.cpp \
      common/gluvi.cpp              common/wallclocktime.cpp \
      common/marching_triangles.cpp  newsparse/sparse_matrix.cpp \
      newsparse/krylov_solvers.cpp   

CSRC = triangle/triangle.c 

# object files
C_RELEASE_OBJ = $(patsubst %.c,cobj/%.o,$(notdir $(CSRC)))
RELEASE_OBJ = $(patsubst %.cpp,obj/%.o,$(notdir $(SRC)))
C_DEBUG_OBJ = $(patsubst %.c,cobj_debug/%.o,$(notdir $(CSRC)))
DEBUG_OBJ = $(patsubst %.cpp,obj_debug/%.o,$(notdir $(SRC)))

# how to make the main target (debug mode)
$(MAIN_PROGRAM): $(C_DEBUG_OBJ) $(DEBUG_OBJ)
	$(LINK) $(DEBUG_LINKFLAGS) -o $@ $^ $(LINK_LIBS) 

# how to make the main target (release mode)
$(MAIN_PROGRAM)_release: $(C_RELEASE_OBJ) $(RELEASE_OBJ)
	$(LINK) $(RELEASE_LINKFLAGS) -o $@ $^ $(LINK_LIBS) 

.PHONY: release
release: $(MAIN_PROGRAM)_release

.PHONY: debug
debug: $(MAIN_PROGRAM)

# how to compile each file
.SUFFIXES:
obj/%.o:
	$(CCPP) -c $(RELEASE_FLAGS) $(INCLUDE_PATH) -o $@ $<
cobj/%.o:
	$(CC) -c $(RELEASE_FLAGS) $(INCLUDE_PATH) -o $@ $< 
obj_debug/%.o:
	$(CCPP) -c $(DEBUG_FLAGS) $(INCLUDE_PATH) -o $@ $<
cobj_debug/%.o:
	$(CC) -c $(DEBUG_FLAGS) $(INCLUDE_PATH) -o $@ $< 

# cleaning up
.PHONY: clean
clean:
	-rm -f cobj/*.o cobj/depend obj/*.o obj/depend $(MAIN_PROGRAM) cobj_debug/*.o cobj_debug/depend obj_debug/*.o obj_debug/depend $(MAIN_PROGRAM)_release *core
	
# dependencies are automatically generated
.PHONY: depend
depend:
	-mkdir obj
	-rm -f obj/depend
	$(foreach srcfile,$(SRC),$(DEPEND) -MM $(srcfile) -MT $(patsubst %.cpp,obj/%.o,$(notdir $(srcfile))) $(INCLUDE_PATH) >> obj/depend;)
	-mkdir cobj
	-rm -f cobj/depend
	$(foreach srcfile,$(CSRC),$(CDEPEND) -MM $(srcfile) -MT $(patsubst %.c,cobj/%.o,$(notdir $(srcfile))) $(INCLUDE_PATH) >> cobj/depend;)
	-mkdir obj_debug
	-rm -f obj_debug/depend
	$(foreach srcfile,$(SRC),$(DEPEND) -MM $(srcfile) -MT $(patsubst %.cpp,obj_debug/%.o,$(notdir $(srcfile))) $(INCLUDE_PATH) >> obj_debug/depend;)
	-mkdir cobj_debug
	-rm -f cobj_debug/depend
	$(foreach srcfile,$(CSRC),$(CDEPEND) -MM $(srcfile) -MT $(patsubst %.c,cobj_debug/%.o,$(notdir $(srcfile))) $(INCLUDE_PATH) >> cobj_debug/depend;)

-include obj/depend
-include cobj/depend
-include obj_debug/depend
-include cobj_debug/depend
