CC          = c++ 

#-----------------------------------------
#Optimization ----------------------------
OPT   =  -O3 -g

#-----------------------------------------

TARGETS = springmass explicitfem implicitfem fluid rigidbody partviewer

OBJECTS = jsoncpp.o 

#-----------------------------------------

CCOPTS = $(OPT) -I/usr/local/include/eigen3/
LDOPTS = 

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS)

clean:
	/bin/rm -f *.o $(TARGETS)

#-----------------------------------------
#-----------------------------------------

springmass: $(OBJECTS) springmass.o
	$(CC) $(OBJECTS) $(LDOPTS) springmass.o -o springmass

explicitfem: $(OBJECTS) explicitfem.o
	$(CC) $(OBJECTS) $(LDOPTS) explicitfem.o -o explicitfem

implicitfem: $(OBJECTS) implicitfem.o
	$(CC) $(OBJECTS) $(LDOPTS) implicitfem.o -o implicitfem

fluid: $(OBJECTS) fluid.o
	$(CC) $(OBJECTS) $(LDOPTS) fluid.o -o fluid

rigidbody: $(OBJECTS) rigidbody.o
	$(CC) $(OBJECTS) $(LDOPTS) rigidbody.o -o rigidbody

partviewer: $(OBJECTS) partviewer.o
	$(CC) $(OBJECTS) $(LDOPTS) -framework OpenGL -framework GLUT -framework foundation partviewer.o -o partviewer

#-----------------------------------------
#-----------------------------------------

.cpp.o: 
	$(CC) $(CCOPTS) -c $< 

#-----------------------------------------
#-----------------------------------------















