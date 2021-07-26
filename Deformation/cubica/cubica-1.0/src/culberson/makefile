
LDLIBS = -lm 
CC = gcc
CFLAGS =  -O2  -Wall -Dlinux

default: doc
	
BINPATH = .


GREEDY = gmain.o greedy.o graph.o colorrtns.o 
greedy: $(GREEDY) $(BINPATH)/greedy
	
$(BINPATH)/greedy: $(GREEDY)
	$(CC) $(GREEDY) -o $(BINPATH)/greedy $(LDLIBS)

$(GREEDY) : mysys.h greedy.h colorrtns.h graph.h makefile

DSATUR = dmain.o dsatur.o graph.o colorrtns.o
dsatur: $(DSATUR) $(BINPATH)/dsatur
	
$(BINPATH)/dsatur: $(DSATUR)
	$(CC) $(DSATUR) -o $(BINPATH)/dsatur $(LDLIBS)

$(DSATUR) : mysys.h dsatur.h colorrtns.h graph.h makefile

MAXIS = mmain.o maxis.o graph.o colorrtns.o
maxis: $(MAXIS) $(BINPATH)/maxis
	
$(BINPATH)/maxis: $(MAXIS)
	$(CC) $(MAXIS) -o $(BINPATH)/maxis $(LDLIB)

$(MAXIS) : mysys.h maxis.h colorrtns.h graph.h makefile

BKTDSAT = bkdmain.o bktdsat.o graph.o colorrtns.o maxclique.o
bktdsat: $(BKTDSAT) $(BINPATH)/bktdsat
	
$(BINPATH)/bktdsat: $(BKTDSAT)
	$(CC) $(BKTDSAT) -o $(BINPATH)/bktdsat $(LDLIBS)

$(BKTDSAT) : mysys.h bktdsat.h colorrtns.h graph.h maxclique.h makefile

ITRGRDY = igmain.o greedy.o graph.o colorrtns.o itrgrdy.o
itrgrdy: $(ITRGRDY) $(BINPATH)/itrgrdy
	
$(BINPATH)/itrgrdy: $(ITRGRDY)
	$(CC) $(ITRGRDY) -o $(BINPATH)/itrgrdy $(LDLIBS)

$(ITRGRDY) : mysys.h itrgrdy.h greedy.h colorrtns.h graph.h makefile


TABU = tbumain.o greedy.o graph.o colorrtns.o itrgrdy.o tabu.o partition.o brute.o
tabu: $(TABU) $(BINPATH)/tabu
	
$(BINPATH)/tabu: $(TABU)
	$(CC) $(TABU) -o $(BINPATH)/tabu $(LDLIBS)

$(TABU) : mysys.h itrgrdy.h greedy.h colorrtns.h graph.h tabu.h\
	partition.h brute.h makefile


ALLH = mysys.h greedy.h colorrtns.h graph.h dsatur.h maxis.h bktdsat.h\
	maxclique.h tabu.h brute.h partition.h itrgrdy.h makefile

all: $(ALLH) greedy dsatur maxis bktdsat itrgrdy tabu 

doc: README
	less README

clean:
	\rm *.o 
