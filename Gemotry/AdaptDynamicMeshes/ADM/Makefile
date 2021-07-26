DIRS = src

all:
	@( echo "Usage: 'make install'" )

install: 
	@( \
	for dir in $(DIRS); do \
		cd $$dir; make install; cd ..; \
	done )

clean:
	@( \
	rm -f *~; \
	for dir in $(DIRS); do \
		cd $$dir; make clean; cd ..; \
	done ); \

uninstall: clean
	cd include; rm -f *.h; cd ..; \
	cd lib; rm -f *.a; cd ..
