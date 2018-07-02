export OPENGL=1
include ../../src/Makefile.defs
LIB+=-L/usr/X11/lib
OPT+=-I/usr/X11/include

MYSOURCES=../../src_spring/spring.c ../../src_spring/kepcart.c ../../src_spring/heat.c ../../src_spring/orb.c ../../src_spring/m_output.c

all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I../../src/ -I./ -I../../src_spring/ -Wl,-rpath,./ $(OPT) $(PREDEF) $(MYSOURCES) problem.c  -L. -lrebound $(LIB) -o rebound_spring
	@echo ""
	@echo "REBOUND-Spring compiled successfully."

librebound: 
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C ../../src/ 
	@-rm -f librebound.so
	@ln -s ../../src/librebound.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C ../../src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
