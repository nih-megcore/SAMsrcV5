include ../../config/Makefile.local

inst = $(LIBDIR)
bin = $(BDIR)

INSTFILES = nolteFwd.so

CFLAGS = -O3 -fPIC

targets = $(INSTFILES)

all: $(targets)
	@mkdir -p $(inst)
	cp $(INSTFILES) $(inst)

nolteFwd.so: nolteFwd.c
	env PYMODNAME=nolteFwd PYMODCFLAGS="-I../../include -fopenmp" PYMODLIBS="-L../../lib -lgsl -lgslcblas -llfu -lsam -lsig -lsamutil -lm -lgomp -lfftw3 -lgiftiio" $(MAKE) -f Makefile.pymod

clean:
	rm -f *~ *.so *.o *.pyc \#*
