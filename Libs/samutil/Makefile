include ../../config/Makefile.local

LIB = libsamutil
SRC = $(shell ls *.c)
OBJ = $(SRC:.c=.o)

include ../Makefile.lib

clean:
	rm -f *~ *.o *.so *.il *.a core $(LIBDIR)/$(LIB).*
