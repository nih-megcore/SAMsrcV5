include ../../config/Makefile.local

all: 1dstats scripts

1dstats: 1dstats.c
	gcc -O2 -s -o 1dstats 1dstats.c -lm
	cp 1dstats $(BDIR)

scripts:
	cp 3dNormalize $(BDIR)

clean:
	rm -f *~ *.pyc 1dstats
