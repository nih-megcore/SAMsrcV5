include ../../config/Makefile.local

CFLAGS = -O2 -Wno-unused-result
LDFLAGS = -s

all: meshnorm scripts

OBJS = meshnorm.o err.o mem.o
AOBJS = ascnorm.o err.o mem.o

meshnorm: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) -lm
	cp meshnorm $(BDIR)/

scripts:
	sed s.@@libdir@@.$(LIBDIR). < orthohull > $(BDIR)/orthohull
	sed s.@@libdir@@.$(LIBDIR). < ply2fid > $(BDIR)/ply2fid
	sed s.@@libdir@@.$(LIBDIR). < orthohull.py > $(BDIR)/orthohull.py
	sed s.@@libdir@@.$(LIBDIR). < fiddist.py > $(BDIR)/fiddist.py
	chmod +x $(BDIR)/orthohull $(BDIR)/ply2fid $(BDIR)/fiddist.py $(BDIR)/orthohull.py
	cp plothull.py $(BDIR)
	cp fid.py samutil.py thd_atr.py null.tag master+orig.* $(LIBDIR)/

ascnorm: $(AOBJS)
	$(CC) $(LDFLAGS) -o $@ $(AOBJS) -lm

clean:
	rm -f *~ *.o meshnorm ascnorm
