include ../../config/Makefile.local

all: scripts

scripts:
	sed s.@@libdir@@.$(LIBDIR). < FSnormals.py > $(BDIR)/FSnormals.py
	sed s.@@libdir@@.$(LIBDIR). < mkGiiAtlas.py > $(BDIR)/mkGiiAtlas.py
	chmod +x $(BDIR)/FSnormals.py $(BDIR)/mkGiiAtlas.py

clean:
	rm -f *~ *.pyc
