include ../../config/Makefile.local

all: scripts

scripts:
	sed 's|@@libdir@@|$(LIBDIR)|g' < FSnormals.py > $(BDIR)/FSnormals.py
	sed 's|@@libdir@@|$(LIBDIR)|g' < mkGiiAtlas.py > $(BDIR)/mkGiiAtlas.py
	chmod +x $(BDIR)/FSnormals.py $(BDIR)/mkGiiAtlas.py

clean:
	rm -f *~ *.pyc
