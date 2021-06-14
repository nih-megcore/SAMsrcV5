include ../../config/Makefile.local

all: scripts

scripts:
	sed s.@@libdir@@.$(LIBDIR). < ROIbuilder > $(BDIR)/ROIbuilder
	chmod +x $(BDIR)/ROIbuilder
	cp ROIbuilder.ui SessionAtlases.niml $(LIBDIR)/

clean:
	rm -f *~
	rm -rf __pycache__
