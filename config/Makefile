include Makefile.local

DEST=$(HOME)/bin

make:
	@echo "Usage: 'make symlinks [DEST=/path/to/bin]' creates symlinks in" \
	"DEST that point to the executables in $(BDIR). DEST defaults" \
	"to $(HOME)/bin" | fmt

symlinks:
	mkdir -p $(DEST)
	@for x in $(BDIR)/*; do \
		y=`basename $$x` ; \
		echo ln -s -f $$x $(DEST)/$$y ; \
		ln -s -f $$x $(DEST)/$$y ; \
	done

clean:
	rm -f *~ Makefile.local
