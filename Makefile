CONF := $(shell cd config; ./configure)

include config/Makefile.local

all:
	for x in Libs Mains; do \
		$(MAKE) -C $$x $@ || exit ; \
	done
ifndef AFNI
	@echo
	@echo "Warning! AFNI is not installed. Some functionality will be unavailable."
endif

clean:
	for x in Libs Mains config; do \
		$(MAKE) -C $$x clean || exit ; \
	done
	rm -f *~
