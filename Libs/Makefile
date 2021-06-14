libs = LFULIB SAMLIB SIGLIB samutil

all clean:
	for x in $(libs) ; do \
		$(MAKE) -C $$x $@ || exit ;\
	done
