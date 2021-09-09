include ../config/Makefile.local

#sprogs = SAMNwts SAMhfo SAMersc MUTimg NIFTIPeak SAMvs fastica SAMsim KurtSensor

sprogs = SAMcov SAMwts SAM3d SAM4d SAMpower OPMsim #SAMepi SAMers SAM3dc SAM4dc

mprogs = orthohull

scripts = 3dNormalize FSnormals ROIbuilder

all clean:
	# Use generic Makefile.prog for simple progs
	for x in $(sprogs) ; do \
		env PROG=$$x $(MAKE) -C $$x $@ -f ../Makefile.prog || exit ; \
	done
	# Use Makefile in prog dir
	for x in $(mprogs) $(scripts) ; do \
		env PROG=$$x $(MAKE) -C $$x $@ || exit ; \
	done
	rm -f *~
