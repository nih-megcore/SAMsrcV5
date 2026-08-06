CONF := $(shell cd config; ./configure)

include config/Makefile.local

.PHONY: all clean test test-unit test-integration test-slow fetch-test-data

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

test:
	$(MAKE) test-unit
	$(MAKE) test-integration

test-unit: all
	$(MAKE) -C test unit

fetch-test-data:
	./test/fetch-test-data.sh

test-integration:
	./test/run-container-tests.sh

test-slow:
	SAM_PYTEST_MARKER=slow SAM_TEST_REPORT=slow-junit.xml ./test/run-container-tests.sh
