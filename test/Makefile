include ../config/Makefile.local

TEST_BINDIR = $(CURDIR)/bin
TEST_LIBS = -L$(LIBDIR) $(SAMLIBS) $(SAMLIBS) -lgsl $(BLAS) -lfftw3 -lm -lgomp
TEST_CFLAGS = $(IFLAGS) -std=gnu17 -O2 -Wall -Wextra
PYTEST_ARGS ?=

.PHONY: unit integration-build clean

unit: $(TEST_BINDIR)/test_core
	$(TEST_BINDIR)/test_core
	python3 -m pytest -m unit -q $(PYTEST_ARGS) unit

integration-build: $(TEST_BINDIR)/inspect_ctf

$(TEST_BINDIR):
	mkdir -p $@

$(TEST_BINDIR)/test_core: unit/test_core.c | $(TEST_BINDIR)
	$(CC) $(TEST_CFLAGS) $< -o $@ $(TEST_LIBS)

$(TEST_BINDIR)/inspect_ctf: integration/inspect_ctf.c | $(TEST_BINDIR)
	$(CC) $(TEST_CFLAGS) $< -o $@ $(TEST_LIBS)

clean:
	rm -f $(TEST_BINDIR)/test_core $(TEST_BINDIR)/inspect_ctf
