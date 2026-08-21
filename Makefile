CONF := $(shell cd config; ./configure)

include config/Makefile.local

.PHONY: all clean test test-unit test-integration test-slow fetch-test-data rebuild_install_wheel

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

rebuild_install_wheel:
	@set -euo pipefail; \
	version="$$(awk -F '"' '/^version[[:space:]]*=/ { print $$2; exit }' pyproject.toml)"; \
	if [[ -z "$$version" ]]; then \
		echo "Could not determine the samsrcv5 version from pyproject.toml" >&2; \
		exit 1; \
	fi; \
	python_executable="$$(python -c 'import sys; print(sys.executable)')"; \
	echo "Using Python: $$python_executable"; \
	echo "Current samsrcv5 version: $$version"; \
	shopt -s nullglob; \
	wheels=(dist/samsrcv5-"$$version"-*.whl); \
	if (( $${#wheels[@]} > 0 )); then \
		echo "Current-version wheel(s) already exist:"; \
		printf '  %s\n' "$${wheels[@]}"; \
		printf 'Remove, rebuild, and reinstall? [y/N] '; \
		if ! IFS= read -r answer; then \
			echo; \
			echo "No confirmation received; aborting." >&2; \
			exit 1; \
		fi; \
		case "$$answer" in \
			y|Y|yes|YES|Yes) ;; \
			*) echo "Cancelled; existing wheel(s) were not changed."; exit 0 ;; \
		esac; \
		echo "Removing current-version wheel(s)..."; \
		rm -f -- "$${wheels[@]}"; \
	fi; \
	echo "Building samsrcv5 $$version wheel..."; \
	python -m build --wheel; \
	wheels=(dist/samsrcv5-"$$version"-*.whl); \
	if (( $${#wheels[@]} != 1 )); then \
		echo "Expected exactly one wheel for samsrcv5 $$version, found $${#wheels[@]}." >&2; \
		printf '  %s\n' "$${wheels[@]}" >&2; \
		exit 1; \
	fi; \
	echo "Installing $${wheels[0]} into $$python_executable..."; \
	python -m pip install --force-reinstall "$${wheels[0]}"
