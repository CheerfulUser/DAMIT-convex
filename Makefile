# Build the binaries used by the TESSELLATE shape pipeline: convexinv and period_scan (convexinv/,
# BLAS-linked, see convexinv/Makefile) and minkowski (fortran/minkowski.f, shipped but never built
# upstream; the canonical facet-to-vertex reconstructor). Run `make` here.
UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
  MINK_LDFLAGS = -L$(shell xcrun --show-sdk-path 2>/dev/null)/usr/lib
endif

all: convexinv minkowski

convexinv:
	$(MAKE) -C convexinv convexinv period_scan

minkowski: fortran/minkowski.f
	gfortran -O2 $(MINK_LDFLAGS) -o minkowski fortran/minkowski.f

clean:
	rm -f convexinv/*.o convexinv/convexinv convexinv/period_scan minkowski

.PHONY: all convexinv clean
