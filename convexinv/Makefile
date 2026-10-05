MAKE=make
CC=gcc
LD=ld
AR=ar

OPTFLAGS=-O3
# BLAS for the batched shape derivatives (mrqcof.c): Accelerate on macOS, OpenBLAS elsewhere
UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
  BLASDEF=-DUSE_ACCELERATE
  BLASLIB=-framework Accelerate
else
  BLASDEF=
  BLASLIB=-lopenblas
endif
# gnu99: current compilers default to C23, which rejects the code's unprototyped function pointers
CFLAGS=$(OPTFLAGS) -Wall -std=gnu99 $(BLASDEF)
LDFLAGS=-lm $(BLASLIB)
MY_LIBS=trifac.o areanorm.o sphfunc.o ellfit.o covsrt.o ludcmp.o lubksb.o mrqmin.o mrqcof.o\
        curv.o blmatrix.o conv.o gauss.o phasec.o matrix.o bright.o memory.o\
	dot_product.o

all: convexinv period_scan

libs: $(MY_LIBS)

convexinv: convexinv.c $(MY_LIBS)
	$(CC) $(CFLAGS) -o $@ $< $(MY_LIBS) $(LDFLAGS)

period_scan: period_scan.c $(MY_LIBS)
	$(CC) $(CFLAGS) -o $@ $< $(MY_LIBS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o
	
