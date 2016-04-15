FORT = gfortran
CC = gcc
OPTS = -O3 -march=native -ffast-math -funroll-loops -fstrict-aliasing -cpp -Wunused

#BLASINC = -I/opt/OpenBLAS/include
#BLASLIB = -L/opt/OpenBLAS/lib/ -lopenblas

#BLASINC = -I/usr/include/atlas/
#BLASLIB = -L/usr/lib/atlas-base/ -lcblas

MKLROOT = /opt/intel/mkl
BLASINC = -m64 -I${MKLROOT}/include
BLASLIB = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
          ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a \
          -Wl,--end-group -lpthread -lm -ldl

SRCDIR = ./src
OBJDIR = $(SRCDIR)/obj
BINDIR = ./test
LIBDIR = ./lib

LIBFILE = $(LIBDIR)/libemdee.a

OBJ = $(OBJDIR)/EmDee.o $(OBJDIR)/mEmDee.o

.PHONY: all test lib

all: test

clean:
	rm -rf $(OBJDIR)
	rm -rf $(LIBDIR)
	rm -rf $(BINDIR)

test: $(BINDIR)/testfortran

lib: $(LIBFILE)

$(BINDIR)/testfortran: $(OBJDIR)/testfortran.o
	mkdir -p $(BINDIR)
	$(FORT) $(OPTS) -o $@ -J$(LIBDIR) $< $(OBJDIR)/mRandom.o -L$(LIBDIR) -lemdee $(BLASLIB)

$(OBJDIR)/testfortran.o: $(SRCDIR)/testfortran.f90 $(OBJDIR)/mRandom.o $(LIBFILE)
	$(FORT) $(OPTS) -c -o $@ -J$(LIBDIR) $<

$(LIBFILE): $(OBJ)
	ar -cr $(LIBFILE) $(OBJ)

$(OBJDIR)/mRandom.o: $(SRCDIR)/mRandom.f90
	mkdir -p $(LIBDIR)
	$(FORT) $(FLAGS) -c -o $@ $< -J$(LIBDIR)

$(OBJDIR)/mEmDee.o: $(SRCDIR)/mEmDee.f90
	mkdir -p $(LIBDIR)
	$(FORT) $(FLAGS) -c -o $@ $< -J$(LIBDIR)

$(OBJDIR)/EmDee.o: $(SRCDIR)/EmDee.c $(SRCDIR)/EmDee.h
	mkdir -p $(OBJDIR)
	$(CC) $(OPTS) -c -o $@ $< $(BLASINC)

