FORT = gfortran
CC   = gcc
OPTS = -march=native -ffast-math -fstrict-aliasing -Ofast -fPIC -fopenmp -Wunused -cpp

SRCDIR = ./src
OBJDIR = $(SRCDIR)/obj
BINDIR = ./test
LIBDIR = ./lib

LIBFILE = $(LIBDIR)/libemdee.a

LIBS = -L$(LIBDIR) -lemdee -lgfortran -lm

OBJ = $(OBJDIR)/EmDee.o $(OBJDIR)/c_binding_extra.o

.PHONY: all test lib testc testfortran

all: test

clean:
	rm -rf $(OBJDIR)
	rm -rf $(LIBDIR)
	rm -f $(BINDIR)/testfortran $(BINDIR)/testc

test: testfortran testc

testc: $(BINDIR)/testc

testfortran: $(BINDIR)/testfortran

lib: $(LIBFILE)

$(BINDIR)/testfortran: $(OBJDIR)/testfortran.o
	$(FORT) $(OPTS) -o $@ -J$(LIBDIR) $< $(LIBS)

$(OBJDIR)/testfortran.o: $(SRCDIR)/testfortran.f90 $(LIBFILE)
	$(FORT) $(OPTS) -c -o $@ -J$(LIBDIR) $<

$(BINDIR)/testc: $(OBJDIR)/testc.o
	$(CC) $(OPTS) -o $@ $< $(LIBS)

$(OBJDIR)/testc.o: $(SRCDIR)/testc.c $(SRCDIR)/EmDee.h $(LIBFILE)
	$(CC) $(OPTS) -c -o $@ $<

$(LIBFILE): $(OBJ)
	ar -cr $(LIBFILE) $(OBJ)

$(OBJDIR)/EmDee.o: $(SRCDIR)/EmDee.f90 $(OBJDIR)/c_binding_extra.o \
                   $(SRCDIR)/pair_compute.f90 $(SRCDIR)/pair_setup.f90
	$(FORT) $(OPTS) -J$(LIBDIR) -c -o $@ $<

$(OBJDIR)/c_binding_extra.o: $(SRCDIR)/c_binding_extra.f90
	mkdir -p $(OBJDIR)
	mkdir -p $(LIBDIR)
	$(FORT) $(OPTS) -J$(LIBDIR) -c -o $@ $<

