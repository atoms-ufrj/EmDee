FORT = gfortran
CC   = gcc
OPTS = -march=native -ffast-math -fstrict-aliasing -Ofast -fPIC -m64 -fopenmp -Wunused -cpp

SRCDIR = ./src
OBJDIR = $(SRCDIR)/obj
BINDIR = ./test
LIBDIR = ./lib

LIBFILE = $(LIBDIR)/libemdee.so

LIBS = -L$(LIBDIR) -lemdee -lgfortran -lm

OBJ = $(OBJDIR)/EmDee_code.o $(OBJDIR)/bonded_structs.o $(OBJDIR)/models.o \
      $(OBJDIR)/lists.o $(OBJDIR)/c_binding.o

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
	$(FORT) -shared -fPIC -o $(LIBFILE) $^

$(OBJDIR)/EmDee_code.o: $(SRCDIR)/EmDee_code.f90 $(SRCDIR)/compute_pair.f90 \
                        $(SRCDIR)/compute_bond.f90 \
                        $(OBJDIR)/bonded_structs.o $(OBJDIR)/models.o $(OBJDIR)/lists.o
	$(FORT) $(OPTS) -J$(LIBDIR) -c -o $@ $<

$(OBJDIR)/bonded_structs.o: $(SRCDIR)/bonded_structs.f90 $(OBJDIR)/models.o
	$(FORT) $(OPTS) -J$(LIBDIR) -c -o $@ $<

$(OBJDIR)/models.o: $(SRCDIR)/models.f90 $(OBJDIR)/c_binding.o
	$(FORT) $(OPTS) -J$(LIBDIR) -c -o $@ $<

$(OBJDIR)/lists.o: $(SRCDIR)/lists.f90 $(OBJDIR)/c_binding.o
	$(FORT) $(OPTS) -J$(LIBDIR) -c -o $@ $<

$(OBJDIR)/c_binding.o: $(SRCDIR)/c_binding.f90
	mkdir -p $(OBJDIR)
	mkdir -p $(LIBDIR)
	$(FORT) $(OPTS) -J$(LIBDIR) -c -o $@ $<

