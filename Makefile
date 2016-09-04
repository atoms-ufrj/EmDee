FORT = gfortran
CC   = gcc
OPTS = -march=native -ffast-math -fstrict-aliasing -Ofast -fPIC -m64 -fopenmp -Wunused -cpp

SRCDIR = ./src
OBJDIR = $(SRCDIR)/obj
BINDIR = ./test
LIBDIR = ./lib
INCDIR = ./include

LIBFILE = $(LIBDIR)/libemdee.so

LIBS = -L$(LIBDIR) -lemdee -lgfortran -lm

OBJ = $(OBJDIR)/EmDeeCode.o $(OBJDIR)/ArBee.o $(OBJDIR)/math.o $(OBJDIR)/structs.o  \
      $(OBJDIR)/models.o $(OBJDIR)/lists.o $(OBJDIR)/global.o

.PHONY: all test clean install lib

all: test

clean:
	rm -rf $(OBJDIR)
	rm -rf $(LIBDIR)
	rm -rf $(BINDIR)

install:
	cp $(LIBFILE) /usr/local/lib/
	cp $(INCDIR)/* /usr/local/include/

test: $(BINDIR)/testfortran $(BINDIR)/testc $(BINDIR)/testjulia

lib: $(LIBFILE)

$(BINDIR)/testfortran: $(OBJDIR)/testfortran.o
	mkdir -p $(BINDIR)
	$(FORT) $(OPTS) -o $@ -I$(INCDIR) -J$(OBJDIR) $< $(LIBS)

$(OBJDIR)/testfortran.o: $(SRCDIR)/testfortran.f90 $(INCDIR)/emdee.f03 $(LIBFILE)
	$(FORT) $(OPTS) -c -o $@ -J$(OBJDIR) $<

$(BINDIR)/testc: $(OBJDIR)/testc.o
	mkdir -p $(BINDIR)
	$(CC) $(OPTS) -o $@ $< $(LIBS)

$(OBJDIR)/testc.o: $(SRCDIR)/testc.c $(INCDIR)/emdee.h $(LIBFILE)
	$(CC) $(OPTS) -c -o $@ $< -I$(INCDIR)

$(BINDIR)/testjulia: $(SRCDIR)/testjulia.jl
	mkdir -p $(BINDIR)
	cp $< $@ && chmod +x $@

$(LIBFILE): $(OBJ)
	mkdir -p $(LIBDIR)
	$(FORT) -shared -fPIC -o $(LIBFILE) $^

$(OBJDIR)/EmDeeCode.o: $(SRCDIR)/EmDeeCode.f90 $(OBJDIR)/ArBee.o $(SRCDIR)/compute_pair.f90 \
                        $(SRCDIR)/compute_bond.f90 $(SRCDIR)/compute_angle.f90    \
                        $(SRCDIR)/compute_dihedral.f90 $(OBJDIR)/structs.o \
                        $(OBJDIR)/models.o $(OBJDIR)/lists.o $(OBJDIR)/global.o
	$(FORT) $(OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/ArBee.o: $(SRCDIR)/ArBee.f90 $(OBJDIR)/math.o $(OBJDIR)/global.o
	$(FORT) $(OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/math.o: $(SRCDIR)/math.f90 $(OBJDIR)/global.o
	$(FORT) $(OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/structs.o: $(SRCDIR)/structs.f90 $(OBJDIR)/models.o
	$(FORT) $(OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/models.o: $(SRCDIR)/models.f90 $(OBJDIR)/global.o
	$(FORT) $(OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/lists.o: $(SRCDIR)/lists.f90
	$(FORT) $(OPTS) -J$(OBJDIR) -c -o $@ $<

: $(SRCDIR)/c_binding.f90 $(OBJDIR)/global.o
	$(FORT) $(OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/global.o: $(SRCDIR)/global.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(OPTS) -J$(OBJDIR) -c -o $@ $<

