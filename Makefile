DEBUG?=0
#make the "fast" version with: `make` or  #`make DEBUG=0`
#make the "debug" version with: #`make DEBUG=1`

FORT = gfortran
CC   = gcc

BASIC_F_OPTS = -march=native -m64 -fPIC -fopenmp -cpp -fmax-errors=1 -Wunused
BASIC_C_OPTS = -march=native -m64 -fPIC -fopenmp -cpp -fmax-errors=1 -Wunused

#Personal preference -Wall
#BASIC_F_OPTS += -Wall
#BASIC_C_OPTS += -Wall

#option FAST, (default)
FAST_F_OPTS = -Ofast
FAST_C_OPTS = -Ofast
# Ofast enables all -O3 optimizations.
# It turns on -ffast-math and
#   the Fortran-specific -fno-protect-parens and -fstack-arrays
# The -fstrict-aliasing option is enabled at level -O3

#option DEBUG
DEBUG_F_OPTS = -g -Og -fopenmp -fcheck=all -Ddebug
DEBUG_C_OPTS = -g -Og -fopenmp -fstack-check -fsanitize=null -fbounds-check -Ddebug

#checks chosen option

ifeq ($(DEBUG), 1)
F_OPTS = $(BASIC_F_OPTS) $(DEBUG_F_OPTS)
C_OPTS = $(BASIC_C_OPTS) $(DEBUG_C_OPTS)
else
F_OPTS = $(BASIC_F_OPTS) $(FAST_F_OPTS)
C_OPTS = $(BASIC_C_OPTS) $(FAST_C_OPTS)
endif

LIBS = -lgfortran -lm -lgomp

SRCDIR = ./src
OBJDIR = $(SRCDIR)/obj
BINDIR = ./test
LIBDIR = ./lib
INCDIR = ./include

EMDEELIB = -L$(LIBDIR) -lemdee

obj = $(addprefix $(OBJDIR)/, $(addsuffix .o, $(1)))
src = $(addprefix $(SRCDIR)/, $(addsuffix .f90, $(1)))

OBJECTS = $(call obj,EmDeeCode ArBee math structs models lists global)

.PHONY: all test clean install uninstall lib

.DEFAULT_GOAL := all

all: lib

clean:
	rm -rf $(OBJDIR)
	rm -rf $(LIBDIR)
	rm -rf $(BINDIR)

install:
	cp $(LIBDIR)/libemdee.* /usr/local/lib/
	cp $(INCDIR)/emdee.* /usr/local/include/
	ldconfig

uninstall:
	rm -f /usr/local/lib/libemdee.a /usr/local/lib/libemdee.so
	rm -f /usr/local/include/emdee.h /usr/local/include/emdee.f03
	ldconfig

# Executables:

test: $(addprefix $(BINDIR)/,testfortran testc testjulia)

$(BINDIR)/testfortran: $(SRCDIR)/testfortran.f90 $(INCDIR)/emdee.f03 $(LIBDIR)/libemdee.so
	mkdir -p $(BINDIR)
	$(FORT) $(F_OPTS) -static-libgfortran -o $@ -J$(OBJDIR) $< $(EMDEELIB)

$(BINDIR)/testc: $(SRCDIR)/testc.c $(INCDIR)/emdee.h $(LIBDIR)/libemdee.so
	mkdir -p $(BINDIR)
	$(CC) -static-libgfortran $(C_OPTS) -o $@ $< $(EMDEELIB) -lm

$(BINDIR)/testjulia: $(SRCDIR)/testjulia.jl
	mkdir -p $(BINDIR)
	cp $< $@ && chmod +x $@

# Static and shared libraries:

lib: $(LIBDIR)/libemdee.so $(LIBDIR)/libemdee.a

$(LIBDIR)/libemdee.so: $(OBJECTS)
	mkdir -p $(LIBDIR)
	$(FORT) -shared -fPIC -o $@ $^ $(LIBS)

$(LIBDIR)/libemdee.a: $(OBJECTS)
	mkdir -p $(LIBDIR)
	ar cr $@ $^

# Object files:

$(OBJDIR)/EmDeeCode.o: $(call src,EmDeeCode compute_pair compute_bond compute_angle compute_dihedral) \
                       $(call obj,ArBee structs models lists global)
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/ArBee.o: $(SRCDIR)/ArBee.f90 $(call obj,math global)
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/math.o: $(SRCDIR)/math.f90 $(OBJDIR)/global.o
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/structs.o: $(SRCDIR)/structs.f90 $(OBJDIR)/models.o
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/models.o: $(SRCDIR)/models.f90 $(OBJDIR)/global.o
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/lists.o: $(SRCDIR)/lists.f90
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/global.o: $(SRCDIR)/global.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<
